"""Replay the paired finite-panel fit using only the supplied input tree."""
import argparse
import csv
import hashlib
import json
import platform
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import scipy
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation
from scipy.stats import qmc

from structures import coords, frames, in_core, native, write_snapshot


BRANCHES = ("B-on", "B-near", "B-off", "ref-B-on")
GATES = {"primary": (2.6, 2.2, 2.6), "strict": (3., 2.6, 3.)}


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


class Panel:
    def __init__(self, dna, sugar, seed=6280):
        draws = np.vstack([np.zeros(6), 2*qmc.Sobol(6, scramble=True, seed=seed).random_base2(13)-1])
        self.rot = Rotation.from_rotvec(np.deg2rad(draws[:,3:]*8)).as_matrix()
        self.shift = dna.mean(0) - dna.mean(0) @ self.rot + 2.5*draws[:,:3]
        self.dna = dna @ self.rot + self.shift[:,None,:]
        delta = self.dna-dna
        self.rms = np.sqrt(np.square(delta[:,sugar]).sum(2).mean(1))
        self.max_move = np.sqrt(np.square(delta).sum(2).max(0))
        self.sugar = sugar
        ids = np.arange(len(self.rms))
        self.parts = dict(union=ids>=0, first=ids<=4096, second=(ids==0)|(ids>4096))

    def contacts(self, frame):
        dist = cdist(frame["ref_core"], self.dna[0,self.sugar])
        ia, js = np.nonzero(dist <= 10)
        fraction = np.zeros(len(self.rms))
        idx = np.flatnonzero(self.rms<=3)
        for chunk in np.array_split(idx, max(1,(len(idx)+255)//256)):
            delta = self.dna[chunk[:,None], self.sugar[js]] - frame["core"][ia]
            fraction[chunk] = (np.square(delta).sum(2) <= (dist[ia,js]+2)**2).mean(1)
        return fraction, len(ia)

    def clear(self, xyz, eligible, threshold):
        """Exact threshold predicate; bound candidates by maximal DNA displacement."""
        keep = eligible.copy()
        initial = self.dna[0]
        neighbors = cKDTree(xyz).query_ball_point(initial, self.max_move+threshold+1e-9)
        dna_ix = np.repeat(np.arange(len(initial)), [len(group) for group in neighbors])
        if not len(dna_ix):
            return keep
        mol_ix = np.concatenate([np.asarray(group,dtype=int) for group in neighbors])
        order = np.square(initial[dna_ix]-xyz[mol_ix]).sum(1).argsort()
        for block in range(0, len(order), 32):
            live = np.flatnonzero(keep)
            if not len(live):
                break
            pairs = order[block:block+32]
            delta = self.dna[live[:,None],dna_ix[pairs]] - xyz[mol_ix[pairs]]
            keep[live] = (np.einsum("ijk,ijk->ij",delta,delta) >= threshold**2).all(1)
        return keep

    def evaluate(self, frame):
        recovery, n_pairs = self.contacts(frame)
        eligible = (recovery>=.7) & (self.rms<=3.)
        result = {}
        for name,(t,l,p) in GATES.items():
            aok = self.clear(frame["A"], eligible, t)
            base = self.clear(frame["L"], aok, l)
            full = self.clear(frame["P"], base, p)
            assert np.all(~full | base)
            if name == "strict":
                for new, old in zip((base,full), result["primary"], strict=True):
                    assert np.all(~new | old)
            result[name] = base, full
        return result, recovery, n_pairs

    def best(self, mask):
        ids = np.flatnonzero(mask)
        return int(ids[np.argmin(self.rms[ids])]) if len(ids) else None


def gap(component, dna):
    return float(cKDTree(component).query(dna)[0].min())


def control(panel, ref):
    keys = [k for k in ref if in_core(k)]
    core = coords([ref[k] for k in keys])
    fractions, pairs = panel.contacts(dict(ref_core=core,core=core))
    eligible = (panel.rms<=3) & (fractions>=.7)
    heavy = coords([v for k,v in ref.items() if k[0]=="A" and v.startswith("ATOM  ")])
    rows = []
    for name, (threshold,_,_) in GATES.items():
        valid = panel.clear(heavy, eligible, threshold)
        best = panel.best(valid)
        # Independent KD-tree confirmation on a deterministic held-out subset.
        sample = np.unique(np.r_[np.arange(0,len(eligible),67), best if best is not None else 0])
        for i in sample:
            assert bool(valid[i]) == bool(eligible[i] and gap(heavy,panel.dna[i])>=threshold)
        rows.append(dict(gates=name,proposals=len(valid),receptor_accepted=int(valid.sum()),
                         native_contact_pairs=pairs,identity_gap_A=gap(heavy,panel.dna[0]),
                         identity_contacts=fractions[0],identity_accepted=int(valid[0]),
                         nearest_accepted_id=best,nearest_accepted_rmsd_A=panel.rms[best] if best is not None else ""))
    assert rows[0]["receptor_accepted"]>0 and rows[0]["identity_contacts"]==1
    return rows


def aggregate(rows):
    counts = defaultdict(lambda: np.zeros(5,int))
    for row in rows:
        base,full = row["baseline_accepted"]>0,row["full_accepted"]>0
        for rep in (row["replica"],"all"):
            for subset in (["all","dual_head"] if row["dual_head"] else ["all"]):
                k=(row["branch"],rep,subset,row["gates"],row["panel"])
                counts[k] += np.array([1,base,full,base and not full,not base],int)
    return [dict(branch=k[0],replica=k[1],subset=k[2],gates=k[3],panel=k[4],
                 **dict(zip(("frames","baseline_pass","full_pass","panel_loss","baseline_failure"),v)))
            for k,v in counts.items()]


def run(root, out):
    started=time.perf_counter()
    out.mkdir(exist_ok=True,parents=True)
    structures=out/"representatives"
    structures.mkdir(exist_ok=True)
    ref,dna_lines,dna,sugar=native(root)
    panel=Panel(dna,sugar)
    crystal_control=control(panel,ref)
    write_csv(out/"native-control.csv",crystal_control)
    np.savez_compressed(out/"transforms.npz",R=panel.rot,t=panel.shift,interior_C1_RMSD_A=panel.rms)
    with (root/"independent-identity-gaps.csv").open() as stream:
        identity={(r["start"].removeprefix("ref-") if r["family"]=="monomer" else "ref-B-on",
                   int(r["replica"]),float(r["time_ps"])):r for r in csv.DictReader(stream)}
    rows, pose_rows, examples = [], [], {}
    checked=0
    for branch in BRANCHES:
        for f in frames(root,branch,ref):
            key = dict(branch=branch,replica=f["replica"],time_ps=f["time_ps"])
            x,y=f["heads"]
            dual=int(max(x,y)<=2)
            result,recovery,n_pairs=panel.evaluate(f)
            trees=[cKDTree(f[c]) for c in ("A","L","P")]
            audits=identity[tuple(key.values())]
            for tree, label in zip(trees,("TBXT","probe","CypA")):
                gap_identity=tree.query(dna)[0].min()
                assert abs(gap_identity-float(audits[label+"_whole_min_A"]))<3e-5, (key,label,gap_identity)
            checked += 1
            for name,(base,full) in result.items():
                # Verify all classifications in an independent, deterministic sample using exact distances.
                for i in sorted({0,(checked*113)%8193,panel.best(base),panel.best(full)} - {None}):
                    ga,gl,gp=[t.query(panel.dna[i])[0].min() for t in trees]
                    a,l,p=GATES[name]
                    truth=recovery[i]>=.7 and panel.rms[i]<=3 and ga>=a and gl>=l
                    assert bool(base[i])==truth and bool(full[i])==(truth and gp>=p)
                for part,mask in panel.parts.items():
                    bmask, fmask=base&mask,full&mask
                    nb,nf=int(bmask.sum()),int(fmask.sum())
                    status="retained_fit" if nf else ("panel_loss" if nb else "baseline_failure")
                    ib,ifull=panel.best(bmask),panel.best(fmask)
                    rows.append(dict(**key,dual_head=dual,TBXT_head_A=x,CypA_head_A=y,
                                     gates=name,panel=part,proposals=int(mask.sum()),native_contact_pairs=n_pairs,
                                     baseline_accepted=nb,full_accepted=nf,status=status,
                                     best_baseline_id=ib,best_full_id=ifull))
                    if part != "union":
                        continue
                    for condition,i in [("reference",0),("closest_baseline",ib),("closest_full",ifull)]:
                        if i is None:
                            continue
                        ga,gl,gp = [float(t.query(panel.dna[i])[0].min()) for t in trees]
                        pose_rows.append(dict(**key,gates=name,condition=condition,proposal=i,
                                     sugar_RMSD_A=panel.rms[i],contact_recovery=recovery[i],
                                     TBXT_gap_A=ga,probe_gap_A=gl,CypA_gap_A=gp))
                    if name=="primary":
                        preference=(1-dual,f["replica"],f["time_ps"])
                        slot=(branch,status)
                        if slot not in examples or preference<examples[slot][0]:
                            i=ifull if ifull is not None else ib if ib is not None else 0
                            filename=f"{branch}-{status}.pdb"
                            notes=[f"{branch} rep {f['replica']} {f['time_ps']} ps; DNA proposal {i}",
                                   f"primary union: {status}; both crystal head RMSDs <=2 A: {dual}",
                                   "DNA is a fitted overlay; assembly atoms are from the simulated frame."]
                            write_snapshot(structures/filename,f,dna_lines,panel.dna[i],notes)
                            examples[slot]=(preference,dict(file="representatives/"+filename,**key,
                                               proposal=i,status=status,dual_head=dual))
            for base,full in result.values():
                for m in (base,full):
                    assert m.sum()==(m&panel.parts["first"]).sum()+(m&panel.parts["second"]).sum()-int(m[0])
        print(f"{branch} complete ({time.perf_counter()-started:.1f} s)",flush=True)
    assert checked==200 and len(identity)==200 and len(rows)==1200
    write_csv(out/"frames.csv",rows)
    write_csv(out/"selected-poses.csv",pose_rows)
    write_csv(out/"summary.csv",aggregate(rows))
    write_csv(out/"representatives.csv",[v[1] for v in examples.values()])
    manifest=dict(frames=checked,elapsed_s=time.perf_counter()-started,
                  python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                  inputs={p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in sorted(root.rglob("*")) if p.is_file()})
    (out/"replay.json").write_text(json.dumps(manifest,indent=2)+"\n")


if __name__=="__main__":
    directory=Path(__file__).resolve().parent
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs",type=Path,default=directory/"inputs")
    parser.add_argument("--out",type=Path,default=directory)
    args=parser.parse_args()
    run(args.inputs,args.out)
