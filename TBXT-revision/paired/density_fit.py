"""Evaluate the predeclared 32,769-pose sensitivity and compare with the primary panel."""
import argparse
import csv
import json
import time
from collections import Counter
from pathlib import Path

import numpy as np

from paired_fit import BRANCHES, Panel, aggregate, control, gap, write_csv
from structures import frames, native, write_snapshot


SEEDS = range(6280, 6284)


def status(base, full):
    return "retained_fit" if full else "panel_loss" if base else "baseline_failure"


def read_csv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def run(inputs, primary, out):
    started = time.perf_counter()
    out.mkdir(exist_ok=True,parents=True)
    (out/"representatives").mkdir(exist_ok=True)
    primary_rows = {(r["branch"],int(r["replica"]),float(r["time_ps"])):r
                    for r in read_csv(primary/"frames.csv")
                    if r["panel"]=="union" and r["gates"]=="primary"}
    source,lines,dna,sugar = native(inputs)
    per_seed, native_rows, totals, best, transforms = [], [], {}, {}, {}
    for block, seed in enumerate(SEEDS):
        panel = Panel(dna,sugar,seed)
        # Store the identity once under global ID 0 and 8,192 new rows per block.
        transforms[seed] = panel.rot.copy(), panel.shift.copy()
        test = control(panel,source)[0]
        native_rows.append(dict(seed=seed,**test))
        for branch in BRANCHES:
            for f in frames(inputs, branch, source):
                key = branch,f["replica"],f["time_ps"]
                info = dict(branch=key[0],replica=key[1],time_ps=key[2])
                row = primary_rows[key]
                recovery,n_pairs = panel.contacts(f)
                target_ok = panel.clear(f["A"], (recovery>=.7)&(panel.rms<=3.), 2.6)
                base = panel.clear(f["L"], target_ok, 2.2)
                full = panel.clear(f["P"], base, 2.6)
                assert np.all(~full | base) and n_pairs==int(row["native_contact_pairs"])
                counts = np.array([base.sum(),full.sum()],dtype=int)
                id_count = np.array([base[0],full[0]],dtype=int)
                if block==0:
                    assert counts.tolist()==[int(row["baseline_accepted"]),int(row["full_accepted"])]
                    totals[key]=counts.copy()
                else:
                    totals[key] += counts-id_count
                assert (totals[key]>=counts).all()
                per_seed.append(dict(**info,seed=seed,proposals=len(base),baseline_accepted=counts[0],
                                     full_accepted=counts[1],identity_baseline=id_count[0],identity_full=id_count[1]))
                for name,mask in (("closest_baseline",base),("closest_full",full)):
                    i=panel.best(mask)
                    if i is None:
                        continue
                    gid=0 if i==0 else block*8192+i
                    rank=float(panel.rms[i]),gid
                    if (key,name) in best and best[key,name][0]<=rank:
                        continue
                    ga,gl,gp=[gap(f[c],panel.dna[i]) for c in ("A","L","P")]
                    assert ga>=2.6 and gl>=2.2 and (name!="closest_full" or gp>=2.6)
                    record=dict(**info,condition=name,proposal=gid,seed=seed,local_proposal=i,
                                sugar_RMSD_A=rank[0],contact_recovery=recovery[i],
                                TBXT_gap_A=ga,probe_gap_A=gl,CypA_gap_A=gp)
                    best[key,name]=rank,record
        print(f"seed {seed} completed ({time.perf_counter()-started:.1f} s)",flush=True)
    rows, transitions, reps = [], Counter(), {}
    for key,(b,f) in totals.items():
        ref = primary_rows[key]
        before,after=ref["status"],status(b,f)
        old_b,old_f=int(ref["baseline_accepted"]),int(ref["full_accepted"])
        assert b>=old_b and f>=old_f and f<=b
        dual=int(ref["dual_head"])
        record=dict(branch=key[0],replica=key[1],time_ps=key[2],dual_head=dual,gates="primary",panel="union",
                    proposals=32769,baseline_accepted=b,full_accepted=f,status=after,
                    original_baseline_accepted=old_b,original_full_accepted=old_f,original_status=before,
                    status_changed=int(before!=after))
        for name in ("closest_baseline","closest_full"):
            record[name+"_id"]=best[key,name][1]["proposal"] if (key,name) in best else ""
        rows.append(record)
        for subset in (["all","dual_head"] if dual else ["all"]):
            transitions[(key[0],subset,before,after)]+=1
        if before != after:
            slot=(key[0],before,after)
            priority=(1-dual,key[1],key[2])
            if slot not in reps or priority<reps[slot][0]:
                reps[slot]=priority,key
    exports=[]
    for slot,(_,key) in reps.items():
        branch,before,after=slot
        f=next(f for f in frames(inputs,branch,source) if (branch,f["replica"],f["time_ps"])==key)
        name="closest_full" if after=="retained_fit" else "closest_baseline"
        pose=best[key,name][1]
        r,t=transforms[pose["seed"]]
        i=pose["local_proposal"]
        filename=f"representatives/{branch}-{before}-to-{after}.pdb"
        notes=[f"{branch} rep {key[1]} {key[2]} ps; overlaid DNA global ID {pose['proposal']}",
               f"32,769-pose sensitivity: {before} to {after}"]
        write_snapshot(out/filename,f,lines,dna @ r[i]+t[i],notes)
        exports.append(dict(file=filename,previous_status=before,dense_status=after,**pose))
    control_union=sum(r["receptor_accepted"] for r in native_rows)-3*native_rows[0]["identity_accepted"]
    control_best=min((r for r in native_rows if r["nearest_accepted_id"]!=""),
                     key=lambda r:r["nearest_accepted_rmsd_A"])
    united=dict(control_best,seed="union",proposals=32769,receptor_accepted=control_union,
                nearest_accepted_id=(control_best["seed"]-6280)*8192+control_best["nearest_accepted_id"])
    native_rows.append(united)
    for name,records in (("frames",rows),("seed-counts",per_seed),("summary",aggregate(rows)),
                         ("transitions",[dict(branch=k[0],subset=k[1],original_status=k[2],
                                             dense_status=k[3],frames=v) for k,v in transitions.items()]),
                         ("selected-poses",[v[1] for v in best.values()]),
                         ("representatives",exports),("native-control",native_rows)):
        if records:
            write_csv(out/(name+".csv"),records)
    dense_r=np.concatenate([transforms[i][0][0 if i==6280 else 1:] for i in SEEDS])
    dense_t=np.concatenate([transforms[i][1][0 if i==6280 else 1:] for i in SEEDS])
    assert len(dense_r)==len(dense_t)==32769 and len(rows)==200 and len(per_seed)==800
    assert len(np.unique(dense_t,axis=0))==32769
    np.savez_compressed(out/"transforms.npz",R=dense_r,t=dense_t)
    (out/"replay.json").write_text(json.dumps(dict(elapsed_s=time.perf_counter()-started,
                        seeds=list(SEEDS),frames=len(rows),poses=len(dense_r)),indent=2)+"\n")


if __name__=="__main__":
    directory=Path(__file__).resolve().parent
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs",type=Path,default=directory/"inputs")
    parser.add_argument("--primary",type=Path,default=directory)
    parser.add_argument("--out",type=Path,default=directory/"density")
    args=parser.parse_args()
    run(args.inputs,args.primary,args.out)
