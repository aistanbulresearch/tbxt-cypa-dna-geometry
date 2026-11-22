"""Rebuild the rigid-helix sensitivity counts and export matched illustrative models."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from rdkit import Chem
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist

HERE = Path(__file__).resolve().parent


def read_atoms(path, record="ATOM  "):
    selected = {}
    for line in path.read_text().splitlines():
        if not line.startswith(record) or line[76:78].strip() in ("H", "D"):
            continue
        key = (line[21], line[22:27], line[17:20], line[12:16])
        if key not in selected or float(line[54:60]) > float(selected[key][54:60]):
            selected[key] = line
    return list(selected.values())


def coords(lines):
    return np.array([[float(line[i:i+8]) for i in (30, 38, 46)] for line in lines])


def alignment(donor, fragment):
    def core(lines):
        return {int(l[22:26]): l for l in lines if l[12:16].strip() == "CA"
                and l[21] == "A" and not l[26].strip() and 45 <= int(l[22:26]) <= 208
                and not 110 <= int(l[22:26]) <= 121}
    s, t = core(donor), core(fragment)
    common = sorted(s.keys() & t.keys())
    assert len(common) == 152
    x, y = [coords([d[i] for i in common]) for d in (s, t)]
    u, _, v = np.linalg.svd((x - x.mean(0)).T @ (y - y.mean(0)))
    r = u @ np.diag([1, 1, np.linalg.det(u @ v)]) @ v
    shift = y.mean(0) - x.mean(0) @ r
    return r, shift, float(np.sqrt(np.sum((x @ r + shift - y)**2, axis=1).mean()))


def write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)


def name(line, chain=None):
    return f"{chain or line[21]}/{line[17:20].strip()}{line[22:26].strip()}/{line[12:16].strip()}"


def move_atom(line, xyz, chain, serial):
    return (line[:6] + f"{serial:5d}" + line[11:16] + " " + line[17:21] + chain
            + line[22:30] + "".join(f"{x:8.3f}" for x in xyz) + line[54:])


def geometry(case, inp):
    donor = read_atoms(inp / "8CDN.pdb")
    folder = inp / case
    frag, cyp = [read_atoms(folder / f) for f in ("tbxt.pdb", "cypa.pdb")]
    lig = read_atoms(folder / "probe.pdb", "HETATM")
    dna = [l for l in donor if l[21] in "CD" and l[17:20].strip() in ("DA", "DC", "DG", "DT")]
    hel = [l for l in donor if l[21] == "A" and 212 <= int(l[22:26]) <= 221]
    numbers = {ch: sorted({int(l[22:26]) for l in dna if l[21] == ch}) for ch in "CD"}
    inner = np.array([int(l[22:26]) in numbers[l[21]][1:-1] for l in dna])
    assert (len(dna), sum(inner), len(hel)) == (486, 410, 74)
    rot, trans, rms = alignment(donor, frag)
    with np.load(folder / "duplex.npz") as reg:
        d = reg["full"]
        assert np.allclose(d[inner], reg["inner"])
    return dict(lines={"A": frag, "H": hel, "P": cyp, "L": lig, "DNA": dna},
                xyz={"A": coords(frag), "H": coords(hel) @ rot + trans, "P": coords(cyp), "DNA": d},
                inner=inner, rms=rms)


def measure(bank, geom):
    x = geom["xyz"]
    dt = cKDTree(x["DNA"])
    tt, ht = [cKDTree(x[k]) for k in ("A", "H")]
    nts = np.array([l[21:27].strip() for l in geom["lines"]["DNA"]])[geom["inner"]]
    rows = []
    for i, (seed, draw) in enumerate(bank["keys"]):
        l = bank["xyz"][i]
        p = x["P"] @ bank["R"][i] + bank["t"][i]
        pt = cKDTree(p)
        pd = pt.query(x["DNA"])[0]
        hits = pd[geom["inner"]] < 2.5
        n_atoms, n_nt = int(hits.sum()), len(set(nts[hits]))
        g = [tt.query(l)[0].min(), pt.query(l)[0].min(), tt.query(p)[0].min(),
             dt.query(l)[0].min(), pd.min()]
        stored = bank["source_scores"][i]
        finite = np.isfinite(stored[:5])
        assert np.allclose(np.asarray(g)[finite], stored[:5][finite], atol=.0001)
        assert np.all(np.asarray(g)[~finite] >= 3.0099)
        hit = n_atoms >= 3 and n_nt >= 2
        assert hit == bool(stored[5])
        rows.append(dict(seed=int(seed), draw=int(draw), bank=bank["bank"][i],
                         T_L=g[0], P_L=g[1], T_P=g[2], L_D=g[3], P_D=g[4],
                         H_L=ht.query(l)[0].min(), H_P=ht.query(p)[0].min(),
                         interior_atoms=n_atoms, interior_nt=n_nt, footprint=hit))
    return rows


def classify(row, a, b):
    fit = min(row["T_L"], row["P_L"]) >= a and row["T_P"] >= b
    allowed = fit and row["H_L"] >= a and row["H_P"] >= b
    clear = row["L_D"] >= a and row["P_D"] >= b
    return [fit, fit and clear, fit and row["footprint"], allowed, allowed and clear,
            allowed and row["footprint"]]


def export(case, label, i, geom, bank, folder, out):
    model = case + "-" + label
    seed, draw = map(int, bank["keys"][i])
    lines, x = geom["lines"], geom["xyz"].copy()
    x["L"] = bank["xyz"][i].astype(float)
    x["P"] = x["P"] @ bank["R"][i] + bank["t"][i]
    serial = 1
    contents = [f"REMARK 900 {model} SEED {seed} DRAW {draw}. HYPOTHETICAL TERNARY MODEL.",
                "REMARK 900 H IS A RIGID HELIX OVERLAY. C AND D ARE VIRTUAL DNA.",
                "REMARK 900 LIGAND CHEMISTRY AND STEREOCHEMISTRY ARE IN THE MATCHING SDF."]
    mol = Chem.MolFromMolFile(str(folder / "probe.sdf"))
    atom_map = list(csv.DictReader((folder / "head_atoms.csv").open()))
    assert len(atom_map) == len(lines["L"]) == mol.GetNumAtoms()
    for chain in ("A", "H", "P", "L", "DNA"):
        if chain == "L":
            first = serial
        for j, (line, point) in enumerate(zip(lines[chain], x[chain], strict=True)):
            contents.append(move_atom(line, point, line[21] if chain == "DNA" else chain, serial))
            if chain == "L":
                m = atom_map[j]
                assert int(m["full_index"]) == j and line[76:78].strip() == m["element"]
                m.update(pdb_atom=line[12:16].strip(), pdb_serial=serial, sdf_atom=j+1,
                         role={"0": "linker", "1": "TBXT head", "2": "CypA head"}[m["head"]])
                mol.GetConformer().SetAtomPosition(j, point.tolist())
            serial += 1
            if chain == "DNA" and j+1 < len(lines[chain]) and line[21] != lines[chain][j+1][21]:
                contents.append("TER")
        contents.append("TER")
    contents.extend(f"CONECT{first+b.GetBeginAtomIdx():5d}{first+b.GetEndAtomIdx():5d}" for b in mol.GetBonds())
    (out / f"{model}.pdb").write_text("\n".join(contents + ["END", ""]), encoding="ascii")
    mol.SetProp("_Name", f"{model} seed {seed} draw {draw}")
    with Chem.SDWriter(str(out / f"{model}.sdf")) as writer:
        writer.write(mol)
    reread = Chem.MolFromMolFile(str(out / f"{model}.sdf"))
    Chem.RemoveStereochemistry(reread)
    Chem.AssignStereochemistryFrom3D(reread)
    prescribed = json.loads((folder / "provenance.json").read_text())["input_chemistry"]
    assert Chem.MolToSmiles(reread) == Chem.MolToSmiles(Chem.MolFromSmiles(prescribed))
    write_csv(out / f"{case}-atom-map.csv", atom_map)
    pairs = []
    for c1, c2 in (("A","P"),("A","L"),("P","L"),("H","L"),("H","P"),
                   ("H","DNA"),("A","DNA"),("P","DNA"),("L","DNA")):
        ds, idx = cKDTree(x[c1]).query(x[c2])
        j = int(ds.argmin())
        pairs.append(dict(model=model, component_1=c1, atom_1=name(lines[c1][idx[j]], c1),
                          component_2=c2, atom_2=name(lines[c2][j], c2 if c2 != "DNA" else None),
                          distance_A=round(ds[j], 6)))
    distances, indices = cKDTree(x["P"]).query(x["DNA"])
    interior = np.flatnonzero(geom["inner"])
    j = interior[np.argmin(distances[interior])]
    pairs.append(dict(model=model, component_1="P", atom_1=name(lines["P"][indices[j]],"P"),
                      component_2="interior DNA", atom_2=name(lines["DNA"][j]), distance_A=round(distances[j],6)))
    hits = [dict(model=model, cypa_atom=name(lines["P"][indices[j]],"P"), dna_atom=name(lines["DNA"][j]),
                 distance_A=round(distances[j],6)) for j in interior if distances[j] < 2.5]
    head_rows = []
    for h, target, attachment in [(1,"A",{"B":"C4","C":"C6"}[case]),(2,"P","NBR")]:
        amap = sorted((r for r in atom_map if int(r["head"]) == h), key=lambda r:int(r["head_index"]))
        ids = [int(r["full_index"]) for r in amap]
        ref = Chem.MolFromMolFile(str(folder / f"head{h}.sdf")).GetConformer().GetPositions()
        if h == 2:
            ref = ref @ bank["R"][i] + bank["t"][i]
        local = x["L"][ids]
        dev = np.linalg.norm(local - ref, axis=1)
        j = next(k for k,r in enumerate(amap) if r["atom_name"] == attachment)
        dm, dm_ref = cdist(local,x[target]), cdist(ref,x[target])
        atom, target_atom = np.unravel_index(dm.argmin(), dm.shape)
        native_pairs = (dm_ref <= 4)
        head_rows.append(dict(model=model, head=target, attachment=attachment, attachment_shift_A=dev[j],
            head_rmsd_A=np.sqrt((dev**2).mean()), max_shift_A=dev.max(),
            closest_head_atom=amap[atom]["atom_name"], closest_pocket_atom=name(lines[target][target_atom],target),
            closest_A=dm.min(), native_pairs_4A=native_pairs.sum(), retained_pairs_4A=(native_pairs&(dm<=4)).sum()))
    return pairs, head_rows, hits


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--inputs", type=Path, default=HERE/"inputs")
    p.add_argument("--out", type=Path, default=HERE)
    opt = p.parse_args()
    opt.out.mkdir(parents=True, exist_ok=True)
    all_rows, counts, examples, close, heads, fits, hits = [], [], [], [], [], [], []
    expected = {"B": [[110,78,15,110,78,15],[22,15,1,22,15,1]],
                "C": [[481,213,26,121,96,1],[213,99,5,56,47,0]]}
    for case in "BC":
        folder = opt.inputs / case
        geom = geometry(case, opt.inputs)
        with np.load(folder / "accepted.npz") as source:
            bank = dict(source)
        rows = measure(bank, geom)
        all_rows.extend(dict(case=case, **r) for r in rows)
        x = geom["xyz"]
        helix_dna = cKDTree(x["H"]).query(x["DNA"])[0].min()
        junction = {}
        for chain, num, atom in [("A",211,"C"),("H",212,"N")]:
            j = next(k for k,line in enumerate(geom["lines"][chain])
                     if int(line[22:26]) == num and line[12:16].strip() == atom)
            junction[chain] = x[chain][j]
        fits.append(dict(case=case, aligned_CA=152, core_RMSD_A=geom["rms"], helix_atoms=len(x["H"]),
                         junction_C211_N212_A=np.linalg.norm(junction["A"]-junction["H"]),
                         helix_DNA_A=helix_dna))
        for k,(a,b) in enumerate([(2.2,2.6),(2.6,3.0)]):
            masks = np.array([classify(r,a,b) for r in rows])
            nums = masks.sum(0).tolist()
            assert nums == expected[case][k], (case, a, b, nums)
            counts.append(dict(case=case, ligand_cutoff_A=a, protein_cutoff_A=b,
                **dict(zip(("fit","whole_clear","footprint","helix_fit","helix_whole_clear","helix_footprint"),nums)),
                helix_DNA_A=helix_dna))
            print(case, a, b, nums, flush=True)
            if (case == "B") != (k == 1):
                continue
            cover, = np.flatnonzero(masks[:,5])
            clear = min(np.flatnonzero(masks[:,4]), key=lambda j:
                        (rows[j]["seed"] != rows[cover]["seed"], rows[j]["seed"], rows[j]["draw"]))
            for label,i in [("cover",cover),("clear",clear)]:
                model = case + "-" + label
                pairs, metrics, interior_hits = export(case, label, i, geom, bank, folder, opt.out)
                close.extend(pairs)
                heads.extend(metrics)
                hits.extend(interior_hits)
                examples.append(dict(model=model, ligand_cutoff_A=a, protein_cutoff_A=b, **rows[i]))
                print(model, rows[i], flush=True)
    for name, rows in [("draw-metrics",all_rows),("counts",counts),("examples",examples),
                       ("closest-pairs",close),("head-contacts",heads),("alignment",fits),
                       ("interior-contacts",hits)]:
        write_csv(opt.out / f"{name}.csv", rows)


if __name__ == "__main__":
    main()
