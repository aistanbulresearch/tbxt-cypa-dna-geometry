"""Read mapped solute coordinates and reproduce the two-step 8CDN registration."""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


def records(path):
    entries = {}
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or line[76:78].strip() in ("H", "D"):
            continue
        key = line[21], int(line[22:26]), line[12:16].strip()
        if key not in entries or float(line[54:60]) > float(entries[key][54:60]):
            entries[key] = line
    return entries


def coords(lines):
    return np.array([[float(line[i:i+8]) for i in (30, 38, 46)] for line in lines])


def fit(mobile, target):
    u, _, v = np.linalg.svd((mobile-mobile.mean(0)).T @ (target-target.mean(0)))
    r = u @ np.diag([1, 1, np.sign(np.linalg.det(u @ v))]) @ v
    return r, target.mean(0)-mobile.mean(0) @ r


def in_core(k):
    return k[0] == "A" and k[2] == "CA" and 45 <= k[1] <= 208 and not 110 <= k[1] <= 121


def native(root):
    ref = records(root / "work/8CDN.pdb")
    lines = [line for key, line in ref.items() if key[0] in ("C", "D")
             and line[17:20].strip() in ("DA", "DC", "DG", "DT")]
    ids = [(line[21], int(line[22:26]), line[12:16].strip()) for line in lines]
    residues = {c: sorted({n for chain, n, _ in ids if chain == c}) for c in ("C", "D")}
    sugar = [i for i, (c, n, k) in enumerate(ids) if k == "C1'" and n in residues[c][1:-1]]
    assert len(lines) == 486 and len(sugar) == 20
    return ref, lines, coords(lines), np.array(sugar)


def protein_lines(mapping, source):
    """Include atoms added during system construction, including terminal OXT."""
    next_serial = max(int(line[6:11]) for line in source.values()) + 1
    out = []
    for c in ("A", "P"):
        for _, n, name in mapping[c]:
            line = source.get((c, n, name))
            if line is None:
                line = source[c, n, "CA"]
                line = (line[:6] + f"{next_serial:5d} {name:>4s}" + line[16:76]
                        + f"{name[0]:>2s}  ")
                next_serial += 1
            out.append(line)
    probe = {(c, name): line for (c, _, name), line in source.items() if c in ("X","Y","Z")}
    out += [probe[c, name] for _, c, name in mapping["probe"]]
    return out


def frames(root, branch, ref):
    alt = branch == "ref-B-on"
    start = branch.removeprefix("ref-")
    base = root / ("work/b-linkers/cases/ref" if alt else "work/monomer")
    src = base / "inputs" / f"{start}.pdb"
    mapping = json.loads((base / "systems" / start / "map.json").read_text())
    assert hashlib.sha256(src.read_bytes()).hexdigest() == mapping["input_sha256"]
    source = records(src)
    pdb_lines = protein_lines(mapping, source)
    conect = [s for s in src.read_text().splitlines() if s.startswith("CONECT")]
    heads = [r for r in csv.DictReader((root/"outputs/B-occupancy/crystal-heads/frames.csv").open())
             if r["phase"] == ("new-ref" if alt else "new") and r["start"] == start]
    heads = {(int(r["rep"]), float(r["time_ps"])): tuple(float(r[k+"_crystal_head_A"])
             for k in ("TBXT","CypA")) for r in heads}
    core = [(c,n,k) for c in ("A",) for _,n,k in mapping[c] if in_core((c,n,k)) and (c,n,k) in source]
    source_core = coords([source[k] for k in core])
    shared = [k for k in core if k in ref]
    crystal_core = coords([ref[k] for k in shared])
    rc, tc = fit(crystal_core, coords([source[k] for k in shared]))
    for rep in range(1, 3 if alt else 4):
        filename = "traj.npz" if alt else "traj-10ps.npz"
        path = root/"work/b-occupancy/runs"/branch/f"rep-{rep}"/filename
        with np.load(path) as traj:
            ix = {int(k): j for j, k in enumerate(traj["ids"])}
            at = {(c,n,k): ix[i] for c in ("A","P") for i,n,k in mapping[c]}
            a,p,l = [[ix[i] for i, _, _ in mapping[c]] for c in ("A","P","probe")]
            xo,yo = [[ix[i] for i,c,_ in mapping["probe"] if c == tag] for tag in ("X","Y")]
            pocket = [at["A",n,"CA"] for n in (119,120,159)]
            pc = [v for (c,n,k), v in at.items() if c=="P" and 2<=n<=163 and k=="CA"]
            ci = [at[k] for k in core]
            shared_in_A = [[(n,k) for _,n,k in mapping["A"]].index((n,k)) for _,n,k in shared]
            assert np.all(traj["stage"] == "unrestrained")
            for raw, box, t in zip(traj["xyz"], traj["boxes"], traj["meta"][:,0], strict=True):
                if t == 0:
                    continue
                x = raw.astype(float)
                inv = np.linalg.inv(box)
                x[l] -= np.rint((x[xo].mean(0)-x[pocket].mean(0)) @ inv) @ box
                x[p] -= np.rint((x[pc].mean(0)-x[yo].mean(0)) @ inv) @ box
                r, shift = fit(x[ci], source_core)
                x = ((x @ r + shift) - tc) @ rc.T
                yield dict(branch=branch, replica=rep, time_ps=float(t),
                           A=x[a], P=x[p], L=x[l], core=x[a][shared_in_A],
                           ref_core=crystal_core, core_res=[k[1] for k in shared],
                           heads=heads[rep,t], lines=pdb_lines, conect=conect)


def write_snapshot(path, frame, dna_lines, dna, notes):
    full = np.vstack((frame["A"], frame["P"], frame["L"], dna))
    lines = list(frame["lines"])
    serial = max(int(s[6:11]) for s in lines) + 1
    lines += [s[:6] + f"{serial+j:5d}" + s[11:] for j,s in enumerate(dna_lines)]
    out = ["REMARK 950 " + s for s in notes]
    for s, xyz in zip(lines, full, strict=True):
        if out[-1].startswith(("ATOM  ", "HETATM")) and out[-1][21] != s[21]:
            out.append("TER")
        out.append(s[:30] + "".join(f"{v:8.3f}" for v in xyz) + s[54:])
    path.write_text("\n".join(out + ["TER"] + frame["conect"] + ["END",""]), encoding="ascii")
