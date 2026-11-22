"""Replay fixed-duplex endpoint sensitivity from the deposited static coordinate banks."""

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


HERE = Path(__file__).resolve().parent
WITNESSES = {("B", 62008, 3069), ("C", 62040, 1261)}
GATES = {"primary": (2.2, 2.6), "strict": (2.6, 3.0)}


def heavy_atoms(path):
    """Preserve PDB order, using the highest-occupancy alternate for each atom."""
    atoms = {}
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM  ") or line[76:78].strip() in {"H", "D"}:
            continue
        key = (line[21], int(line[22:26]), line[26], line[17:20].strip(), line[12:16].strip())
        occupation = float(line[54:60])
        if key not in atoms or occupation > atoms[key][0]:
            atoms[key] = occupation, [float(line[i:i + 8]) for i in (30, 38, 46)]
    return {k: np.array(v[1]) for k, v in atoms.items()}


def xyz(atoms):
    return np.array(list(atoms.values()))


def transferred_helix(donor, target):
    # The donor and target can differ at an amino acid. Match positions for these sites.
    def core(atoms):
        return {k[1]: v for k, v in atoms.items() if k[0] == "A" and k[2] == " "
                and k[4] == "CA" and 45 <= k[1] <= 208 and not 110 <= k[1] <= 121}
    dcore, tcore = core(donor), core(target)
    ids = sorted(dcore.keys() & tcore.keys())
    assert len(ids) == 152
    d, t = [np.array([table[i] for i in ids]) for table in (dcore, tcore)]
    u, _, vt = np.linalg.svd((d - d.mean(0)).T @ (t - t.mean(0)))
    rotation = u @ np.diag([1, 1, np.linalg.det(u @ vt)]) @ vt
    helix = np.array([v for k, v in donor.items() if k[0] == "A" and 212 <= k[1] <= 221])
    assert len(helix) == 74
    return (helix - d.mean(0)) @ rotation + t.mean(0)


def atom_label(key, chain=None):
    return f"{chain or key[0]}/{key[3]}{key[1]}/{key[4]}"


def register(dna):
    keys = list(dna)
    pairs = np.array([k[1] + 2 if k[0] == "C" else 17 - k[1] for k in keys])
    assert sorted(set(pairs)) == list(range(1, 13)) and len(keys) == 486
    rows = []
    for bp in range(1, 13):
        row = {"bp": bp}
        anchors = []
        for strand in "CD":
            members = [k for i, k in enumerate(keys) if pairs[i] == bp and k[0] == strand]
            [(number, res)] = sorted({(k[1], k[3]) for k in members})
            atom = "N1" if res in {"DG", "DA"} else "N3"
            [anchor] = [k for k in members if k[4] == atom]
            row.update({f"{strand}_residue": f"{res}{number}", f"{strand}_anchor": atom})
            anchors.append(anchor)
        assert (anchors[0][3][1] + anchors[1][3][1]) in {"AT", "TA", "GC", "CG"}
        row["anchor_distance_A"] = round(np.linalg.norm(dna[anchors[0]] - dna[anchors[1]]), 6)
        rows.append(row)
    return pairs, rows


def footprint(mask, hits, keys, pairs):
    ids = np.flatnonzero(mask & hits)
    return len(ids), len({keys[i][:2] for i in ids}), len(set(pairs[ids]))


def describe_witness(key, dna_keys, p_keys, pairs, nearest, distances, exported_path):
    pocket, seed, draw = key
    identity = dict(pocket=pocket, seed=seed, draw=draw)
    atom_rows, nt_rows, min_rows = [], [], []
    exported = heavy_atoms(exported_path)
    exported_cypa = np.array([v for k, v in exported.items() if k[0] == "P"])
    exported_dna = np.array([exported[k] for k in dna_keys])
    pdb_distances = cKDTree(exported_cypa).query(exported_dna)[0]
    # Each endpoint of a distance is rounded to 0.001 A independently on three axes.
    np.testing.assert_allclose(distances, pdb_distances, atol=np.sqrt(3) * 0.001, rtol=0)
    np.testing.assert_array_equal(distances < 2.5, pdb_distances < 2.5)
    for bp in range(1, 13):
        for strand in "CD":
            group = [i for i, k in enumerate(dna_keys) if k[0] == strand and pairs[i] == bp]
            i_min = min(group, key=distances.__getitem__)
            residue = dna_keys[i_min]
            common = dict(**identity, bp=bp, strand=strand, residue_number=residue[1],
                          residue_name=residue[3])
            hit_indices = [i for i in group if distances[i] < 2.5]
            nt_rows.append(dict(**common, hit_atoms=len(hit_indices),
                                min_A=round(distances[i_min], 6)))
            for i in hit_indices:
                atom_rows.append(dict(**common, dna_atom=dna_keys[i][4],
                                      cypa_atom=atom_label(p_keys[nearest[i]], "P"),
                                      distance_A=round(distances[i], 6)))
    for name, lo, hi in [("whole", 1, 12), ("inner_10", 2, 11), ("central_8", 3, 10)]:
        mask = (pairs >= lo) & (pairs <= hi)
        ids = np.flatnonzero(mask)
        i = ids[distances[mask].argmin()]
        stats = footprint(mask, distances < 2.5, dna_keys, pairs)
        min_rows.append(dict(**identity, mask=name, hit_atoms=stats[0], hit_nt=stats[1],
                             hit_bp=stats[2], closest_bp=pairs[i], dna_atom=atom_label(dna_keys[i]),
                             cypa_atom=atom_label(p_keys[nearest[i]], "P"), min_A=round(distances[i], 6),
                             exported_PDB_min_A=round(pdb_distances[mask].min(), 6)))
    return atom_rows, nt_rows, min_rows


def write_table(out, name, rows):
    with (out / name).open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def replay(structures, out):
    inputs = structures / "inputs"
    donor = heavy_atoms(inputs / "8CDN.pdb")
    dna = {k: v for k, v in donor.items() if k[0] in {"C", "D"} and k[3] in {"DA", "DC", "DG", "DT"}}
    dna_keys = list(dna)
    pairs, register_rows = register(dna)
    inner = (pairs >= 2) & (pairs <= 11)
    central = (pairs >= 3) & (pairs <= 10)
    assert int(inner.sum()) == 410
    with (structures / "draw-metrics.csv").open() as f:
        reference = {(r["case"], int(r["seed"]), int(r["draw"])): r for r in csv.DictReader(f)}
    draws, summary, atom_rows, nt_rows, min_rows = [], [], [], [], []
    seen = set()
    for pocket in "BC":
        folder = inputs / pocket
        tbxt = heavy_atoms(folder / "tbxt.pdb")
        cypa = heavy_atoms(folder / "cypa.pdb")
        cypa_xyz = xyz(cypa)
        p_keys = list(cypa)
        target = cKDTree(xyz(tbxt))
        helix = cKDTree(transferred_helix(donor, tbxt))
        with np.load(folder / "duplex.npz") as archive:
            duplex = archive["full"]
            np.testing.assert_allclose(duplex[inner], archive["inner"])
        d_tree = cKDTree(duplex)
        assert min(target.query(duplex)[0].min(), helix.query(duplex)[0].min()) >= 2.6
        with np.load(folder / "accepted.npz") as source:
            bank = dict(source)
        for index, (seed, draw) in enumerate(bank["keys"]):
            key = pocket, int(seed), int(draw)
            assert key not in seen
            seen.add(key)
            lig = bank["xyz"][index]
            moved = cypa_xyz @ bank["R"][index] + bank["t"][index]
            p_tree = cKDTree(moved)
            distances, nearest = p_tree.query(duplex)
            hits = distances < 2.5
            inner_stats, central_stats = [footprint(mask, hits, dna_keys, pairs) for mask in (inner, central)]
            gaps = dict(T_L=target.query(lig)[0].min(), P_L=p_tree.query(lig)[0].min(),
                        T_P=target.query(moved)[0].min(), H_L=helix.query(lig)[0].min(),
                        H_P=helix.query(moved)[0].min(), L_D=d_tree.query(lig)[0].min(),
                        P_D=distances.min())
            original = reference[key]
            np.testing.assert_allclose(list(gaps.values()), [float(original[k]) for k in gaps], atol=1e-5, rtol=0)
            assert inner_stats[:2] == (int(original["interior_atoms"]), int(original["interior_nt"]))
            published = inner_stats[0] >= 3 and inner_stats[1] >= 2
            assert published == bool(bank["source_scores"][index, 5]) == (original["footprint"] == "True")
            row = dict(pocket=pocket, seed=key[1], draw=key[2], bank=str(bank["bank"][index]),
                       inner_atoms=inner_stats[0], inner_nt=inner_stats[1], inner_bp=inner_stats[2],
                       central_atoms=central_stats[0], central_nt=central_stats[1], central_bp=central_stats[2],
                       published_inner_10=int(published),
                       two_bp_inner_10=int(inner_stats[0] >= 3 and inner_stats[2] >= 2),
                       two_nt_central_8=int(central_stats[0] >= 3 and central_stats[1] >= 2))
            for gates, (a, b) in GATES.items():
                admitted = min(gaps["T_L"], gaps["P_L"], gaps["H_L"]) >= a and min(gaps["T_P"], gaps["H_P"]) >= b
                row[f"{gates}_admitted"] = int(admitted)
                row[f"{gates}_whole_clear"] = int(admitted and gaps["L_D"] >= a and gaps["P_D"] >= b)
            draws.append(row)
            if key in WITNESSES:
                atoms, nts, minima = describe_witness(key, dna_keys, p_keys, pairs, nearest, distances,
                                                      structures / f"{pocket}-cover.pdb")
                atom_rows.extend(atoms)
                nt_rows.extend(nts)
                min_rows.extend(minima)
        for gates in GATES:
            selected = [r for r in draws if r["pocket"] == pocket and r[f"{gates}_admitted"]]
            summary.append(dict(pocket=pocket, gates=gates, admitted=len(selected),
                                whole_clear=sum(r[f"{gates}_whole_clear"] for r in selected),
                                **{k: sum(r[k] for r in selected) for k in
                                   ("published_inner_10", "two_bp_inner_10", "two_nt_central_8")}))
    assert seen == reference.keys() and len(seen) == 591 and len(min_rows) == 6
    with (structures / "counts.csv").open() as f:
        for row, old in zip(summary, csv.DictReader(f), strict=True):
            assert (row["pocket"], *GATES[row["gates"]]) == (old["case"], float(old["ligand_cutoff_A"]), float(old["protein_cutoff_A"]))
            assert [row[k] for k in ("admitted", "whole_clear", "published_inner_10")] == [
                int(old[k]) for k in ("helix_fit", "helix_whole_clear", "helix_footprint")]
    out.mkdir(parents=True, exist_ok=True)
    for name, rows in [("register.csv", register_rows), ("endpoint-sensitivity.csv", summary),
                       ("draw-endpoints.csv", draws), ("witness-atoms.csv", atom_rows),
                       ("witness-nucleotides.csv", nt_rows), ("witness-minima.csv", min_rows)]:
        write_table(out, name, rows)
    print(f"Verified {len(seen)} deposited draws from two pockets and four gate rows.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structures", type=Path, default=HERE.parent / "structures")
    parser.add_argument("--out", type=Path, default=HERE)
    options = parser.parse_args()
    replay(options.structures, options.out)
