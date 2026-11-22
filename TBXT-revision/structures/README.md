# Static structures and rigid-helix replay

The two B models show the same 64-heavy-atom probe in the same fragment-receptor frame tested against the same duplex. One passes whole-duplex clearance. The other places CypA across the duplex interior under the strict assembly cutoffs. The two C models illustrate the weaker primary criterion for the 65-heavy-atom C probe.

| Model | Seed/draw | Probe/protein cutoff (Å) | Closest CypA–DNA (Å) | Closest probe–DNA (Å) | Interior DNA atoms / nucleotides contacted |
|---|---|---|---:|---:|---:|
| B-cover | 62008/3069 | 2.6/3.0 | 0.705 | 7.168 | 8 / 2 |
| B-clear | 62006/1306 | 2.6/3.0 | 5.930 | 10.016 | 0 / 0 |
| C-cover | 62040/1261 | 2.2/2.6 | 0.480 | 10.274 | 6 / 2 |
| C-clear | 62003/2937 | 2.2/2.6 | 6.896 | 6.880 | 0 / 0 |

Each covering example is the sole qualifying interior hit in its indicated class after the helix test. We selected a clearing draw at the same cutoffs, preferring the covering example's seed if available, then the smallest seed and draw. Neither covering seed supplied a qualifying clear example. Each B/C pair therefore matches chemistry, templates, duplex and gates across two torsion draws. The selection gives existence examples, with no ranking by energy or structural likelihood.

## Inspect a model

Open the corresponding PDB to inspect all components together. Chain A is fragment-template TBXT (5QS1 for B, 5QSE with the Asp177→Gly back-mutation for C). Chain P is the 9GHY-derived CypA receptor. Chain L is the covalently intact probe. Chain H carries the 74 observed WT 8CDN heavy atoms at TBXT residues 212–221, aligned as an independent rigid body. Chains C and D are the virtual 8CDN DNA strands after registration. Hiding chains C/D reveals the modeled ternary assembly. Superimposed CypA–DNA atoms in the covering examples illustrate steric exclusion.

The companion SDF supplies bond order and stereochemistry at the exported chain-L coordinates. B-atom-map.csv and C-atom-map.csv relate PDB serial/atom names to the one-based SDF atom index, zero-based full and head indices, and original head/linker names. Each SDF was reopened and its coordinate-derived isomeric SMILES compared with the prescribed chemistry in inputs/B or inputs/C/provenance.json. PDB positions are rounded to 0.001 Å and SDF positions to 0.0001 Å. The calculations use the archived bank precision.

## Geometry and interpretation

Each fragment ends at residue 211. A least-squares fit of 152 common 8CDN chain-A Cα atoms at positions 45–208, excluding 110–121, transfers the helix into the fragment frame. The B/C core RMSDs are 1.010/1.144 Å. Chain H supplies an exclusion volume. The fragment retains its terminal OXT. Its C211 to transferred N212 distance is 1.567 Å for B and 2.306 Å for C. These files contain no repaired peptide junction or relaxed helix. The registered DNA is the separately selected closest admissible displaced duplex, draw 815 for B (1.498 Å from the aligned reference) and draw 623 for C (2.353 Å). Both retain 22/23 reference Cα–sugar contacts under the register generator's distance allowance.

All distances use heavy atoms. At probe cutoff a and protein cutoff b, an assembly requires probe–TBXT and probe–CypA minima ≥a, and CypA–TBXT ≥b. Adding the helix requires probe–H ≥a and CypA–H ≥b. Whole-duplex clearance additionally requires probe–DNA ≥a and CypA–DNA ≥b. Both DNA registers were prequalified against fragment TBXT at 2.6 Å. The added helix clears those registers by 3.118 Å (B) and 2.660 Å (C). An interior hit requires ≥3 distinct DNA atoms across ≥2 nucleotides within strictly 2.5 Å of CypA. The interior omits both terminal nucleotides on each strand, leaving 410 of the 486 DNA atoms. These two endpoints leave some assembly draws in neither category.

The B-cover closest full-duplex pair includes terminal DT10 on chain C. The separate interior-contact list identifies the atoms supporting the footprint. Within the assembly, B-cover passes the coarse gates with only 2.671 Å between TBXT Met159 CE and linker C059, and 2.685 Å between CypA Thr73 CG2 and probe N045. C-cover has only 2.255 Å between CypA Thr73 O and linker O065, and 2.798 Å between CypA Gly80 O and helix Ala216 N. These short distances limit the physical interpretation of the examples.

The torsion bank fixes each head in its cognate binary frame. In all four models the TBXT attachment atom (C4 for B or C6 for C) and CypA attachment atom (NBR) are within 0.0008 Å of their binary reference coordinates after the appropriate protein transform. The models retain all parent head–protein heavy-atom pairs initially within 4 Å (46 B-head pairs, 13 C-head pairs and 73 CypA-head pairs). These diagnostics check atom mapping and anchor placement by construction. They do not measure linker strain, binding affinity or occupancy. No energy minimization was performed for this comparison.

## Reproduce the tables and models

Run from this directory with Python, NumPy, SciPy and RDKit. The checked runtime used NumPy 2.4.6, SciPy 1.18.0 and RDKit 2026.09.2. To write to a fresh comparison directory, use

```console
python replay.py --out replay-output
```

The included inputs contain the donor structure, per-pocket receptor coordinates, probe and binary head chemistry, atom maps, the chosen full/interior duplex positions, and accepted.npz. Each NPZ stores every original primary accepted (seed, draw), its bank name, full ligand coordinates, row-vector CypA rotation R and translation t, and the archived gap/footprint row. Reconstruct moved CypA as source coordinates @ R + t. The script recomputes distances and interior counts from coordinates, checks finite archived gaps and the original footprint flag for every draw, and asserts all four count rows. The set has 110 B and 481 C draws.

The source hashes in inputs/source-sha256.json identify the original banks and upstream files. The distributed inputs support replay conditional on the original 591 primary fits. Reconstructing the full torsion prior and testing rejection of the remaining draws requires the complete banks or their generator and inputs.

## Tables

| File | Contents |
|---|---|
| counts.csv | Original and helix-conditioned ternary, whole-duplex and interior counts at both thresholds |
| draw-metrics.csv | Exact finite minima for all 591 draws. T, P, L, H and D denote fragment TBXT, CypA, ligand, helix and DNA. Distances are in Å |
| examples.csv | Four selected seed/draw identities with bank names, thresholds and unrounded metrics |
| closest-pairs.csv | Closest heavy-atom pair for each reported component pair, including separate full and interior CypA–DNA minima |
| interior-contacts.csv | Every interior DNA atom within 2.5 Å of CypA in each displayed model and the nearest CypA atom |
| head-contacts.csv | Per-model and per-head displacement from binary coordinates, closest contact and ≤4 Å pair retention |
| alignment.csv | Core fit, helix size, apparent junction distance and helix–duplex distance |
