# Vina scoring of the isolated binary pockets

This pilot tests whether two intact TBXT–CypA probes retain a known head pose during local optimization in one isolated protein at a time. It uses AutoDock Vina 1.2.7, its standard score and a rigid receptor. The Vina interface reports empirical scores in kcal/mol. The scores describe the sampled poses and provide no measured affinity, thermodynamic cycle closure, DNA displacement estimate, or equilibrium population.

## Design and run order

1. **Initial test.** One independently embedded full probe from each reference B/C design was rigidly fitted by its mapped head to each cognate parent structure. This fit can leave substantial head RMSD or linker and distal-head clashes. Those initial inputs and their optimized outputs are retained in `vina-results.tsv`.
2. **Prefiltered test.** In the base B and C conformer banks from the binary insertion study, for each probe–receptor pair we took the first conformer clearing the primary 2.6 Å receptor–probe minimum heavy-atom cutoff in each of the first five qualifying seed streams ordered numerically. All 20 selected rows came from streams 62001–62005. No empirical score entered inclusion. Both banks supply internally self-compatible conformations. We used their direct TBXT frame and inverse stored CypA transforms and checked that isomeric SMILES reassigned from each 3D pose matched the input. `vina-bank-results.tsv` preserves the exact stream, draw, row and distance cutoff evidence along with each raw result.
3. **Receptor audit.** Meeko's original template preparation omitted incompletely modeled protein residues. The receptor audit measured the distances to the known head, the complete input and optimized bank probes, and the independently embedded full probe. 5QS1 Arg43 stays more than 26 Å from all five selected full B probes before and after optimization. The 5QSE-derived WT input lost no residue. On CypA, missing Lys82 is 6.493 Å from the crystallographic head and about 3.34–3.78 Å from the nearest locally optimized full probe. This can affect scores. The original CypA outputs are labeled by this preparation history.
4. **CypA completion sensitivity.** PDBFixer added 13 heavy atoms while preserving every heavy atom coordinate from the CypA input. Meeko retained all 167 input residues. The same bank rows and comparator types were used to compare the completed receptor with the initial preparation. The modeled side chains provide one structural completion.

## Full-probe results

| Probe | Rigid receptor | Retained head | Local score range | Head RMSD (Å) |
| --- | --- | ---: | ---: | ---: |
| B | TBXT 5QS1 | 5/5 | −4.427 to −3.047 | 0.528–0.622 |
| C | TBXT 5QSE-derived WT | 5/5 | −4.281 to −3.141 | 0.683–1.054 |
| B | Completed CypA 9GHY | 5/5 | −7.578 to −6.874 | 0.345–0.510 |
| C | Completed CypA 9GHY | 5/5 | −7.075 to −6.042 | <0.001–0.413 |

These short optimizations retained the prepositioned head in every selected conformer. The initially embedded B probe, by contrast, lost its head pose on both proteins after local optimization. The placement sensitivity rules out using that single input to dismiss B. The selected conformers establish local pose compatibility under this rigid empirical score. The experiment does not estimate how often such poses form from solution or how the B/C probes rank by affinity.

The independently embedded C probe scored +183.493 against the original CypA preparation and its local-only run exited without a pose. This failed attempt remains in the results table and contributes no RMSD, score, or retention count after optimization.

The B TBXT receptor and head are from [5QS1/NZ1](https://www.rcsb.org/structure/5QS1). The C head is [5QSE/O0V](https://www.rcsb.org/structure/5QSE). Its TBXT receptor is a model derived from 5QSE by changing Asp177 to Gly, so O0V pose recovery is a cross-variant check. The CypA and macrocycle use [9GHY/A1ILH](https://www.rcsb.org/structure/9GHY). The crystal growth solution was reported at pH 8.8. The macrocycle parent 19c was prepared as both neutral and terminal primary ammonium; this is a protonation sensitivity check. For 19f, we added the experimentally specified N-acetyl cap to the 19c graph and minimized the cap and hydrogens with RDKit MMFF while fixing the original 19c heavy atoms. A 19f pose is therefore model-derived, measured only at the level of experimental affinity. Table 1 reports CypA FP Ki of 5.2 ± 1.5 μM for 19c and 0.43 ± 0.02 μM for 19f. [Dreizler et al., Table 1 and Fig. 2](https://doi.org/10.1021/acsmedchemlett.4c00427).

## Scoring and pose checks

Meeko adds hydrogens and AutoDock atom types to the isolated, water-free receptors and explicit-hydrogen ligands. Macrocycle internal geometry is fixed. All scores use a box centered on the reference head. The score-only and local-only boxes are at least 22 Å on each axis and enclose the starting molecule with margin. All local placements are free to move during optimization. The retained-head criterion is symmetry-aware heavy-atom RMSD of at most 2 Å relative to the deposited head coordinates, measured directly in the receptor frame. The setup excludes the second protein, DNA, head restraints, and explicit solvent. The Vina empirical terms also omit a physical electrostatic charging free energy.

Parent site-directed docking uses a 22 Å cube, seeds 104729, 130363 and 15485863, a single CPU thread, up to nine poses and a six-score-unit reporting range. The TBXT/NZ1 and TBXT/O0V searches use exhaustiveness 24 with the default evaluation budget. The neutral 19c, protonated 19c and modeled 19f searches on both CypA receptor preparations each use exhaustiveness 4 and `--max_evals 750` per Vina invocation. The higher-effort CypA sessions at exhaustiveness 24 had no completed ranked poses before cancellation. From the first saved local result to the decision to cancel, approximately 10 minutes elapsed for the initial receptor and 3 minutes for the completed receptor. These elapsed intervals cover the broader CypA stages with overlapping jobs and are only wall-time estimates. The three higher-effort NZ1 jobs finished and remain included. An unsuccessful CypA redock is inconclusive at the capped search budget.

The docking search is blinded to the experimental orientation given the pocket and the supplied starting macrocycle internal geometry. For 19f, RMSD is measured over the common 19c heavy atoms. Independent full-probe docking outside the prespecified placements was out of scope for this local test.

All three rank-1 TBXT parent poses were within 2 Å of their crystal ligand. NZ1 on 5QS1 gave 0.740–0.787 Å. O0V on the G177 WT model derived from 5QSE gave 1.298–1.353 Å against the crystallographic D177 complex.

| Ligand on completed CypA | Fixed-geometry score | Local score | Capped rank-1 median score | Rank-1 head RMSD ≤2 Å |
| --- | ---: | ---: | ---: | ---: |
| 19c, neutral | −8.396 | −9.000 | −5.667 | 0/3 |
| 19c, protonated | −8.396 | −9.000 | −6.327 | 1/3 |
| Modeled 19f | −7.649 | −9.569 | −6.019 | 1/3 |

The fixed-geometry scores favor 19c while the FP measurements favor 19f by about 12.1-fold. Local optimization reverses that score ordering, and the median rank-1 scores change the ordering with 19c protonation. These controls support treating Vina as a local compatibility test here. The 19f redocking count compares its common core with the 19c crystal structure and gives no direct measure of 19f pose accuracy.

## Software

The calculations used AutoDock Vina 1.2.7, Python 3.12, NumPy 2.5.5, RDKit 2026.09.2, OpenMM 8.6.1, PDBFixer 1.12, Meeko 0.8.0 and Gemmi 0.7.5. The Windows Vina executable had SHA-256 `e0c4b2715e0c1a74f6e92d0f3be0328ac97542eafbc111e6b1efad897a73cce5`. [Vina release](https://github.com/forlilab/AutoDock-Vina/releases/tag/v1.2.7) · [Meeko documentation](https://meeko.readthedocs.io/en/develop/).
