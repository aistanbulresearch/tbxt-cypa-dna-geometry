# Checking both ends of the B probe against their binary crystal poses

For each saved frame, the TBXT-binding end is compared with NZ1 from 5QS1 and the CypA-binding end with the 19c ligand from 9GHY. These are two experimentally observed isolated head poses. Each protein has its own coordinate fit. The calculation uses the same binary reference structures used for the docking benchmark.

For TBXT, superpose corresponding Cα atoms in residues 45–208, omitting the modeled 110–121 loop. For CypA, use residues 2–163. Move each ligand head to the nearest periodic image of its own pocket. Calculate its heavy-atom RMSD in the protein frame, allowing graph symmetries that preserve the linker-attachment atom and chirality. The protein fit is the only spatial superposition. The 2 Å cutoff applies to each head individually.

The [replicate table](replicates.csv) uses one denominator convention for old and new segments, excluding time zero from both. The crystal-head count gives frames where both ends meet the criterion, and the joint count also requires CypA to cover single-T-box DNA in the matched footprint analysis. The final source-head count repeats the 2 Å test relative to the chosen modeled starting coordinates. [Per-frame crystal distances](frames.csv)

The crystal references establish the locations of the parent ligands. The intact B probe's binary affinity, its equilibrium ternary occupancy and its residence time require other evidence. [TBXT crystal reference](https://www.rcsb.org/structure/5QS1) · [CypA crystal reference](https://www.rcsb.org/structure/9GHY)
