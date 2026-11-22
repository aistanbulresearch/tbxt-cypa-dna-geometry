# Can the fitted DNA duplex clear the entire MD snapshot?

The primary panel found a complete-assembly fit in 14 near-branch frames and two far-branch frames. All 14 near frames met both crystal-head RMSD cutoffs. In ten original B-on frames, the same fixed TBXT and intact probe admitted at least one proposal after deleting CypA, and every such proposal clashed with CypA when it was restored. The other 174 frames lacked a baseline fit. The 16 positive full-fit frames and ten paired losses came from only three of the eleven continuations.

## Primary 8,193-pose panel

Counts refer to post-start frames passing the specified tests. Entries in parentheses count the subset with both crystal-head RMSDs ≤2 Å. The baseline retains the entire probe and frame TBXT.

| Start | Frames | Baseline has fit | Complete assembly has fit | Paired panel loss | Baseline failure |
|---|---:|---:|---:|---:|---:|
| Original B-on | 60 (15) | 10 (5) | 0 (0) | 10 (5) | 50 (10) |
| B-near | 60 (50) | 14 (14) | 14 (14) | 0 (0) | 46 (36) |
| B-off | 60 (7) | 2 (0) | 2 (0) | 0 (0) | 58 (7) |
| Alternative B-on | 20 (20) | 0 (0) | 0 (0) | 0 (0) | 20 (20) |

The ten original B-on losses all came from replica 2 and involved one to six baseline placements per frame. The 14 B-near fits and the two B-off fits all came from their respective replica 3. Every baseline placement in those 16 frames also cleared CypA. Each trajectory inherits its starting geometry, and successive frames are serially correlated.

For the first and second halves of the Sobol bank, the numbers of baseline-positive frames were 4/10 for original B-on, 9/13 for near, 1/2 for far, and 0/0 for alternative B-on. Each half found complete fits in near and far and lost all baseline fits in original B-on. The two halves sampled different sets of individual frames. At the strict 3.0/2.6/3.0 Å TBXT/probe/CypA thresholds even the baseline failed in every frame.

The independent native control used the observed 8CDN protein chain A and its own two DNA strands. The observed pose recovered all 23 core Cα–interior sugar contacts. Its shortest heavy-atom distance, 2.320 Å, fell below both clearance thresholds. The primary panel recovered 21 passing native receptor placements, the closest at 0.631 Å interior sugar RMSD. The strict panel recovered zero even for the native receptor. The strict negative result thus has little power to diagnose incompatibility.

## Inspectable coordinates

The nearest retained near example is selected within its frame. In near replica 3 at 3.0 ps, proposal 6251 had an interior sugar RMSD of 2.874 Å, recovered 17/23 native contacts, and cleared TBXT, probe and CypA by 2.728, 11.995 and 6.312 Å respectively. Both ligand ends met their crystal-pose cutoff. In original B-on replica 2 at 0.5 ps, four baseline proposals passed. The nearest, 4951, recovered 19/23 contacts at 2.633 Å RMSD and had gaps of 2.623, 10.305 and 1.114 Å. None of the four cleared CypA. Coordinate files contain the full mapped heavy-atom solute and the placed duplex in native 8CDN axes. The intact ligand retains source atom names and CONECT records. Chain A is TBXT, P is CypA, X/Y/Z encode probe heads/linker, and C/D are the duplex.

## Reconstruction and replay

Use Python with NumPy and SciPy. The recorded reference run used Python 3.14.3, NumPy 2.4.6 and SciPy 1.18.0 and took about 41 seconds. From this directory run

```sh
python -B paired_fit.py --out replay-primary
python -B density_fit.py --primary replay-primary --out replay-density
```

The included inputs contain 11 solute-coordinate trajectories with atom IDs and cell boxes, four source PDBs and maps, native 8CDN, the archived partner-aligned head scores and the independent identity-overlay distance table. All coordinates are in Å. The script omits copied 0 ps frames. It shifts the ligand to the image nearest the centroid of TBXT pocket Cα atoms 119, 120 and 159 using its TBXT head, then shifts CypA relative to the ligand's CypA head using Cα atoms of residues 2–163. Each frame core is fitted to its source core. An inverse crystal-to-source fit places the entire frame in native 8CDN axes. Cores use chain-A Cα 45–208 except residues 110–121. The chain-A all-heavy-atom clearance includes the mapped helix and terminal atoms. The paired baseline omits CypA alone from the clearance test and keeps all coordinates identical.

The rotation transforms use row vectors: placed DNA = native DNA multiplied by R plus t. The rotation-vector components and translational increments are expressed along native 8CDN axes. Each pose first requires ≤3 Å unaligned RMSD on the 20 interior C1′ atoms and retention of at least 17/23 native core Cα/sugar distance pairs. All 486 duplex heavy atoms then undergo clearance tests. For each DNA atom, exact threshold checks can be restricted to component atoms initially no farther than its maximum panel displacement plus the threshold, by the triangle inequality. Direct KD-tree distance queries independently check the identity, the closest passing poses and a deterministic sample. For all 200 frames the three identity component minima reproduce the independent coordinate audit within 0.00003 Å.

The independently archived head scores are joined by family, start, replica and time. These RMSDs describe the pre-existing simulated ligand heads. Rigid DNA displacement leaves them fixed. The geometries use uniform heavy-atom clearance cutoffs and a reduced set of core–sugar distance constraints. The search samples a finite local rigid panel without induced protein or DNA relaxation. Accepted proposals and correlated frames measure existence under these rules. Panel losses are conditional on the searched register, cutoffs and observed snapshot.

## Files

| File | Contents |
|---|---|
| `protocol.md` | Criteria recorded before scoring the primary panel |
| `frames.csv` | Counts, classifications and closest IDs for all 200 frames, two cutoff choices and both halves plus their union |
| `summary.csv` | Counts by branch, replica, crystal-head subset, cutoff and panel |
| `selected-poses.csv` | Identity/nearest baseline/nearest full geometry and component distances by frame and gate |
| `transforms.npz` | Rotation and translation arrays indexed by proposal ID, plus interior sugar RMSD |
| `native-control.csv` | Receptor-only crystal self-test at both cutoffs |
| `representatives.csv`, `representatives/` | Coordinate-file keys and exported snapshot structures |
| `replay.json` | Input SHA-256 values, runtime and package versions |
| `density/` | Prespecified four-scramble expansion with framewise results and changes |
