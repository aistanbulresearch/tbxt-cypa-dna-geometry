# Four-scramble density sensitivity

Expanding the original 8,193 poses to 32,769 rescued baselines in eight more frames. Six of the eight lost every baseline proposal when CypA was included. Two gained a full fit. All ten original framewise losses and all 16 original full-fit classifications persisted. The expanded search left 166 of 200 frames without a matched baseline fit.

| Start | Frames | Baseline fits, original → expanded | Full fits, original → expanded | Expanded panel losses | Expanded baseline failures |
|---|---:|---:|---:|---:|---:|
| Original B-on | 60 | 10 → 15 | 0 → 0 | 15 | 45 |
| B-near | 60 | 14 → 14 | 14 → 14 | 0 | 46 |
| B-off | 60 | 2 → 3 | 2 → 3 | 0 | 57 |
| Alternative B-on | 20 | 0 → 2 | 0 → 1 | 1 | 18 |

Among the 92 frames with both partner-aligned crystal-head RMSDs ≤2 Å, the expanded baseline/full frame counts were 8/0 of 15 original B-on, 14/14 of 50 B-near, 0/0 of seven B-off, and 2/1 of 20 alternative B-on. Each branch's baseline-positive frames came from a single replica (original B-on 2, B-near 3, B-off 3, alternative B-on 2). These are descriptive counts over correlated frames.

The clearest new qualification concerns the alternative B-on geometry. In replica 2 at 1.5 ps, its single accepted baseline proposal (ID 14785) approached CypA to 2.128 Å and failed the 2.6 Å gate. At 4.0 ps of the same replica, only one of two accepted baseline proposals cleared CypA. This full fit (ID 30702) displaced the interior sugars by 1.737 Å and recovered 22/23 native contacts. Its TBXT/probe/CypA gaps were 2.60146, 10.97335 and 2.65351 Å. The TBXT margin above cutoff was 0.00146 Å. Both frames met the crystal-head criteria. The observed examples are sensitive to search density and to small changes in clearance or coordinates.

The native 8CDN receptor-only control gave 21, 18, 22 and 31 accepted proposals in seeds 6280–6283, respectively. Identity failed its primary clearance, so 92 unique poses passed in the union. The original seed reproduced all 200 per-frame primary baseline counts and all 200 full counts.

The original search and this expansion share all predicates and the same 200 mapped snapshots. In `transforms.npz`, global ID 0 is identity. For positive within-seed proposal number j (1–8192), global ID = 8192 × (seed − 6280) + j. Scale and orientation follow the primary protocol. `seed-counts.csv` includes identity in each seed. `frames.csv` counts identity only once, gives both primary and expanded results and records whether the frame classification changed. `summary.csv` counts frames by branch, replica and head subset; `transitions.csv` compares the two classifications. `selected-poses.csv` contains exact gap and recovery data for the closest available passing baseline/full pose in each frame. `representatives.csv` identifies the four PDB files representing changed branch/status combinations. All representative coordinate files were read back and their three whole-duplex minima checked to within 0.002 Å of their full-precision values.

From the parent directory run `python -B density_fit.py --primary . --out replay-density`. The four-seed calculation took about 106 seconds in the recorded environment. The expanded panel remains a finite set of rigid placements. Its failure in an individual frame does not establish exclusion throughout the continuous space of near-native DNA conformations.
