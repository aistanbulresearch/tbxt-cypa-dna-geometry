# Base-pair sensitivity of the fixed-duplex footprint

The two exceptional footprints each involve one penultimate base pair. In the strict B example (seed 62008, draw 3069), seven contacted DNA atoms belong to C-chain cytosine 9 and one to D-chain guanine 6. All eight lie at position 11. The whole-duplex minimum of 0.705 Å is at terminal C-chain thymine 10 (position 12). The minimum within positions 2–11 is 1.470 Å. In the primary C example (seed 62040, draw 1261), one contacted atom belongs to C-chain guanine 0 and five to D-chain cytosine 15. All six lie at position 2. The whole-duplex minimum of 0.480 Å also occurs at that pair. Central-eight minima from the deposited banks are 3.566 Å for B and 3.973 Å for C.

| Pocket | Helix-inclusive assembly gates | Admitted | Whole duplex clears | Published inner 10, ≥2 nucleotides | Inner 10, ≥2 base pairs | Central 8, ≥2 nucleotides |
|---|---|---:|---:|---:|---:|---:|
| B | Primary | 110 | 78 | 15 | 14 | 12 |
| B | Strict | 22 | 15 | 1 | 0 | 0 |
| C | Primary | 121 | 96 | 1 | 0 | 0 |
| C | Strict | 56 | 47 | 0 | 0 | 0 |

## Definitions

Number the 12 base pairs along 8CDN author chain C. A residue numbered *n* belongs to position *n* + 2 on chain C and 17 − *n* on chain D. The measured purine N1–pyrimidine N3 distances across all complementary pairs appear in `register.csv`. The inner ten comprise positions 2–11 and the central eight comprise positions 3–10. Both masks describe position in the duplex. Each footprint requires at least three distinct DNA heavy atoms with nearest-CypA distance strictly below 2.5 Å. The published endpoint counts contacted nucleotides on opposite strands separately. The alternative inner-ten endpoint counts each pair once. The central-eight endpoint retains the published nucleotide rule. These two sensitivity endpoints were defined after inspecting the original results.

At probe/protein cutoffs of 2.2/2.6 Å (primary) or 2.6/3.0 Å (strict), each admitted draw clears the probe from both proteins and the two proteins from each other, including the transferred TBXT 212–221 helix. Whole-duplex clearance additionally requires the probe and CypA to clear DNA at their respective cutoffs. Fragment and helix clearance from the fixed DNA uses 2.6 Å in both rows. The deposited banks contain the 110 B and 481 C draws admitted by the original primary fragment-only screen. These counts describe that finite sample, with a fixed displaced 8CDN duplex for each pocket.

## Replay

Keep this directory beside `structures/` in the extracted supplement. From this directory run the following with Python 3.10+ and NumPy and SciPy installed.

```console
python replay_endpoints.py --out reproduced
```

The script reads the deposited PDB, duplex and accepted-bank coordinates. It recomputes the rigid helix alignment, all seven draw-dependent distance minima, admission and footprint endpoints. For all 591 seed/draw identifiers it checks the seven minima and the original interior atom/nucleotide counts against `structures/draw-metrics.csv`, the footprint against both archived flags, and the four published aggregate rows against `structures/counts.csv`. In each exceptional model it also checks every DNA atom's nearest-CypA distance and hit status against the separately exported PDB. The latter stores each coordinate to 0.001 Å precision. In the C example its central-eight minimum is 3.972 Å, compared with 3.973 Å from the bank.

| File | Contents |
|---|---|
| `endpoint-sensitivity.csv` | Four helix-inclusive count rows shown above |
| `draw-endpoints.csv` | All 591 original pocket/seed/draw identifiers, contact counts, endpoint indicators and admission flags. Footprint indicators apply before gating. Whole-clear flags already include the indicated admission gate. |
| `register.csv` | All 12 complementary pairs, anchor labels and distances in the 8CDN donor |
| `witness-atoms.csv` | One row per DNA atom strictly within 2.5 Å of CypA, including terminal positions, its nearest CypA atom and distance |
| `witness-nucleotides.csv` | All 24 nucleotides per witness with contact counts and minimum distances |
| `witness-minima.csv` | Atoms, nucleotides and pairs hit for each mask, closest atom identities, and minima from the bank and exported PDB |
