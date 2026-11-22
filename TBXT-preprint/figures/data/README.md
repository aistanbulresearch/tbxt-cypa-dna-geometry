# Plotting inputs

These coordinate snapshots and tables let the figure renderer run from this folder without a project checkout. The table identifies the corresponding analyses in this release.

| Bundled file | Source |
|---|---|
| `B-cover.pdb`, `B-clear.pdb` | [Static structures](../../../TBXT-revision/structures/README.md) |
| `crystal-heads-replicates.csv` | [Crystal-head counts](../../source-data/occupancy/crystal-heads/replicates.csv) |
| `overlay-audit.csv` | [Duplex coordinate audit](../../source-data/occupancy/coordinate-audit.csv) |
| `helix-sensitivity.csv` | [Helix sensitivity table](helix-sensitivity.csv) |
| `paired-frames.csv` | [Paired placements](../../../TBXT-revision/paired/frames.csv) |
| `paired-density-summary.csv` | [Density sensitivity](../../../TBXT-revision/paired/density/summary.csv) (density sensitivity in Figure 4 legend) |

The coordinate snapshots preserve the processed chain identifiers. Chain A is TBXT core, H is the restored 8CDN helix, P is CypA, L is probe, and C/D are the same displaced duplex. Probe bond records come from the exported PDB connectivity. The renderer uses heavy atoms and projects both complexes using the same coordinates and scale.

The 12 base-pair positions are numbered along chain C. Author residue number `n` maps to position `n + 2` on chain C and `17 − n` on chain D. Figure 1A highlights the two residues of base pair 11, C-chain cytosine 9 and D-chain guanine 6. Seven and one heavy atoms, respectively, fall within 2.5 Å of CypA. The diamond marks the atom on terminal C-chain thymine 10 that gives the whole-duplex CypA minimum.
