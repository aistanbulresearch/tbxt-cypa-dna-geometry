# Distinct cyclophilin A orientations around a modeled brachyury probe permit different fits of a single T-box duplex

Ozge A. Cavus and Aysegul Kuskucu

[Read the manuscript (PDF)](TBXT-preprint/TBXT-preprint.pdf) · [Word](TBXT-preprint/TBXT-preprint.docx) · [Markdown](TBXT-preprint/manuscript.md) · [Evidence index](TBXT-preprint/evidence.md)

Brachyury (TBXT) is a transcription factor that recognizes a T-box DNA sequence and is a dependency in chordoma. We asked whether a molecule with one end bound to TBXT and the other to cyclophilin A (CypA) could position CypA in the space occupied by DNA. The question was inspired by daraxonrasib, which recruits CypA to active RAS and restricts access by RAS effectors. Crystal structures of TBXT supplied two possible starting points near its DNA interface, pockets B and C. [Background and references](TBXT-preprint/manuscript.md#introduction)

**The main finding is that the same pocket-B probe can place CypA in two different positions while preserving both modeled binding poses.** In one model CypA intersects the penultimate base pair of a 12-base-pair duplex containing a single T-box element. In another, the entire duplex fits alongside TBXT, the probe and CypA. The models establish a geometric possibility. How often either arrangement forms in solution, and whether recruitment weakens DNA binding, remain open questions. [Overlapping model](TBXT-revision/structures/B-cover.pdb) · [Whole-duplex fit](TBXT-revision/structures/B-clear.pdb)

## How we tested the geometry

1. **Construct two intact probes.** We linked the crystallographic TBXT fragments from pocket B (PDB 5QS1) and pocket C (5QSE) to a CypA-binding head derived from 9GHY. For each probe we sampled 65,536 linker conformations, placed each head in its own protein pocket and checked for collisions. We repeated the geometry test after adding the resolved volume of TBXT residues 212–221 from the DNA-bound structure 8CDN. This helix enters the DNA minor groove and is missing from both fragment structures. [Static models and methods](TBXT-revision/structures/README.md)
2. **Distinguish a local overlap from a whole-duplex fit.** A local overlap required at least three DNA heavy atoms from at least two interior nucleotides within 2.5 Å of CypA. Two nucleotides can be partners in a single base pair. For a whole-duplex fit, every part of the DNA, including both ends, had to clear TBXT, the intact probe and CypA. [Base-pair definitions and sensitivity tests](TBXT-revision/base-pairs/README.md)
3. **Test the effect of CypA in the same simulation frame.** We ran eleven short, DNA-free, explicit-solvent OpenMM continuations from four prepared pocket-B geometries and analyzed 200 post-start frames. In each frame we tested 8,193 rigid placements of the 8CDN duplex for near-native contacts and steric clearance. We tested exactly the same placements again with CypA removed, keeping TBXT and the probe fixed. A denser panel used 32,769 placements per frame. [Paired test and replay instructions](TBXT-revision/paired/README.md) · [Denser panel](TBXT-revision/paired/density/README.md)

## What we found

### The DNA-binding helix changes the pocket-C result

Each row counts modeled assemblies after adding the helix. “Admitted” means the probe and proteins passed the assembly clearance filters. The last two columns test different properties of the same chosen duplex placement, so their counts need not sum to the admitted total. [Figure 3 source table](TBXT-preprint/figures/data/helix-sensitivity.csv)

| Pocket | Clearance filters | Admitted assemblies | Whole duplex fits | CypA overlaps interior DNA |
|---|---|---:|---:|---:|
| B | Primary | 110 | 78 | 15 |
| B | Stricter | 22 | 15 | 1 |
| C | Primary | 121 | 96 | 1 |
| C | Stricter | 56 | 47 | 0 |

Adding the helix preserved the B counts at both filter settings. For C, it reduced the primary overlap count from 26 of 481 admitted assemblies to 1 of 121. The lone strict B overlap contacted both strands at **base pair 11 of 12**. The surviving primary C overlap contacted **base pair 2 of 12**. Both positions are one pair in from a duplex end. At the stricter settings, neither pocket produced an overlap meeting the same criterion within the central eight base pairs. Direct obstruction of the bases responsible for T-box sequence recognition remains unestablished. [Base-pair evidence](TBXT-revision/base-pairs/README.md)

### The pocket-B outcome depends on where CypA starts

With the primary panel, **14 of 60 near-start frames** accommodated DNA with CypA present and both probe heads within 2 Å of their crystal poses. In **10 of 60 original covering-start frames**, a near-native DNA placement fit only after CypA was removed. Each result came from a single replica of its respective start. The denser panel preserved this distinction and found examples of both outcomes in the alternative covering start. In 166 of all 200 frames, even the denser panel found no fit after CypA removal, so those frames cannot establish CypA-specific exclusion. These are counts of correlated snapshots from chosen starts. They cannot estimate how frequently the probe covers DNA in solution. [Primary frame data](TBXT-revision/paired/frames.csv) · [Denser-panel results by replica](TBXT-revision/paired/density/summary.csv)

### Independent docking leaves ternary binding unresolved

DeepTernary received the intact probes and independently positioned protein inputs, without DNA or the proposed ternary pose. Of 40 predictions per probe, 6 B models and 4 C models placed both heads near their reference poses and overlapped the DNA footprint when checked afterward. All ten had a protein–protein clash. The same independent-input workflow also struggled to recover a known ternary complex (PDB 5T35), which limits the interpretation of this negative result. Local AutoDock Vina optimization retained prepositioned heads, but scoring controls gave inconsistent ligand rankings. Binding of each intact probe to its two partners, ternary complex formation and a CypA-dependent change in DNA binding all await experimental measurement. [Docking, scoring and calibration records](TBXT-preprint/evidence.md#local-scoring-and-independent-searches)

## Inspect or reproduce the results

| To check | Start here |
|---|---|
| Compare the two pocket-B orientations in Figure 1 | [Overlapping coordinates](TBXT-revision/structures/B-cover.pdb) · [Clearing coordinates](TBXT-revision/structures/B-clear.pdb) · [Chain identifiers and replay instructions](TBXT-revision/structures/README.md) |
| Recalculate the helix and base-pair results | [Rescore the 591 previously admitted assemblies](TBXT-revision/structures/README.md) · [Recount the affected base pairs](TBXT-revision/base-pairs/README.md) |
| Repeat the with/without-CypA comparison | [Saved frames, duplex placements and paired analysis](TBXT-revision/paired/README.md) |
| Recreate the four manuscript figures | [Figure captions](TBXT-preprint/figures/captions.md) · [Plotting inputs](TBXT-preprint/figures/data/README.md) · [Rendering script](TBXT-preprint/figures/make_figures.py) |
| Trace a manuscript claim to its underlying data | [Evidence index](TBXT-preprint/evidence.md) · [Historical data manifest with SHA-256 checksums](TBXT-preprint/source-data/manifest.csv) |

The deposited coordinates and scripts support replay of the structural comparisons above. The historical docking and simulation tables support recounting frames and reranking models. Regenerating the original trajectories, docking predictions or full conformer banks requires additional upstream inputs and software, as detailed in the [evidence index](TBXT-preprint/evidence.md).
