"""Render the four revised TBXT figures from the adjacent data directory."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from scipy.interpolate import make_interp_spline
from scipy.spatial import cKDTree


HERE = Path(__file__).resolve().parent
INK = "#18363b"
MUTED = "#486268"
GRID = "#d6e1e2"
PAIR_COLOR = "#59267d"
HIT_COLOR = "#9b3425"
COLORS = {
    "A": "#087266",
    "H": "#d49316",
    "P": "#b75231",
    "L": "#a3196f",
    "DNA": "#4b62a8",
}
GROUPS = [
    ("monomer", "B-on", "Original B-on", "B-on"),
    ("monomer", "B-near", "B-near", "B-near"),
    ("monomer", "B-off", "B-far (B-off)", "B-off"),
    ("ref", "B-on", "Alternative B-on", "ref-B-on"),
]
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "text.color": INK,
    "axes.labelcolor": INK,
    "axes.edgecolor": GRID,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "savefig.facecolor": "white",
    "pdf.fonttype": 42,
})


def save(fig, dest, name):
    for extension in ("png", "pdf"):
        fig.savefig(dest / f"{name}.{extension}", dpi=320, facecolor="white")
    plt.close(fig)


def read_pose(path):
    atoms, partner_serials = {}, []
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            if line[76:78].strip() not in {"H", "D"}:
                index = int(line[6:11])
                xyz = np.array([float(line[i:i + 8]) for i in (30, 38, 46)])
                atoms[index] = (line[21], int(line[22:26]), line[12:16].strip(), xyz)
        elif line.startswith("CONECT"):
            ids = [int(line[i:i + 5]) for i in range(6, len(line.rstrip()), 5)]
            partner_serials.extend((ids[0], neighbor) for neighbor in ids[1:])
    chains = {
        chain: np.array([xyz for c, _, _, xyz in atoms.values() if c == chain])
        for chain in "AHPLCD"
    }
    curves = {
        chain: np.array([
            xyz for c, _, atom, xyz in atoms.values() if c == chain and atom == marker
        ])
        for chain, marker in [("A", "CA"), ("H", "CA"), ("P", "CA"), ("C", "C4'"), ("D", "C4'")]
    }
    probe_bonds = {
        tuple(sorted((i, j))) for i, j in partner_serials
        if i in atoms and j in atoms and atoms[i][0] == atoms[j][0] == "L"
    }
    if not probe_bonds or any(array.size == 0 for array in curves.values()):
        raise ValueError(f"Missing chain markers or ligand bonds in {path}")
    ligand = [np.array([atoms[i][3], atoms[j][3]]) for i, j in sorted(probe_bonds)]
    chains["DNA"] = np.vstack([chains["C"], chains["D"]])
    dna_atoms = {
        (chain, residue, atom): xyz
        for chain, residue, atom, xyz in atoms.values()
        if chain in "CD"
    }
    interior = {}
    for chain in "CD":
        residues = sorted({residue for c, residue, _ in dna_atoms if c == chain})
        interior.update(
            (key, xyz) for key, xyz in dna_atoms.items()
            if key[0] == chain and key[1] in residues[1:-1]
        )
    distances = cKDTree(chains["P"]).query(list(interior.values()))[0]
    contacts = {key: xyz for (key, xyz), distance in zip(interior.items(), distances) if distance < 2.5}
    return chains, curves, ligand, dna_atoms, contacts


def add_pose(ax, pose, basis, bounds):
    _, curves, ligand, dna_atoms, contacts = pose
    hit_residues = {key[:2] for key in contacts}
    segments = []

    def trace(coordinates, color, width, smooth=False):
        if smooth and len(coordinates) > 3:
            original = np.arange(len(coordinates))
            dense = np.linspace(0, original[-1], 4 * len(original))
            coordinates = make_interp_spline(original, coordinates)(dense)
        projected = coordinates @ basis
        for first, second in zip(projected, projected[1:]):
            depth = (first[2] + second[2]) / 2
            fraction = (depth - bounds[0][2]) / (bounds[1][2] - bounds[0][2])
            color_rgb = np.array(to_rgb(color))
            tint = color_rgb + (1 - color_rgb) * 0.32 * (1 - fraction)
            segments.append((depth, [first[:2], second[:2]], tint, width))

    for chain, coordinates in curves.items():
        color = COLORS[chain if chain in COLORS else "DNA"]
        max_step = 5 if chain in "AHP" else 9
        breaks = np.flatnonzero(np.linalg.norm(np.diff(coordinates, axis=0), axis=1) > max_step) + 1
        for section in np.split(coordinates, breaks):
            trace(section, color, 2.15, smooth=True)
    for chain in "CD":
        keys = [key for key in dna_atoms if key[0] == chain]
        coordinates = np.array([dna_atoms[key] for key in keys])
        for i, j in sorted(cKDTree(coordinates).query_pairs(1.93)):
            if keys[i][:2] == keys[j][:2] and keys[i][:2] in hit_residues:
                trace(coordinates[[i, j]], PAIR_COLOR, 1.3)
            else:
                trace(coordinates[[i, j]], COLORS["DNA"], 0.72)
    for coordinates in ligand:
        trace(coordinates, COLORS["L"], 2.7)
    segments.sort(key=lambda item: item[0])
    ax.add_collection(LineCollection(
        [item[1] for item in segments],
        colors=[item[2] for item in segments],
        linewidths=[item[3] for item in segments],
        capstyle="round",
    ))
    if len(contacts):
        to_plot = np.array(list(contacts.values())) @ basis
        ax.scatter(to_plot[:, 0], to_plot[:, 1], facecolors="white", edgecolors=HIT_COLOR,
                   s=16, linewidth=.8, zorder=4)
    ax.set(
        xlim=(bounds[0][0] - 2, bounds[1][0] + 2),
        ylim=(bounds[0][1] - 2, bounds[1][1] + 2),
        aspect="equal",
    )
    ax.set_axis_off()


def annotate_cover(ax, pose, basis):
    chains, _, _, dna_atoms, contacts = pose
    c_hits = sum(key[:2] == ("C", 9) for key in contacts)
    d_hits = sum(key[:2] == ("D", 6) for key in contacts)
    if (c_hits, d_hits, len(contacts)) != (7, 1, 8):
        raise ValueError("Expected all eight interior hits on C9 and G6 at base pair 11")
    pairing_atoms = np.array([dna_atoms["C", 9, "N3"], dna_atoms["D", 6, "N1"]])
    midpoint = pairing_atoms.mean(0) @ basis
    ax.annotate(
        f"Base pair 11 of 12\nC:C9 ({c_hits} hits)\nD:G6 ({d_hits} hit)",
        xy=midpoint[:2], xytext=(.015, .94), textcoords="axes fraction",
        ha="left", va="top", fontsize=8.2, color=PAIR_COLOR, linespacing=1.25,
        arrowprops={"arrowstyle": "-", "color": PAIR_COLOR, "lw": .8, "relpos": (.65, 0)},
        zorder=6,
    )
    dna_keys = list(dna_atoms)
    dna_coordinates = np.array(list(dna_atoms.values()))
    cypa_distances = cKDTree(chains["P"]).query(dna_coordinates)[0]
    nearest = np.argmin(cypa_distances)
    if dna_keys[nearest][:2] != ("C", 10):
        raise ValueError("Expected the closest DNA atom to CypA on terminal C:T10")
    x, y, _ = dna_coordinates[nearest] @ basis
    ax.scatter(x, y, marker="D", s=20, color=INK, edgecolor="white", linewidth=.7, zorder=5)


def molecular_figure(data, dest):
    poses = [read_pose(data / f"{pose}.pdb") for pose in ("B-cover", "B-clear")]
    cover, clear = (pose[0] for pose in poses)
    for chain in ("A", "H", "DNA"):
        np.testing.assert_array_equal(cover[chain], clear[chain])
    vertical = np.linalg.svd(cover["DNA"] - cover["DNA"].mean(0))[2][0]
    right = clear["P"].mean(0) - cover["P"].mean(0)
    right -= vertical * (right @ vertical)
    right /= np.linalg.norm(right)
    basis = np.column_stack([right, vertical, np.cross(right, vertical)])
    coordinates = np.vstack([pose[key] for pose in (cover, clear) for key in ("A", "H", "DNA", "L", "P")])
    projected = coordinates @ basis
    bounds = (projected.min(0), projected.max(0))

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.4))
    fig.subplots_adjust(left=.025, right=.975, bottom=.21, top=.89, wspace=.1)
    titles = ["CypA overlaps one interior base pair", "Whole duplex clears"]
    for letter, ax, pose, title in zip("AB", axes, poses, titles):
        add_pose(ax, pose, basis, bounds)
        ax.set_title(f"{letter}   {title}", loc="left", fontweight="bold", fontsize=11, pad=7)
        chains = pose[0]
        distances = {key: cKDTree(chains["DNA"]).query(chains[key])[0].min() for key in "AHLP"}
        cypa_label = f"Minimum to CypA   {distances['P']:.3f} Å"
        if letter == "A":
            annotate_cover(ax, pose, basis)
            cypa_label += "   ◆ Terminal C:T10"
        ax.text(.5, -.035, cypa_label, transform=ax.transAxes,
                ha="center", va="top", fontsize=9.2, fontweight="bold", color=INK)
        ax.text(
            .5, -.115,
            f"Core {distances['A']:.2f} Å   ·   Helix {distances['H']:.2f} Å   ·   Probe {distances['L']:.2f} Å",
            transform=ax.transAxes, ha="center", va="top", fontsize=8.3, color=MUTED,
        )
    handles = [
        Line2D([], [], color=COLORS[key], lw=3.4, label=label)
        for key, label in [("A", "TBXT core"), ("H", "8CDN helix"), ("P", "CypA"),
                           ("L", "Probe B"), ("DNA", "Registered duplex")]
    ]
    handles.append(Line2D([], [], marker="o", mfc="white", mec=HIT_COLOR, ls="none",
                          markersize=4, label="Interior hit atom"))
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, .013), ncol=6,
               frameon=False, columnspacing=1.05, fontsize=8.2)
    save(fig, dest, "figure1-structures")


def select_group(table, group):
    family, start, _, _ = group
    return table[(table.family == family) & (table.start == start)]


def dynamics(data, dest):
    counts = pd.read_csv(data / "crystal-heads-replicates.csv")
    distances = pd.read_csv(data / "overlay-audit.csv")
    fig = plt.figure(figsize=(8.6, 7.05))
    grid = fig.add_gridspec(2, 6, height_ratios=[3, 1.35], hspace=.52,
                           left=.205, right=.98, top=.89, bottom=.09, wspace=.97)
    cmap = LinearSegmentedColormap.from_list("frequency", ["#f1f5f4", "#aed5cd", "#15685f"])
    for i, title in enumerate(["A   Parent trajectories", "B   Continuations"]):
        ax = fig.add_subplot(grid[0, 3*i:3*(i+1)])
        subset = counts[(counts.phase == "parent") == (i == 0)]
        ordered = pd.concat([select_group(subset, group).sort_values("rep") for group in GROUPS])
        den = ordered.post_start.to_numpy(int)
        values = ordered[["occupied", "two_crystal_heads", "covered_and_two_crystal_heads"]].to_numpy(int)
        fraction = values / den[:, None]
        x_edges = np.arange(fraction.shape[1] + 1) - .5
        y_edges = np.arange(fraction.shape[0] + 1) - .5
        ax.pcolormesh(x_edges, y_edges, fraction, cmap=cmap, vmin=0, vmax=1,
                      shading="flat", rasterized=False)
        ax.set_ylim(y_edges[-1], y_edges[0])
        for row, column in np.ndindex(values.shape):
            ax.text(column, row, f"{values[row, column]}/{den[row]}",
                    color="white" if fraction[row, column] > .71 else INK, va="center", ha="center", fontsize=8.5)
        labels = [
            f"{name}  ·  R{rep}"
            for group in GROUPS for name in [group[2]]
            for rep in select_group(subset, group).sort_values("rep").rep
        ]
        ax.set_xticks(range(3), ["Interior\nhit", "Both\nheads", "Joint"])
        ax.set_yticks(range(len(labels)), labels if i == 0 else [])
        ax.xaxis.tick_top()
        ax.tick_params(length=0, pad=4)
        ax.set_title(title, loc="left", fontweight="bold", pad=10, fontsize=10)
        for boundary in np.cumsum([len(select_group(subset, group)) for group in GROUPS])[:-1]:
            ax.axhline(boundary - .5, color="white", lw=1.8)
        for spine in ax.spines.values():
            spine.set_visible(False)

    rng = np.random.default_rng(20)
    ymax = 3.55
    names = ["TBXT", "Probe", "CypA"]
    fields = ["TBXT_whole_min_A", "probe_whole_min_A", "CypA_whole_min_A"]
    fig.text(.205, .357, "C   Whole-duplex distances in the continuations", fontweight="bold", fontsize=10)
    for j, (component, field, color) in enumerate(zip(names, fields, [COLORS["A"], COLORS["L"], COLORS["P"]])):
        ax = fig.add_subplot(grid[1, 2*j:2*(j+1)])
        for index, group in enumerate(GROUPS):
            frame_values = select_group(distances, group)[field].to_numpy()
            offset = rng.uniform(-.18, .18, len(frame_values))
            ax.scatter(frame_values, index + offset, s=9, color=color, alpha=.49, edgecolor="none")
        ax.set_title(component, fontsize=9, pad=3)
        ax.set_ylim(ymax, -.55)
        ax.set_xlim(0, max(2.125, distances[fields].max().max() * 1.04))
        ax.set_xticks([0, 5, 10])
        ax.set_yticks(range(4), [group[2] for group in GROUPS] if j == 0 else [])
        ax.set_axisbelow(True)
        ax.grid(axis="x", linewidth=.6, color=GRID)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0, labelsize=8)
        if j == 1:
            ax.set_xlabel("Minimum heavy-atom distance (Å)", fontsize=8.5)
    save(fig, dest, "figure2-dynamics")


def designs(data, dest):
    source = pd.read_csv(data / "helix-sensitivity.csv")
    rows = [
        (getattr(row, f"{variant}_dual_fit"), getattr(row, f"{variant}_whole_duplex_clear"),
         getattr(row, f"{variant}_CypA_interior_hit"))
        for row in source.itertuples() for variant in ("original", "helix")
    ]
    fits, clears, hits = np.asarray(rows).T
    positions = np.asarray([0, 1, 2.5, 3.5, 5, 6, 7.5, 8.5])
    labels = [
        f"{row.probe}  {'primary' if row.probe_clearance_A == 2.2 else 'strict'}  ·  {variant}"
        for row in source.itertuples() for variant in ("Core", "+ helix")
    ]
    palette = ["#738b91", "#127a6a"] * (len(rows) // 2)
    fig, axes = plt.subplots(1, 3, figsize=(8.6, 4.64), sharey=True,
                             gridspec_kw={"width_ratios": [1, 1, 1]})
    fig.subplots_adjust(left=.216, right=.98, top=.84, bottom=.19, wspace=.27)
    panels = [
        ("A   Admitted", "Admitted models", fits, np.full_like(fits, 65536), [0, 200, 400], 700),
        ("B   Clears probe and CypA", "% of admitted", 100 * clears / fits, fits, [0, 25, 50, 75, 100], 125),
        ("C   CypA hits interior", "% of admitted", 100 * hits / fits, fits, [0, 5, 10, 15], 21),
    ]
    for j, (ax, (title, xlabel, plotted, denominators, ticks, xmax), numerators) in enumerate(
        zip(axes, panels, (fits, clears, hits))
    ):
        ax.barh(positions, plotted, color=palette, height=.64, zorder=2)
        for ypos, value, numerator, denominator in zip(positions, plotted, numerators, denominators):
            label = f"{numerator:,}/{denominator:,}"
            ax.text(value + xmax * .023, ypos, label, fontsize=7.7, va="center")
        ax.set_title(title, loc="left", fontsize=9.3, fontweight="bold", pad=13)
        ax.set(xlim=(0, xmax), ylim=(9.22, -.72), xticks=ticks, xlabel=xlabel)
        ax.xaxis.label.set_size(8.3)
        ax.set_yticks(positions, labels)
        ax.tick_params(axis="y", labelleft=j == 0, length=0, labelsize=8.4, pad=6)
        ax.tick_params(axis="x", length=0, pad=5, labelsize=8.1)
        ax.spines["left"].set_visible(False)
        ax.grid(axis="x", color=GRID, linewidth=.65)
        ax.set_axisbelow(True)
        for position in (1.75, 4.25, 6.75):
            ax.axhline(position, color=GRID, linewidth=.65)
    fig.text(.216, .04, "Primary 2.2/2.6 Å   ·   Strict 2.6/3.0 Å   ·   Probe–protein / protein–protein",
             fontsize=8.4, color=MUTED)
    save(fig, dest, "figure3-designs")


def paired_fit(data, dest):
    frames = pd.read_csv(data / "paired-frames.csv").query("gates == 'primary' and panel == 'union'")
    maximum = int(frames[["baseline_accepted", "full_accepted"]].max().max())
    upper = max(1, maximum)
    fig, axes = plt.subplots(2, 2, figsize=(8.6, 6.25), sharex=True, sharey=True)
    fig.subplots_adjust(left=.124, right=.973, bottom=.16, top=.9, hspace=.38, wspace=.23)
    if len(frames) != 200 or not (frames.full_accepted <= frames.baseline_accepted).all():
        raise ValueError("Expected 200 valid same-frame baseline/complete pairs")
    for ax, group, letter in zip(axes.flat, GROUPS, "ABCD"):
        subset = frames[frames.branch == group[3]]
        total = len(subset)
        baseline_failures = (subset.baseline_accepted == 0).sum()
        complete_passes = (subset.full_accepted > 0).sum()
        losses = total - baseline_failures - complete_passes
        ax.set_title(f"{letter}   {group[2]}", loc="left", fontsize=10, fontweight="bold", pad=8)
        ax.plot([0, upper], [0, upper], color="#7f979c", lw=1, linestyle=(0, (3, 2)), zorder=1)
        for flag, color, marker in [(True, "#b95131", "^"), (False, "#66738b", "o")]:
            points = (subset[subset.dual_head.astype(int) == int(flag)]
                      .groupby(["baseline_accepted", "full_accepted"]).size().reset_index(name="frames"))
            ax.scatter(points.baseline_accepted, points.full_accepted,
                       s=15 + (8 if flag else 22) * np.sqrt(points.frames), marker=marker,
                       facecolors=color if flag else "none", linewidth=.85,
                       edgecolors="white" if flag else color, alpha=.89, zorder=3)
        ax.text(.04, .95, f"{complete_passes}/{total} retained fit\n"
                f"{losses}/{total} panel loss\n{baseline_failures}/{total} baseline failure",
                ha="left", va="top", transform=ax.transAxes, fontsize=8.1,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": .91, "pad": 3})
        ax.set(xlim=(-.055*upper, 1.09*upper), ylim=(-.055*upper, 1.09*upper))
        ax.xaxis.set_major_locator(MaxNLocator(6, integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(6, integer=True))
        ax.tick_params(length=2, labelsize=8)
        ax.set_axisbelow(True)
        ax.grid(axis="both", color=GRID, linewidth=.6)
    fig.supxlabel("Accepted duplex placements without CypA", x=.55, y=.095, fontsize=9.4)
    fig.supylabel("Accepted duplex placements with CypA", x=.037, fontsize=9.4)
    handles = [
        Line2D([], [], marker=marker, markerfacecolor=fill, markeredgecolor=color,
               label=label, ls="none", markersize=6)
        for marker, fill, color, label in [
            ("^", "#b95131", "#b95131", "Both heads ≤2 Å"),
            ("o", "none", "#66738b", "Head criterion unmet"),
        ]
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.55, .013), ncol=2,
               frameon=False, fontsize=8.5)
    save(fig, dest, "figure4-paired")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=HERE / "data")
    parser.add_argument("--dest", type=Path, default=HERE)
    parser.add_argument("--figures", nargs="+", type=int, choices=range(1, 5), default=range(1, 5))
    options = parser.parse_args()
    options.dest.mkdir(parents=True, exist_ok=True)
    functions = [molecular_figure, dynamics, designs, paired_fit]
    for figure in options.figures:
        functions[figure - 1](options.data, options.dest)


if __name__ == "__main__":
    main()
