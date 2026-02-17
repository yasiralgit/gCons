"""
gCons_fonctionsTK.py
--------------------
Core logic for the gCons genomic consensus pipeline.

This module handles:
  - Invocation of external tools (redoak, gkampi) via subprocess
  - Common k-mer extraction and positional indexing
  - Common zone detection across multiple genomes
  - FASTA consensus output generation
  - Interactive matplotlib visualisation embedded in a Tkinter window

Dependencies:
  - redoak and gkampi binaries must be present in the working directory
  - biopython, customtkinter, matplotlib, numpy

Author: AL-YOUSSFI Yasir
"""

import os
import csv
import time
from subprocess import call

import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from customtkinter import *
from Bio import SeqIO


# ---------------------------------------------------------------------------
# Global state
# ---------------------------------------------------------------------------

_submit_count: int = 0  # Tracks the number of submitted analysis runs


# ---------------------------------------------------------------------------
# Debug helpers  (not called in production — safe to remove)
# ---------------------------------------------------------------------------

def _debug_print_matrix(matrix: list) -> None:
    """Pretty-print every element of a genome matrix."""
    print("[")
    for i, genome in enumerate(matrix):
        print(f"\t genome {i + 1}")
        print("\t [")
        for element in genome:
            print(f"\t\t {element}")
        print("\t ]")
    print("]")


def _debug_print_zone(matrix: list, zone_index: int) -> None:
    """Pretty-print a single zone across all genomes."""
    print("[")
    for i, genome in enumerate(matrix):
        print(f"\t genome {i + 1}  zone {zone_index}")
        print("\t [")
        print(f"\t\t {genome[zone_index]}")
        print(f"\t\t length: {len(genome[zone_index])}")
        print("\t ]")
    print("]")


def _debug_print_all_zones(matrix: list) -> None:
    """Pretty-print every zone across all genomes simultaneously."""
    for i in range(len(matrix[0])):
        _debug_print_zone(matrix, i)


def _debug_print_reference(matrix: list) -> None:
    """Pretty-print the reference genome (index 0) zone by zone."""
    print("[")
    print("\t genome 1")
    print("\t [")
    for j, zone in enumerate(matrix[0]):
        print(f"\n\t zone {j}")
        print(f"\t\t {zone}")
        print(f"\t\t length: {len(zone)}")
    print("\t ]")
    print("]")


# ---------------------------------------------------------------------------
# External tool wrappers
# ---------------------------------------------------------------------------

def _run_gkampi(k: int, fasta_file: str, output_index: int) -> None:
    """
    Run gkampi on a single FASTA file.

    Args:
        k: k-mer length.
        fasta_file: Path to the input FASTA file.
        output_index: Integer suffix used to name the output CSV (result<n>.csv).
    """
    cmd = [
        "./gkampi", fasta_file,
        "--output", f"result{output_index}.csv",
        "--pos", "--column",
        "--kmer-len", str(k),
        "--quiet",
    ]
    call(cmd)


def _run_redoak(k: int, fasta_files: list) -> None:
    """
    Run redoak on a set of FASTA files to produce the shared k-mer output.

    Args:
        k: k-mer length.
        fasta_files: List of paths to input FASTA files.
    """
    cmd = ["./redoak"]
    for f in fasta_files:
        cmd += ["--genome", f]
    cmd += ["--output", "test.txt", "--kmer", str(k)]
    call(cmd)


def _grep_common_kmers(alpha: int, fasta_files: list) -> None:
    """
    Filter redoak output to retain only k-mers present in at least *alpha* genomes.

    Writes matching lines to commun.txt.

    Args:
        alpha: Minimum number of genomes that must share a k-mer.
        fasta_files: Used solely to determine the upper bound of the grep pattern.
    """
    pattern_parts = [f"\\({v}\\)" for v in range(alpha, len(fasta_files) + 1)]
    pattern = "|".join(pattern_parts)
    os.system(f'grep -E "{pattern}" test.txt > commun.txt')


# ---------------------------------------------------------------------------
# Data structure builders
# ---------------------------------------------------------------------------

def _load_common_kmers() -> list:
    """
    Parse commun.txt and return the list of k-mer sequences shared across genomes.

    Returns:
        List of k-mer strings (ACGT characters only, up to the first space).
    """
    kmers = []
    with open("commun.txt", "r") as fh:
        for line in fh:
            kmer = ""
            for ch in line:
                if ch == " " or ch not in "ACGT":
                    break
                kmer += ch
            kmers.append(kmer)
    print(f"Common k-mers found: {len(kmers)}")
    return kmers


def _build_position_table(k: int, fasta_files: list) -> tuple:
    """
    Build the position table mapping each common k-mer to its positions in every genome.

    For each genome, gkampi is called and its CSV output is parsed. The reference genome
    (index 0) also populates a reverse-lookup dictionary (position → k-mer sequence).

    Args:
        k: k-mer length.
        fasta_files: Ordered list of FASTA file paths; index 0 is the reference.

    Returns:
        Tposition: 3-D list [genome][kmer_index][occurrence] = genomic position.
        ref_kmer_map: Dict mapping each reference position to its k-mer sequence.
    """
    common_kmers = _load_common_kmers()
    Tposition = []
    ref_kmer_map = {}

    for idx, fasta_file in enumerate(fasta_files):
        genome_num = idx + 1
        _run_gkampi(k, fasta_file, genome_num)

        # One sub-list per common k-mer for this genome
        Tposition.append([[] for _ in common_kmers])

        t0 = time.time()
        with open(f"result{genome_num}.csv", newline="") as fh:
            rows = list(csv.reader(fh, delimiter=" "))

        if genome_num == 1:
            # Build the reference reverse-lookup: position (int) → k-mer string
            for row in rows:
                ref_kmer_map[int(row[1][1:])] = row[0]

        # Single-pass match: rows and common_kmers are both k-mer-ordered
        row_idx = 0
        for kmer_idx, kmer in enumerate(common_kmers):
            while row_idx < len(rows) and rows[row_idx][0] != kmer:
                row_idx += 1
            while row_idx < len(rows) and rows[row_idx][0] == kmer:
                Tposition[-1][kmer_idx].append(int(rows[row_idx][1][1:]))
                row_idx += 1

        print(f"  genome {genome_num} processed in {time.time() - t0:.2f}s")

    return Tposition, ref_kmer_map


# ---------------------------------------------------------------------------
# Common-zone detection
# ---------------------------------------------------------------------------

def _is_extensible(pos: int, kmer_idx: int, Tposition: list, beta: int) -> tuple:
    """
    Determine whether a k-mer at *pos* in the reference can be extended (chained).

    An extension is valid when a consecutive k-mer exists at pos+1 in the reference
    and that same chain appears in at least *beta* genomes.

    Args:
        pos: Current position in the reference genome.
        kmer_idx: Index of the current k-mer in Tposition.
        Tposition: 3-D position table.
        beta: Minimum number of genomes required to validate an extension.

    Returns:
        (True, next_kmer_idx, positions_per_genome) if extensible, else (False, -1, -1).
    """
    for next_kmer in range(len(Tposition[0])):
        if pos + 1 not in Tposition[0][next_kmer]:
            continue

        chain_count = 1
        positions_per_genome = [[pos + 1]]
        seen_genomes = []

        for genome in range(1, len(Tposition)):
            positions_per_genome.append([])
            for p in Tposition[genome][kmer_idx]:
                if (p + 1) in Tposition[genome][next_kmer] and genome not in seen_genomes:
                    positions_per_genome[-1].append(p + 1)
                    chain_count += 1
                    seen_genomes.append(genome)

        if chain_count >= beta:
            return True, next_kmer, positions_per_genome
        return False, -1, -1

    return False, -1, -1


def _find_common_zones(Tposition: list, beta: int) -> tuple:
    """
    Identify all common zones across genomes by chaining extensible k-mers.

    Args:
        Tposition: 3-D position table.
        beta: Minimum genomes required to validate a k-mer chain.

    Returns:
        final_zones: Sorted list of zones (each zone = list of reference positions).
        zone_edges: 4-D list [genome][zone][edge] = [start, end] of each chain link.
    """
    final_zones = []
    visited = []
    zone_edges = [[] for _ in Tposition]
    raw_zones = []

    for kmer_idx in range(len(Tposition[0])):
        for pos in Tposition[0][kmer_idx]:
            if pos in visited:
                continue

            visited.append(pos)
            raw_zones.append([pos])
            for genome in range(len(Tposition)):
                zone_edges[genome].append([])

            extensible, next_kmer, positions = _is_extensible(pos, kmer_idx, Tposition, beta)
            while extensible:
                for genome_idx, gpos in enumerate(positions):
                    if gpos:
                        zone_edges[genome_idx][-1].append([gpos[0] - 1, gpos[0]])
                raw_zones[-1].append(positions[0][0])
                visited.append(positions[0][0])
                extensible, next_kmer, positions = _is_extensible(positions[0][0], next_kmer, Tposition, beta)

    # Deduplicate: keep the largest zone when overlap is detected
    if len(raw_zones) > 1:
        final_zones.append(raw_zones[0])
        for zone in raw_zones[1:]:
            i, j = 0, 0
            while i < len(final_zones) and zone[j] not in final_zones[i]:
                j += 1
                if j == len(zone):
                    i += 1
                    j = 0
            if i == len(final_zones):
                final_zones.append(zone)
            else:
                final_zones[i] = zone  # Larger zone supersedes the stored one
    else:
        final_zones = raw_zones

    final_zones.sort()
    return final_zones, zone_edges


def _deduplicate_zone_edges(zone_edges: list) -> list:
    """
    Remove duplicate zones from the zone_edges structure, mirroring _find_common_zones dedup logic.

    Args:
        zone_edges: 4-D list produced by _find_common_zones.

    Returns:
        Deduplicated 4-D list.
    """
    result = [[] for _ in zone_edges]
    if len(zone_edges[0]) <= 1:
        return zone_edges

    for genome in range(len(zone_edges)):
        result[genome].append(zone_edges[genome][0])

    for z_idx, zone in enumerate(zone_edges[0]):
        if not zone:
            continue
        i, j = 0, 0
        while i < len(result[0]) and zone[j] not in result[0][i]:
            j += 1
            if j == len(zone):
                i += 1
                j = 0
        if i == len(result[0]):
            for genome in range(len(result)):
                result[genome].append(zone_edges[genome][z_idx])
        else:
            for genome in range(len(result)):
                result[genome][i] = zone_edges[genome][z_idx]

    return result


def _sort_zone_edges(zone_edges: list) -> list:
    """
    Sort zone_edges by the starting position of each zone in the reference genome (index 0).

    Args:
        zone_edges: 4-D list produced by _find_common_zones.

    Returns:
        Sorted 4-D list.
    """
    sorted_edges = [[] for _ in zone_edges]
    start_positions = {}

    for z_idx, zone in enumerate(zone_edges[0]):
        if zone:
            start_positions[zone[0][0]] = z_idx

    for start in sorted(start_positions):
        z_idx = start_positions[start]
        for genome in range(len(zone_edges)):
            sorted_edges[genome].append(zone_edges[genome][z_idx])

    return sorted_edges


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def _positions_to_sequences(final_zones: list, ref_kmer_map: dict) -> list:
    """
    Convert reference genome positions to nucleotide sequences using ref_kmer_map.

    Args:
        final_zones: List of zones, each a list of reference positions.
        ref_kmer_map: Dict mapping position → k-mer sequence string.

    Returns:
        List of nucleotide sequence strings, one per zone.
    """
    sequences = []
    for zone in final_zones:
        seq = ref_kmer_map[zone[0]]
        for pos in zone[1:]:
            seq += ref_kmer_map[pos][-1]  # Append only the last nucleotide
        sequences.append(seq)
    return sequences


def _write_consensus_fasta(final_zones: list, sequences: list, k: int, fasta_files: list) -> None:
    """
    Write the consensus genome to ResulFasta.fasta.

    Each entry is annotated with its gkampi index range and scaffold membership.

    Args:
        final_zones: List of zones (reference positions).
        sequences: Corresponding nucleotide sequences.
        k: k-mer length (used to compute the true end position).
        fasta_files: List of FASTA paths; index 0 is the reference.
    """
    scaffold_ranges = []
    scaffold_names = []

    reference = SeqIO.to_dict(SeqIO.parse(fasta_files[0], "fasta"))
    cumulative = 0
    for seq_id, record in reference.items():
        scaffold_ranges.append([cumulative, cumulative + len(record)])
        scaffold_names.append(seq_id)
        cumulative += len(record) + 1

    with open("ResulFasta.fasta", "w") as fh:
        for zone_id, (zone, seq) in enumerate(zip(final_zones, sequences)):
            start, end = zone[0], zone[-1] + k - 1
            scaffold_note = ""
            for sc_idx, (sc_start, sc_end) in enumerate(scaffold_ranges):
                if sc_start <= start and end <= sc_end:
                    rel_start = start - sc_start
                    scaffold_note = (
                        f" | scaffold: {scaffold_names[sc_idx]}"
                        f" | relative start: {rel_start}"
                    )
            header = (
                f">zone_{zone_id + 1}"
                f" | ref_pos: {start}..{end}"
                f" | source: {fasta_files[0]}"
                f"{scaffold_note}\n"
            )
            fh.write(header + seq + "\n")


def _write_zones_csv(final_zones: list) -> None:
    """
    Write final_zones to resultat.csv for debugging and reproducibility.

    Args:
        final_zones: List of zones (reference positions).
    """
    with open("resultat.csv", "w", newline="") as fh:
        csv.writer(fh, delimiter=";").writerows(final_zones)


# ---------------------------------------------------------------------------
# Visualisation
# ---------------------------------------------------------------------------

def _rgba_to_hex(rgba: tuple) -> str:
    """Convert a matplotlib RGBA tuple to a hex colour string."""
    r, g, b, _ = rgba
    return f"#{int(r * 255):02x}{int(g * 255):02x}{int(b * 255):02x}"


def display_common_zones(
    zone_edges: list,
    results_frame,
    app,
    fasta_files: list,
    k: int,
    alpha: int,
    beta: int,
) -> None:
    """
    Render an interactive matplotlib figure of common zones inside the Tkinter UI.

    Each genome is drawn on a separate horizontal track. Zones are colour-coded by
    group. Clicking a segment displays its position range in a status label.

    Args:
        zone_edges: Sorted, deduplicated 4-D zone structure.
        results_frame: CTkScrollableFrame widget that hosts the plot.
        app: Root CTk application (used to bind the close protocol).
        fasta_files: List of FASTA paths (used for axis labels).
        k: k-mer length.
        alpha: Alpha threshold (for display only).
        beta: Beta threshold (for display only).
    """
    global _submit_count

    # -- Colour map ----------------------------------------------------------
    num_groups = max(len(genome) for genome in zone_edges)
    cmap = plt.cm.get_cmap("viridis", num_groups)

    all_positions = [pos for genome in zone_edges for zone in genome for edge in zone for pos in edge]
    min_val, max_val = min(all_positions), max(all_positions)
    max_zones = max(len(genome) for genome in zone_edges)

    # -- Build figure --------------------------------------------------------
    fig, ax = plt.subplots(figsize=(5, 2))
    legend_handles = []
    segment_info = {}
    y_ticks = []
    y_labels = []

    for genome_idx, genome_zones in enumerate(zone_edges):
        y_base = 3 * max_zones * genome_idx
        y_ticks.append(y_base + len(genome_zones) / 2)
        y_labels.append(f"Genome {genome_idx + 1}  ({fasta_files[genome_idx]})")

        for group_idx, group in enumerate(genome_zones):
            y_base += 1
            color = cmap(group_idx / num_groups)
            hex_color = _rgba_to_hex(color)

            if genome_idx == 0 and group:
                lo, hi = group[0][0], group[-1][1]
                label = f"{lo}–{hi + k - 1}"
                legend_handles.append((hex_color, label))
            else:
                label = legend_handles[group_idx][1] if group_idx < len(legend_handles) else ""

            for edge in group:
                line, = ax.plot(edge, [y_base, y_base], linewidth=10, color=color)
                segment_info[line] = label

    ax.set_xlim(min_val - 1, max_val + 1)
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels)
    fig.tight_layout()

    # -- Embed in Tkinter ----------------------------------------------------
    def _on_close():
        plt.close("all")
        for task in app.tk.eval("after info").split():
            app.after_cancel(task)
        app.destroy()

    def _on_segment_click(event):
        if event.xdata is None or event.ydata is None:
            return
        for line, label in segment_info.items():
            xd, yd = line.get_xdata(), line.get_ydata()
            if xd[0] <= event.xdata <= xd[1] and abs(event.ydata - yd[0]) <= 1:
                status_label.config(text=f"Selected zone: {label}")
                break

    def _delete_last_plot():
        global _submit_count
        results_frame.nametowidget(f"plot_frame_{_submit_count}").destroy()
        if _submit_count == 1:
            results_frame.nametowidget("btn_delete").destroy()
        _submit_count -= 1

    # Container frame for this run
    plot_frame = ttk.Frame(results_frame, padding="0.05i", name=f"plot_frame_{_submit_count}")

    # Parameter summary labels
    ref_label = f"Reference: {fasta_files[0]}  |  Compared: {', '.join(fasta_files[1:])}"
    param_label = f"k={k}  |  α={alpha}  |  β={beta}"
    ttk.Label(plot_frame, text=ref_label, background="white").pack(fill=tk.BOTH, expand=True)
    ttk.Label(plot_frame, text=param_label, background="white").pack(fill=tk.BOTH, expand=True)

    # Status bar
    status_label = tk.Label(
        plot_frame, text="Click a segment to inspect its position range.",
        bg="white", fg="black", font=("Arial", 11),
    )
    status_label.pack(side=tk.BOTTOM, fill=tk.X)

    # Scrollable canvas for the plot
    plot_canvas_widget = tk.Canvas(plot_frame, background="#1a1a1a")
    h_scroll = ttk.Scrollbar(plot_frame, orient=tk.HORIZONTAL, command=plot_canvas_widget.xview)
    v_scroll = ttk.Scrollbar(plot_frame, orient=tk.VERTICAL, command=plot_canvas_widget.yview)
    h_scroll.pack(side=tk.BOTTOM, fill=tk.X)
    v_scroll.pack(side=tk.RIGHT, fill=tk.Y)
    plot_frame.pack(side=tk.TOP, expand=True, fill=tk.BOTH)
    plot_canvas_widget.pack(side=tk.LEFT, expand=True, fill=tk.BOTH)
    plot_canvas_widget.configure(xscrollcommand=h_scroll.set, yscrollcommand=v_scroll.set)

    inner_frame = ttk.Frame(plot_canvas_widget)
    plot_canvas_widget.create_window((0, 0), window=inner_frame, anchor="center", height=300)
    inner_frame.pack(expand=True, fill=tk.BOTH)

    mpl_canvas = FigureCanvasTkAgg(fig, master=inner_frame)
    mpl_canvas.draw()
    mpl_canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

    inner_frame.update_idletasks()
    plot_canvas_widget.config(scrollregion=plot_canvas_widget.bbox("all"))

    NavigationToolbar2Tk(mpl_canvas, plot_frame)

    # Legend panel
    legend_frame = ttk.Frame(plot_frame)
    legend_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
    legend_scroll = ttk.Scrollbar(legend_frame, orient=tk.VERTICAL)
    legend_scroll.pack(side=tk.RIGHT, fill=tk.Y)
    legend_canvas = tk.Canvas(legend_frame)
    legend_inner = ttk.Frame(legend_canvas)
    legend_scroll.config(command=legend_canvas.yview)
    legend_canvas.create_window((0, 0), window=legend_inner, anchor="nw")
    legend_canvas.config(yscrollcommand=legend_scroll.set)
    legend_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

    for color, text in legend_handles:
        row = ttk.Frame(legend_inner)
        tk.Label(row, bg=color, width=2).pack(side=tk.LEFT, fill=tk.Y)
        tk.Label(row, text=text, anchor="w").pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        row.pack(fill=tk.X, pady=1)

    legend_inner.update_idletasks()
    legend_canvas.config(scrollregion=legend_canvas.bbox("all"))
    plot_canvas_widget.update_idletasks()
    plot_canvas_widget.config(scrollregion=plot_canvas_widget.bbox("all"))

    mpl_canvas.mpl_connect("button_press_event", _on_segment_click)
    app.protocol("WM_DELETE_WINDOW", _on_close)

    if _submit_count == 1:
        ttk.Button(
            results_frame, text="Remove last plot",
            command=_delete_last_plot, name="btn_delete",
        ).pack(pady=10)


# ---------------------------------------------------------------------------
# Pipeline entry point
# ---------------------------------------------------------------------------

def run_pipeline(
    fasta_files: list,
    kmer_size: str,
    alpha_pct: str,
    beta_pct: str,
    results_frame,
    app,
) -> None:
    """
    Execute the full gCons pipeline for a given set of genomes and parameters.

    Steps:
      1. Run redoak to find shared k-mers.
      2. Filter by alpha threshold (grep).
      3. Run gkampi per genome to get positional data.
      4. Detect common zones via k-mer chaining.
      5. Write FASTA and CSV outputs.
      6. Render the interactive visualisation.

    Args:
        fasta_files: List of FASTA file paths; index 0 is the reference genome.
        kmer_size: k-mer length as a string (converted internally to int).
        alpha_pct: Percentage of genomes that must share a k-mer (0–100).
        beta_pct: Percentage of genomes that must share a k-mer chain link (0–100).
        results_frame: CTkScrollableFrame widget that will host the result plot.
        app: Root CTk application instance.
    """
    global _submit_count
    _submit_count += 1

    k = int(kmer_size)
    n = len(fasta_files)
    alpha = int(np.round(int(alpha_pct) * n / 100))
    beta = int(np.round(int(beta_pct) * n / 100))

    # Step 1 & 2 — shared k-mer discovery
    _run_redoak(k, fasta_files)
    _grep_common_kmers(alpha, fasta_files)

    # Step 3 — positional indexing
    Tposition, ref_kmer_map = _build_position_table(k, fasta_files)

    # Step 4 — zone detection
    final_zones, zone_edges = _find_common_zones(Tposition, beta)
    sorted_edges = _sort_zone_edges(_deduplicate_zone_edges(zone_edges))
    _debug_print_all_zones(_sort_zone_edges(_deduplicate_zone_edges(zone_edges)))

    # Step 5 — outputs
    sequences = _positions_to_sequences(final_zones, ref_kmer_map)
    _write_consensus_fasta(final_zones, sequences, k, fasta_files)
    _write_zones_csv(final_zones)

    # Step 6 — visualisation
    display_common_zones(sorted_edges, results_frame, app, fasta_files, k, alpha, beta)