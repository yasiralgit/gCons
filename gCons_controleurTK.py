"""
gCons_controleurTK.py
---------------------
Entry point and UI controller for the gCons genomic consensus tool.

Builds a CustomTkinter interface with two tabs:
  - Parameters: user inputs (FASTA files, k-mer size, alpha/beta thresholds)
  - Results:    interactive matplotlib visualisation of detected common zones

Usage:
    python gCons_controleurTK.py

Dependencies:
    customtkinter, gCons_fonctionsTK (same directory)
    redoak and gkampi binaries must be present in the working directory.

Author: AL-YOUSSFI Yasir
"""

from customtkinter import *
from gCons_fonctionsTK import run_pipeline


# ---------------------------------------------------------------------------
# Widgets — module-level references populated by build_ui()
# ---------------------------------------------------------------------------

root: CTk = None
results_frame: CTkScrollableFrame = None
entry_files: CTkEntry = None
entry_kmer: CTkEntry = None
entry_alpha: CTkEntry = None
entry_beta: CTkEntry = None


# ---------------------------------------------------------------------------
# Event handler
# ---------------------------------------------------------------------------

def _on_submit() -> None:
    """
    Validate inputs and launch the gCons pipeline when the user clicks Submit.

    Reads the four input fields and delegates execution to run_pipeline().
    All type conversion and threshold scaling is handled inside run_pipeline().
    """
    fasta_files = entry_files.get().split()
    kmer_size = entry_kmer.get()
    alpha_pct = entry_alpha.get()
    beta_pct = entry_beta.get()

    if not fasta_files or not kmer_size or not alpha_pct or not beta_pct:
        print("[gCons] Error: all fields must be filled before submitting.")
        return

    print(f"[gCons] Starting pipeline — files: {fasta_files}  k={kmer_size}  α={alpha_pct}%  β={beta_pct}%")
    run_pipeline(fasta_files, kmer_size, alpha_pct, beta_pct, results_frame, root)


# ---------------------------------------------------------------------------
# UI builder
# ---------------------------------------------------------------------------

def build_ui() -> None:
    """
    Construct the full application window and initialise all widget references.

    Layout:
        ┌────────────────────────────────────┐
        │  Parameters tab  │  Results tab    │
        │  ─────────────── │  ───────────── │
        │  [FASTA files  ] │  (plots appear  │
        │  [k-mer size   ] │   here after    │
        │  [Alpha %      ] │   each submit)  │
        │  [Beta  %      ] │                 │
        │  [ Submit ]      │                 │
        └────────────────────────────────────┘
    """
    global root, results_frame, entry_files, entry_kmer, entry_alpha, entry_beta

    root = CTk()
    root.geometry("560x520")
    root.title("gCons — Genomic Consensus Tool")
    set_default_color_theme("green")

    # Main tab container
    tabview = CTkTabview(master=root)
    tabview.pack(expand=True, anchor="n", padx=20, pady=20, fill="both")
    tabview.add("Parameters")
    tabview.add("Results")

    # ------------------------------------------------------------------ #
    #  Parameters tab                                                      #
    # ------------------------------------------------------------------ #
    params_scroll = CTkScrollableFrame(
        master=tabview.tab("Parameters"), border_width=2, orientation="vertical"
    )
    params_scroll.pack(expand=True, anchor="n", pady=10, fill="both")

    def _labeled_entry(parent, label_text: str, width: int = 300) -> CTkEntry:
        """Helper — add a label + entry pair inside *parent*."""
        frame = CTkFrame(master=parent, border_width=2)
        frame.pack(expand=True, anchor="n", pady=15, padx=10)
        CTkLabel(master=frame, text=label_text).pack(anchor="w", padx=10, pady=8)
        entry = CTkEntry(master=frame, width=width, text_color="#FFCC70")
        entry.pack(anchor="center", pady=8, padx=20)
        return entry

    entry_files = _labeled_entry(
        params_scroll,
        "FASTA file paths — space-separated, first file is the reference:",
    )
    entry_kmer = _labeled_entry(params_scroll, "k-mer size:", width=80)
    entry_alpha = _labeled_entry(
        params_scroll,
        "Alpha — % of genomes that must share a k-mer (0–100):",
        width=80,
    )
    entry_beta = _labeled_entry(
        params_scroll,
        "Beta — % of genomes that must share a k-mer chain link (0–100):",
        width=80,
    )

    CTkButton(
        master=params_scroll, text="Run analysis", command=_on_submit
    ).pack(anchor="s", padx=30, pady=25)

    # ------------------------------------------------------------------ #
    #  Results tab                                                         #
    # ------------------------------------------------------------------ #
    results_frame = CTkScrollableFrame(
        master=tabview.tab("Results"), border_width=2, orientation="vertical"
    )
    results_frame.pack(expand=True, anchor="n", pady=10, fill="both")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    build_ui()
    root.mainloop()