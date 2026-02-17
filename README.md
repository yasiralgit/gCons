# gCons

A desktop tool for genomic consensus analysis across multiple FASTA genomes.

gCons identifies **common zones** — genomic regions shared by a configurable proportion of input genomes — using k-mer based comparison. Results are exported as a consensus FASTA file and visualised interactively.

---

## How it works

gCons chains two external tools:

1. **redoak** finds k-mers shared across genomes and reports how many genomes contain each one.
2. **gkampi** maps each common k-mer to its exact position within every genome.

From this positional data, gCons builds chains of consecutive k-mers that appear consistently across genomes (controlled by the α and β thresholds), which it calls *common zones*. Each zone is exported as a nucleotide sequence anchored to the reference genome.

```
FASTA files
    │
    ▼
redoak ──► shared k-mer list
    │
    ▼
gkampi ──► positional index (per genome)
    │
    ▼
Zone detection (α, β filtering)
    │
    ├──► ResulFasta.fasta   (consensus sequences)
    ├──► resultat.csv       (zone positions)
    └──► Interactive plot   (Tkinter/matplotlib)
```

---

## Requirements

- Python 3.10+
- `redoak` and `gkampi` binaries in the working directory
- tkinter (see note below)

```bash
pip install -r requirements.txt
```

> **macOS (Homebrew):** tkinter requires a separate install.
> ```bash
> brew install python-tk@3.x   # replace 3.x with your Python version
> ```
> **Ubuntu/Debian:**
> ```bash
> sudo apt install python3-tk
> ```

---

## Usage

```bash
python gCons_controleurTK.py
```

Fill in the four fields and click **Run analysis**:

| Field | Description |
|---|---|
| FASTA file paths | Space-separated paths. The **first file is the reference genome**. |
| k-mer size | Length of the k-mers used for comparison (e.g. `21`). |
| Alpha (%) | Minimum percentage of genomes that must share a k-mer for it to be considered common. |
| Beta (%) | Minimum percentage of genomes that must share a consecutive k-mer link for a zone to be extended. |

**Example:**

```
FASTA files  →  genome1.fasta genome2.fasta genome3.fasta
k-mer size   →  21
Alpha        →  80
Beta         →  60
```

This finds k-mers present in at least 80% of genomes, then chains them into zones where at least 60% of genomes share consecutive links.

---

## Outputs

| File | Description |
|---|---|
| `ResulFasta.fasta` | Consensus sequences, one entry per common zone, annotated with reference coordinates and scaffold membership. |
| `resultat.csv` | Raw zone positions in the reference genome (semicolon-separated). |
| `commun.txt` | Intermediate redoak output filtered by alpha. |
| `result<n>.csv` | Raw gkampi positional output for genome *n*. |

---

## Project structure

```
gCons/
├── gCons_controleurTK.py   # UI entry point
├── gCons_fonctionsTK.py    # Pipeline logic and visualisation
├── requirements.txt
├── redoak                  # External binary (not versioned)
├── gkampi                  # External binary (not versioned)
└── data/
    └── example/            # Sample FASTA files for testing
```

---

## Parameters — choosing α and β

Both thresholds are expressed as a percentage of the total number of input genomes and are rounded to the nearest integer.

- **High α** (e.g. 90%): only k-mers nearly universal across all genomes are retained. Yields fewer, more conserved zones.
- **Low α** (e.g. 50%): includes k-mers present in just half the genomes. Yields more zones, potentially less conserved.
- **β ≤ α** is the natural constraint: a chain link cannot be more prevalent than the k-mer itself.

---

## Internship context

This tool was developed during a research internship at [LIRMM](https://www.lirmm.fr/) (Laboratoire d'Informatique, de Robotique et de Microélectronique de Montpellier) as part of a genomic analysis pipeline for inter-genome comparison using k-mer methods (redoak, gkampi).