# streamlit_cross_library_comparisons

Interactive Streamlit app for comparing genetic interaction calls (dLFC, HSA, lfc_combo) across three Genmap library versions (MinLib, Library 1, Library 2).

---

## Data provenance

| Library | Input path | Used for |
|---------|-----------|----------|
| MinLib | `/lustre/scratch126/casm/teams/team215/projects/Genmap/library2_analysis/n_guides_analysis/minlib_guides_only/results/GI_calls.txt` | MinLib-guides-only GI calls |
| Library 1 | `/lustre/scratch126/casm/teams/team215/projects/Genmap/library1_results/GI_calls.txt` | Library 1 GI calls |
| Library 2 | `/lustre/scratch126/casm/teams/team215/projects/Genmap/library2_analysis/results/GI_calls.txt` | Full Library 2 GI calls |

Columns used from each file: `gene_combo`, `lfc_combo`, `dLFC`, `HSA`.

---

## Analyses

### 1. Pairwise library comparison app (`scatterplots_and_stats.py`)

- **Directory:** project root
- **Description:** Streamlit app that compares any two of the three libraries on a selected metric. Displays a scatterplot with identity line and annotated Spearman/Pearson r values, plus a correlation statistics table with Benjamini-Hochberg-adjusted p-values.
- **Inputs:** `GI_calls.txt` files listed in the data provenance table above
- **Outputs:** interactive figures and statistics (browser only; no files written)
- **NFS archive path:** N/A (no output files)
- **Environment:** `library_comparisons_env` (see Environments)

---

## Directory structure

```
streamlit_cross_library_comparisons/
├── README.md
├── environments.md              # environment documentation
├── environments/                # dependency exports
│   └── library_comparisons_env.txt   # pip freeze output
├── .gitignore
├── scatterplots_and_stats.py    # Streamlit app
└── library_comparisons_env/    # pip venv (gitignored)
```

---

## Environments

| Environment | Purpose |
|-------------|---------|
| `library_comparisons_env` | Run Streamlit app |

See [environments.md](environments.md) for full details.

**Quick start:**
```bash
source library_comparisons_env/bin/activate
streamlit run scatterplots_and_stats.py
```

---

## Known issues / future work

- Data paths are hardcoded in `scatterplots_and_stats.py`; if input files move, update `LIBRARY_PATHS` at the top of the script.
- BH correction is applied separately within each correlation type (Spearman p-values corrected together across metrics; Pearson p-values corrected together across metrics).
