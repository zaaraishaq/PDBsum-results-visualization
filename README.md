# PDBsum-results-visualization

A dependency-light, offline replacement for the parts of the [PDBsum](http://www.ebi.ac.uk/thornton-srv/databases/pdbsum/)
web server that generate ligand-plots and interface diagnostics for a docked or
crystallographic protein complex -- built because the PDBsum web interface has been
unreliable / down for extended stretches.

Given any two-chain PDB coordinate file (a peptide bound to a target, an antibody-antigen
complex, a docked antimicrobial peptide against a receptor, etc.), it computes:

- **Hydrogen bonds** (polar-polar atom pairs, 2.4-3.5 A)
- **Salt bridges** (charged side-chain atom pairs, <= 4.0 A)
- **Hydrophobic / non-bonded contacts** (<= 4.0 A)
- **Per-residue interface catalogs** for both chains
- **Secondary structure** (from HELIX/SHEET header records, when present)
- **Ramachandran (phi/psi) backbone geometry**, with a simplified favored/outlier summary

...and renders everything as a text report plus six figures (five static PNGs and one
standalone interactive 3D HTML viewer).

No Biopython required -- just `numpy` and `matplotlib`.

## Quick start

```bash
pip install -r requirements.txt

# Single complex, chains given explicitly (recommended when known)
python src/pdb_interface_analyzer.py path/to/complex.pdb \
    --ligand-chain B --receptor-chain A \
    --ligand-label p53_peptide --receptor-label MDM2

# Single complex, auto-detect which chain is the (smaller) ligand vs. (larger) receptor
python src/pdb_interface_analyzer.py path/to/complex.pdb

# Batch mode: every .pdb file in a folder
python src/pdb_interface_analyzer.py path/to/folder/ --batch
```

Results land in `./interface_results/<complex_name>/` by default (override with
`--output`), each folder containing:

| File | Contents |
|---|---|
| `Complex_Interface_Report.txt` | Full text report: stats, contact lists, residue catalog |
| `interface_summary_bubble.png` | Proportional bubble summary (residue counts, bond counts) |
| `interface_interactions_network.png` | Residue-to-residue Bezier contact network |
| `interface_contact_heatmap.png` | Contact-density heatmap between interface residues |
| `docked_complex_3d_static.png` | Static 3D backbone trace, interface highlighted |
| `docked_complex_3d.html` | Standalone interactive 3Dmol.js viewer (open in any browser) |
| `secondary_structure.png` | Per-chain helix/sheet/coil track along the sequence |
| `ramachandran_plot.png` | Phi/psi scatter plot per chain |
| `Ramachandran_Summary.txt` | Favored-region / outlier residue counts per chain |

In batch mode you additionally get a `Summary_All_Complexes.txt` comparing every
complex processed in that run.

## Google Colab

Prefer not to install anything locally? Open
[`notebooks/PDBsum_Interface_Analyzer_Colab.ipynb`](notebooks/PDBsum_Interface_Analyzer_Colab.ipynb)
in Colab. It covers both single-complex and batch analysis, displays every figure
inline, and lets you download the results as a zip. It can fetch a structure directly
from RCSB by PDB ID, or you can upload your own `.pdb` file(s).

## Chain assignment

If `--ligand-chain`/`--receptor-chain` (single-file mode) or a `--chain-map` JSON file
(batch mode) aren't given, the script auto-detects: whichever protein chain has fewer
residues is treated as the ligand (e.g. a short peptide), and the largest remaining
chain as the receptor. For anything with more than two chains, or where "biggest/smallest"
isn't the right split, specify chains explicitly.

## Important caveats

- **Secondary structure** is read straight from the PDB file's `HELIX`/`SHEET` header
  records. Raw docking output frequently has no header at all -- in that case the plot
  will show the whole chain as coil and say so explicitly. This is expected, not a bug.
- **Ramachandran "favored regions"** shown here are a simplified, illustrative guide
  (rough alpha/beta/left-handed-helix boxes) meant for an at-a-glance sanity check --
  **not** a validated statistical potential. For publication-grade validation, cross-check
  with [MolProbity](http://molprobity.biochem.duke.edu/) or the structure's wwPDB
  validation report.
- Distance thresholds (2.4-3.5 A for H-bonds, 4.0 A for salt bridges/hydrophobic contacts)
  are geometry-only heuristics, the same kind PDBsum/LIGPLOT use -- they don't account for
  hydrogen positions or electrostatics, so treat borderline contacts as candidates worth a
  closer look rather than certainties.

## Repository layout

```
src/pdb_interface_analyzer.py   - the analyzer (single-file + batch CLI)
notebooks/                      - Google Colab notebook
examples/data/                  - example input PDB (RCSB 1YCR: MDM2 / p53 peptide)
examples/output/                - example output from running the script on it
requirements.txt
LICENSE
```

## Example

[`examples/`](examples/) contains a full worked example using PDB entry
[1YCR](https://www.rcsb.org/structure/1YCR) (MDM2 bound to a 15-residue p53
transactivation-domain peptide) -- chosen because, like an AMP-target docking result,
it's a short peptide (chain B) bound to a larger target domain (chain A). See
[`examples/README.md`](examples/README.md) for how it was generated and how to read
the output.

## License

MIT -- see [LICENSE](LICENSE). Free to use, modify, and redistribute for research or
any other purpose.

## Citing PDBsum / LIGPLOT concepts

This tool is an independent reimplementation inspired by the *type* of diagram PDBsum
and LIGPLOT popularized; it is not affiliated with, and does not reuse code from,
either project. If your work also benefits from the original concepts, consider citing:

> Laskowski, R.A. (2001) PDBsum: summaries and analyses of PDB structures.
> *Nucleic Acids Research*, 29, 221-222.
