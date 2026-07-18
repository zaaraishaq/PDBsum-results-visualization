# Worked example: MDM2 / p53 peptide (PDB 1YCR)

## Why this structure

[1YCR](https://www.rcsb.org/structure/1YCR) is MDM2's N-terminal domain (chain A,
109 residues) bound to a 15-residue transactivation-domain peptide from p53 (chain B).
It was picked as the example because that shape -- a short peptide docked onto a larger
target protein -- mirrors an AMP-vs-target-protein docking result, so the output here is
a reasonable preview of what you'll get running this tool on your own docked complexes.

Reference: Kussie et al. (1996) *Science* 274:948-953,
[doi:10.1126/science.274.5289.948](https://doi.org/10.1126/science.274.5289.948).
Coordinates downloaded directly from RCSB
([files.rcsb.org/download/1YCR.pdb](https://files.rcsb.org/download/1YCR.pdb));
PDB coordinate data is not copyrightable and is freely redistributable.

## Files

- `data/1YCR_MDM2_p53.pdb` -- the input structure (header trimmed of refinement-statistics
  boilerplate to keep the file small; all `ATOM`/`HELIX`/`SHEET`/`SEQRES` records are intact).
- `output/1YCR_MDM2_p53/` -- everything the script produces from that input.

## How it was generated

```bash
python src/pdb_interface_analyzer.py examples/data/1YCR_MDM2_p53.pdb \
    --ligand-chain B --receptor-chain A \
    --ligand-label p53_peptide --receptor-label MDM2 \
    --output examples/output
```

## Reading the results

- **`Complex_Interface_Report.txt`** lists 11 hydrogen bonds, 1 salt bridge, and 11
  hydrophobic contacts. The hydrophobic-contact list correctly recovers the interface's
  three headline residues from the original paper -- **Phe19, Trp23, and Leu26** -- as
  the peptide's main contacts burying into MDM2's hydrophobic cleft.
- **`interface_summary_bubble.png`** / **`interface_contact_heatmap.png`** /
  **`interface_interactions_network.png`** visualize the same contacts three different ways.
- **`secondary_structure.png`** shows the p53 peptide folding into a single alpha helix upon
  binding (matching `HELIX 5` in the source file) while MDM2 shows its four native helices.
- **`ramachandran_plot.png`** + **`Ramachandran_Summary.txt`** show backbone phi/psi angles;
  ~86-91% of residues fall in the illustrative favored regions here, with the outliers
  concentrated at loop/turn regions -- as expected for this resolution (2.6 A) structure.
- **`docked_complex_3d.html`** -- open this directly in a browser for an interactive,
  rotatable 3D view with toggleable chains, bonds, and render styles.

Your own docking output will very likely differ in one respect: real docking results
usually have **no HELIX/SHEET header records at all** (that information comes from
crystallographic refinement, which a docking run doesn't do), so `secondary_structure.png`
will show your chains as plain coil and say so on the plot. That's expected, not an error.
