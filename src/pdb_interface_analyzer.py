#!/usr/bin/env python3
"""
pdb_interface_analyzer.py
==========================
A dependency-light, PDBsum/LIGPLOT-style interface analyzer for protein-protein
and peptide-protein docked/crystal complexes in the classic PDB coordinate
format.

WHY THIS EXISTS
---------------
PDBsum's web interface (which used to generate exactly this kind of
ligand-plot / interface diagram) has been unreliable / down for extended
periods. This script reproduces the core analysis -- hydrogen bonds, salt
bridges, hydrophobic contacts, and per-residue interface catalogs -- entirely
offline, plus a set of static and interactive visualizations, so researchers
aren't blocked on that service being available.

FEATURES
--------
- Zero third-party PDB parsing dependency: ships its own small, fixed-width
  PDB ATOM/HETATM reader (no Biopython required).
- Works on a single complex or a whole folder of them (batch mode).
- Auto-detects which chain is the "ligand" (smaller, e.g. a peptide) and
  which is the "receptor" (larger, e.g. the target protein) when not told
  explicitly -- or you can pin exact chain IDs per file.
- Outputs, per complex:
    1. Complex_Interface_Report.txt      - full text report
    2. interface_summary_bubble.png      - proportional bubble summary
    3. interface_interactions_network.png- residue-to-residue Bezier network
    4. interface_contact_heatmap.png     - contact-density heatmap
    5. docked_complex_3d_static.png      - static 3D backbone trace
    6. docked_complex_3d.html            - standalone interactive 3Dmol.js viewer
    7. secondary_structure.png           - per-chain helix/sheet/coil track
    8. ramachandran_plot.png             - phi/psi backbone dihedral plot
    9. Ramachandran_Summary.txt          - favored/outlier residue counts
  plus a combined Summary_All_Complexes.txt when run in batch mode.

Notes on #7/#8/#9: HELIX/SHEET assignment is read directly from the PDB
header. Raw docking output often has no header at all, in which case the
secondary-structure figure will show everything as coil and say so on the
plot -- that's expected, not a bug. The Ramachandran "favored region" shading
is a simplified illustrative guide, not a validated statistical potential;
treat it as a quick sanity check, not a publication-grade validation (use
MolProbity or the wwPDB validation report for that).

INSTALL
-------
    pip install numpy matplotlib

USAGE
-----
Single complex, auto-detecting which chain is the ligand vs. receptor:
    python pdb_interface_analyzer.py path/to/complex.pdb

Single complex, pinning chains explicitly (recommended when you know them):
    python pdb_interface_analyzer.py path/to/complex.pdb \\
        --ligand-chain A --receptor-chain W \\
        --ligand-label Peptide --receptor-label Target

Batch mode -- every .pdb file in a folder, auto-detecting chains per file:
    python pdb_interface_analyzer.py path/to/folder/ --batch

Batch mode with a per-file chain map (JSON: {"filename.pdb": {"ligand": "A",
"receptor": "W"}, ...}):
    python pdb_interface_analyzer.py path/to/folder/ --batch --chain-map chains.json

Results are written to --output (default: ./interface_results/), one
subfolder per complex.
"""

import argparse
import collections
import glob
import json
import math
import os
import sys

try:
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.path as mpath
    import matplotlib.patches as mpatches
    from matplotlib.patches import Ellipse
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3D projection)
except ImportError:
    print("[ERROR] Missing required libraries. Please install them:")
    print("    pip install numpy matplotlib")
    sys.exit(1)

# ------------------------------------------------------------------------
# Standard PDBsum/LIGPLOT-style residue color coding
# ------------------------------------------------------------------------
RESIDUE_COLORS = {
    'ALA': '#8e8e8e', 'VAL': '#5e5e5e', 'LEU': '#5e5e5e', 'ILE': '#5e5e5e', 'MET': '#8e8e8e',  # Grey / aliphatic
    'PHE': '#c58bfc', 'TYR': '#c58bfc', 'TRP': '#c58bfc',                                       # Purple / aromatic
    'SER': '#00e676', 'THR': '#00e676', 'ASN': '#00e676', 'GLN': '#00e676', 'CYS': '#00e676',   # Green / polar
    'ARG': '#29b6f6', 'LYS': '#29b6f6', 'HIS': '#29b6f6',                                       # Blue / positive
    'ASP': '#ff1744', 'GLU': '#ff1744',                                                         # Red / negative
    'PRO': '#ff9100', 'GLY': '#ff9100',                                                         # Orange / special
}

HYDROPHOBIC_RESIDUES = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TYR', 'TRP', 'PRO'}
POSITIVE_RESIDUES = {'ARG', 'LYS', 'HIS'}
NEGATIVE_RESIDUES = {'ASP', 'GLU'}
POSITIVE_ATOMS = {'NZ', 'NH1', 'NH2', 'ND1', 'NE2'}
NEGATIVE_ATOMS = {'OD1', 'OD2', 'OE1', 'OE2'}

HBOND_MIN_DIST = 2.4
HBOND_MAX_DIST = 3.5
CONTACT_MAX_DIST = 4.0


def format_res_label(res_name, res_id):
    """3-letter code + residue number -> e.g. 'Trp23'."""
    camel = res_name[0].upper() + res_name[1:].lower()
    return f"{camel}{res_id}"


def calculate_distance(coord1, coord2):
    return math.sqrt(sum((c1 - c2) ** 2 for c1, c2 in zip(coord1, coord2)))


# ------------------------------------------------------------------------
# Minimal, dependency-free PDB reader
# ------------------------------------------------------------------------
class SimpleAtom:
    __slots__ = ("name", "coord")

    def __init__(self, name, coord):
        self.name = name
        self.coord = coord


class SimpleResidue:
    __slots__ = ("resname", "resid", "is_hetatm", "atoms")

    def __init__(self, resname, resid, is_hetatm):
        self.resname = resname
        self.resid = resid
        self.is_hetatm = is_hetatm
        self.atoms = {}

    def get_resname(self):
        return self.resname

    def __contains__(self, atom_name):
        return atom_name in self.atoms

    def __getitem__(self, atom_name):
        return self.atoms[atom_name]


def parse_pdb_chains(pdb_path):
    """Parse ATOM/HETATM records into {chain_id: [SimpleResidue, ...]} using
    standard fixed-width PDB columns. First-encountered altLoc wins."""
    chains = collections.OrderedDict()

    with open(pdb_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            rectype = line[0:6].strip()
            if rectype not in ("ATOM", "HETATM"):
                continue
            if len(line) < 54:
                continue

            is_hetatm = (rectype == "HETATM")
            atom_name = line[12:16].strip()
            resname = line[17:20].strip()
            chain_id = line[21].strip() or " "
            resid_str = line[22:26].strip()
            try:
                resid = int(resid_str)
            except ValueError:
                continue

            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue

            if chain_id not in chains:
                chains[chain_id] = collections.OrderedDict()

            res_key = (resid, is_hetatm)
            if res_key not in chains[chain_id]:
                chains[chain_id][res_key] = SimpleResidue(resname, resid, is_hetatm)

            residue = chains[chain_id][res_key]
            if atom_name not in residue.atoms:
                residue.atoms[atom_name] = SimpleAtom(atom_name, (x, y, z))

    flat_chains = collections.OrderedDict()
    for chain_id, res_dict in chains.items():
        flat_chains[chain_id] = list(res_dict.values())
    return flat_chains


def guess_ligand_receptor_chains(chains):
    """When chains aren't specified explicitly: the chain with fewer ATOM
    (non-HETATM) residues is treated as the ligand (e.g. a short peptide);
    the largest remaining chain is the receptor. Only the two chains with
    the most protein residues are considered."""
    sizes = []
    for chain_id, residues in chains.items():
        n_protein_res = sum(1 for r in residues if not r.is_hetatm)
        if n_protein_res > 0:
            sizes.append((chain_id, n_protein_res))
    sizes.sort(key=lambda x: x[1])
    if len(sizes) < 2:
        raise ValueError(
            "Could not auto-detect ligand/receptor chains: need at least 2 "
            "protein chains. Found: %s. Specify --ligand-chain/--receptor-chain "
            "explicitly." % sizes
        )
    ligand_chain = sizes[0][0]
    receptor_chain = sizes[-1][0]
    return [ligand_chain], [receptor_chain]


# ------------------------------------------------------------------------
# Core interface analyzer
# ------------------------------------------------------------------------
class DockedInterfaceAnalyzer:
    def __init__(self, pdb_path, receptor_chains, ligand_chains,
                 ligand_label="Ligand", receptor_label="Receptor", complex_label=None):
        self.pdb_path = pdb_path
        self.filename = os.path.basename(pdb_path)
        self.receptor_chains = receptor_chains
        self.ligand_chains = ligand_chains
        self.ligand_label = ligand_label
        self.receptor_label = receptor_label
        self.complex_label = complex_label or os.path.splitext(self.filename)[0]

        self.h_bonds = []
        self.salt_bridges = []
        self.hydrophobic_contacts = []
        self.r_interface_residues = set()
        self.l_interface_residues = set()
        self.ca_coordinates = collections.defaultdict(list)
        self.interface_ca_coords = []
        self.raw_pdb_content = ""
        self.chains = None  # populated by load_and_parse(); {chain_id: [SimpleResidue, ...]}

    def load_and_parse(self):
        if not os.path.exists(self.pdb_path):
            raise FileNotFoundError(f"PDB file not found: {self.pdb_path}")

        with open(self.pdb_path, 'r', encoding='utf-8', errors='ignore') as f:
            self.raw_pdb_content = f.read()

        chains = parse_pdb_chains(self.pdb_path)
        self.chains = chains

        if not self.receptor_chains or not self.ligand_chains:
            self.ligand_chains, self.receptor_chains = guess_ligand_receptor_chains(chains)
            print(f"   Auto-detected ligand chain(s): {self.ligand_chains}  "
                  f"receptor chain(s): {self.receptor_chains}")

        missing = [c for c in (self.receptor_chains + self.ligand_chains) if c not in chains]
        if missing:
            raise ValueError(
                f"Chain(s) {missing} not found in {self.filename}. "
                f"Chains present: {list(chains.keys())}"
            )

        self._calculate_interactions(chains)
        return True

    def _calculate_interactions(self, chains):
        r_atoms, l_atoms = [], []

        for chain_id, residues in chains.items():
            if chain_id not in self.receptor_chains and chain_id not in self.ligand_chains:
                continue
            for residue in residues:
                if residue.is_hetatm:
                    continue
                if 'CA' in residue:
                    self.ca_coordinates[chain_id].append((residue.resid, residue['CA'].coord))
                for atom_name, atom in residue.atoms.items():
                    if chain_id in self.receptor_chains:
                        r_atoms.append((chain_id, residue, atom))
                    elif chain_id in self.ligand_chains:
                        l_atoms.append((chain_id, residue, atom))

        for r_chain, r_res, r_atom in r_atoms:
            r_res_name, r_res_id, r_atom_name = r_res.get_resname(), r_res.resid, r_atom.name

            for l_chain, l_res, l_atom in l_atoms:
                l_res_name, l_res_id, l_atom_name = l_res.get_resname(), l_res.resid, l_atom.name
                dist = calculate_distance(r_atom.coord, l_atom.coord)

                if HBOND_MIN_DIST <= dist <= HBOND_MAX_DIST:
                    r_polar = any(ch in r_atom_name for ch in ('N', 'O', 'S'))
                    l_polar = any(ch in l_atom_name for ch in ('N', 'O', 'S'))
                    if r_polar and l_polar:
                        self.h_bonds.append({
                            'r_res': f"{r_res_name} {r_res_id}", 'r_chain': r_chain, 'r_atom': r_atom_name,
                            'l_res': f"{l_res_name} {l_res_id}", 'l_chain': l_chain, 'l_atom': l_atom_name,
                            'dist': dist,
                        })
                        self.r_interface_residues.add((r_res_name, r_res_id, r_chain))
                        self.l_interface_residues.add((l_res_name, l_res_id, l_chain))
                        if 'CA' in r_res and 'CA' in l_res:
                            self.interface_ca_coords.append(r_res['CA'].coord)
                            self.interface_ca_coords.append(l_res['CA'].coord)

                if dist <= CONTACT_MAX_DIST:
                    r_pos = r_res_name in POSITIVE_RESIDUES and r_atom_name in POSITIVE_ATOMS
                    l_neg = l_res_name in NEGATIVE_RESIDUES and l_atom_name in NEGATIVE_ATOMS
                    r_neg = r_res_name in NEGATIVE_RESIDUES and r_atom_name in NEGATIVE_ATOMS
                    l_pos = l_res_name in POSITIVE_RESIDUES and l_atom_name in POSITIVE_ATOMS
                    if (r_pos and l_neg) or (r_neg and l_pos):
                        self.salt_bridges.append({
                            'r_res': f"{r_res_name} {r_res_id}", 'r_chain': r_chain, 'r_atom': r_atom_name,
                            'l_res': f"{l_res_name} {l_res_id}", 'l_chain': l_chain, 'l_atom': l_atom_name,
                            'dist': dist,
                        })
                        self.r_interface_residues.add((r_res_name, r_res_id, r_chain))
                        self.l_interface_residues.add((l_res_name, l_res_id, l_chain))

                if dist <= CONTACT_MAX_DIST:
                    if (r_res_name in HYDROPHOBIC_RESIDUES and l_res_name in HYDROPHOBIC_RESIDUES
                            and r_atom_name.startswith('C') and l_atom_name.startswith('C')):
                        self.hydrophobic_contacts.append({
                            'r_res': f"{r_res_name} {r_res_id}", 'r_chain': r_chain,
                            'l_res': f"{l_res_name} {l_res_id}", 'l_chain': l_chain,
                            'dist': dist,
                        })
                        self.r_interface_residues.add((r_res_name, r_res_id, r_chain))
                        self.l_interface_residues.add((l_res_name, l_res_id, l_chain))

    def unique_hydrophobic_contacts(self):
        seen, unique = set(), []
        for c in self.hydrophobic_contacts:
            key = (c['r_res'], c['r_chain'], c['l_res'], c['l_chain'])
            if key not in seen:
                seen.add(key)
                unique.append(c)
        return unique

    def write_report(self, report_path):
        unique_hydro = self.unique_hydrophobic_contacts()
        with open(report_path, "w", encoding='utf-8') as f:
            f.write("=" * 72 + "\n")
            f.write("      Docked Complex Interface Diagnostics Report\n")
            f.write("=" * 72 + "\n\n")
            f.write(f"Analyzed Complex   : {self.filename}\n")
            f.write(f"Receptor Chain(s)  : {', '.join(self.receptor_chains)} ({self.receptor_label})\n")
            f.write(f"Ligand Chain(s)    : {', '.join(self.ligand_chains)} ({self.ligand_label})\n\n")

            f.write("1. STRUCTURAL INTERFACE GLOBAL STATISTICS\n")
            f.write("-" * 72 + "\n")
            f.write(f"  - Hydrogen Bonds Detected            : {len(self.h_bonds)}\n")
            f.write(f"  - Salt Bridges Detected              : {len(self.salt_bridges)}\n")
            f.write(f"  - Non-Bonded/Hydrophobic Contacts    : {len(unique_hydro)}\n")
            f.write(f"  - Total Interface Residues ({self.ligand_label:<8}): {len(self.l_interface_residues)}\n")
            f.write(f"  - Total Interface Residues ({self.receptor_label:<8}): {len(self.r_interface_residues)}\n\n")

            f.write("2. INTERFACIAL CONTACTS REGISTRY\n")
            f.write("-" * 72 + "\n")
            f.write("A. HYDROGEN BONDS\n")
            if not self.h_bonds:
                f.write("  No hydrogen bonds detected.\n")
            else:
                for hb in self.h_bonds:
                    f.write(f"  - {self.ligand_label} [Chain {hb['l_chain']}] {hb['l_res']} ({hb['l_atom']}) "
                            f"<--- {hb['dist']:.2f} A ---> "
                            f"{self.receptor_label} [Chain {hb['r_chain']}] {hb['r_res']} ({hb['r_atom']})\n")
            f.write("\n")

            f.write("B. SALT BRIDGES\n")
            if not self.salt_bridges:
                f.write("  No salt bridges detected.\n")
            else:
                for sb in self.salt_bridges:
                    f.write(f"  - {self.ligand_label} [Chain {sb['l_chain']}] {sb['l_res']} ({sb['l_atom']}) "
                            f"<--- {sb['dist']:.2f} A ---> "
                            f"{self.receptor_label} [Chain {sb['r_chain']}] {sb['r_res']} ({sb['r_atom']})\n")
            f.write("\n")

            f.write("C. NON-BONDED CONTACT ENVELOPE (up to 40 shown)\n")
            if not unique_hydro:
                f.write("  No hydrophobic contacts detected.\n")
            else:
                for hc in unique_hydro[:40]:
                    f.write(f"  - {self.ligand_label} [Chain {hc['l_chain']}] {hc['l_res']} ... "
                            f"{hc['dist']:.2f} A ... "
                            f"{self.receptor_label} [Chain {hc['r_chain']}] {hc['r_res']}\n")
                if len(unique_hydro) > 40:
                    f.write(f"  - ... ({len(unique_hydro) - 40} additional contacts omitted)\n")
            f.write("\n")

            f.write("3. BINDING RESIDUE CATALOG\n")
            f.write("-" * 72 + "\n")
            f.write(f"{self.ligand_label} Contacting Residues:\n")
            v_sorted = sorted(self.l_interface_residues, key=lambda x: x[1])
            f.write("  " + ", ".join(format_res_label(r[0], r[1]) for r in v_sorted) + "\n\n")

            f.write(f"{self.receptor_label} Contacting Residues:\n")
            t_sorted = sorted(self.r_interface_residues, key=lambda x: x[1])
            f.write("  " + ", ".join(format_res_label(r[0], r[1]) for r in t_sorted) + "\n")


# ------------------------------------------------------------------------
# Visualizations
# ------------------------------------------------------------------------
def render_summary_bubble_plot(analyzer, output_path):
    fig, ax = plt.subplots(figsize=(6.5, 5.5), dpi=300)
    ax.set_facecolor('#ffffff')

    ligand_res_count = len(analyzer.l_interface_residues)
    receptor_res_count = len(analyzer.r_interface_residues)
    sb_count = len(analyzer.salt_bridges)
    hb_count = len(analyzer.h_bonds)
    nb_count = len(analyzer.unique_hydrophobic_contacts())

    v_circle = Ellipse((1.5, 5.0), width=1.7, height=1.7, facecolor='#ca8aff',
                        edgecolor='#000000', linewidth=1.8, zorder=2)
    ax.add_patch(v_circle)
    ax.text(1.5, 6.1, analyzer.ligand_label, fontsize=15, color='#6b21a8', fontweight='black', ha='center')
    ax.text(1.5, 5.0, f"{ligand_res_count}res", fontsize=12, color='#000000', fontweight='black', ha='center', va='center')

    t_circle = Ellipse((4.5, 5.0), width=2.0, height=2.0, facecolor='#ff1744',
                        edgecolor='#000000', linewidth=1.8, zorder=2)
    ax.add_patch(t_circle)
    ax.text(4.5, 6.2, analyzer.receptor_label, fontsize=15, color='#b91c1c', fontweight='black', ha='center')
    ax.text(4.5, 5.0, f"{receptor_res_count}res", fontsize=12, color='#000000', fontweight='black', ha='center', va='center')

    ax.plot([2.3, 3.5], [5.3, 5.3], color='#ef4444', lw=4.5, zorder=3)
    ax.text(2.9, 5.3, f" {sb_count} ", color='white', fontsize=9.5, fontweight='bold', ha='center', va='center',
            bbox=dict(boxstyle='square,pad=0.2', facecolor='#ef4444', edgecolor='none'))

    ax.plot([2.35, 3.5], [5.0, 5.0], color='#29b6f6', lw=4.5, zorder=3)
    ax.text(2.9, 5.0, f" {hb_count} ", color='white', fontsize=9.5, fontweight='bold', ha='center', va='center',
            bbox=dict(boxstyle='square,pad=0.2', facecolor='#29b6f6', edgecolor='none'))

    ax.plot([2.3, 3.5], [4.7, 4.7], color='#f59e0b', lw=4.5, zorder=3)
    ax.text(2.9, 4.7, f" {nb_count} ", color='white', fontsize=9.5, fontweight='bold', ha='center', va='center',
            bbox=dict(boxstyle='square,pad=0.2', facecolor='#f59e0b', edgecolor='none'))

    ax.set_xlim(0.2, 5.8)
    ax.set_ylim(2.5, 6.8)
    ax.axis('off')

    ly = 3.0
    ax.plot([0.5, 0.9], [ly, ly], color='#ef4444', lw=3.5)
    ax.text(1.0, ly, "Salt\nbridges", fontsize=9, va='center', fontweight='bold', color='#1e293b')
    ax.plot([1.7, 2.1], [ly, ly], color='#29b6f6', lw=3.5)
    ax.text(2.2, ly, "Hydrogen\nbonds", fontsize=9, va='center', fontweight='bold', color='#1e293b')
    ax.plot([3.1, 3.5], [ly, ly], color='#ffa726', lw=3.0, linestyle=':')
    ax.text(3.6, ly, "Non-bonded\ncontacts", fontsize=9, va='center', fontweight='bold', color='#1e293b')
    ax.text(0.4, ly + 0.4, "Key:", fontsize=12, fontweight='bold', color='#000000')
    ax.set_title(analyzer.complex_label, fontsize=11, fontweight='bold', color='#334155', pad=8)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def draw_bezier_curve(ax, x1, y1, x2, y2, color, lw, linestyle='-', alpha=0.8):
    Path = mpath.Path
    dx = abs(x2 - x1) * 0.45
    path_data = [
        (Path.MOVETO, (x1, y1)),
        (Path.CURVE4, (x1 + dx, y1)),
        (Path.CURVE4, (x2 - dx, y2)),
        (Path.CURVE4, (x2, y2)),
    ]
    codes, verts = zip(*path_data)
    patch = mpatches.PathPatch(mpath.Path(verts, codes), edgecolor=color, facecolor='none',
                                lw=lw, linestyle=linestyle, alpha=alpha, zorder=1)
    ax.add_patch(patch)


def render_detailed_interactions_network(analyzer, output_path):
    v_res_list = sorted(analyzer.l_interface_residues, key=lambda x: x[1])
    t_res_list = sorted(analyzer.r_interface_residues, key=lambda x: x[1])
    if not v_res_list or not t_res_list:
        print("[Warning] No interface contacts detected; skipping network plot.")
        return

    v_labels = [format_res_label(r[0], r[1]) for r in v_res_list]
    t_labels = [format_res_label(r[0], r[1]) for r in t_res_list]
    num_v, num_t = len(v_labels), len(t_labels)
    max_rows = max(num_v, num_t)
    fig_height = max(8.5, max_rows * 0.38)

    fig, ax = plt.subplots(figsize=(7.5, fig_height), dpi=300)
    ax.set_facecolor('#ffffff')

    y_v = np.linspace(max_rows - 1, 0, num_v) if num_v > 1 else [max_rows / 2]
    y_t = np.linspace(max_rows - 1, 0, num_t) if num_t > 1 else [max_rows / 2]
    v_pos = {v_labels[i]: y_v[i] for i in range(num_v)}
    t_pos = {t_labels[i]: y_t[i] for i in range(num_t)}

    for name, y in v_pos.items():
        color = RESIDUE_COLORS.get(name[:3].upper(), '#8e8e8e')
        ax.text(1.0, y, f"  {name}  ", ha='right', va='center', fontsize=9.5, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.35', facecolor=color, edgecolor='#2d3748', lw=1.2, alpha=0.9))
    for name, y in t_pos.items():
        color = RESIDUE_COLORS.get(name[:3].upper(), '#8e8e8e')
        ax.text(2.0, y, f"  {name}  ", ha='left', va='center', fontsize=9.5, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.35', facecolor=color, edgecolor='#2d3748', lw=1.2, alpha=0.9))

    seen_hydro = set()
    for hc in analyzer.hydrophobic_contacts:
        hc_r = format_res_label(hc['r_res'].split()[0], hc['r_res'].split()[1])
        hc_l = format_res_label(hc['l_res'].split()[0], hc['l_res'].split()[1])
        if hc_l in v_pos and hc_r in t_pos and (hc_l, hc_r) not in seen_hydro:
            seen_hydro.add((hc_l, hc_r))
            draw_bezier_curve(ax, 1.01, v_pos[hc_l], 1.99, t_pos[hc_r], '#f59e0b', 0.8, linestyle=':', alpha=0.55)

    for hb in analyzer.h_bonds:
        hb_r = format_res_label(hb['r_res'].split()[0], hb['r_res'].split()[1])
        hb_l = format_res_label(hb['l_res'].split()[0], hb['l_res'].split()[1])
        if hb_l in v_pos and hb_r in t_pos:
            draw_bezier_curve(ax, 1.01, v_pos[hb_l], 1.99, t_pos[hb_r], '#0ea5e9', 2.2, linestyle='-', alpha=0.85)

    for sb in analyzer.salt_bridges:
        sb_r = format_res_label(sb['r_res'].split()[0], sb['r_res'].split()[1])
        sb_l = format_res_label(sb['l_res'].split()[0], sb['l_res'].split()[1])
        if sb_l in v_pos and sb_r in t_pos:
            draw_bezier_curve(ax, 1.01, v_pos[sb_l], 1.99, t_pos[sb_r], '#ef4444', 3.0, linestyle='-', alpha=0.95)

    ax.set_xlim(0.3, 2.7)
    ax.set_ylim(-0.8, max_rows)
    ax.axis('off')
    header_y = max_rows - 0.2
    ax.text(0.85, header_y + 0.3, f"{analyzer.ligand_label} (Chain {'/'.join(analyzer.ligand_chains)})",
            fontsize=13, color='#6b21a8', fontweight='black', ha='right')
    ax.text(2.15, header_y + 0.3, f"{analyzer.receptor_label} (Chain {'/'.join(analyzer.receptor_chains)})",
            fontsize=13, color='#b91c1c', fontweight='black', ha='left')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def render_contact_heatmap(analyzer, output_path):
    v_res_list = sorted(analyzer.l_interface_residues, key=lambda x: x[1])
    t_res_list = sorted(analyzer.r_interface_residues, key=lambda x: x[1])
    if not v_res_list or not t_res_list:
        print("[Warning] No interface contacts detected; skipping heatmap.")
        return

    v_labels = [format_res_label(r[0], r[1]) for r in v_res_list]
    t_labels = [format_res_label(r[0], r[1]) for r in t_res_list]
    v_idx = {v_labels[i]: i for i in range(len(v_labels))}
    t_idx = {t_labels[i]: i for i in range(len(t_labels))}
    matrix = np.zeros((len(t_labels), len(v_labels)))

    for hc in analyzer.hydrophobic_contacts:
        hc_r = format_res_label(hc['r_res'].split()[0], hc['r_res'].split()[1])
        hc_l = format_res_label(hc['l_res'].split()[0], hc['l_res'].split()[1])
        if hc_l in v_idx and hc_r in t_idx:
            matrix[t_idx[hc_r], v_idx[hc_l]] += 1.0
    for hb in analyzer.h_bonds:
        hb_r = format_res_label(hb['r_res'].split()[0], hb['r_res'].split()[1])
        hb_l = format_res_label(hb['l_res'].split()[0], hb['l_res'].split()[1])
        if hb_l in v_idx and hb_r in t_idx:
            matrix[t_idx[hb_r], v_idx[hb_l]] += 2.0
    for sb in analyzer.salt_bridges:
        sb_r = format_res_label(sb['r_res'].split()[0], sb['r_res'].split()[1])
        sb_l = format_res_label(sb['l_res'].split()[0], sb['l_res'].split()[1])
        if sb_l in v_idx and sb_r in t_idx:
            matrix[t_idx[sb_r], v_idx[sb_l]] += 3.0

    fig, ax = plt.subplots(figsize=(8.5, 7.5), dpi=300)
    im = ax.imshow(matrix, cmap=plt.cm.YlOrRd, origin='lower', aspect='auto')
    ax.set_xticks(np.arange(len(v_labels)))
    ax.set_yticks(np.arange(len(t_labels)))
    ax.set_xticklabels(v_labels, rotation=90, fontsize=8, fontweight='medium')
    ax.set_yticklabels(t_labels, fontsize=8, fontweight='medium')
    ax.set_xticks(np.arange(len(v_labels) + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(len(t_labels) + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="#e2e8f0", linestyle='-', linewidth=0.5)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_xlabel(f"{analyzer.ligand_label} Interfacial Residues", fontsize=11, fontweight='bold', labelpad=10)
    ax.set_ylabel(f"{analyzer.receptor_label} Interfacial Residues", fontsize=11, fontweight='bold', labelpad=10)
    ax.set_title(f"Intermolecular Contact Density Hotspots\n{analyzer.complex_label}", fontsize=12, fontweight='black', pad=15)
    cbar = fig.colorbar(im, ax=ax, pad=0.03, shrink=0.85)
    cbar.set_label("Relative Binding Contact Weighting", fontsize=10, fontweight='bold')
    cbar.ax.tick_params(labelsize=8)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def render_3d_complex_plot(analyzer, output_path):
    fig = plt.figure(figsize=(7, 7), dpi=300)
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('#ffffff')

    ligand_chain = analyzer.ligand_chains[0]
    receptor_chain = analyzer.receptor_chains[0]

    l_coords = analyzer.ca_coordinates.get(ligand_chain, [])
    if l_coords:
        coords = np.array([pt[1] for pt in l_coords])
        ax.plot(coords[:, 0], coords[:, 1], coords[:, 2], color='#ca8aff', lw=3.0, alpha=0.9,
                label=f'{analyzer.ligand_label} (Chain {ligand_chain})')
        ax.plot(coords[:, 0], coords[:, 1], coords[:, 2], color='#512da8', lw=1.0, alpha=0.4)

    r_coords = analyzer.ca_coordinates.get(receptor_chain, [])
    if r_coords:
        coords = np.array([pt[1] for pt in r_coords])
        ax.plot(coords[:, 0], coords[:, 1], coords[:, 2], color='#ff1744', lw=3.0, alpha=0.9,
                label=f'{analyzer.receptor_label} (Chain {receptor_chain})')
        ax.plot(coords[:, 0], coords[:, 1], coords[:, 2], color='#b71c1c', lw=1.0, alpha=0.4)

    if analyzer.interface_ca_coords:
        int_coords = np.array(analyzer.interface_ca_coords)
        ax.scatter(int_coords[:, 0], int_coords[:, 1], int_coords[:, 2], color='#facc15', s=35,
                   edgecolor='#000000', lw=0.6, alpha=0.85, label='Binding Interface')

    ax.grid(False)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.set_axis_off()
    ax.set_title(f"3D Coordinates Docking Trajectory\n{analyzer.complex_label}", fontsize=12, fontweight='black', pad=10)
    ax.legend(loc='lower center', bbox_to_anchor=(0.5, 0.02), ncol=3, frameon=True, fontsize=8)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>__COMPLEX_LABEL__ - Interactive 3D Interface Viewer</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.6.0/jquery.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js"></script>
<style>
  body { margin: 0; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; background:#f8fafc; }
  header { position: absolute; top: 20px; left: 20px; z-index: 10; background: rgba(255,255,255,0.95);
           padding: 15px 25px; border-radius: 12px; border: 1px solid #e2e8f0; box-shadow: 0 4px 20px rgba(0,0,0,0.08); }
  h1 { margin: 0; font-size: 1.15rem; font-weight: 800; color: #1e293b; }
  p { margin: 4px 0 0 0; font-size: 0.75rem; color: #64748b; }
  #viewer-container { width: 100vw; height: 100vh; position: absolute; top: 0; left: 0; z-index: 1; }
  #controls { position: absolute; bottom: 30px; left: 20px; z-index: 10; background: rgba(255,255,255,0.95);
              padding: 18px; border-radius: 12px; border: 1px solid #e2e8f0; display: flex; flex-direction: column;
              gap: 10px; width: 240px; box-shadow: 0 4px 20px rgba(0,0,0,0.08); }
  .control-label { font-size: 0.7rem; font-weight: 700; text-transform: uppercase; color: #64748b; }
  .checkbox-item { display: flex; align-items: center; gap: 8px; font-size: 0.8rem; color: #334155; cursor: pointer; }
  select, button { background: #fff; color: #334155; border: 1px solid #cbd5e1; padding: 6px 10px; border-radius: 6px;
                   font-size: 0.75rem; font-weight: 600; cursor: pointer; }
  button.primary { background: #10b981; color: #fff; border: none; }
  .legend { position: absolute; top: 20px; right: 20px; z-index: 10; background: rgba(255,255,255,0.95);
            padding: 14px; border-radius: 12px; border: 1px solid #e2e8f0; font-size: 0.75rem; }
  .legend-item { display: flex; align-items: center; gap: 8px; margin-bottom: 6px; }
  .legend-color { width: 12px; height: 12px; border-radius: 3px; }
</style>
</head>
<body>
<header><h1>Interactive 3D Interface Viewer</h1><p>__LIGAND_LABEL__ (purple) vs __RECEPTOR_LABEL__ (red)</p></header>
<div class="legend">
  <div class="legend-item"><div class="legend-color" style="background:#ca8aff;"></div><span>__LIGAND_LABEL__ (Chain __LIGAND_CHAIN__)</span></div>
  <div class="legend-item"><div class="legend-color" style="background:#ff1744;"></div><span>__RECEPTOR_LABEL__ (Chain __RECEPTOR_CHAIN__)</span></div>
  <div class="legend-item"><div class="legend-color" style="background:#facc15;"></div><span>Interface residues</span></div>
  <div class="legend-item"><div class="legend-color" style="background:#0ea5e9;"></div><span>H-bonds</span></div>
</div>
<div id="controls">
  <span class="control-label">Visibility</span>
  <label class="checkbox-item"><input type="checkbox" checked onchange="toggleLigand(this)"> __LIGAND_LABEL__</label>
  <label class="checkbox-item"><input type="checkbox" checked onchange="toggleReceptor(this)"> __RECEPTOR_LABEL__</label>
  <label class="checkbox-item"><input type="checkbox" checked onchange="toggleInterface(this)"> Interface sticks</label>
  <label class="checkbox-item"><input type="checkbox" checked onchange="toggleBonds(this)"> H-bonds</label>
  <span class="control-label">Style</span>
  <select onchange="changeStyle(this.value)">
    <option value="cartoon">Cartoon</option>
    <option value="sphere">Spacefill</option>
    <option value="stick">Stick</option>
    <option value="line">Line</option>
  </select>
  <div style="display:flex; gap:8px;">
    <button onclick="toggleSpin()" style="flex:1;">Spin</button>
    <button onclick="viewer.zoomTo(); viewer.render();" style="flex:1;">Reset</button>
  </div>
  <button class="primary" onclick="saveImage()">Save PNG</button>
</div>
<div id="viewer-container"><div id="mol-viewer" style="width:100%; height:100%;"></div></div>
<script>
  let viewer=null, spinning=false, showL=true, showR=true, showI=true, showB=true, style="cartoon";
  const rawPdb = "PDB_DATA_TOKEN";
  const ligandInterfaceRes = LIGAND_RES_TOKEN;
  const receptorInterfaceRes = RECEPTOR_RES_TOKEN;
  const hbonds = HBOND_TOKEN;
  const LIGAND_CHAIN = "__LIGAND_CHAIN__", RECEPTOR_CHAIN = "__RECEPTOR_CHAIN__";

  $(document).ready(function() {
    viewer = $3Dmol.createViewer("mol-viewer", {backgroundColor: "#ffffff"});
    viewer.addModel(rawPdb, "pdb");
    styleComplex();
    viewer.zoomTo();
    viewer.render();
  });

  function styleComplex() {
    viewer.clear();
    viewer.addModel(rawPdb, "pdb");
    let lStyle = {}, rStyle = {};
    if (style === "cartoon") { lStyle = {cartoon:{color:'#ca8aff'}}; rStyle = {cartoon:{color:'#ff1744'}}; }
    else if (style === "sphere") { lStyle = {sphere:{color:'#ca8aff'}}; rStyle = {sphere:{color:'#ff1744'}}; }
    else if (style === "stick") { lStyle = {stick:{colorscheme:'purpleCarbon'}}; rStyle = {stick:{colorscheme:'redCarbon'}}; }
    else { lStyle = {line:{color:'#ca8aff'}}; rStyle = {line:{color:'#ff1744'}}; }

    if (showL) {
      viewer.setStyle({chain: LIGAND_CHAIN}, lStyle);
      if (showI && style !== "stick" && style !== "sphere") {
        viewer.setStyle({chain: LIGAND_CHAIN, resi: ligandInterfaceRes}, {cartoon:{color:'#ca8aff'}, stick:{colorscheme:'purpleCarbon', radius:0.22}});
      }
    } else { viewer.setStyle({chain: LIGAND_CHAIN}, {}); }

    if (showR) {
      viewer.setStyle({chain: RECEPTOR_CHAIN}, rStyle);
      if (showI && style !== "stick" && style !== "sphere") {
        viewer.setStyle({chain: RECEPTOR_CHAIN, resi: receptorInterfaceRes}, {cartoon:{color:'#ff1744'}, stick:{colorscheme:'redCarbon', radius:0.22}});
      }
    } else { viewer.setStyle({chain: RECEPTOR_CHAIN}, {}); }

    if (showB) {
      hbonds.forEach(function(bond) {
        try {
          viewer.addCylinder({start:{chain: RECEPTOR_CHAIN, resi: bond.r_res, atom:'O'},
                               end:{chain: LIGAND_CHAIN, resi: bond.l_res, atom:'N'},
                               radius:0.08, color:'#0ea5e9', dashed:true, fromCap:true, toCap:true});
        } catch(err) {}
      });
    }
    viewer.render();
  }

  function toggleLigand(el){ showL = el.checked; styleComplex(); }
  function toggleReceptor(el){ showR = el.checked; styleComplex(); }
  function toggleInterface(el){ showI = el.checked; styleComplex(); }
  function toggleBonds(el){ showB = el.checked; styleComplex(); }
  function changeStyle(v){ style = v; styleComplex(); }
  function toggleSpin(){ spinning = !spinning; if (spinning) viewer.spin("y", 1); else viewer.spin(false); }
  function saveImage(){
    viewer.render();
    const canvas = document.querySelector("#mol-viewer canvas");
    if (canvas) {
      const link = document.createElement("a");
      link.download = "docked_complex.png";
      link.href = canvas.toDataURL('image/png');
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
    }
  }
</script>
</body>
</html>
"""


def generate_interactive_3d_html(analyzer, output_path):
    pdb_data_js = analyzer.raw_pdb_content.replace('\r', '').replace('\n', '\\n').replace("'", "\\'")
    l_res_ids = sorted({r[1] for r in analyzer.l_interface_residues})
    r_res_ids = sorted({r[1] for r in analyzer.r_interface_residues})

    hb_connections = []
    for hb in analyzer.h_bonds:
        try:
            hb_connections.append({
                "r_res": int(hb['r_res'].split()[1]),
                "l_res": int(hb['l_res'].split()[1]),
            })
        except (IndexError, ValueError):
            continue

    html = HTML_TEMPLATE
    html = html.replace("PDB_DATA_TOKEN", pdb_data_js)
    html = html.replace("LIGAND_RES_TOKEN", str(l_res_ids))
    html = html.replace("RECEPTOR_RES_TOKEN", str(r_res_ids))
    html = html.replace("HBOND_TOKEN", json.dumps(hb_connections))
    html = html.replace("__COMPLEX_LABEL__", analyzer.complex_label)
    html = html.replace("__LIGAND_LABEL__", analyzer.ligand_label)
    html = html.replace("__RECEPTOR_LABEL__", analyzer.receptor_label)
    html = html.replace("__LIGAND_CHAIN__", analyzer.ligand_chains[0])
    html = html.replace("__RECEPTOR_CHAIN__", analyzer.receptor_chains[0])

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)


# ------------------------------------------------------------------------
# Secondary structure (from PDB HELIX/SHEET header records)
# ------------------------------------------------------------------------
def parse_secondary_structure(raw_pdb_content):
    """Read HELIX/SHEET records straight from the PDB header.
    Returns {chain_id: [(start_resid, end_resid, 'H' or 'S'), ...]}.
    Docking output that has no header (no HELIX/SHEET lines) will simply
    come back empty -- callers should treat that as "not available", not
    as "all coil"."""
    ss_by_chain = collections.defaultdict(list)
    for line in raw_pdb_content.splitlines():
        rectype = line[0:6].strip()
        if rectype == "HELIX" and len(line) >= 38:
            chain_id = line[19].strip()
            try:
                start = int(line[21:25])
                end = int(line[33:37])
            except ValueError:
                continue
            ss_by_chain[chain_id].append((start, end, 'H'))
        elif rectype == "SHEET" and len(line) >= 38:
            chain_id = line[21].strip()
            try:
                start = int(line[22:26])
                end = int(line[33:37])
            except ValueError:
                continue
            ss_by_chain[chain_id].append((start, end, 'S'))
    return ss_by_chain


def summarize_secondary_structure(ss_by_chain, chains, chain_ids):
    """Returns {chain_id: {'helix': n, 'sheet': n, 'coil': n, 'total': n}}
    counting residues (not just record spans) for chains actually present."""
    summary = {}
    for chain_id in chain_ids:
        residues = chains.get(chain_id, [])
        protein_resids = sorted({r.resid for r in residues if not r.is_hetatm})
        total = len(protein_resids)
        helix_set, sheet_set = set(), set()
        for start, end, kind in ss_by_chain.get(chain_id, []):
            for resid in range(start, end + 1):
                if kind == 'H':
                    helix_set.add(resid)
                else:
                    sheet_set.add(resid)
        helix_n = len(helix_set & set(protein_resids))
        sheet_n = len(sheet_set & set(protein_resids))
        coil_n = max(total - helix_n - sheet_n, 0)
        summary[chain_id] = {'helix': helix_n, 'sheet': sheet_n, 'coil': coil_n, 'total': total}
    return summary


def render_secondary_structure_figure(analyzer, output_path):
    """Draws a simple per-chain secondary-structure track along the sequence:
    helices as rounded rectangles, strands as arrows, everything else as a
    thin coil line. Chains with no HELIX/SHEET records in the file are drawn
    as coil and flagged in the title (this happens often for raw docking
    output, which typically has no header)."""
    if analyzer.chains is None:
        print("[Warning] No parsed chains available; skipping secondary structure figure.")
        return

    ss_by_chain = parse_secondary_structure(analyzer.raw_pdb_content)
    chain_ids = analyzer.ligand_chains + analyzer.receptor_chains
    labels = {analyzer.ligand_chains[0]: analyzer.ligand_label,
              analyzer.receptor_chains[0]: analyzer.receptor_label}

    fig, axes = plt.subplots(len(chain_ids), 1, figsize=(9, 1.6 * len(chain_ids) + 1), dpi=300)
    if len(chain_ids) == 1:
        axes = [axes]

    any_ss_found = any(ss_by_chain.get(c) for c in chain_ids)

    for ax, chain_id in zip(axes, chain_ids):
        residues = [r for r in analyzer.chains.get(chain_id, []) if not r.is_hetatm]
        residues.sort(key=lambda r: r.resid)
        if not residues:
            ax.axis('off')
            continue
        resids = [r.resid for r in residues]
        rmin, rmax = min(resids), max(resids)

        ax.plot([rmin, rmax], [0, 0], color='#94a3b8', lw=2.5, zorder=1, solid_capstyle='round')

        for start, end, kind in ss_by_chain.get(chain_id, []):
            start, end = max(start, rmin), min(end, rmax)
            if start > end:
                continue
            if kind == 'H':
                ax.add_patch(mpatches.FancyBboxPatch((start, -0.28), end - start, 0.56,
                                                       boxstyle="round,pad=0,rounding_size=0.3",
                                                       facecolor='#ef4444', edgecolor='#7f1d1d', zorder=2))
            else:
                width = end - start
                ax.add_patch(mpatches.FancyArrow(start, 0, max(width, 0.4), 0, width=0.5,
                                                  head_width=0.9, head_length=min(1.2, max(width * 0.3, 0.4)),
                                                  length_includes_head=True,
                                                  facecolor='#facc15', edgecolor='#78350f', zorder=2))

        ax.set_xlim(rmin - 1, rmax + 1)
        ax.set_ylim(-1, 1)
        ax.set_yticks([])
        ax.set_xlabel(f"Residue number (Chain {chain_id})", fontsize=8)
        title = f"{labels.get(chain_id, chain_id)} (Chain {chain_id})"
        if not ss_by_chain.get(chain_id):
            title += "  [no HELIX/SHEET records in file]"
        ax.set_title(title, fontsize=10, fontweight='bold', loc='left')
        for spine in ('top', 'right', 'left'):
            ax.spines[spine].set_visible(False)

    handles = [
        mpatches.Patch(facecolor='#ef4444', edgecolor='#7f1d1d', label='Helix'),
        mpatches.Patch(facecolor='#facc15', edgecolor='#78350f', label='Sheet/strand'),
        mpatches.Patch(facecolor='#94a3b8', label='Coil/loop'),
    ]
    fig.legend(handles=handles, loc='upper right', ncol=3, fontsize=8, frameon=False)
    suptitle = f"Secondary Structure - {analyzer.complex_label}"
    if not any_ss_found:
        suptitle += "\n(source PDB has no HELIX/SHEET header records - shown as coil)"
    fig.suptitle(suptitle, fontsize=11, fontweight='black', y=1.02)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


# ------------------------------------------------------------------------
# Ramachandran plot (phi/psi backbone dihedrals)
# ------------------------------------------------------------------------
def _dihedral(p0, p1, p2, p3):
    """Dihedral angle (degrees) defined by four 3D points."""
    b0 = np.array(p0) - np.array(p1)
    b1 = np.array(p2) - np.array(p1)
    b2 = np.array(p3) - np.array(p2)
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return math.degrees(math.atan2(y, x))


def compute_phi_psi(chains, chain_ids):
    """Returns a list of dicts: {chain, resid, resname, phi, psi}.
    phi/psi are None at chain termini where the flanking residue/atom is
    missing. Requires contiguous backbone N/CA/C atoms; residues with gaps
    or missing backbone atoms are skipped for the angle(s) they can't form."""
    results = []
    for chain_id in chain_ids:
        residues = [r for r in chains.get(chain_id, []) if not r.is_hetatm]
        residues.sort(key=lambda r: r.resid)
        for i, res in enumerate(residues):
            if not all(a in res for a in ('N', 'CA', 'C')):
                continue
            phi = psi = None
            if i > 0:
                prev = residues[i - 1]
                if 'C' in prev and prev.resid == res.resid - 1:
                    phi = _dihedral(prev['C'].coord, res['N'].coord, res['CA'].coord, res['C'].coord)
            if i < len(residues) - 1:
                nxt = residues[i + 1]
                if 'N' in nxt and nxt.resid == res.resid + 1:
                    psi = _dihedral(res['N'].coord, res['CA'].coord, res['C'].coord, nxt['N'].coord)
            if phi is not None or psi is not None:
                results.append({'chain': chain_id, 'resid': res.resid, 'resname': res.get_resname(),
                                 'phi': phi, 'psi': psi})
    return results


# Rough, illustrative "favored-region" boxes (NOT a statistical potential --
# just enough to sanity-check gross outliers at a glance). Real Ramachandran
# validation should use a proper reference (e.g. MolProbity/wwPDB validation).
_RAMA_FAVORED_BOXES = [
    (-160, -40, 90, 180),    # beta sheet region
    (-160, -40, -180, -90),  # beta sheet region (wrap)
    (-100, -30, -60, 5),     # right-handed alpha helix
    (30, 100, 0, 90),        # left-handed helix (mostly Gly)
]


def _in_favored_region(phi, psi):
    for phi_lo, phi_hi, psi_lo, psi_hi in _RAMA_FAVORED_BOXES:
        if phi_lo <= phi <= phi_hi and psi_lo <= psi <= psi_hi:
            return True
    return False


def render_ramachandran_plot(phi_psi_list, analyzer, output_path):
    complete = [d for d in phi_psi_list if d['phi'] is not None and d['psi'] is not None]
    if not complete:
        print("[Warning] No complete phi/psi pairs computed; skipping Ramachandran plot.")
        return None

    fig, ax = plt.subplots(figsize=(6.5, 6.5), dpi=300)
    for lo1, hi1, lo2, hi2 in _RAMA_FAVORED_BOXES:
        ax.add_patch(mpatches.Rectangle((lo1, lo2), hi1 - lo1, hi2 - lo2,
                                         facecolor='#d1fae5', edgecolor='none', zorder=0))

    chain_colors = {analyzer.ligand_chains[0]: '#ca8aff', analyzer.receptor_chains[0]: '#ff1744'}
    chain_labels = {analyzer.ligand_chains[0]: analyzer.ligand_label, analyzer.receptor_chains[0]: analyzer.receptor_label}
    for chain_id, color in chain_colors.items():
        pts = [(d['phi'], d['psi']) for d in complete if d['chain'] == chain_id]
        if not pts:
            continue
        xs, ys = zip(*pts)
        ax.scatter(xs, ys, s=18, color=color, edgecolor='#1e293b', linewidth=0.3, alpha=0.85,
                   label=f"{chain_labels[chain_id]} (Chain {chain_id})", zorder=2)

    ax.axhline(0, color='#cbd5e1', lw=0.8, zorder=1)
    ax.axvline(0, color='#cbd5e1', lw=0.8, zorder=1)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xticks(range(-180, 181, 60))
    ax.set_yticks(range(-180, 181, 60))
    ax.set_xlabel("Phi (degrees)", fontsize=10, fontweight='bold')
    ax.set_ylabel("Psi (degrees)", fontsize=10, fontweight='bold')
    ax.set_title(f"Ramachandran Plot - {analyzer.complex_label}\n"
                 f"(shaded regions are an illustrative favored-angle guide, not a validated potential)",
                 fontsize=9.5, fontweight='bold')
    ax.legend(loc='upper right', fontsize=8, frameon=True)
    ax.set_aspect('equal')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    return complete


def write_ramachandran_summary(phi_psi_list, analyzer, output_path):
    complete = [d for d in phi_psi_list if d['phi'] is not None and d['psi'] is not None]
    chain_labels = {analyzer.ligand_chains[0]: analyzer.ligand_label, analyzer.receptor_chains[0]: analyzer.receptor_label}

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write("=" * 72 + "\n")
        f.write("      Ramachandran / Backbone Geometry Summary\n")
        f.write("=" * 72 + "\n\n")
        f.write(f"Analyzed Complex   : {analyzer.complex_label}\n\n")
        f.write("NOTE: 'favored region' below is an illustrative simplified guide\n")
        f.write("(rough alpha/beta/left-handed-helix boxes), not a validated\n")
        f.write("statistical potential. For publication-grade Ramachandran\n")
        f.write("validation, cross-check with MolProbity or the wwPDB validation report.\n\n")

        for chain_id in (analyzer.ligand_chains[0], analyzer.receptor_chains[0]):
            chain_pts = [d for d in complete if d['chain'] == chain_id]
            n = len(chain_pts)
            n_favored = sum(1 for d in chain_pts if _in_favored_region(d['phi'], d['psi']))
            n_outlier = n - n_favored
            f.write(f"{chain_labels[chain_id]} (Chain {chain_id}):\n")
            f.write(f"  - Residues with complete phi/psi : {n}\n")
            if n:
                f.write(f"  - In illustrative favored regions : {n_favored} ({100.0 * n_favored / n:.1f}%)\n")
                f.write(f"  - Outside favored regions          : {n_outlier} ({100.0 * n_outlier / n:.1f}%)\n")
                outlier_residues = [format_res_label(d['resname'], d['resid']) for d in chain_pts
                                     if not _in_favored_region(d['phi'], d['psi'])]
                if outlier_residues:
                    f.write("  - Outlier residues: " + ", ".join(outlier_residues) + "\n")
            f.write("\n")


# ------------------------------------------------------------------------
# Driver
# ------------------------------------------------------------------------
def analyze_one_complex(pdb_path, receptor_chains, ligand_chains, output_root,
                         ligand_label="Ligand", receptor_label="Receptor"):
    complex_label = os.path.splitext(os.path.basename(pdb_path))[0]
    print("=" * 72)
    print(f"  Analyzing: {complex_label}")
    print("=" * 72)

    analyzer = DockedInterfaceAnalyzer(pdb_path, receptor_chains, ligand_chains,
                                        ligand_label=ligand_label, receptor_label=receptor_label,
                                        complex_label=complex_label)
    out_dir = os.path.join(output_root, complex_label)
    os.makedirs(out_dir, exist_ok=True)

    try:
        analyzer.load_and_parse()
        analyzer.write_report(os.path.join(out_dir, "Complex_Interface_Report.txt"))
        render_summary_bubble_plot(analyzer, os.path.join(out_dir, "interface_summary_bubble.png"))
        render_detailed_interactions_network(analyzer, os.path.join(out_dir, "interface_interactions_network.png"))
        render_contact_heatmap(analyzer, os.path.join(out_dir, "interface_contact_heatmap.png"))
        render_3d_complex_plot(analyzer, os.path.join(out_dir, "docked_complex_3d_static.png"))
        generate_interactive_3d_html(analyzer, os.path.join(out_dir, "docked_complex_3d.html"))

        render_secondary_structure_figure(analyzer, os.path.join(out_dir, "secondary_structure.png"))
        chain_ids = analyzer.ligand_chains + analyzer.receptor_chains
        phi_psi = compute_phi_psi(analyzer.chains, chain_ids)
        render_ramachandran_plot(phi_psi, analyzer, os.path.join(out_dir, "ramachandran_plot.png"))
        write_ramachandran_summary(phi_psi, analyzer, os.path.join(out_dir, "Ramachandran_Summary.txt"))

        print(f"-> Done. Results in: {out_dir}\n")
        return analyzer
    except Exception as e:
        print(f"[Error] Failed to analyze {complex_label}: {e}")
        return None


def write_summary_across_complexes(analyzers, output_root):
    summary_path = os.path.join(output_root, "Summary_All_Complexes.txt")
    unique_analyzers = [a for a in analyzers if a is not None]
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("=" * 72 + "\n")
        f.write("      Summary Comparison Across All Analyzed Complexes\n")
        f.write("=" * 72 + "\n\n")
        header = f"{'Complex':<38}{'H-Bonds':>10}{'SaltBridges':>14}{'Hydrophobic':>14}{'Ligand-Res':>12}{'Receptor-Res':>14}\n"
        f.write(header)
        f.write("-" * len(header) + "\n")
        for a in unique_analyzers:
            unique_hydro = len(a.unique_hydrophobic_contacts())
            f.write(f"{a.complex_label:<38}{len(a.h_bonds):>10}{len(a.salt_bridges):>14}{unique_hydro:>14}"
                    f"{len(a.l_interface_residues):>12}{len(a.r_interface_residues):>14}\n")
    print(f"\n-> Combined cross-complex summary written: {summary_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Dependency-light PDBsum-style interface analyzer for docked/crystal complexes.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("input", help="Path to a single .pdb file, or a folder (use --batch).")
    parser.add_argument("--batch", action="store_true", help="Treat 'input' as a folder and process every .pdb file in it.")
    parser.add_argument("--chain-map", help="JSON file mapping {filename: {\"ligand\": \"A\", \"receptor\": \"W\"}} for batch mode.")
    parser.add_argument("--ligand-chain", help="Chain ID of the ligand/peptide (single-file mode). Auto-detected if omitted.")
    parser.add_argument("--receptor-chain", help="Chain ID of the receptor/target (single-file mode). Auto-detected if omitted.")
    parser.add_argument("--ligand-label", default="Ligand", help="Display label for the ligand chain (default: Ligand).")
    parser.add_argument("--receptor-label", default="Receptor", help="Display label for the receptor chain (default: Receptor).")
    parser.add_argument("--output", default="interface_results", help="Output directory (default: ./interface_results).")
    args = parser.parse_args()

    output_root = args.output
    os.makedirs(output_root, exist_ok=True)

    if args.batch:
        chain_map = {}
        if args.chain_map:
            with open(args.chain_map) as f:
                chain_map = json.load(f)

        pdb_files = sorted(glob.glob(os.path.join(args.input, "*.pdb")))
        if not pdb_files:
            print(f"[Error] No .pdb files found in {args.input}")
            sys.exit(1)

        results = []
        for pdb_path in pdb_files:
            fname = os.path.basename(pdb_path)
            entry = chain_map.get(fname, {})
            ligand_chains = [entry["ligand"]] if "ligand" in entry else None
            receptor_chains = [entry["receptor"]] if "receptor" in entry else None
            analyzer = analyze_one_complex(pdb_path, receptor_chains, ligand_chains, output_root,
                                            ligand_label=args.ligand_label, receptor_label=args.receptor_label)
            results.append(analyzer)
        write_summary_across_complexes(results, output_root)
    else:
        ligand_chains = [args.ligand_chain] if args.ligand_chain else None
        receptor_chains = [args.receptor_chain] if args.receptor_chain else None
        analyze_one_complex(args.input, receptor_chains, ligand_chains, output_root,
                             ligand_label=args.ligand_label, receptor_label=args.receptor_label)

    print("\n" + "=" * 72)
    print("                 ALL COMPLEXES PROCESSED!")
    print("=" * 72)
    print(f"Results directory: {os.path.abspath(output_root)}")


if __name__ == "__main__":
    main()
