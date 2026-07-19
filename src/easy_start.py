#!/usr/bin/env python3
"""
easy_start.py
=============
A no-flags, question-and-answer version of pdb_interface_analyzer.py for
people who don't want to learn command-line options.

Just run:

    python easy_start.py

...and answer the plain-English questions it asks. It calls the exact same
analysis code as pdb_interface_analyzer.py --  this is only a friendlier
front door, not a different tool. If you're comfortable with the command
line, pdb_interface_analyzer.py directly (see the --help text, or
GETTING_STARTED.md) gives you more control (batch mode, custom chain maps).

New to all this? Read GETTING_STARTED.md in this repo first -- it explains,
in plain English, what a "chain" and a "PDB file" even are.
"""

import os
import sys

# Let this script be run either from inside src/ or from the repo root.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    import pdb_interface_analyzer as pia
except ImportError:
    print("Could not find pdb_interface_analyzer.py. Make sure easy_start.py")
    print("is sitting in the same folder as pdb_interface_analyzer.py.")
    sys.exit(1)


def ask(prompt, default=None):
    suffix = f" [{default}]" if default is not None else ""
    answer = input(f"{prompt}{suffix}: ").strip()
    return answer if answer else default


def ask_yes_no(prompt, default_yes=True):
    default_label = "Y/n" if default_yes else "y/N"
    answer = input(f"{prompt} ({default_label}): ").strip().lower()
    if not answer:
        return default_yes
    return answer.startswith("y")


def main():
    print("=" * 72)
    print("  PDB Interface Analyzer -- easy start")
    print("=" * 72)
    print("This will ask a few questions, then run the same analysis as the")
    print("main script and tell you exactly where to find your results.\n")

    # ---- 1. Get the PDB file ----
    while True:
        pdb_path = ask("Path to your .pdb file (drag-and-drop the file into this window works too)")
        if not pdb_path:
            print("  Please enter a file path.\n")
            continue
        pdb_path = pdb_path.strip('"').strip("'")
        if not os.path.exists(pdb_path):
            print(f"  Couldn't find a file at: {pdb_path}\n  Please check the path and try again.\n")
            continue
        break

    # ---- 2. Figure out chains ----
    print("\nA PDB file is usually made of a few 'chains' -- separate molecules in")
    print("the same file, each labeled with a single letter (A, B, W, etc). We need")
    print("to know which chain is your peptide/ligand, and which is the target/receptor.\n")

    knows_chains = ask_yes_no("Do you already know the chain letters?", default_yes=False)

    ligand_chain = receptor_chain = None
    if knows_chains:
        ligand_chain = ask("Which chain letter is the ligand/peptide? (e.g. B)")
        receptor_chain = ask("Which chain letter is the receptor/target? (e.g. A)")

    ligand_label = ask("What should we call the ligand/peptide in the report?", default="Ligand")
    receptor_label = ask("What should we call the receptor/target in the report?", default="Receptor")

    output_root = ask("Folder to save results in", default="interface_results")

    # ---- 3. Run it ----
    print("\nRunning the analysis now...\n")
    analyzer = pia.analyze_one_complex(
        pdb_path,
        receptor_chains=[receptor_chain] if receptor_chain else None,
        ligand_chains=[ligand_chain] if ligand_chain else None,
        output_root=output_root,
        ligand_label=ligand_label,
        receptor_label=receptor_label,
    )

    if analyzer is None:
        print("\nSomething went wrong -- scroll up for the error message.")
        print("Common fixes:")
        print("  - Double check the chain letters actually exist in your file")
        print("    (open the .pdb file in a text editor and look at column 22")
        print("    of any line starting with 'ATOM' -- that's the chain letter).")
        print("  - Make sure the file is a standard PDB coordinate file, not a")
        print("    PDF, image, or Word document with a similar name.")
        sys.exit(1)

    result_dir = os.path.join(output_root, analyzer.complex_label)

    # ---- 4. Explain the results in plain English ----
    print("\n" + "=" * 72)
    print("  Done! Here's what was created and what each file means:")
    print("=" * 72)
    print(f"\nEverything is saved in: {os.path.abspath(result_dir)}\n")
    print(f"  Complex_Interface_Report.txt")
    print(f"      Plain-text summary: every hydrogen bond, salt bridge, and")
    print(f"      close contact found between {ligand_label} and {receptor_label},")
    print(f"      plus the full list of residues involved on each side.")
    print(f"\n  interface_summary_bubble.png")
    print(f"      A quick-glance picture: two bubbles (one per chain) with")
    print(f"      lines between them counting each type of contact.")
    print(f"\n  interface_interactions_network.png")
    print(f"      Every contacting residue pair drawn as a connected diagram,")
    print(f"      similar to a PDBsum 'ligplot'.")
    print(f"\n  interface_contact_heatmap.png")
    print(f"      A grid showing which residues touch which -- darker squares")
    print(f"      mean stronger/more numerous contacts at that residue pair.")
    print(f"\n  docked_complex_3d_static.png")
    print(f"      A snapshot of the actual 3D shape of the complex.")
    print(f"\n  docked_complex_3d.html")
    print(f"      Double-click this file to open an interactive, rotatable")
    print(f"      3D view in your web browser -- no installation needed.")
    print(f"\n  secondary_structure.png")
    print(f"      Shows which parts of each chain are folded into a helix,")
    print(f"      a flat strand, or neither (coil). If your file came from a")
    print(f"      docking program rather than a real crystal structure, this")
    print(f"      may show everything as coil -- that's normal, not an error.")
    print(f"\n  ramachandran_plot.png  and  Ramachandran_Summary.txt")
    print(f"      A structural sanity-check plot of backbone angles, plus a")
    print(f"      plain-text count of how many residues look geometrically")
    print(f"      normal vs. unusual.")
    print("\nOpen the .png files like any image, and the .txt files in any text")
    print("editor (Notepad, TextEdit, VS Code, etc).")


if __name__ == "__main__":
    main()
