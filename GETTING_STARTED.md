# Getting Started (no coding experience required)

This guide assumes you've never used the command line before. It walks through
everything from scratch. If you get stuck, re-read the step slowly -- these
tools are picky about exact spelling and spacing, which trips everyone up at
first.

## What this tool actually does, in plain English

You give it a **PDB file** -- a plain text file describing the 3D positions of
every atom in a protein structure (a crystal structure you downloaded, or the
output of a docking program). The file usually contains two or more molecules
glued together, called **chains**, each labeled with a single letter (A, B, W...).

The tool looks at the gap between two chains -- say, a small peptide (chain B)
sitting on a larger target protein (chain A) -- and works out exactly which
atoms are close enough to be touching, what kind of contact each one is
(hydrogen bond, salt bridge, or just a loose hydrophobic contact), and draws
you several pictures of that.

## Step 1 -- Install Python

Skip this if `python --version` in a terminal already prints something like
`Python 3.10.x` or higher.

- **Windows**: go to [python.org/downloads](https://www.python.org/downloads/),
  download the installer, run it, and **tick the box that says "Add Python to
  PATH"** before clicking Install. This checkbox is the single most common
  thing people forget.
- **Mac**: go to the same link, or if you have [Homebrew](https://brew.sh)
  installed, run `brew install python` in Terminal.
- **Linux**: Python is almost always already installed. Check with
  `python3 --version`.

## Step 2 -- Download this project

You don't need `git` for this. On the project's GitHub page:

1. Click the green **"Code"** button
2. Click **"Download ZIP"**
3. Find the downloaded `.zip` file (usually in your Downloads folder) and
   double-click it to unzip it into a regular folder

## Step 3 -- Open a terminal in that folder

- **Windows**: open the unzipped folder in File Explorer, click the address
  bar at the top, type `cmd`, and press Enter. A black terminal window opens
  already in the right folder.
- **Mac**: open the folder in Finder, right-click inside it, and choose
  "New Terminal at Folder" (or open Terminal and type `cd ` then drag the
  folder into the window and press Enter).

## Step 4 -- Install the two libraries this tool needs

Type this exactly and press Enter:

```
pip install -r requirements.txt
```

Wait for it to finish (you'll see a bunch of text scroll by, ending in
something like `Successfully installed ...`). If `pip` isn't recognized, try
`pip3` instead, or `python -m pip install -r requirements.txt`.

## Step 5 -- Get a PDB file to try

If you already have your own `.pdb` file (from a docking program, or
downloaded from [rcsb.org](https://www.rcsb.org)), you can use that. Otherwise,
use the one already included at `examples/data/1YCR_MDM2_p53.pdb` to try the
tool out first.

## Step 6 -- Run the easy version

In the same terminal window, type:

```
python src/easy_start.py
```

It will ask you a few plain questions one at a time:

- **"Path to your .pdb file"** -- type the location of your file, e.g.
  `examples/data/1YCR_MDM2_p53.pdb`, or drag-and-drop the file into the
  terminal window (this types the full path for you automatically), then
  press Enter.
- **"Do you already know the chain letters?"** -- see the box below if you're
  not sure. Type `y` or `n` and press Enter.
- If yes, it asks which letter is your **ligand/peptide** and which is your
  **receptor/target**.
- It asks what to call each one in the report (just a label, type anything
  you like, or press Enter to accept the default shown in brackets).
- It asks where to save the results (press Enter to accept the default).

Then it runs, and tells you exactly what got created and what each file means.

### How do I find the chain letters in my file?

Open your `.pdb` file in any plain text editor (Notepad, TextEdit, VS Code --
**not** Microsoft Word). Scroll until you find a line starting with `ATOM`,
which looks like this:

```
ATOM      1  N   GLU A  25      10.801 -12.147  -5.180  1.00 49.08           N
```

Counting from the start of the line, the **chain letter** is the single
character right after the 3-letter residue name (`GLU` here) -- in this
example, it's `A`. Scroll further down and look for where that letter changes
to a different one (e.g. `B`) -- that's your second chain. If your file has
more than two chains, note all their letters, then decide which one is your
ligand and which is your target.

## Step 7 -- Look at your results

Open the folder the script told you about (by default,
`interface_results/<your file name>/`). Double-click any `.png` file to view
it like a normal picture, and any `.txt` file to read it like a normal
document. Double-click `docked_complex_3d.html` to open an interactive,
spinnable 3D view in your web browser.

## Common problems

- **"python is not recognized"** -- Python wasn't added to PATH during
  install. Re-run the Python installer and make sure to tick that checkbox
  (Windows), or use `python3` instead of `python` (Mac/Linux).
- **"No such file or directory"** when giving your PDB path -- you're
  probably not in the right folder, or there's a typo. Try dragging the file
  into the terminal window instead of typing the path by hand.
- **"Chain(s) [...] not found in [file]"** -- the letter you typed doesn't
  exist in your file. Re-check using the steps above (chain letters are
  case-sensitive).
- **Everything shows as "coil" in the secondary structure picture** -- this
  is normal for files produced by docking software; that information only
  exists in real crystal-structure files. Not a bug.

## Ready for more control?

Once you're comfortable, `src/pdb_interface_analyzer.py` (run with
`python src/pdb_interface_analyzer.py --help`) gives you batch mode (analyze
a whole folder of complexes at once) and finer control via command-line
flags. See the main [README.md](README.md) for details, or use the
[Google Colab notebook](notebooks/PDBsum_Interface_Analyzer_Colab.ipynb) to
skip local installation entirely.
