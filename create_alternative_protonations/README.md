# CREATE ALTERNATIVE PROTONATION STATES OF A RECEPTOR

This is a Python script to create all alternative protonation state combinations of a protein given a ligand and a specified radius. It could be useful in case you have a receptor but you are not sure about the protonation states of some residues in the binding site and you want to do docking of MD using all of them.

The script must be executed with UCSF Chimera, which searches for protonatable standard residues (ASP, GLU, HIS) around the ligand, finds all possible combinations of protonation states, and writes a pdb file for each combination. Since the number of combinations of more than 4 protonatable residues becomes very large, the user can fix some residues to a give protonated/unprotonated state. See the examples below. You can also get the same info by typing `protonate_receptor.py -h`.

### EXAMPLE 1: list all protonatable residues within 8 Angstroms from the ligand.

`chimera --nogui --nostatus --script "$(which protonate_receptor.py) -rec 1a30_protein.pdb -lig 1a30_ligand.sdf -r 8.0 -list"`

### EXAMPLE 2: keep GLU_34.B, ASP_30.B, ASP_29.B, ASP_29.A, and ASP_30.A fixed, and create alternative protonations for all the rest (namely ASP_25.B and ASP_25.A).
`chimera --nogui --nostatus --script "$(which protonate_receptor.py) -rec 1a30_protein.pdb -lig 1a30_ligand.sdf -r 4.0 -fix GLU_34.B -fix ASP_30.B -fix ASP_29.B -fix ASP_29.A -fix ASP_30.A"`

### EXAMPLE 3: protonate all protein residues.
`chimera --nogui --nostatus --script "$(which protonate_receptor.py) -rec 1a30_protein.pdb -lig 1a30_ligand.sdf"`


