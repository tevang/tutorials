#!/usr/bin/env python

__author__="Thomas Evangelidis"
__email__="tevang3@gmail.com"



import os
from argparse import ArgumentParser, RawDescriptionHelpFormatter


## Parse command line arguments
def cmdlineparse():
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter, description="""
DESCRIPTION:

This is a Python script to prepare the receptor-ligand complex for scoring. If you encounter problems with the input pdb file then try correcting it using:
1) pdb4amber from AmberTools (https://github.com/Amber-MD/pdb4amber)
2) pdbfixer (https://github.com/pandegroup/pdbfixer)

TODO: add optional support for LYS and CYS protonated forms.
https://www.cgl.ucsf.edu/chimera/docs/ContributedSoftware/addh/addh.html

                         """,
                            epilog="""
### EXAMPLE 1: 
pychimera $(which dockprep.py) -complex 3K5C-BACE_150_complex.pdb -cmethod gas -neut


    """)
    parser.add_argument("-complex", dest="COMPLEX", required=False, default=None, type=str,
                        help="pdb file with the holo form of the receptor.")
    parser.add_argument("-cmethod", dest="CHARGE_METHOD", required=False, default='gas', type=str, choices=['gas', 'am1'],
                        help="Method to calculate charges fo the ligand. Default: %(default)s")
    parser.add_argument("-neut", dest="NEUTRALIZE", required=False, default=False, action='store_true',
                        help="Neutralize the system by adding coutner ions.")
    parser.add_argument("-rec", dest="RECEPTOR", required=False, default=None, type=str,
                        help="Instead of -complex give the pdb file with the apo form of the receptor.")
    parser.add_argument("-lig", dest="LIGAND", required=False, default=None, type=str,
                        help="Instead of -complex give an sdf or mol2 file with optimized ligand structure from which to find the "
                             "binding site residues.")
    args = parser.parse_args()
    return args


if __name__ == "__main__":

    args = cmdlineparse()
    import chimera
    from chimera import runCommand as rc
    from Addions import initiateAddions
    from DockPrep import prep
    from AddCharge import estimateFormalCharge
    from chimera.selection import currentResidues

    if args.COMPLEX:
        rc("open %s" % args.COMPLEX)  # load the protein-ligand complex
        rc("split #0 ligands")
        rc("sel #0.2")  # select the ligand
        ligres = currentResidues()[0]
        ligres.type = 'LIG'  # change the resname of the ligand to 'LIG'
        rc("combine #0.1 modelId 1")  # create a new molecule containing just the receptor
        rc("combine #0.2 modelId 2")  # create a new molecule containing just the ligand
        # Now that we calculated the charges of the protein and the ligand, we just need the complex
        rc("combine #1,2 modelId 3")  # create a new molecule containing the protein-ligand complex
        rc("del #0-2")
        pdb = args.COMPLEX.replace(".pdb", "_prep.pdb")
    elif args.RECEPTOR and args.LIGAND:
        rc("open %s" % args.RECEPTOR)  # load the receptor
        rc("open %s" % args.LIGAND)  # load the ligand
        rc("sel #1")  # select the ligand
        ligres = currentResidues()[0]
        ligres.type = 'LIG'  # change the resname of the ligand to 'LIG'
        rc("combine #0,1 modelId 2")  # create a new molecule containing the protein-ligand complex
        rc("combine #2 modelId 3")  # create a new molecule containing the protein-ligand complex
        rc("del #0-2")
        pdb = os.path.splitext(os.path.basename(args.RECEPTOR))[0] + "_" + os.path.splitext(os.path.basename(args.LIGAND))[0] + "_prep.pdb"

    print "Preparing receptor for docking and calculating ligand '%s' charges (may be slow)." % args.CHARGE_METHOD
    models = chimera.openModels.list(modelTypes=[chimera.Molecule])
    prep(models, nogui=True, method=args.CHARGE_METHOD)
    net_charge = estimateFormalCharge(models[0].atoms)
    # Neutralize system
    if args.NEUTRALIZE:
        if net_charge < 0:
            initiateAddions(models, "Na+", "neutralize", chimera.replyobj.status)
        elif net_charge > 0:
            initiateAddions(models, "Cl-", "neutralize", chimera.replyobj.status)
        if net_charge != 0:
            # change the resids of the ions, which by default they are all 1
            rc("sel ~ions")
            existing_resids = [int(str(r.id).split('.')[0]) for r in currentResidues()]
            start = max(existing_resids) + 2
            rc("resrenumber %i ions" % start)   # renumber the resids of the added ions
    # change the resid of the ligand
    rc('sel #3 & ~ #3:LIG')
    existing_resids = [int(str(r.id).split('.')[0]) for r in currentResidues()]
    start = max(existing_resids) + 2
    rc("resrenumber %i #3:LIG" % start)
    rc("combine #3 modelId 4")  # create a new molecule to split it into receptor and ligand
    rc("split #4 atoms ~#4:LIG")
    rc("combine #4.1 modelId 5")  # create a new molecule containing just the receptor
    rc("combine #4.2 modelId 6")  # create a new molecule containing just the ligand
    models = chimera.openModels.list(modelTypes=[chimera.Molecule])
    # for m in models: print len(m.atoms), estimateFormalCharge(m.atoms)    # DEBUGGING
    rec_charge = estimateFormalCharge(models[3].atoms)
    lig_charge = estimateFormalCharge(models[4].atoms)
    rc("del #4-6")

    # Finally, write the complex pdb file with headers
    rc("write format pdb #3 %s" % pdb)
    with open(pdb, "r+") as f:
        s = f.read()
        f.seek(0)
        f.write("# receptor net charge = %i\n# ligand net charge = %i\n" % (rec_charge, lig_charge))  # after system neutralization
        f.write(s)
