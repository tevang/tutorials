#!/usr/bin/env python

__author__="Thomas Evangelidis"
__email__="tevang3@gmail.com"



import os, sys, traceback
import random
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from collections import defaultdict
from operator import itemgetter


def cmdlineparse():
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter, description="""
DESCRIPTION:

This is a Python script to prepare the receptor-ligand complex for scoring. If you encounter problems with the input pdb file 
then try correcting it using:
1) pdb4amber from AmberTools (https://github.com/Amber-MD/pdb4amber)
2) pdbfixer (https://github.com/pandegroup/pdbfixer)

TODO: add optional support for LYS and CYS protonated forms.
https://www.cgl.ucsf.edu/chimera/docs/ContributedSoftware/addh/addh.html

                         """,
                            epilog="""
### EXAMPLE 1: 
pychimera $(which dockprep.py) -complex 3K5C-BACE_150_complex.pdb -cmethod gas -neut -lignetcharge -2


    """)
    parser.add_argument("-complex", dest="COMPLEX", required=False, default=None, type=str,
                        help="pdb file with the holo form of the receptor.")
    parser.add_argument("--charge-method", dest="CHARGE_METHOD", required=False, default='gas', type=str, choices=['gas', 'am1'],
                        help="Method to calculate charges fo the ligand. Default: %(default)s")
    parser.add_argument("-neut", dest="NEUTRALIZE", required=False, default=False, action='store_true',
                        help="Neutralize the system by adding counter ions.")
    parser.add_argument("-stripions", dest="STRIP_IONS", required=False, default=False, action='store_true',
                        help="Strip out all ions.")
    parser.add_argument("-keepchains", dest="KEEP_CHAINIDS", required=False, default=False, action='store_true',
                        help="Keep the original chain IDs. Default is False, ligand and protein will be chain A for homology modeling.")
    parser.add_argument("-keephydrogens", dest="KEEP_PROTEIN_HYDROGENS", required=False, default=False, action='store_true',
                        help="Keep the protein's hydrogens (default is to strip them).")
    parser.add_argument("-rec", dest="RECEPTOR", required=False, default=None, type=str,
                        help="Instead of -complex give the pdb file with the apo form of the receptor.")
    parser.add_argument("-lig", dest="LIGAND", required=False, default=None, type=str,
                        help="Instead of -complex give an sdf or mol2 file with optimized ligand structure from which to find the "
                             "binding site residues.")
    parser.add_argument("-o", dest="OUT_PDB", required=False, default=None, type=str,
                        help="Output PDB file name of the prepared protein-ligand complex.")
    parser.add_argument("-lignetcharge", dest="LIG_NET_CHARGE", required=False, default=None, type=int,
                        help="Optionaly (but RECOMMENDED) give the net charge of the ligand, otherwise it will be estimated by Chimera.")
    args = parser.parse_args()
    return args


################################################### FUNCTION DEFINITIONS #########################################

def standardize_terminal_protein_residues(receptor_pdb, molID="#0"):
    """
    To prevent errors like 'ValueError: Cannot determine GAFF type for :11.A@HD14 (etc.)' raised by
    initiateAddions(), originating from termini capped by Shrodinger's Maestro or incomplete or missing protein
    residues, mutate the terminal residues to their original aa type (applies to every protein chain).
    """
    print("Standardizing protein's terminal residues.")
    rc("sel %s & protein" % molID)
    chaindID_resids = defaultdict(list)
    for r in currentResidues():
        chaindID_resids[r.id.chainId].append((r.id.position, r.type))
    for chainID in chaindID_resids.keys():
        chaindID_resids[chainID].sort(key=itemgetter(0))
        # NOTE: the side-chain mutations (swapaa) were not necessary in the proteins tested so far.
        # rc("swapaa %s %s:%i.%s" %
        #    (chaindID_resids[chainID][0][1], molID, chaindID_resids[chainID][0][0], str(chainID)))     # N-term
        rc("del element.H & %s:%i.%s" %
           (molID, chaindID_resids[chainID][0][0], str(chainID)))
        # rc("swapaa %s %s:%i.%s" %
        #    (chaindID_resids[chainID][-1][1], molID, chaindID_resids[chainID][-1][0], str(chainID)))     # C-term
        rc("del element.H & %s:%i.%s" %
           (molID, chaindID_resids[chainID][-1][0], str(chainID)))
    # Only if you write and load the PDB then Chimera will reset the valence of the N-terminal amide and
    # thus will not again the H1,H2,H3 which cause the error.
    rc("write format pdb #0 " + receptor_pdb.replace(".pdb", "_tmp.pdb"))
    rc("del #0")
    rc("open %s %s" % (molID, receptor_pdb.replace(".pdb", "_tmp.pdb")))
    os.remove(receptor_pdb.replace(".pdb", "_tmp.pdb"))

##################################################################################################################

if __name__ == "__main__":

    try:

        args = cmdlineparse()
        import chimera
        from chimera import runCommand as rc
        from Addions import initiateAddions
        from DockPrep import prep
        import AddH
        from AddCharge import estimateFormalCharge, addNonstandardResCharges
        from chimera.selection import currentResidues

        if args.COMPLEX:
            rc("open %s" % args.COMPLEX)  # load the protein-ligand complex
            if args.KEEP_PROTEIN_HYDROGENS:
                rc("delete element.H")
            standardize_terminal_protein_residues("#0")  # TODO: UNTESTED
            if args.STRIP_IONS:
                rc("delete ions")
            rc("split #0 ligands")
            rc("sel #0.2")  # select the ligand
            ligres = currentResidues()[0]
            ligres.type = 'LIG'  # change the resname of the ligand to 'LIG'
            rc("combine #0.1 modelId 1")  # create a new molecule containing just the receptor
            rc("combine #0.2 modelId 2")  # create a new molecule containing just the ligand
            rc("del #0")
            # We will estimate the receptor's net charge. For this we need to DockPrep the receptor (is fast).
            models = chimera.openModels.list(modelTypes=[chimera.Molecule])
            # For a full list of DockPrep options, look into file Chimera-alpha_py2.7/share/DockPrep/__init__.py
            prep([models[0]], nogui=True, method=args.CHARGE_METHOD, addHFunc=AddH.simpleAddHydrogens)
            rec_charge = estimateFormalCharge(models[0].atoms)  # DockPred does not assign charges to receptor atoms, only to ligand atoms
            # Now that we calculated the charges of the protein and the ligand, we just need the complex
            rc("combine #1,2 modelId 3")  # create a new molecule containing the protein-ligand complex
            rc("del #1-2")
            pdb = args.COMPLEX.replace(".pdb", "_prep.pdb")
        elif args.RECEPTOR and args.LIGAND:
            rc("open %s" % args.RECEPTOR)  # load the receptor
            if args.KEEP_PROTEIN_HYDROGENS:
                rc("delete element.H")
            standardize_terminal_protein_residues("#0")  # read function's definition to understand why is here
            rc("open %s" % args.LIGAND)  # load the ligand
            if args.STRIP_IONS:
                rc("delete ions")
            # We will estimate the receptor's net charge. For this we need to DockPrep the receptor (is fast).
            models = chimera.openModels.list(modelTypes=[chimera.Molecule])
            # For a full list of DockPrep options, look into file Chimera-alpha_py2.7/share/DockPrep/__init__.py
            prep([models[0]], nogui=True, method=args.CHARGE_METHOD, addHFunc=AddH.simpleAddHydrogens)
            rec_charge = estimateFormalCharge(models[0].atoms)  # DockPred does not assign charges to receptor atoms, only to ligand atoms
            rc("sel #1")  # select the ligand
            ligres = currentResidues()[0]
            ligres.type = 'LIG'  # change the resname of the ligand to 'LIG'
            rc("combine #0,1 modelId 2")  # create a new molecule containing the protein-ligand complex
            rc("combine #2 modelId 3")  # create a new molecule containing the protein-ligand complex
            rc("del #0-2")
            pdb = os.path.splitext(os.path.basename(args.RECEPTOR))[0] + "_" + os.path.splitext(os.path.basename(args.LIGAND))[0] + "_prep.pdb"
        elif args.RECEPTOR:
            rc("open %s" % args.RECEPTOR)  # load the receptor
            if args.KEEP_PROTEIN_HYDROGENS:
                rc("delete element.H")
            standardize_terminal_protein_residues("#0")  # TODO: UNTESTED
            if args.STRIP_IONS:
                rc("delete ions")
            # We will estimate the receptor's net charge. For this we need to DockPrep the receptor (is fast).
            models = chimera.openModels.list(modelTypes=[chimera.Molecule])
            # For a full list of DockPrep options, look into file Chimera-alpha_py2.7/share/DockPrep/__init__.py
            prep([models[0]], nogui=True, method=args.CHARGE_METHOD, addHFunc=AddH.simpleAddHydrogens)
            rec_charge = estimateFormalCharge(models[0].atoms)  # DockPred does not assign charges to receptor atoms, only to ligand atoms
            pdb = os.path.splitext(os.path.basename(args.RECEPTOR))[0] + "_prep.pdb"

        print("Preparing receptor for docking and calculating ligand '%s' charges (may be slow)." % args.CHARGE_METHOD)
        models = chimera.openModels.list(modelTypes=[chimera.Molecule]) # actually only one model is left
        # For a full list of DockPrep options, look into file Chimera-alpha_py2.7/share/DockPrep/__init__.py
        prep(models, nogui=True, method=args.CHARGE_METHOD, addHFunc=AddH.simpleAddHydrogens)
        # NOTE: the default option addHFunc=AddH.hbondAddHydrogens raised an Error in Carbonic Unhydrase with the Zn+2 ion.
        if args.LIGAND != None and args.LIG_NET_CHARGE != None:
            net_charge = args.LIG_NET_CHARGE + rec_charge
        elif args.LIGAND != None and args.LIG_NET_CHARGE == None:
            net_charge = estimateFormalCharge(models[0].atoms)
        elif args.LIGAND == None and args.LIG_NET_CHARGE == None:
            net_charge = estimateFormalCharge(models[0].atoms)
        # Neutralize system
        if args.NEUTRALIZE:
            if net_charge < 0:
                initiateAddions(models, "Na+", str(abs(net_charge)), chimera.replyobj.status)
            elif net_charge > 0:
                initiateAddions(models, "Cl-", str(net_charge), chimera.replyobj.status)

            if net_charge != 0:
                # change the resids of the ions, which by default they are all 1
                rc("sel ~ions")
                existing_resids = [r.id.position for r in currentResidues()]
                start = max(existing_resids) + 2
                rc("resrenumber %i ions" % start)   # renumber the resids of the added ions

        if args.COMPLEX or args.LIGAND:
            # change the resid of the ligand
            rc('sel #3 & ~ #3:LIG')
            existing_resids = [r.id.position for r in currentResidues()]
            start = max(existing_resids) + 2
            rc("resrenumber %i #3:LIG" % start)
            rc("combine #3 modelId 4")  # create a new molecule to split it into receptor and ligand
            rc("split #4 atoms ~#4:LIG")
            rc("combine #4.1 modelId 5")  # create a new molecule containing just the receptor
            rc("combine #4.2 modelId 6")  # create a new molecule containing just the ligand
            models = chimera.openModels.list(modelTypes=[chimera.Molecule])
            # for m in models: print(len(m.atoms), estimateFormalCharge(m.atoms)    # DEBUGGING
            if args.LIG_NET_CHARGE:
                lig_charge = args.LIG_NET_CHARGE
            else:
                lig_charge = estimateFormalCharge(models[3].atoms)
            rc("del #4-6")

        # Finally, write the complex pdb file with headers
        if args.KEEP_CHAINIDS == False:
            rc("changechains B A all")  # <== OPTIONAL (ligand and protein will be chain A for homology modeling)
        if args.OUT_PDB:
            pdb = args.OUT_PDB
        if args.COMPLEX or args.LIGAND:
            rc("write format pdb #3 %s" % pdb)
        else:
            rc("write format pdb #0 %s" % pdb)
        with open(pdb, "r+") as f:
            s = f.read()
            f.seek(0)
            if args.COMPLEX or args.LIGAND:
                f.write("# receptor net charge = %i\n# ligand net charge = %i\n" % (rec_charge, lig_charge))  # after system neutralization
            else:
                f.write("# receptor net charge = %i\n" % (rec_charge))  # after system neutralization
            f.write(s)

    except:
        type, value, tb = sys.exc_info()
        lines = traceback.format_exception(type, value, tb)
        print((''.join(lines)))
        raise