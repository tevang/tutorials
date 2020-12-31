#!/usr/bin/env python

__author__="Thomas Evangelidis"
__email__="tevang3@gmail.com"

SMILES="CCC1=CC=C(C=C1)C=C2C(=O)NC(=S)S2"   # 10058-F4 c-Myc inhibitor
SMILES='CC(C)c5cc(CNC[C@@H](O)[C@@H]4C[C@H](C)CCCCCN([C@H](C)c1ccccc1)C(=O)c2cc(cc(c2)C3=NC=CO3)C(=O)N4)ccc5'   # macrocycle

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from rdkit.Chem.Draw.IPythonConsole import *
from rdkit.Chem.Draw import MolToFile
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdmolfiles import MolFromMol2File
import os, shutil

## Parse command line arguments
def cmdlineparse():
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter, description="""
DESCRIPTION:

This is a Python script to create fragments as in ECFP fingerprints.

                         """,
                            epilog="""
### EXAMPLE 1:  

    """)
    parser.add_argument("-smiles", dest="SMILES", required=False, default=None,
                        help="The molecule to be fragmented in SMILES format.")
    parser.add_argument("-mol2", dest="MOL2", required=False, default=None,
                        help="The molecule to be fragmented in MOL2 format.")
    parser.add_argument("-outfolder", dest="OUT_FOLDER", required=False, default="fragments", type=str,
                        help="The folder name which will be created (or erased if it already exists) where the "
                             "PNG images of the fragments will be saved.")
    parser.add_argument("-fpradius", dest="FP_RADIUS", required=False, default=2, type=int,
                        help="The ECFP radius parameter value (distance in number of bonds). Default: %(default)s")

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = cmdlineparse()
    if args.SMILES:
        mol = Chem.MolFromSmiles(args.SMILES)
    if args.MOL2:
        mol = MolFromMol2File(args.MOL2, sanitize=False, removeHs=False)
    if os.path.exists(args.OUT_FOLDER):
        shutil.rmtree(args.OUT_FOLDER)
    os.mkdir(args.OUT_FOLDER)

    MolToFile(mol, "original_molecule.png")
    shutil.move("original_molecule.png", args.OUT_FOLDER + "/original_molecule.png")
    bi = {}
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=args.FP_RADIUS, bitInfo=bi)
    for k in bi.keys(): 
        mfp2_svg = DrawMorganBit(mol, k, bi)
        mfp2_svg.save(fp="%s/%i_frag.png" % (args.OUT_FOLDER, k), format="PNG")

    # TODO: show all fragments in one figure
    # https://stackoverflow.com/questions/37365824/pandas-ipython-notebook-include-and-display-an-image-in-a-dataframe
