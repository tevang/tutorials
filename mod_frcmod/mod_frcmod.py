#!/usr/bin/env python
import os, shutil
import traceback, sys, re
from subprocess import call

import pandas as pd
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from collections import defaultdict
from collections import OrderedDict
import parmed as pmd
from collections import defaultdict
import numpy as np

class tree(OrderedDict):
    def __missing__(self, key):
        self[key] = type(self)()
        return self[key]

## Parse command line arguments
def cmdlineparse():
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter, description="""
DESCRIPTION:
    This script reads an optimized ligand structure, measures the bond bonds, dihedrals and bond lengths and writes them to
    a provided frcmod file by replacing the existing ones.
    
    Some examples input files:
    * correct ligand frcmod file:
    /home2/thomas/Documents/Consensus_Scoring_Project/D3R_2018/BACE/MD_FEset/BACE_from_3dv5_apo/BACE68_frcmod.ligand
    * wrong ligand frcmod file:
    /home2/thomas/Documents/Consensus_Scoring_Project/D3R_2018/BACE/MD_FEset/BACE_from_3dv5_apo/BACE68_wrong/frcmod.ligand
    * optimized ligand geometry file:
    /home2/thomas/Documents/Consensus_Scoring_Project/D3R_2018/BACE/MD_FEset/BACE_from_3dv5_apo/ligands/bcc/bace68.bcc.mol2
    
    
    
                            """,
                            epilog="""
EXAMPLE:
        
        mod_frcmod.py -ligfile bace68.bcc.mol2 -frcmod BACE68_frcmod.ligand -ofrcmod BACE68_frcmod.ligand_corrected

    """)
    parser.add_argument("-ligfile", dest="LIGFILE", required=False, default=None, type=str,
                        help="sdf or mol2 file with optimized ligand structure from which to measure the equilibrium "
                             "bond angles, dihedrals and bond lengths.")
    parser.add_argument("-frcmod", dest="FRCMOD", required=True, default=None,
                        help="the frcmod parameter file of the ligand.")
    parser.add_argument("-ofrcmod", dest="OUT_FRCMOD", required=False, default=None,
                        help="the name of the modified frcmod parameter file of the ligand, namely the output.")
    parser.add_argument("-ff", dest="FF", required=False, default="gaff2",
                        help="the ligand force field.")
    parser.add_argument("-verbose", dest="VERBOSE", required=False, default=False, action='store_true',
                        help="Print more details.")

    args = parser.parse_args()
    return args


#################################################### FUNCTION DEFINITIONS ################################################


# THE FOLLOWING CODE IS USELESS SINCE PARMED CAN READ AND WRITE FRCMOD FILES
# ##~~~~~~~~~~~~~~~~~~~~`` DataFrames to store the force field parameters ``~~~~~~~~~~~~~~~~~~##
# mass_cols = ["KNDSYM", "AMASS", "ATPOL", "comment"]
# # NOTE: by defining the dtype you will be able to retrieve the value of a column by simply doing row[colname]
# mass_df = pd.DataFrame([], columns=mass_cols)
# mass_format = "%2s %-6.3f%13.3f\t%s\n"    # (A2,2X,F10.2x,f10.2)
# mass_pattern = "^([ 0-9a-z]{2})\s+([0-9.-]+)\s+([0-9.-]+)[\s$]+(.*)"  # 4 groups
#
# bond_cols = ["IBT", "JBT", "RK", "REQ", "comment"]
# bond_df = pd.DataFrame([], columns=bond_cols)
# bond_format = "%2s-%2s%8.2f%8.3f\t%s\n"    # A2,1X,A2,2F10.2
# bond_pattern = "^([ 0-9a-z]{2})-([ 0-9a-z]{2})\s+([0-9.-]+)\s+([0-9.-]+)[\s$]+(.*)"   # 5 groups
#
# angl_cols = ["ITT" , "JTT" , "KTT" , "TK" , "TEQ", "comment"]
# angl_df = pd.DataFrame([], columns=angl_cols)
# angl_format = "%2s-%2s-%2s%9.3f%12.3f\t%s\n"    # A2,1X,A2,1X,A2,2F10.2
# angl_pattern = "^([ 0-9a-z]{2})-([ 0-9a-z]{2})-([ 0-9a-z]{2})\s+([0-9.-]+)\s+([0-9.-]+)[\s$]+(.*)"   # 6 groups
#
# dihe_cols = ["IPT" , "JPT" , "KPT" , "LPT" , "IDIVF" , "PK" , "PHASE" , "PN", "comment"]
# dihe_df = pd.DataFrame([], columns=dihe_cols)
# dihe_format = "%2s-%2s-%2s-%2s%4i%9.3f%14.3f%16.3f\t%s\n"    # A2,1X,A2,1X,A2,1X,A2,I4,3F15.2
# dihe_pattern = "^([ 0-9a-z]{2})-([ 0-9a-z]{2})-([ 0-9a-z]{2})-([ 0-9a-z]{2})\s+([0-9]+)\s+([0-9.-]+)\s+([0-9.-]+)\s+([0-9.-]+)[\s$]+(.*)"   # 9 groups
#
# impr_cols = ["IPT" , "JPT" , "KPT" , "LPT"  , "PK" , "PHASE" , "PN", "comment"]
# impr_df = pd.DataFrame([], columns=impr_cols)
# impr_format = "%2s-%2s-%2s-%2s%12.1f%15.1f%12.1f\t%s\n"    # A2,1X,A2,1X,A2,1X,A2,I4,3F15.2
# impr_pattern = "^([ 0-9a-z]{2})-([ 0-9a-z]{2})-([ 0-9a-z]{2})-([ 0-9a-z]{2})\s+([0-9.-]+)\s+([0-9.-]+)\s+([0-9.-]+)[\s$]+(.*)"   # 8 groups
#
# # H-BOND 10-12 POTENTIAL PARAMETERS
# hbon_cols = ["KT1" , "KT2" , "A" , "B", "comment"]
# hbon_df = pd.DataFrame([], columns=hbon_cols)
# hbon_df = hbon_df.astype({"KT1": 'str', "KT2": 'str', "A": 'float', "B": 'float', "comment": 'str'})
# hbon_format = ""    # 2X,A2,2X,A2,2x,5F10.2,I2
#
# # ONLY IF KINDNB .EQ. 'RE' ???
# nonb_cols = ["LTYNB" , "R" , "EDEP", "comment"]
# nonb_df = pd.DataFrame([], columns=nonb_cols)
# nonb_format = "%4s%16.4f%8.4f\t%s\n"    # A2,1X,A2,1X,A2,1X,A2,I4,3F15.2
# nonb_pattern = "^\s*([ 0-9a-z]{2})\s+([0-9.-]+)\s+([0-9.-]+)[\s$]+(.*)"  # 4 groups
#
# # Put all dataframes together into a dict
# forcefield = {"MASS": mass_df,
#               "BOND": bond_df,
#               "ANGLE": angl_df,
#               "DIHE": dihe_df,
#               "IMPROPER": impr_df,
#               "NONBON": nonb_df}
# fields = ["MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"]
# columns = [mass_cols, bond_cols, angl_cols, dihe_cols, impr_cols, nonb_cols]
# patterns = [mass_pattern, bond_pattern, angl_pattern, dihe_pattern, impr_pattern, nonb_pattern]
# formats = [mass_format, bond_format, angl_format, dihe_format, impr_format, nonb_format]
#
# def update_forcefield_dtypes():
#
#     global forcefield
#
#     forcefield["MASS"] = forcefield["MASS"].astype({'KNDSYM': 'str', 'AMASS': 'float', 'ATPOL': 'float', "comment": 'str'})
#     forcefield["BOND"] = forcefield["BOND"].astype({"IBT": 'str', "JBT": 'str', "RK": 'float', "REQ": 'float', "comment": 'str'})
#     forcefield["ANGLE"] = forcefield["ANGLE"].astype({"ITT": 'str', "JTT": 'str', "KTT": 'str', "TK": 'float', "TEQ": 'float', "comment": 'str'})
#     forcefield["DIHE"] = forcefield["DIHE"].astype({"IPT": str, "JPT": 'str', "KPT": 'str', "LPT": 'str', "IDIVF": 'int', "PK": 'float',
#                           "PHASE": 'float', "PN": 'float', "comment": 'str'})
#     forcefield["IMPROPER"] = forcefield["IMPROPER"].astype({"IPT": 'str', "JPT": 'str', "KPT": 'str', "LPT": 'str',
#                           "PK": 'float', "PHASE": 'float', "PN": 'float', "comment": 'str'})
#     forcefield["NONBON"] = forcefield["NONBON"].astype({"LTYNB": 'str', "R": 'float', "EDEP": 'float', "comment": 'str'})
#
# def load_frcmod(fname):
#     """
#     For the format of frcmod file look at:
#     http://ambermd.org/formats.html#frcmod
#
#     :param fname:
#     :return:
#     """
#     global forcefield
#     with open(fname, 'r') as f:
#         contents = f.readlines()
#
#     starts = [contents.index(f+"\n") for f in fields]
#     ends = [s-1 for s in starts[1:]]
#     ends.append(len(contents)-1)
#     for i in range(len(fields)):
#         field, cols, start, end, pattern = fields[i], columns[i], starts[i], ends[i], patterns[i]
#         for line in contents[start+1:end+1]:
#             m = re.search(pattern, line)
#             if not m:
#                 continue
#             row_dict = {c:w for w,c in zip(m.groups(), cols)}
#             forcefield[field] = forcefield[field].append(row_dict, ignore_index=True)  # save this line to the dataframe
#     update_forcefield_dtypes()
#
# def write_frcmod(outfname):
#
#     global forcefield, args
#
#     out = open(outfname, 'w')
#     for i in range(len(fields)):
#         field, format = fields[i], formats[i]
#         out.write(field + "\n")
#         for i, row in forcefield[field].iterrows():
#             out.write(format % tuple(row.values))
#         out.write("\n")


def run_commandline(commandline, logname="log", append=False, return_out=False, error_keywords=[], skip_fail=False,
                    verbose=True):
    """
        FUNCTION to run a single command on the UNIX shell. The worker will only receive an index from network.
    """
    if append:
        fout = open(logname, 'a')
    else:
        fout = open(logname, 'w')
    if verbose:
        print("Running commandline:", commandline)
    return_code = call(commandline, stdout=fout, stderr=fout, shell=True, executable='/bin/bash')

    if (return_code != 0):
        print("ERROR, THE FOLLOWING COMMAND FAILED TO RUN:", "FAIL")
        print(commandline)
        print("return_code=", return_code)
        fout.close()
        print("Output:")
        with open(logname, 'r') as f:
            contents = f.readlines()
            for line in contents:
                print(line)
        if not skip_fail:
            raise Exception()
    fout.close()

    if len(error_keywords) > 0:
        with open(logname, 'r') as f:
            contents = f.readlines()
            for line in contents:
                for word in error_keywords:
                    if word in line:
                        print("ERROR, THE FOLLOWING COMMAND FAILED TO RUN:")
                        print(commandline)
                        print("COMMAND OUTPUT:")
                        for line in contents:
                            print(line)
                        raise Exception()

    if return_out:
        with open(logname, 'r') as f:
            contents = f.readlines()
            return contents

def create_prmtop(frcmod, ligfile):

    if os.path.exists("tmp/"):
        shutil.rmtree("tmp/")
    os.mkdir("tmp/")
    run_commandline("ln -s %s %s/frcmod.ligand" % (os.path.abspath(frcmod), os.path.abspath("tmp/")))

    # convert with antechamber to mol2 with GAFF2 atom types
    # NOTE: -at gaff2 writes some unknown atom types that are not in the frcmod file (e.g. nh->nu, n->ns, n3->n7).
    run_commandline("antechamber -i %s -fi %s -o tmp/ligand.gaff2.mol2 -fo mol2 -rn LIG -at gaff2 -dr n"
                        % (ligfile, ligfile.split('.')[-1]))

    ligand_leap = """
source leaprc.gaff2
loadAmberParams tmp/frcmod.ligand
LIG = loadMol2 tmp/ligand.gaff2.mol2
saveAmberParm LIG tmp/ligand.prmtop tmp/ligand.inpcrd
quit
    """

    with open("tmp/ligand_leap.in", 'w') as f:
        f.write(ligand_leap)
    leap_out = run_commandline("tleap -s -f tmp/ligand_leap.in", return_out=True, error_keywords=['FATAL:'])


def write_corrected_frcmod(ligfile, frcmod, out_frcmod, verbose=False):
    """
    This method takes the equilibrium bond lengths and angles from the ligfile and writes a new
    frcmod file with corrected GAFF2 ligand parameters for MD.

    :param ligfile: mol2 or sdf file with optimized ligand geometry from where to copy bond lengths and angles.
    :param frcmod: the frcmod file that needs corrections.
    :param out_frcmod: the name of the output frcmod file that carries the corrections.
    :return:
    """
    global args

    # create the prmtop and inpcrd file within a 'tmp/' folder
    create_prmtop(frcmod, ligfile)

    # load them to PARMED
    mol = pmd.load_file("tmp/ligand.prmtop", xyz="tmp/ligand.inpcrd", structure=True)
    bond_dict = defaultdict(list)
    for bond in mol.bonds:
        # print("%s-%s XXX %f" % (bond.atom1.type, bond.atom2.type, bond.measure()))
        bond_dict["%s-%s" % (bond.atom1.type, bond.atom2.type)].append(bond.measure())
        bond_dict["%s-%s" % (bond.atom2.type, bond.atom1.type)].append(bond.measure())    # add the reverse bond, too

    if verbose:
        print("\nBond = mean value += stdev, min-max")
        for bondname, distlist in bond_dict.items():
            print("%s = %f +- %f, %f" % (bondname, np.mean(distlist), np.std(distlist), np.ptp(distlist)))

    angle_dict = defaultdict(list)
    for angle in mol.angles:
        angle_dict["%s-%s-%s" % (angle.atom1.type, angle.atom2.type, angle.atom3.type)].append(angle.measure())
        angle_dict["%s-%s-%s" % (angle.atom3.type, angle.atom2.type, angle.atom1.type)].append(angle.measure())   # add the reverse angle, too

    if verbose:
        print("\nAngle = mean value += stdev, min-max")
        for anglename, anglelist in angle_dict.items():
            print("%s = %f +- %f, %f" % (anglename, np.mean(anglelist), np.std(anglelist), np.ptp(anglelist)))

    # par = pmd.load_file(frcmod)

    for bond in mol.bonds:
        bondname = "%s-%s" % (bond.atom1.type, bond.atom2.type)
        assert bondname in bond_dict.keys(), "ERROR: bond %s does not exist in the mol2 file with " \
                                               "the optimized geometry!" % bondname
        idx = bond.type.idx
        bond.type.req = round(np.mean(bond_dict[bondname]), 3)    # replace with the mean bond value
        mol.bond_types[idx].req = round(np.mean(bond_dict[bondname]), 3)    # replace with the mean bond value

    for angle in mol.angles:
        anglename = "%s-%s-%s" % (angle.atom1.type, angle.atom2.type, angle.atom3.type)
        assert anglename in angle_dict.keys(), "ERROR: angle %s does not exist in the mol2 file with " \
                                               "the optimized geometry!" % anglename
        idx = angle.type.idx
        angle.type.theteq = round(np.mean(angle_dict[anglename]), 3)    # replace with the mean angle value
        mol.angle_types[idx].theteq = round(np.mean(angle_dict[anglename]), 3)    # replace with the mean angle value

    # par.write('edited_'+frcmod, title="Created by mod_frcmod.py script.", style='frcmod')
    pmd.tools.writeFrcmod(mol, out_frcmod).execute()

    # clean intermediate files
    shutil.rmtree("tmp/")


################################################### END OF FUNCTION DEFINITIONS ##########################################

if __name__ == "__main__":

    try:
        args = cmdlineparse()
        if args.OUT_FRCMOD == None:
            args.OUT_FRCMOD = "mod_%s" % args.FRCMOD
        write_corrected_frcmod(args.LIGFILE, args.FRCMOD, args.OUT_FRCMOD)

    except:
        type, value, tb = sys.exc_info()
        lines = traceback.format_exception(type, value, tb)
        print(''.join(lines))
        raise