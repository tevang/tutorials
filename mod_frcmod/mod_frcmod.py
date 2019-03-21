#!/usr/bin/env python

import traceback, sys, re
import pandas as pd
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from collections import defaultdict
from collections import OrderedDict

class tree(OrderedDict):
    def __missing__(self, key):
        self[key] = type(self)()
        return self[key]

## Parse command line arguments
def cmdlineparse():
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter, description="""
DESCRIPTION:
    This script reads an optimized ligand structure, measures the bond angles, dihedrals and bond lengths and writes them to
    a provided frcmod file by replacing the existing ones.
                            """,
                            epilog="""
EXAMPLE:


    """)
    parser.add_argument("-ligfile", dest="LIGAND", required=False, default=None, type=str,
                        help="sdf or mol2 file with optimized ligand structure from which to measure the equilibrium "
                             "bond angles, dihedrals and bond lengths.")
    parser.add_argument("-frcmod", dest="FRCMOD", required=True, default=None,
                        help="the frcmod parameter file of the ligand.")
    parser.add_argument("-ofrcmod", dest="OUT_FRCMOD", required=False, default=None,
                        help="the name of the modified frcmod parameter file of the ligand, namely the output.")
    parser.add_argument("-ff", dest="FF", required=False, default="gaff2",
                        help="the ligand force field.")

    args = parser.parse_args()
    return args


#################################################### FUNCTION DEFINITIONS ################################################

##~~~~~~~~~~~~~~~~~~~~`` DataFrames to store the force field parameters ``~~~~~~~~~~~~~~~~~~##
mass_cols = ["KNDSYM", "AMASS", "ATPOL", "comment"]
# NOTE: by defining the dtype you will be able to retrieve the value of a column by simply doing row[colname]
mass_df = pd.DataFrame([], columns=mass_cols)
mass_format = "%2s %-6.3f%13.3f\t%s\n"    # (A2,2X,F10.2x,f10.2)
mass_pattern = "^([ 0-9a-z]{2})\s+([0-9.-]+)\s+([0-9.-]+)[\s$]+(.*)"  # 4 groups

bond_cols = ["IBT", "JBT", "RK", "REQ", "comment"]
bond_df = pd.DataFrame([], columns=bond_cols)
bond_format = "%2s-%2s%8.2f%8.3f\t%s\n"    # A2,1X,A2,2F10.2
bond_pattern = "^([ 0-9a-z]{2})-([ 0-9a-z]{2})\s+([0-9.-]+)\s+([0-9.-]+)[\s$]+(.*)"   # 5 groups

angl_cols = ["ITT" , "JTT" , "KTT" , "TK" , "TEQ", "comment"]
angl_df = pd.DataFrame([], columns=angl_cols)
angl_format = "%2s-%2s-%2s%9.3f%12.3f\t%s\n"    # A2,1X,A2,1X,A2,2F10.2
angl_pattern = "^([ 0-9a-z]{2})-([ 0-9a-z]{2})-([ 0-9a-z]{2})\s+([0-9.-]+)\s+([0-9.-]+)[\s$]+(.*)"   # 6 groups

dihe_cols = ["IPT" , "JPT" , "KPT" , "LPT" , "IDIVF" , "PK" , "PHASE" , "PN", "comment"]
dihe_df = pd.DataFrame([], columns=dihe_cols)
dihe_format = "%2s-%2s-%2s-%2s%4i%9.3f%14.3f%16.3f\t%s\n"    # A2,1X,A2,1X,A2,1X,A2,I4,3F15.2
dihe_pattern = "^([ 0-9a-z]{2})-([ 0-9a-z]{2})-([ 0-9a-z]{2})-([ 0-9a-z]{2})\s+([0-9]+)\s+([0-9.-]+)\s+([0-9.-]+)\s+([0-9.-]+)[\s$]+(.*)"   # 9 groups

impr_cols = ["IPT" , "JPT" , "KPT" , "LPT"  , "PK" , "PHASE" , "PN", "comment"]
impr_df = pd.DataFrame([], columns=impr_cols)
impr_format = "%2s-%2s-%2s-%2s%12.1f%15.1f%12.1f\t%s\n"    # A2,1X,A2,1X,A2,1X,A2,I4,3F15.2
impr_pattern = "^([ 0-9a-z]{2})-([ 0-9a-z]{2})-([ 0-9a-z]{2})-([ 0-9a-z]{2})\s+([0-9.-]+)\s+([0-9.-]+)\s+([0-9.-]+)[\s$]+(.*)"   # 8 groups

# H-BOND 10-12 POTENTIAL PARAMETERS
hbon_cols = ["KT1" , "KT2" , "A" , "B", "comment"]
hbon_df = pd.DataFrame([], columns=hbon_cols)
hbon_df = hbon_df.astype({"KT1": 'str', "KT2": 'str', "A": 'float', "B": 'float', "comment": 'str'})
hbon_format = ""    # 2X,A2,2X,A2,2x,5F10.2,I2

# ONLY IF KINDNB .EQ. 'RE' ???
nonb_cols = ["LTYNB" , "R" , "EDEP", "comment"]
nonb_df = pd.DataFrame([], columns=nonb_cols)
nonb_format = "%4s%16.4f%8.4f\t%s\n"    # A2,1X,A2,1X,A2,1X,A2,I4,3F15.2
nonb_pattern = "^\s*([ 0-9a-z]{2})\s+([0-9.-]+)\s+([0-9.-]+)[\s$]+(.*)"  # 4 groups

# Put all dataframes together into a dict
forcefield = {"MASS": mass_df,
              "BOND": bond_df,
              "ANGLE": angl_df,
              "DIHE": dihe_df,
              "IMPROPER": impr_df,
              "NONBON": nonb_df}
fields = ["MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"]
columns = [mass_cols, bond_cols, angl_cols, dihe_cols, impr_cols, nonb_cols]
patterns = [mass_pattern, bond_pattern, angl_pattern, dihe_pattern, impr_pattern, nonb_pattern]
formats = [mass_format, bond_format, angl_format, dihe_format, impr_format, nonb_format]

def update_forcefield_dtypes():

    global forcefield

    forcefield["MASS"] = forcefield["MASS"].astype({'KNDSYM': 'str', 'AMASS': 'float', 'ATPOL': 'float', "comment": 'str'})
    forcefield["BOND"] = forcefield["BOND"].astype({"IBT": 'str', "JBT": 'str', "RK": 'float', "REQ": 'float', "comment": 'str'})
    forcefield["ANGLE"] = forcefield["ANGLE"].astype({"ITT": 'str', "JTT": 'str', "KTT": 'str', "TK": 'float', "TEQ": 'float', "comment": 'str'})
    forcefield["DIHE"] = forcefield["DIHE"].astype({"IPT": str, "JPT": 'str', "KPT": 'str', "LPT": 'str', "IDIVF": 'int', "PK": 'float',
                          "PHASE": 'float', "PN": 'float', "comment": 'str'})
    forcefield["IMPROPER"] = forcefield["IMPROPER"].astype({"IPT": 'str', "JPT": 'str', "KPT": 'str', "LPT": 'str',
                          "PK": 'float', "PHASE": 'float', "PN": 'float', "comment": 'str'})
    forcefield["NONBON"] = forcefield["NONBON"].astype({"LTYNB": 'str', "R": 'float', "EDEP": 'float', "comment": 'str'})

def load_frcmod(fname):
    """
    For the format of frcmod file look at:
    http://ambermd.org/formats.html#frcmod

    :param fname:
    :return:
    """
    global forcefield
    with open(fname, 'r') as f:
        contents = f.readlines()

    starts = [contents.index(f+"\n") for f in fields]
    ends = [s-1 for s in starts[1:]]
    ends.append(len(contents)-1)
    for i in range(len(fields)):
        field, cols, start, end, pattern = fields[i], columns[i], starts[i], ends[i], patterns[i]
        for line in contents[start+1:end+1]:
            m = re.search(pattern, line)
            if not m:
                continue
            row_dict = {c:w for w,c in zip(m.groups(), cols)}
            forcefield[field] = forcefield[field].append(row_dict, ignore_index=True)  # save this line to the dataframe
    update_forcefield_dtypes()

def write_frcmod(outfname):

    global forcefield, args

    out = open(outfname, 'w')
    for i in range(len(fields)):
        field, format = fields[i], formats[i]
        out.write(field + "\n")
        for i, row in forcefield[field].iterrows():
            out.write(format % tuple(row.values))
        out.write("\n")


################################################### END OF FUNCTION DEFINITIONS ##########################################

if __name__ == "__main__":

    try:
        args = cmdlineparse()
        if args.OUT_FRCMOD == None:
            args.OUT_FRCMOD = "mod_%s" % args.FRCMOD

        load_frcmod(args.FRCMOD)
        write_frcmod(args.OUT_FRCMOD)

    except:
        type, value, tb = sys.exc_info()
        lines = traceback.format_exception(type, value, tb)
        print(''.join(lines))
        raise