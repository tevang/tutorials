#!/usr/bin/env python

__author__="Thomas Evangelidis"
__email__="tevang3@gmail.com"



from argparse import ArgumentParser, RawDescriptionHelpFormatter
from itertools import combinations, permutations
import sys, gc, os

## Parse command line arguments
def cmdlineparse():
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter, description="""
DESCRIPTION:

This is a Python script to create all alternative protonation state combinations of a protein given a ligand and a specified radius. 
It could be useful in case you have a receptor but you are not sure about the protonation states of some residues in the binding site 
and you want to do docking of MD using all of them.

The script must be executed with UCSF Chimera, which searches for protonatable standard residues (ASP, GLU, HIS) around the ligand, 
finds all possible combinations of protonation states, and writes a pdb file for each combination. Since the number of combinations 
of more than 4 protonatable residues becomes very large, the user can fix some residues to a give protonated/unprotonated state. See 
the examples below.

### EXAMPLE 1: list all protonatable residues within 8 Angstroms from the ligand.

`chimera --nogui --nostatus --script "$(which protonate_receptor.py) -rec 1a30_protein.pdb -lig 1a30_ligand.sdf -r 8.0 -list"`

### EXAMPLE 2: keep GLU_34.B, ASP_30.B, ASP_29.B, ASP_29.A, and ASP_30.A fixed, and create alternative protonations for all the rest 
(namely ASP_25.B and ASP_25.A).
`chimera --nogui --nostatus --script "$(which protonate_receptor.py) -rec 1a30_protein.pdb -lig 1a30_ligand.sdf -r 4.0 -fix GLU_34.B -fix ASP_30.B -fix ASP_29.B -fix ASP_29.A -fix ASP_30.A"`

### EXAMPLE 3: protonate all protein residues.
`chimera --nogui --nostatus --script "$(which protonate_receptor.py) -rec 1a30_protein.pdb -lig 1a30_ligand.sdf"`
                            """,
                            epilog="""
EXAMPLE:

chimera --nogui --nostatus --script "$(which protonate_receptor.py) -rec 1a30_protein.pdb -lig 1a30_ligand.sdf -r 8.0"

    """)
    parser.add_argument("-list", dest="LIST_PROTONATABLE", required=False, default=False, action='store_true',
                        help="List the protonatable residues within the binding site and exit.")
    parser.add_argument("-rec", dest="RECEPTOR", required=False, default=None, type=str,
                        help="pdb file with the apo form of the receptor.")
    parser.add_argument("-lig", dest="LIGAND", required=False, default=None, type=str,
                        help="sdf or mol2 file with optimized ligand structure from which to find the binding site residues.")
    parser.add_argument("-r", dest="RADIUS", required=False, default=8.0, type=float,
                        help="The distance around the ligand within which residues will be protonated. Use '-r 0' if you "
                             "want to protonate the whole protein. Default: %(default)s.")
    # parser.add_argument("-p", dest="PYTHONPATH", required=False, default=None, type=str,
    #                     help="the PYTHONPATH environment variable.")
    parser.add_argument("-fix", dest="FIXED_STATES", required=False, default=[], type=str, action='append',
                        help="the residue to fixed to one state. E.g. '-fix ASP_30.A GLH_24.B' will NOT produce any structure with"
                             "ASH_30.A or GLU_24.B. This is useful when you have >4 protonatable residues within the binding site"
                             " and you want to reduce the number of combinations.")

    args = parser.parse_args()
    return args

########################################################## FUNCTION DEFINITIONS ####################################################

def change_protonation(remaining_residues):
    if not remaining_residues:
        global pdb_num
        rc("del H")
        rc("addh")
        ### For Debugging
        # print "List of current models:"
        # rc("list models")
        # print "List of current selections:"
        # rc("list selection level residue")
        ###
        rc("write format pdb #0 prot%d.pdb" % (pdb_num))
        pdb_num += 1
        return
    r = remaining_residues[0]
    if r.type == "GLU":
        states = ["GLU", "GLH"]
    elif r.type == "ASP":
        states = ["ASP", "ASH"]
    elif r.type == "HIS":
        states = ["HIE", "HID", "HIP"]
    else:
        states = [r.type]
    for state in states:
        r.type = state
        change_protonation(remaining_residues[1:])

def protonation(remaining_resids):
    global residue_states, state_dict
    if not remaining_resids:
        return
    resid = remaining_resids[0]
    states = residue_states[resid]
    for state in states:
        state_dict[resid] = state
        print state_dict
        protonation(remaining_resids[1:])


def populate_leaves(Peptide_Tree, resid, residue_states):
    """
        FUNCTION that adds new branches to the leaves of the Tree.
        ARGUMENTS:
        Peptide_Tree:    The Tree structure with connectivities
        RETURNS:
        (Peptide_Tree, BOOLEAN):    A tuple with elements the input Tree structure with new branches (if applicable), and a BOOLEAN value which is True if the function added
                                       new leaves to the Tree, or False otherwise
    """

    number_of_new_leaves = 0
    for leaf, prev_resid in zip(Peptide_Tree.iter_leaves(), Peptide_Tree.iter_leaf_names()):
        try:
            for state in residue_states[resid]:
                new_child = leaf.add_child(name=resid)  # add a new brach to the current TOCSY add index (leaf) with length the respective probability
                new_child.add_features(state=state)
                number_of_new_leaves += 1
                # print "DEBUG: adding connection: ",name,"-->",NOESYaaindex
        except(KeyError, IndexError):
            continue

    # print Peptide_Tree.get_ascii(show_internal=True, compact=False, attributes=["name", "dist", "occupancy", "numOfResonances"])
    # print Peptide_Tree.get_ascii(show_internal=True, compact=False)
    if number_of_new_leaves > 0:
        return (Peptide_Tree, True)
    else:
        return (Peptide_Tree, False)


def build_Protonation_Tree(peptide, residue_states):

    print "Building Protonation Trees from peptide %s" % peptide
    expand_tree = True
    Peptide_Tree = Tree()
    Root = Peptide_Tree.get_tree_root()
    Root.add_feature("name", "root")
    Root.add_feature("state", "delete")
    level = 0
    sys.stdout.write("Expanding tree from level ")
    while level < len(peptide):
        sys.stdout.write(str(level) + " ")
        sys.stdout.flush()
        Peptide_Tree, expand_tree = populate_leaves(Peptide_Tree, peptide[level], residue_states)
        level += 1
    # Print the Tree
    # print Peptide_Tree.get_ascii(show_internal=True, compact=False)
    # print Peptide_Tree.get_ascii(show_internal=True, compact=False, attributes=["name", "dist", "occupancy", "numOfResonances"])

    print "\nSaving protonations from Tree..."

    all_protonations_set = set()
    for leaf in Peptide_Tree.iter_leaves():
        protonations = []
        protonations.append("%s_%s" % (leaf.state, leaf.name))
        for ancestor in leaf.get_ancestors()[:-1]:  # skip the root
            protonations.append("%s_%s" % (ancestor.state, ancestor.name))
        all_protonations_set.add(tuple(reversed(protonations)))
        del protonations
        del ancestor
        del leaf
        # Peptide_Tree = None
    del Peptide_Tree
    gc.collect()
    return all_protonations_set

######################################################################################################################################


if __name__ == "__main__":
    args = cmdlineparse()

    from chimera import runCommand as rc
    from chimera.selection import currentResidues

    pdb_num = 1
    rc("open %s" % args.RECEPTOR)   # load the receptor
    rc("open %s" % args.LIGAND) # load the ligand
    if args.RADIUS > 0:
        rc("sel #1 z<%f & ~ #1" % args.RADIUS)
    elif args.RADIUS == 0:
        rc("sel #0")
    residues = currentResidues()    # get the residue of the pocket
    residue_states = {}
    protonatable_resids = []
    protonatable_resnames = []
    for r in residues:
        if r.type == "GLU":
            states = ["GLU", "GLH"]
            protonatable_resids.append(str(r.id))
            protonatable_resnames.append(r.type)
        elif r.type == "ASP":
            states = ["ASP", "ASH"]
            protonatable_resids.append(str(r.id))
            protonatable_resnames.append(r.type)
        elif r.type == "HIS":
            states = ["HIE", "HID", "HIP"]
            protonatable_resids.append(str(r.id))
            protonatable_resnames.append(r.type)
        else:
            states = [r.type]
        residue_states[r.id] = states

    if args.LIST_PROTONATABLE:
        protonatable_rstates = ["%s_%s" % (name,id) for name,id in zip(protonatable_resnames, protonatable_resids)]
        print "\n~~~ The protonatable residues within %.3f Angstroms from the ligand are: %s\n" % (args.RADIUS, " ".join(protonatable_rstates))
        sys.exit(0)


    # for p in args.PYTHONPATH.split(':'):
    #     sys.path.insert(0, p)
    # os.environ['PYTHONPATH'] = args.PYTHONPATH
    # print sys.path
    from ete3 import Tree

    for rstate in args.FIXED_STATES:
        state, resid = rstate.split('_')
        residue_states[resid] = [state]
        protonatable_resids.remove(resid)

    all_protonations = set()
    for peptide in permutations(protonatable_resids, len(protonatable_resids)):
        all_protonations = all_protonations.union(build_Protonation_Tree(peptide, residue_states))

    print all_protonations
