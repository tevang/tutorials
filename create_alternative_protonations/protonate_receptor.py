#!/usr/bin/env python

__author__="Thomas Evangelidis"
__email__="tevang3@gmail.com"



from argparse import ArgumentParser, RawDescriptionHelpFormatter
from itertools import combinations, permutations
import sys, gc, os
from operator import itemgetter
from ete3 import Tree


## Parse command line arguments
def cmdlineparse():
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter, description="""
DESCRIPTION:

This is a Python script to create all alternative protonation state combinations of a protein given a ligand and a specified radius. 
It could be useful in case you have a receptor but you are not sure about the protonation states of some residues in the binding site 
and you want to do docking, MD, or any other structure-based drug design method using all of alternative receptor protonations.

The script must be executed with PyChimera, a Python wrappen of UCSF Chimera, which searches for protonatable standard residues (ASP, GLU, HIS) 
around the ligand, finds all possible combinations of protonation states, and writes a pdb file for each combination. Since the number 
of combinations of more than 6 protonatable residues becomes very large, the user can fix some residues to a give protonated/unprotonated state. 
See the examples below. You can also get the same info by typing `protonate_receptor.py -h`.

                         """,
                            epilog="""
### EXAMPLE 1: list all protonatable residues within 4 Angstroms from the ligand.
pychimera $(which protonate_receptor.py) -rec 1a30_protein.pdb -lig 1a30_ligand.sdf -r 4.0 -list

### EXAMPLE 2: keep ASP_29.A ASP_30.A fixed to the unprotonated state and create alternative protonations for all the rest.
pychimera $(which protonate_receptor.py) -rec 1a30_protein.pdb -lig 1a30_ligand.sdf -r 4.0 -fix ASP_29.A -fix ASP_30.A

### EXAMPLE 3: protonate all protein residues.
pychimera $(which protonate_receptor.py) -rec 1a30_protein.pdb -lig 1a30_ligand.sdf   

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

def write_protonated_structure(protonations):

    global residues, args

    id2state = {}
    pdb = args.RECEPTOR.replace(".pdb", "")
    for rstate in protonations:
        state, resid = rstate.split('_')
        id2state[resid] = state
        pdb += "_%s%s" % (state , resid.replace(".",""))
    pdb += ".pdb"
    # Alter the protonation states
    for r in residues:
        try:
            r.type = id2state[str(r.id)]
        except KeyError:
            continue
    # Write the structure
    rc("del H")
    rc("addh")
    rc("write format pdb #0 %s" % pdb)

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

    print "Building Protonation Trees from peptide %s" % list(peptide)
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
        resid, chain = leaf.name.split(".")
        protonations.append((leaf.state, resid, chain))
        for ancestor in leaf.get_ancestors()[:-1]:  # skip the root
            resid, chain = ancestor.name.split(".")
            protonations.append((ancestor.state, resid, chain))
        protonations.sort(key=itemgetter(2, 1))     # sort by chain and resid to avoid permutations of the same combination
        protonations = tuple(["%s_%s.%s" % (t[0], t[1], t[2]) for t in protonations])
        all_protonations_set.add(protonations)
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
        residue_states[str(r.id)] = states

    if args.LIST_PROTONATABLE:
        protonatable_rstates = ["%s_%s" % (name,id) for name,id in zip(protonatable_resnames, protonatable_resids)]
        print "\n~~~ The protonatable residues within %.3f Angstroms from the ligand are: %s\n" % (args.RADIUS, " ".join(protonatable_rstates))
        sys.exit(0)


    for rstate in args.FIXED_STATES:
        state, resid = rstate.split('_')
        residue_states[resid] = [state]
        try:
            protonatable_resids.remove(resid)
            print "Fixed resid %s to %s state." % (resid, state)
        except ValueError:
            print "Warning: residue %s is not within the specified distance from the ligand or is not a valid residue, " \
                  "therefore it will be ignored." % rstate

    all_protonations = set()
    for peptide in permutations(protonatable_resids, len(protonatable_resids)):
        all_protonations = all_protonations.union(build_Protonation_Tree(peptide, residue_states))

    # Finally create and write the protonated structures
    all_protonations = list(all_protonations)
    all_protonations.sort(key=lambda x: x.count)
    for protonations in all_protonations:
        print "Writing structure with the following protonation states: ", protonations
        write_protonated_structure(protonations)
