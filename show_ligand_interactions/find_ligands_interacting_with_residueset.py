"""
This script was writen to find docking poses of ihibitors that interact with residues 399+402+403+404 of c-Myc.
"""

from pymol import cmd, util
import os, re

################################################### FUNCTION DEFINITIONS ####################################################
def list_files(folder, pattern, full_path=False):
    """
        FUNCTION to list the files in 'folder' that match the 'pattern'.
    """
    if not folder:
        folder = "."
    folder = os.path.abspath(folder)
    fpaths = os.listdir(folder)
    fpattern = re.compile(pattern)
    file_list = list(filter(fpattern.search, fpaths))
    if full_path:
        file_list = [folder + "/" + f for f in file_list]
    return file_list

from collections import OrderedDict
class tree(OrderedDict):
    def __missing__(self, key):
        self[key] = type(self)()
        return self[key]

#############################################################################################################################
for inhibitor in ["IF4", "M19", "M19a"]:
    WORKDIR="/home/thomas/Dropbox/Myc/Manuscript_Myc_structure/Docking/vsw_%s"  % inhibitor
    resids="399+402+403+404"
    receptor_pattern="vsw_%s-SP_OUT_[0-9]_pv_receptor[0-9]+\.pdb"   % inhibitor
    contacts_mdict = tree()
    for receptor_pdb in list_files(WORKDIR, receptor_pattern, full_path=False):
        receptor = receptor_pdb.replace(".pdb", "")
        cmd.load(WORKDIR + "/" + receptor_pdb)
        ligname_pattern = receptor.replace("receptor1", "ligand[0-9]+") + ".pdb"
        ligpdb_list = list_files(WORKDIR, ligname_pattern, full_path=False)
        ligmol_list = [m.replace(".pdb", "") for m in ligpdb_list]
        cmd.set('h_bond_max_angle', 30)
        cmd.set('h_bond_cutoff_center', 3.6)
        cmd.set('h_bond_cutoff_edge', 3.2)
        for ligpdb in ligpdb_list:
            cmd.load(WORKDIR + "/" + ligpdb)
            ligmol = ligpdb.replace(".pdb", "")
            # contacts = cmd.dist('contacts', ligmol, receptor + " resid " + resids, 3.5, mode=2)
            cont399 = cmd.dist('cont399', ligmol, receptor + " and resid 399", 3.5, mode=2)
            cont402 = cmd.dist('cont402', ligmol, receptor + " and resid 402", 3.5, mode=2)
            cont403 = cmd.dist('cont403', ligmol, receptor + " and resid 403", 3.5, mode=2)
            cont404 = cmd.dist('cont404', ligmol, receptor + " and resid 404", 3.5, mode=2)
            # Save the contacts only if they exist
            if cont399+cont402+cont403+cont404 > 0:
                contacts_mdict[receptor_pdb][ligpdb] = (cont399, cont402, cont403, cont404)
            cmd.delete("cont*")
            cmd.delete(ligmol)
        cmd.delete(receptor)

    print("\nCONTACT RESULTS FOR INHIBITOR %s:"   % inhibitor)
    print("receptor_pdb\tligand_pdb\tcontact_399\tcontact_402\tcontact_403\tcontact_404\n")
    struct_files = ""
    for receptor_pdb in contacts_mdict.keys():
        struct_files += " " + receptor_pdb
        for ligpdb in contacts_mdict[receptor_pdb].keys():
            struct_files += " " + ligpdb
            print(receptor_pdb, ligpdb, contacts_mdict[receptor_pdb][ligpdb])
    print("To load the poses:")
    print("pymol " + struct_files)


