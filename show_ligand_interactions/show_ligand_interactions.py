#show_ligand_interactions v.1.0
# author: Thomas Evangelidis, 2019
# License: BSD-2-Clause

from pymol import cmd, util
import show_bumps

def show_ligand_interactions(recsel="not hetatm", ligsel="hetatm", cutoff=5):
    """
DESCRIPTION

    Visualize interactions between receptor and ligand.

ARGUMENTS

    recsel = string: atom selection of the receptor {default: "not hetatm"}

    ligsel = string: atom selections of the ligand {default: "hetatm"}

    cutoff = float: show as sticks all receptor residues within this distance from the ligand {default: 5.0}
    """
    cmd.select('ligand', ligsel)
    cmd.select('receptor', recsel)

    cmd.bg_color('white')
    cmd.show_as('cartoon')
    cmd.show_as('sticks', 'hetatm or ligand')
    cmd.show_as('nonbonded', "resn HOH+T3P+WAT within %s of ligand" % cutoff)
    cmd.set('cartoon_transparency', 0.2)
    cmd.spectrum(selection=recsel+" or "+ligsel,byres=1)
    util.cbag('not name C*')
    cmd.set('cartoon_fancy_helices', 1);
    cmd.show("sticks", "(hydro)");
    cmd.select("pocket", "byres (receptor within %s of ligand)" % cutoff);
    cmd.show("sticks", "pocket")
    cmd.hide('(h. and (e. c extend 1))')
    cmd.set('h_bond_max_angle', 30)
    cmd.set('h_bond_cutoff_center', 3.6)
    cmd.set('h_bond_cutoff_edge', 3.2)
    cmd.dist('ligand_Hbonds', 'ligand', 'receptor', 3.5, mode=2)
    cmd.set('dash_radius', 0.15)
    # now set the label options
    cmd.set('label_size', 20)
    cmd.set('label_position', [0,0,10])
    cmd.orient("ligand")

cmd.extend('show_ligand_interactions', show_ligand_interactions)
