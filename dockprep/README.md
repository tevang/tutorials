This is a Python script to prepare a receptor-ligand complex for scoring (e.g. for PM6/COSMO scoring) using [UCSF Chimera](https://www.cgl.ucsf.edu/chimera/). It takes as input either (i) the protein-ligand complex in one pdb file, or (ii) the receptor in a pdb file and the ligand in a separate file of any format (pdb, mol, mol2, sdf). The output is always a pdb file of the protein-ligand complex in which the residue name of the ligand is "LIG".

If you encounter problems with the input pdb file then try correcting it using one of the following programs:
1) pdb4amber from AmberTools (https://github.com/Amber-MD/pdb4amber)
2) pdbfixer (https://github.com/pandegroup/pdbfixer)
Sometimes you may need to rename atoms manually in Chimera (e.g. the N- and C-terminal caps), but this is out of the scope of this script.


For a full list of options run `dockprep.py -h`.


You can launch the script either using [PyChimera](https://pychimera.readthedocs.io/en/latest/) (to install it read this [tutorial](https://github.com/tevang/tutorials/tree/master/create_alternative_protonations)) 
```
pychimera $(which dockprep.py) -rec example_files/3K5C-BACE.pdb -lig $mol -cmethod gas -neut
```
or the vanilla UCSF Chimera executable.
```
chimera --nogui --nostatus --script "$(which dockprep.py) -rec example_files/3K5C-BACE.pdb -lig example_files/3K5C-BACE_1.mol -cmethod gas -neut"
```
To prepare multiple protein-ligand complexes with dockprep.py **in parallel** on the Unix Shel, first create a file with all the individual commands like this:
```
for mol in $(ls example_files/*mol)
do
echo "chimera --nogui --nostatus --script \"$(which dockprep.py) -rec example_files/3K5C-BACE.pdb -lig $mol -cmethod gas -neut\""
done > commands.txt
```

Then launch it using [GNU parallel script](https://www.gnu.org/software/parallel/)
``` 
parallel -j3 < commands.txt
```
