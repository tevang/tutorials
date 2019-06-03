# UCSF Chimera Tutorial: visualize ligand properties (charges, bond lengths & angles, etc.) 

Load the two ligand files into Chimera and label atom names and charges
```
labelopt info "%(name)s %(charge)+.3f"
label
```
or if you want the atom type as well
```
labelopt info "%(name)s(%(idatmType)s) %(charge)+.3f"
label
```
or to show only one molecule and label by atom type
```
~display #1
labelopt info "%(idatmType)s"
label
```
to change the label font size go Favorites->Preferences->"Category:"Background.
to change the background color (althouhg black is good to display labels)
```
background solid white
```
change representation (only the "wire" works well with .mol2 files)
```
represent wire
```
play with scale command to zoom in the image as much as you want
```
scale 1.4
```


You can label the bonds with their length with this Python:

```python
from chimera import openModels, Molecule
for mol in openModels.list(modelTypes=[Molecule]):
	for b in mol.bonds:
		b.label = "%.2f" % b.length()
```

You can label atoms, bonds, and residues in this fashion, but there is no provision to label bond angles per se.  If you meant torsion angles and wanted to add that to the bond label, you can get the floating-point dihedral value of four atoms with this code:

```python
import chimera
dihed_val = chimera.dihedral(a1.coord(), a2.coord(), a3.coord(), a4.coord())
```
