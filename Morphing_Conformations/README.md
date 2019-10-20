Morphing protein-ligand conformations with UCSF Chimera
==========================

In this tutorial we will do simple morphing between 2 protein-ligand complex conformations, but using the same steps
we can use even more structures.

```python
open data/state1.pdb
open data/state2.pdb
del :Na+
del :WAT
ribspline cardinal smoothing both stiffness 0.8
ribrepr edged
set silhouette
set silhouette_width 2
set dcstart 0.3
# under Tools->Depiction->Ribbon Style Editor, set Coil width and height to 0.15. The rest look fine.
sel ~:LIG
namesel rec
sel :LIG
namesel lig
color pink rec
color green lig
zonesel lig 3.5 rec
namesel pocket
~display rec
display pocket
background solid white

# then Tools->Structure Comparison->Morph Conformations, add the two models, increase the number of intermediate conformations (e.g. 60),
# do NOT select minimize, click create. Alternatively use the following command lines:
morph start #0 name holo2apo frames 60
morph interpolate #1 name holo2apo
morph movie name holo2apo nogui true minimize false

```

