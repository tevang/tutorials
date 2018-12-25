Requirements:
* Perl
* GROMACS Tools
* Sparta+

Usage: calc_shifts_from_xtc.pl <pdb> <xtc> <start time in ps> <end time in ps> <timestep in ps>

Input: you find the start and end time of your trajectory by: gmx check -f <xtc>.

Output: 'SHIFTS' file contains chemical shifts of all frames, 'consensus_chemical_shifts.tab' file contains the average chemical shifts of each atom.

