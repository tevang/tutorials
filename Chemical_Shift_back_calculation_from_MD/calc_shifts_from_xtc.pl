#!/usr/bin/perl -w
 
(@ARGV == 5) or die "Usage: calc_shifts_from_xtc.pl <pdb> <xtc> <start time in ps> <end time in ps> <timestep in ps>\n
Input: as described above. You can find the start and end time of your trajectory by doing: gmx check -f <xtc>.\n
Output: 'SHIFTS' file contains chemical shifts of all frames, 'consensus_chemical_shifts.tab' file contains the average chemical shifts of each atom.\n";

$pdbfile=$ARGV[0];
$trj=$ARGV[1];
$start=$ARGV[2];
$end=$ARGV[3];
$timestep=$ARGV[4];

#
# How many shifts we will be collecting ?
#

`echo 0 | gmx trjconv -s $pdbfile -f $trj -dump 1 -o ps0.pdb`;
`sparta+ -in ps0.pdb > /dev/null 2>&1`;
`/bin/rm -rf ps*.pdb`;
 
# Find the line where the stucture description and CS predictions start
open ( IN, "pred.tab" ) or die "Can not open pred.tab. Usage: calc_shifts <dcd> <psf>\n";
while ( $line = <IN> )
  {
    if ( $line =~ /^FORMAT/ )
      {
        last;
      }
  }
 
$line = <IN>;
$tot = 0;	# the total number of atoms with CS prediction
while ( $line = <IN> )
  {
    $ids[ $tot ] = substr( $line, 0, 14 );
    $tot++;
  }
 
close( IN );

`/bin/rm -rf *.tab`;
 
if ( $tot < 1 )
  {
    print "Too few atoms for calculating shifts. Something is wrong. Bye.\n";
    exit;
  }
 
 
print "Will be collecting data for $tot atoms. Starting ...\n";


#
# Will do it in sets of 400 structures ...
#
 
$first = int($start);

`/bin/rm SHIFTS`;

while( $first < int($end) ) {
	 
	$last = $first + 400 * $timestep -1;
	printf("Now processing set starting at picosecond %8d and ending at picosecond %8d\n", $first, $last );

	#for ($ps = $first; $ps <= $last; $ps += $timestep) { # not all ps are present in the trajectory (timestep for saving coordinates in MD was 2ps)
	#	#print "DEBUG: echo 0 | gmx trjconv -dump $ps -s $pdbfile -f $trj -o ps${ps}.pdb > /dev/null 2>&1";
	#	`echo 0 | gmx trjconv -dump $ps -s $pdbfile -f $trj -o ps${ps}.pdb > /dev/null 2>&1`;
	#}
	
	# real	17m26.892s
	# user	17m7.694s
	# sys	0m15.413s
	`echo 0 | gmx trjconv -b $first -e $last -sep -s $pdbfile -f $trj -o ps.pdb > /dev/null 2>&1`;
	
	`sparta+ -in ps*.pdb > /dev/null 2>&1`;
	 
	`/bin/rm -rf ps*.pdb *_struct.tab`;
	
	@files = glob("ps*.tab");
	
	
	if ( @files == 0 )
	  {
		last;
	  }
	
	
	@files = map {s/_pred\.tab//g; $_; } @files;
	@files = map {s/^ps//g; $_; } @files;
	my @ordered_files = sort {$a <=> $b} @files;
	
	
	foreach $file ( @ordered_files )
	{
	  #print "DEBUG: appending shifts from file ps${file}_pred.tab \n";
	  `tail -$tot ps${file}_pred.tab | awk '{printf "%8.3f ", \$5}' >> SHIFTS`; # append all chemical shifts of the current frame to a single line of SHIFTS file
	  `echo >> SHIFTS`;	# change line
	}
	 
	`/bin/rm -rf ps*.tab`; # remove the frames that have been processed
	 
	$first += 400 * $timestep;
 
}
 
print "\n\n";
 
#
# Calculate means + sigmas using SHIFTS file
#
open ( IN, "SHIFTS" ) or die "Can not open SHIFTS ??? How did this happen ???\n";
open ( OUT, '>>consensus_chemical_shifts.tab');

for ( $i=0 ; $i < $tot ; $i++ )
{
	$mean= 0.0;
	$nof_lines = 0;
	$std = 0.0;
	while ( $line = <IN> )
	  {
	    @data = split( ' ', $line );
 
	    $nof_lines++;
	    $delta = $data[ $i ] - $mean;
	    $mean += $delta / $nof_lines;
	    $std += $delta * ($data[ $i ] - $mean);
	  }

        printf "%s    %8.4f %8.4f\n", $ids[ $i ], $mean, sqrt( $std / ($nof_lines -1));
		printf OUT "%s    %8.4f %8.4f\n", $ids[ $i ], $mean, sqrt( $std / ($nof_lines -1));
        seek( IN, 0, 0 );        
}
 
close( IN );
close( OUT );
 
print "\nAll done.\n\n";
