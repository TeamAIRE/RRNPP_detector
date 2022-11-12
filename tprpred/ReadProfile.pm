package ReadProfile;

use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
use Exporter;
our $VERSION=1.00;
our @ISA          = qw(Exporter);
our @EXPORT       = qw(&ReadProfile);


#############################################################################
####             Subroutine  ReadMatrix                                  ####
#############################################################################
sub ReadProfile {
    # AA freqs: A    C    D    E    F    G    H    I    K     L    M    N    P    Q    R    S    T    V    W    Y
    my @pbg=(7.68,1.89,5.41,5.99,3.50,7.56,3.69,5.06,5.97,10.01,2.20,4.02,4.54,3.27,5.14,4.67,7.12,7.28,1.25,3.95);
    my @i2i=(0,2,3,4,5,6,7,8,10,11,12,13,15,16,17,18,19,21,22,24);
    my $line;                                    # input line
    my @row = ();                                # One row of pssm file
    my $a;                                       # amino acid of current row of pssm matrix
    my $i=0;                                     # position in profile
    my $pssm=$_[0];
    my $v=$main::v;
    my $W;

    open (PSSM, "<$pssm") || die "ERROR: Couldn't open $pssm: $!\n";
    $line=<PSSM>;
    if ($line=~/^\s*(\S*\d+\.?\S*)\s+(\d+\.?\S*)\s*$/) {  # profile already calibrated?
	$main::mu = $1; $main::sigma = $2;             # read mu and lamda
	$line=<PSSM>;
    }
    if ($line=~/NAME/) {
	while ($line=<PSSM>) {if ($line=~/^HMM/) {<PSSM>; <PSSM>; last;}} 
	while ($line=<PSSM>) { 
	    if (!($line=~/^\w \d+\s/)) {next;}
	    @row = split (/\s+/, $line); 
	    my $aa = shift (@row);                   # read character at beginning of line
	    if (($i+1)!=shift(@row)) {
		printf(STDERR "WARNING: inconsistent column numbering will be ingnored\n");
	    }
	    if (scalar(@row)>21 || scalar(@row)<20) {
		die ("ERROR: state $i of file $pssm contains ".scalar(@row)." instead of 20 aa frequencies \n$line\n");
	    }
	    for ($a=0; $a<20; $a++) {
		$main::Profile[$i][$i2i[$a]] = $row[$a]*0.001;
	    }
	    $i++;
	}
	$W=$i;
    } else {
	while ($line) { 
	    @row = split (/\s+/, $line); 
	    $a=ord(shift (@row))-65;                 # read character at beginning of line
	    if ($W && scalar(@row)!=$W) {
		die ("ERROR: pssm file contains rows of unqual lengths");}
	    else {$W=scalar(@row);}
	    $i=0;
	    while (@row) {
		$main::Profile[$i++][$a] = log( (shift(@row)+$main::pc)/(1+$main::pc) ) / log(2);
	    }
	    $line=<PSSM>;
	}
    }
    close (PSSM);
    
    # Set all non-standard one-letter codes 
    foreach ($i=0; $i<$W; $i++) {
	my $ScoreAny=0;
	for ($a=0; $a<20; $a++) {
	    $ScoreAny += 0.01*$pbg[$a]*$main::Profile[$i][$i2i[$a]];
	}
	$main::Profile[$i][ord("X")-65] = $ScoreAny;
	$main::Profile[$i][ord("J")-65] = $ScoreAny;
	$main::Profile[$i][ord("O")-65] = $ScoreAny;
	$main::Profile[$i][ord("B")-65] = $main::Profile[$i][ord("D")-65];
	$main::Profile[$i][ord("Z")-65] = $main::Profile[$i][ord("E")-65];
	$main::Profile[$i][ord("U")-65] = $main::Profile[$i][ord("C")-65];
    }

    if ($v>=3) {
	printf(STDERR "Repeat profile:\n     ");
	for ($a=0; $a<26; $a++) {
		printf(STDERR "    %1s ",chr($a+65));
	    }
	printf(STDERR "\n");
	for ($i=0; $i<$W; $i++) {
	    printf(STDERR "%3i: ",$i);
	    for ($a=0; $a<26; $a++) {	
		printf(STDERR "%5.2f ",$main::Profile[$i][$a]);
	    }
	    printf(STDERR "\n");
	}
	printf(STDERR "\n");
    }
    $main::W = $W;
    return;
}

