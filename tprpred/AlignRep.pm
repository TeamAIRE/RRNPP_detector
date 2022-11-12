#############################################################################
####             Subroutine AlignRep                                     ####
#############################################################################

# Usage: &AlignRep($seq,$v,$W,$K,$p0,$p1,$g)
# Finds the best alignment of sequence $seq (in ASCII) with repeat profile $Profile[$i][$a]
# Output:  return value:    whole-protein raw score
#          $PvalProt:       whole-protein p-value
#          @ScoreRep:       single-repeat raw scores
#          @PvalRep:        single-repeat p-values
#          @StartRep:       first residues of repeats
package AlignRep;

use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);
use Exporter;
our $VERSION=1.00;
our @ISA          = qw(Exporter);
our @EXPORT       = qw(&AlignRep &Align @StartRep @ResRep @ScoreRep @PvalRep $ScoreProt $PvalProt);
use Erf;


# Result variables to be exported
our @StartRep;             # start positions of single repeats
our @ResRep;               # Residues of single repeats
our @ScoreRep;             # single-repeat scores
our @PvalRep;              # single-repeat p-value
our $ScoreProt;            # whole-protein score
our $PvalProt;             # whole-protein p-value

sub AlignRep {                             
 
    # Local variables
    my $SeqLen = 0;                          # seq length l
    my $ScoreWin;                            # bit score over the window
    my $j;                                   # start position of window in sequence (first=0)
    my $jj;                                  # position in sequence: jj=j+i (first=0)
    my $i;                                   # position in profile (first=0)
    my @score=();                            # stores the window scores S[j] 

    my $seq = shift (@_);                    # @Residues contains the residues of sequence
    my @Residues = split (//,$seq);          # @Residues contains the residues of sequence
    my $L=scalar(@Residues);                 # length of the sequence ;
    my ($v,$W,$K,$p0,$p1,$g)=@_;
    my @Profile=@main::Profile;
    my @AlignPar=@_;

    @StartRep=();           #clean single-repeat hits
    @ScoreRep=();           #clean single-repeat hits
    @PvalRep=();            #clean single-repeat hitsm
    @ResRep=();             #clean single-repeat hits
    $ScoreProt=-1000;
    $PvalProt=1;
    if ($L<$W) {return -1000;}

    # Convert residues from ASCII to integers 
    for ($j=0; $j<$L; $j++) {
	$Residues[$j]=ord($Residues[$j])-65; # 65=ord("A")
	if ($Residues[$j]<0 || $Residues[$j]>25) {$Residues[$j]=23;}  #X=23
    }
	

    if ($W!=34) {
	# Calculate the score for each window starting at $j in sequence
	for ($j=0; $j<=$L-$W; $j++) 
	{
	    $ScoreWin = 0; $jj=$j;
	    for ($i=0; $i<$W;) {
		$ScoreWin += $Profile[$i++][$Residues[$jj++]];
	    } 
	    push (@score , $ScoreWin);    
	} 
    } else {
	# Window has exactly length 34 residues => use fast version
	my @profcol0=@{$Profile[0]}[0..33];
	my @profcol1=@{$Profile[1]}[0..33];
	my @profcol2=@{$Profile[2]}[0..33];
	my @profcol3=@{$Profile[3]}[0..33];
	my @profcol4=@{$Profile[4]}[0..33];
	my @profcol5=@{$Profile[5]}[0..33];
	my @profcol6=@{$Profile[6]}[0..33];
	my @profcol7=@{$Profile[7]}[0..33];
	my @profcol8=@{$Profile[8]}[0..33];
	my @profcol9=@{$Profile[9]}[0..33];
	my @profcol10=@{$Profile[10]}[0..33];
	my @profcol11=@{$Profile[11]}[0..33];
	my @profcol12=@{$Profile[12]}[0..33];
	my @profcol13=@{$Profile[13]}[0..33];
	my @profcol14=@{$Profile[14]}[0..33];
	my @profcol15=@{$Profile[15]}[0..33];
	my @profcol16=@{$Profile[16]}[0..33];
	my @profcol17=@{$Profile[17]}[0..33];
	my @profcol18=@{$Profile[18]}[0..33];
	my @profcol19=@{$Profile[19]}[0..33];
	my @profcol20=@{$Profile[20]}[0..33];
	my @profcol21=@{$Profile[21]}[0..33];
	my @profcol22=@{$Profile[22]}[0..33];
	my @profcol23=@{$Profile[23]}[0..33];
	my @profcol24=@{$Profile[24]}[0..33];
	my @profcol25=@{$Profile[25]}[0..33];
	my @profcol26=@{$Profile[26]}[0..33];
	my @profcol27=@{$Profile[27]}[0..33];
	my @profcol28=@{$Profile[28]}[0..33];
	my @profcol29=@{$Profile[29]}[0..33];
	my @profcol30=@{$Profile[30]}[0..33];
	my @profcol31=@{$Profile[31]}[0..33];
	my @profcol32=@{$Profile[32]}[0..33];
	my @profcol33=@{$Profile[33]}[0..33];
	
	# Calculate the score for each window starting at $j in sequence
	for ($j=0; $j<=$L-$W; $j++) 
	{
	    $jj=$j;
	    $ScoreWin = $profcol0[$Residues[$jj++]]  + $profcol1[$Residues[$jj++]]
		+ $profcol2[$Residues[$jj++]]  + $profcol3[$Residues[$jj++]]
		+ $profcol4[$Residues[$jj++]]  + $profcol5[$Residues[$jj++]]
		+ $profcol6[$Residues[$jj++]]  + $profcol7[$Residues[$jj++]]
		+ $profcol8[$Residues[$jj++]]  + $profcol9[$Residues[$jj++]]
		+ $profcol10[$Residues[$jj++]] + $profcol11[$Residues[$jj++]]
		+ $profcol12[$Residues[$jj++]] + $profcol13[$Residues[$jj++]]
		+ $profcol14[$Residues[$jj++]] + $profcol15[$Residues[$jj++]]
		+ $profcol16[$Residues[$jj++]] + $profcol17[$Residues[$jj++]]
		+ $profcol18[$Residues[$jj++]] + $profcol19[$Residues[$jj++]]
		+ $profcol20[$Residues[$jj++]] + $profcol21[$Residues[$jj++]]
		+ $profcol22[$Residues[$jj++]] + $profcol23[$Residues[$jj++]]
		+ $profcol24[$Residues[$jj++]] + $profcol25[$Residues[$jj++]]
		+ $profcol26[$Residues[$jj++]] + $profcol27[$Residues[$jj++]]
		+ $profcol28[$Residues[$jj++]] + $profcol29[$Residues[$jj++]]
		+ $profcol30[$Residues[$jj++]] + $profcol31[$Residues[$jj++]]
		+ $profcol32[$Residues[$jj++]] + $profcol33[$Residues[$jj++]];
	    push (@score , $ScoreWin);    
	} 
    }


    # Call algorithm to calculate whole-protein score and p-value, and single-repeat p-values
    unshift(@score,0);                       # &Score() ignores zero'th vector element
    &Align(\@score,@AlignPar);                         # generate subsequences 

    # Store residues of single-repeat hits in @ResRep
    foreach (@StartRep) {
	push(@ResRep,substr($seq,$_-1,$W));  # "-1" because $seq starts with residue 0 (not 1)
    }
} # End of the ScanSeq subroutine




#############################################################################
####             Subroutine Align                                       ####
#############################################################################

# &Score(@S,$v,$W,$K,$p0,$p1,$g) calculates the best alignment according to the vector of scores S(j)
# Input   @S:           vector of scores of length l=L-W+1
# Output  return value: whole-protein score as return value,
#         @ScoreRep:    single-repeat scores
#         @StartRep:    start positions of single repeats 
#
# Call function with &Score(\@S);
# Call function for subarray of @S with   
# my @subS = @S[$j1..$j2]; 
# &Score(\@subS);

sub Align
{
    my @S=@{shift @_};  # $S[$j] = score of profile aligned with subsequence x(j)..x(j+W-1)
    my $l=scalar(@S)-1; # $l = $L-$W+1 ("-1" because index 0 of @S is not useD!)
    my ($v,$W,$K,$p0,$p1,$g)=@_;
    my $L=$l+$W-1;      # length of sequence
    my $b=&ScoreOffset($p0,$l); # whole-score offset (bits)
    my $r=&ScoreOffset($p1,$K); # snug-fit reward (including score offset) (bits)
    my $j;              # current position in sequence
    my $k;              # length of gap between position j and last residue of last aligned repeat
    my @M;              # dynamic programming vectors
    my @B;              # backtracing vectors
    my $max;            # maximum score found so far
    my $jmax;           # position in sequence of maximum score
    my $bt;             # backtracing variable set by &maxbt: which argument was largest? (first=0)

    #clean single-repeat hits (if Align is called without AlignRep)
    @StartRep=();
    @ScoreRep=();
    @PvalRep=();
    @ResRep=(); 

    if ($v>=4) {
	printf ("l=%-3i   p0=%-6.2f => b=%6.2f      p1=%-6.2f => r=%6.2f\n",$l,$p0,$b,$p1,$r);
	printf ("p-val(%-6.2f,%2i)=%-6.2f=p0?   p-val(%-6.2f,%2i)=%-6.2f=p1?\n",$b,$l,&Pvalue(-$b,$l),$r,$K,&Pvalue(-$r,$K));
    }

    # @M[$j][$k] = max score of aligning seq up to residue j with any number of repeating profiles
    #              and a gap of length k after last repeat (of length AT LEAST K for k=K)
    
    # Initialize @M
    for ($j=0; $j<$W; $j++) {
	for ($k=0; $k<=$K; $k++) {
	    $M[$j][$k]= -1000;
	    $B[$j][$k]=-1;  # -1 means STOP when backtracing 
	}
    }

    # Dynamic programming 
    for ($j=$W; $j<=$L; $j++) {
    	$M[$j][0] = $r + $S[$j-$W+1] + &maxbt(@{$M[$j-$W]}[0..$K], 0, \$bt); # array slice
	$B[$j][0] = $bt;
	for ($k=1; $k<$K; $k++)	{
	    $M[$j][$k] = $M[$j-1][$k-1];
	}
	$M[$j][$K] = &max2bt($M[$j-1][$K-1]-$r+$b, $M[$j-1][$K]-$g, \$bt);
	$B[$j][$K] = $bt;
    }
    
    # Find maximum score
    $max = $M[$W][0];
    $jmax=$W;
    for ($j=$W+1; $j<=$L; $j++) {
	if($M[$j][0]>$max) {
	    $max = $M[$j][0];
	    $jmax = $j;
	}
    }


    # DEBUG: print @S
    if ($v>=5) {
	for ($j=$W; $j<=$l; $j++) {
	    printf("S[%3.3i]=%6.2f  ",$j-$W+1,$S[$j-$W+1]);
	    $k=0;
	    printf("M[%3i][%i]=%-5.1f ",$j,$k,$M[$j][$k]);
	    for ($k=1; $k<=$K; $k++) {
		printf("M[%i]=%-5.1f ",$k,$M[$j][$k]);
	    }
	    printf("\n");
	}
    }
    
    # Backtracing: find all single-repeat hits and their start positions
    $k=0;                   #initial state is M0 (k=0)
    $j=$jmax;               #initial position is j=jmax
    while ($k>=0) {         # k=-1 means STOP backtracing
	if ($k==0) {        #current state is M0
	    $bt = $B[$j][$k];
	    $j+=-$W+1;
	    if ($v>=4) {printf("Hit: S[%-3i] = %6.3f\n",$j,$S[$j]);}
	    if ($bt>=$K) {
		#Repeat has no predecessor 
		if (exists $StartRep[0] && $StartRep[0]<$j+$W+$K) {
		    #Repeat has successor   => Push p-value of new repeat based on $K 
		    if ($v>=4) {printf("Hit has no previous but next with %2i possible positions\n",$K);}
		    unshift (@PvalRep,&Pvalue($S[$j],$K));           
		} else {
		    #Repeat is single with no other nearby => Push p-value of new repeat based on $l
		    if ($v>=4) {printf("Hit has no previous nad no next with %2i possible positions\n",$l);}
		    unshift (@PvalRep,&Pvalue($S[$j],$l)); #Single repeat
		}
	    } else {
		#Repeat has predecessor
		if (exists $StartRep[0] && $StartRep[0]<$j+$W+$K) {
		    #Repeat has successor   => Push p-value of new repeat based on space between repeats
		    unshift (@PvalRep,&Pvalue($S[$j],$K));           
		    if ($v>=4) {printf("Hit has previous and next with %2i possible positions\n",$K);}
		} else {
		    #Repeat has no successor => Push p-value of new repeat based on $K
		    if ($v>=4) {printf("Hit has previous but no next with %2i possible positions\n",$K);}
		    unshift (@PvalRep,&Pvalue($S[$j],$K));           
		}
	    }
	    unshift (@StartRep,$j);                #store start position of repeat
	    unshift (@ScoreRep,$S[$j]);            #store score of repeat
	    if ($bt==$K+1 ) {last;}                #previous state was START state (0)
	    $k=$bt;
	    $j--;
	}
	elsif ($k==$K) {    #current state is MK
	    $bt = $B[$j][$k];
	    if ($bt==0) {$k=$K-1;}  else {$k=$K;}
	    $j--;
	}
	else {$k--; $j--;}  #current state is Mk (0<k<K)
    }
    
    $PvalProt = 1;
    for ($k=0; $k<scalar(@ScoreRep); $k++) {
	$PvalProt *= &Pvalue($S[$StartRep[$k]],$l-$k*$W)
    }
    $ScoreProt=$max-$r;
    return;

}

# max2bt(val1,val2,\$bt) finds maximum of values and puts index of maximum into $bt
sub max2bt {
    if ($_[0] > $_[1]) {
	${$_[2]}=0;
	return $_[0];
    } else {
	${$_[2]}=1;
	return $_[1];
    }
}

# maxbt(val1,...,valx,\$bt) finds maximum of values and puts index of maximum into $bt
sub maxbt {
    my $rbt=pop @_; # last element of @_ is address of $bt
    my $max = shift;
    my $i=0;
    $$rbt = 0;
    foreach $_ (@_) {
	$i++;
	if ($_>$max) {$max=$_; $$rbt=$i;} 
    }
    return $max;
}

# &Pvalue($score,$l) is the single-repeat p-value of a score in a protein of length L=l+W-1
sub Pvalue {
    my $score=$_[0];
    my $pvalue = 0.5*&erfc(($score-$main::mu)/$main::sigma);
    return 1-(1-$pvalue)**$_[1];
}    

# &ScoreOffset($prob,$l) is the score offset needed to make P(S+offset>0) = prob
sub ScoreOffset {
    my $p = 1-(1-$_[0])**(1/$_[1]);
    return -$main::mu - $main::sigma*&erfcinv(2*$p);
}

1;
