#!/usr/bin/perl -w
use strict;
use lib './';  
use AlignRep;
use ReadProfile;

#############################################################################
####                    Main                                             ####
#############################################################################

no warnings 'shadow';
# Default parameters
my $K=10;                 # snug-fit reward = r*(K-k) for k<=K  
my $g=0;                  # gap penalty = g*(k-K) for k>=K
my $p0=0.01;              # prob of at least one window score above 0 (no matter what sequence length)
my $p1=0.01;              # prob of at least one window score above 0 in the $K residues following a match
my $EvalProtThr=1;        # Maximum whole-protein e-value for output to alifile
my $PvalRepThr=1;         # Maximum single-repeat p-value for output to alifile  
my $nhitsThr=0;           # Minimum number of single-repeat hits for output to alifile
my $lmax=10000;           # Maximum number of lines in hit list
my $EvalList=10000;       # Maximum whole-protein e-value for output to hit list
my $nhitsList=0;          # Minimum number of single-repeat for output to hit list
my ($sigma, $mu);         # Parameters of single-residue score distribution for repeat profile
my $v;                    # verbose mode 0:only important warnings  2:DEBUG
my $pc;                   # pseudocounts for simple pssm format
our $pc=0.05;  
our $v=0;

# Variable declarations 
my @Profile = ();         # store the profile value in 2D matrix
my $W=0;                  # length of repeat profile ('window') (will be determined when loading pssm)
my $Identifier = '';      # seq IDs
my $Description;          # sequence description
my $SeqNum = 0;           # counts sequences in inputfile
my @HitList=();           # Hit list with result summary, one hit per array element
my @AliList=();           # Array of pointers to arrays with sequences of single-hit repeats 
my @AlignRepPar;          # parameters to hand over to &AlignRep
my $pssm="pssm3_34";      # repeat pssm file
my $infile="";            # name  of input file to be scanned
my $alifile="";           # output alignment file with repeats in fasta format
my $append=0;
my $line;
my @lines;


our $sigma =0;
our $mu    =0;
our @Profile;
our $W=0;

my $usage = "
Repeat identification: 
Scan (multiple/single) sequence fasta file for occurence of repeat proteins. 
Return sorted list of sequences with their whole-protein p-values, repeat positions and repeat p-values
Usage: repid.pl infile [options] 
Options:
 -r  pssm-file pssm of repeat   
 -v  int       verbose mode  (def=$v)
 -o  alifile   Print out the aligned repeats in fasta format into alifile
 -a  alifile   Same, but append to alifile
 -e  ]0,inf[   Maximum whole-protein e-value for output of repeat residues to alifile (def=$EvalProtThr)
 -p  ]0,inf[   Maximum single-repeat p-value for output of repeat residues to alifile (def=$PvalRepThr)
 -n  ]0,inf[   Minimum number of single-repeat hits for output to alifile
 -l  int       Maximum number of lines in hit list (def=$lmax)
 -E  ]0,inf[   Maximum whole-protein e-value for output in hit list (def=$EvalList)
 -N  ]0,inf[   Minimum number of single-repeat hits for output to hit list
Search options:
 -p0 [0,1]     Probability for single-repeat chance score >0 in whole protein def=$p0)
 -p1 [0,1]     Probability for single-repeat chance score >0 in $K residues following repeat (def=$p1) 
 -K  [1,inf[   When a repeat follows another by <=K residues a reward is paid (def=$K)
 -g  [0,inf[   Gap penalty (in bits) for gaps between repeat hits of > K residues

Example: repid.pl lal.fas -v -o lal.afas\n\n";

# Processing command line options
if (@ARGV<1) {die $usage;}

# if input file is given without -i option: insert -i before pssm name
if ($ARGV[0]!~/^-/) {unshift(@ARGV,"-i");}
my $options="";
for (my $i=0; $i<@ARGV; $i++) {$options.="$ARGV[$i] ";}

# Set verbose mode?
if ($options=~s/-v\s+(\d+)//g) {$v=$1;}
elsif ($options=~s/-v\s+//g) {$v=1;}

# Set options
if ($options=~s/-i\s+(\S+)//g) {$infile=$1;}
if ($options=~s/-r\s+(\S+)//g) {$pssm=$1;}
if ($options=~s/-o\s+(\S+)//g) {$alifile=$1; $append=0;}
if ($options=~s/-a\s+(\S+)//g) {$alifile=$1; $append=1;}
if ($options=~s/-e\s+(\S+)//g) {$EvalProtThr=$1;}
if ($options=~s/-p\s+(\S+)//g) {$PvalRepThr=$1;}
if ($options=~s/-n\s+(\S+)//g) {$nhitsThr=$1;}
if ($options=~s/-N\s+(\S+)//g) {$nhitsList=$1;}
if ($options=~s/-E\s+(\S+)//g) {$EvalList=$1;}
if ($options=~s/-l\s+(\S+)//g)  {$lmax=$1;}

# Set expert options: p0, p1, K, g, nmax
if ($options=~s/-p0\s+(\S+)//g) {$p0=$1; if ($p0>=1 || $p0<=0) {$p0=0.99;}}
if ($options=~s/-p1\s+(\S+)//g) {$p1=$1; if ($p0>=1 || $p1<=0) {$p1=0.99;}}
if ($options=~s/-K\s+(\S+)//g)  {$K=$1;}
if ($options=~s/-g\s+(\S+)//g)  {$g=$1;}

# Warn if unknown options found
if ($options!~/^\s*$/) {
    $options=~s/^\s*(.*?)\s*$/$1/g; 
    print("WARNING: unknown options '$options'\n");
}

if ($infile eq "") {die $usage;}


# Read the pssm file 
&ReadProfile($pssm);  
@AlignRepPar=($v,$W,$K,$p0,$p1,$g);

my $Seq="";
my $Nseqs=0;
    
# Count number of sequences in infile
open (INFILE, "<$infile") || die "ERROR: Couldn't open $infile: $!\n";
while ($line=<INFILE>) {
    if( ($line=~/^>/) ) {$Nseqs++;}
}
close (INFILE);

# Go through infile sequence by sequence and call AlignRep
open (INFILE, "<$infile") || die "ERROR: Couldn't open $infile: $!\n";
while (1)
{
    $line=<INFILE>;
    if (!$line) {
	if ($Seq ne "") {&AlignAndPushResults;} 
	last;
    } else { 
	@lines = split (/
/,$line);
    }
    foreach $line (@lines) {
	if($line=~/^>/)
	{
	    if ($Seq ne "") {&AlignAndPushResults;} 
	    $line=~/^>\s*(\S+)\s*(.*)/;
	    $Identifier = $1; 
	    $Description = $2;
	    $Seq="";
	} else {
	    chomp($line);
	    $Seq.=$line;
	}
    }
}
close (INFILE);



# Print sorted HitList
my @SortedHitList = sort CompareHitList @HitList;
if ($v>=1) {printf(STDERR " $SeqNum\n");}
print("Nr.\tSequence\tLen\tScore\tP-value\tE-value\tProbab\t#\tDescription\tpos:score ...\n");
for (my $i=0; $i<scalar(@HitList) && $i<$lmax; $i++) 
{
    printf("%-3i\t$SortedHitList[$i]\n",$i);
}
printf("\n");

# Print residues of single-hit repeats in sorted order
my @SortedAliList = sort CompareAliList @AliList;
if ($alifile) {
    if ($append) {
	open (ALIFILE, ">>$alifile")  or die ("ERROR: cannot open $alifile: $!\n");
    } else {
	open (ALIFILE, ">$alifile") or die ("ERROR: cannot open $alifile: $!\n");
    }
    for (my $i=0; $i<scalar(@SortedAliList); $i++) 
    {
	my ($Identifier,$L,$Pval,$Eval,$ResRep,$StartRep,$PvalRep) = @{$SortedAliList[$i]};
	if ($Eval<=$EvalProtThr && scalar(@{$ResRep})>=$nhitsThr) {
	    for (my $k=0; $k<scalar(@{$ResRep}); $k++) {
		if (@{$PvalRep}[$k]<=$PvalRepThr) {
		    my $IdPos = sprintf("%s(%i-%i)",$Identifier,@{$StartRep}[$k],@{$StartRep}[$k]+$W-1);
		    printf(ALIFILE ">%-37s Len=%-4i E-prot=%-9.2e P-rep=%-8.1e\n",
			   $IdPos,$L,$Eval,@{$PvalRep}[$k]);
		    printf(ALIFILE "%s\n",@{$ResRep}[$k]);
		}
	    }
	}
    }
    printf(ALIFILE "\n");
    close (ALIFILE);
}

exit;

#############################################################################
####        Prepare sequence $Seq, call AlignRep and record results      ####
#############################################################################
sub AlignAndPushResults {
    my $line; 
    $SeqNum++;
    $Seq =~ tr/a-z/A-Z/;                         # convert to upper case
    if ($v>=4) {printf(STDERR ">$Identifier\n$Seq\n");}
    &AlignRep($Seq,@AlignRepPar);                        # calling the main subroutine
    if ($v ==1 && !($SeqNum % 100)) {printf(STDERR ".");}
    if ($v ==1 && !($SeqNum % 10000)) {printf(STDERR " $SeqNum\n");}
    if ($v >1 && !($SeqNum % 10)) {printf(STDERR ".");}
    if ($v >1 && !($SeqNum % 1000)) {printf(STDERR " $SeqNum\n");}
    
    # Push HitList for sequence onto @HitList ?
    if ($PvalProt*$Nseqs<=$EvalList && scalar(@ScoreRep)>=$nhitsList) {
	my $logp;
	if ($PvalProt<=0) {$logp=999;} else {$logp=-log($PvalProt)/0.693147;}
	my $probpos=&ProbPos($logp);
	$line=sprintf ("%s\t%4i\t%6.1f\t%8.1E\t%8.1E \t%6.2f%%%%\t%2i\t%s\t",$Identifier,length($Seq),$logp,$PvalProt,$PvalProt*$Nseqs,100*$probpos,scalar(@ScoreRep),$Description);
	for (my $i=0; $i<scalar(@ScoreRep); $i++) {
#	    $line.=sprintf(" %4i:%5.1f %7.1e",$StartRep[$i],$ScoreRep[$i],$PvalRep[$i]);
#	    $line.=sprintf(" %4i:%5.1f",$StartRep[$i],-log($PvalRep[$i])/log(2));
	    $line.=sprintf(" %4i: %7.1e",$StartRep[$i],$PvalRep[$i]);
	}
	push(@HitList,$line);
    }	

    if ($v>=4) {printf(STDERR "%4i $line\n",$SeqNum);}
    elsif ($v>=3 && $PvalProt*$Nseqs<=$EvalProtThr) {printf(STDERR "%4i $line\n",$SeqNum);}
    
    if ($alifile) {
	my @AliListEl=();
	my @locResRep  =@ResRep;
	my @locStartRep=@StartRep;
	my @locPvalRep =@PvalRep;
	push(@AliListEl,$Identifier);
	push(@AliListEl,length($Seq));
	push(@AliListEl,$PvalProt);
	push(@AliListEl,$PvalProt*$Nseqs);
	push(@AliListEl,\@locResRep);
	push(@AliListEl,\@locStartRep);
	push(@AliListEl,\@locPvalRep);
	push(@AliList,\@AliListEl);
    }
    return;
}


#############################################################################
####             Compare functions for sorting                           ####
#############################################################################
sub CompareHitList
{
    $a=~/^\s*\S+\s+\S+\s+\S+\s+(\S+)/;
    my $Evala=$1;
    $b=~/^\s*\S+\s+\S+\s+\S+\s+(\S+)/;
    my $Evalb=$1;
    return ($Evala <=> $Evalb);
}

sub CompareAliList
{
    my @a=@{$a};
    my @b=@{$b};
    return ($a[3] <=> $b[3]);
}


#############################################################################
####         Calculate probability that protein is positive              ####
#############################################################################
sub ProbPos {
    my $Nnegs=652000;
    my $log2=0.693147;
    my $logp=shift @_;
    my $dbin=5;
    my $a=0.20666;
    my $b=7.6649;
    my $c=1.99086E-02;
    my $A=592.5;          # correction for negative distr 
    my $C=10;             # correction for negative distr
    my $D=0.5;            # correction for negative distr

    my $negcount= $Nnegs*$log2*(exp(-$logp*$log2) + exp(-$D*($logp-$C)*$log2)/$A);
    my $poscount  = 2**(&min($a*$logp,$b-$c*$logp))/$dbin;
    return $poscount/($negcount+$poscount);
}

# minimum
sub min {
    my $min = shift @_;
    foreach (@_) {
	if ($_<$min) {$min=$_} 
    }
    return $min;
}
