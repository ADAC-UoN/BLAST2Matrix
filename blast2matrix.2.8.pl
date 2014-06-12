#!/usr/bin/perl 
use strict;
use warnings;



#########################################################
# RD Emes JAN 2005  BLAST2MATRIX.2.0.pl
# Run a series of blastall searches using a fasta file supplied as a query
# against a series of fasta files which will be formatted using Formatdb
# outputs a matrix of presence/absence from a database depending on the eval cut off defined in blast parsing stage. 
#
# INPUT = db.txt a list of all fasta files to be used as subjects
#         fasta.txt fasta of all fasta formatted sequences
# BLAST Y or N (if N this will create a new matrix by re-parsing blast results with given cutoff value
# SAFE MODE Y or N (see README for definition of cutoff values).
########################################################


#################################################
my $path = "/usr/bin/"; # CHANGE THIS TO PATH OF BLAST EXECS
#################################################

# USERS SHOULD NOT CHANGE ANYTHING BELOW THIS LINE
##################################################

my $usage = "\nRD Emes Jan 2005\nRun a series of blastall searches using a fasta file supplied as a query against a series of fasta files which will be formatted using Formatdb. outputs a matrix of presence/absence from a database depending on the eval cut off defined in blast parsing stage.\n~~~~~~~~~~~~~~~~~~~~~~~~~~~\nINPUT\n~~~~~~~~~~~~~~~~~~~~~~~~~~~\nDB.TXT a list of all fasta files to be used as subjects\nDB TYPE P = peptide N = nucleotide\nQUERY FASTA fasta formatted query seqences\nQUERY TYPE P = peptide N = nucleotide\nEval cut (1e-3 = 0.001, 1e3 = 1000) \n~~~~~~~~~~~~~~~~~~~~~~~~~~~\nEXAMPLE: ./blast2matrix.pl db.txt N query.pep P 1e-3\n~~~~~~~~~~~~~~~~~~~~~~~\n";

unless (exists $ARGV[4])
{print "$usage\n";
exit;
}




my $dbs = shift or die "$usage";
my $dbtype = shift or die "$usage";
my $queryfasta = shift or die "$usage";
my $seqtype = shift or die "$usage";
my $cutoff = shift or die "$usage";
my $pidcut = 0;
my $lengthcut = 0;
my $cutprint = $cutoff;

my %queryfasta;
my @queryids = ();


unless (-e $dbs && -r $dbs)
{print "Error: $dbs does not exist or is not readable\n"; exit;}

unless ($dbtype eq "P" | $dbtype eq "N")
{print "Error: Db type ($dbtype) not recognised (N or P only)\n"; exit;}

unless (-e $queryfasta && -r $queryfasta)
{print "Error: $queryfasta does not exist or is not readable\n"; exit;}

unless ($seqtype eq "P" | $seqtype eq "N")
{print "Error: Sequence type ($seqtype) not recognised (N or P only)\n"; exit;}

unless ($cutoff =~ /\d+e-*\d+/)
{print "Error: Eval cutoff ($cutoff) not in the form 1e-3\n"; exit;}


$cutoff = EXPCONVERT($cutoff);

my $formatcommand = $path."formatdb";
my $blastcommand = $path."blastall";

print "Rename query sequences sequentially? [Y/N]\n";
my $rename = <STDIN>;
chomp $rename;

while ($rename !~ "Y|N"){
	print "Rename query sequences sequentially? [Y/N]\n";
	$rename = <STDIN>;
	chomp $rename;
	}


print "Run BLAST searches? [Y/N]\n";
my $runblast = <STDIN>;
chomp $runblast;

while ($runblast !~ "Y|N"){
	print "Run BLAST searches? [Y/N]\n";
	$runblast = <STDIN>;
	chomp $runblast;
	}

	
print "Run SAFE MODE? [Y/N]\n(hsp pecent identity and percent length cutoffs defined)\n";
my $safe = <STDIN>;
chomp $safe;

while ($safe !~ "Y|N"){
	print "Run SAFE MODE? [Y/N]\n(hsp pecent identity and percent length cutoffs defined)\n";
	$safe = <STDIN>;
	chomp $safe;
	}
	
if ($safe eq "Y")
	{
	print "hsp Percent identity cutoff [0-100]\n";
	$pidcut = <STDIN>;
	chomp $pidcut;

	while ($pidcut !~ '\d+')
		{
		print "hsp Percent identity cutoff [0-100]\n";
		$pidcut = <STDIN>;
		chomp $pidcut;
		}
	}

if ($safe eq "Y")
	{
	print "hsp length cutoff (percentage of query sequence length) [0-100]\n";
	$lengthcut = <STDIN>;
	chomp $lengthcut;

	while ($pidcut !~ '\d+')
		{
		print "hsp length cutoff (percentage of query sequence length) [0-100]\n";
		$lengthcut = <STDIN>;
		chomp $lengthcut;
		}
	}



#####################################
# determine the blast search to use.
my $blastflava = 0;

if ($dbtype eq "P" && $seqtype eq "P"){
    $blastflava = "blastp";
}
elsif ($dbtype eq "N" && $seqtype eq "P"){
    $blastflava = "tblastn";
    }
elsif ($dbtype eq "N" && $seqtype eq "N"){
    $blastflava = "blastn";
}
elsif ($dbtype eq "P" && $seqtype eq "N"){
    $blastflava = "blastx";
}
else {print "what are $dbtype and $seqtype? (UPPER CASE)\n";}



############################################################################
#MAKE HASH OF FASTA FILE OF QUERY SEQUENCES
if ($rename eq "Y")
{
open FASTA, $queryfasta;
open README, ">$queryfasta\.README";

{                           ### { }'S NEEDED TO DEFINE A SCOPE FOR LOCAL COMMAND
    local $/ = '>'; 
    <FASTA>;                                  # throw away the first line 'cos will only contain ">"
my $num = "a";
    while (<FASTA>) {	
	chomp $_;
	my ($queryseq_id, @querysequence) = split "\n";      # split the fasta input into Id and sequence
	my $queryfasta_sequence = join '',@querysequence;    # reassembles the sequence
	my $num_id = $num;
	print README "$num_id\t$queryseq_id\n";
	$queryseq_id = $num_id;
	$queryfasta{$queryseq_id} = $queryfasta_sequence;
	push (@queryids, $queryseq_id);
	$num++;
    }
}
close README;
close FASTA;
}

elsif ($rename eq "N")
{

open FASTA, $queryfasta;

{                           ### { }'S NEEDED TO DEFINE A SCOPE FOR LOCAL COMMAND
    local $/ = '>'; 
    <FASTA>;                                  # throw away the first line 'cos will only contain ">"
    while (<FASTA>) {	
	chomp $_;
	my ($queryseq_id, @querysequence) = split "\n";      # split the fasta input into Id and sequence
	my $queryfasta_sequence = join '',@querysequence;    # reassembles the sequence
	$queryfasta{$queryseq_id} = $queryfasta_sequence;
	push (@queryids, $queryseq_id);
    }
}
close FASTA;
}


my $querycount = scalar @queryids;

open FAIN, ">$queryfasta\.IN\.fa"; # file generated that is query fasta file

foreach  (@queryids){
    if (exists $queryfasta{$_}){
	print FAIN "\>$_\n", $queryfasta{$_}, "\n";
    }
}
close FAIN;

#####################################
# READ IN DB NAMES
my @dbnames = ();
open TXT, "$dbs";

while (<TXT>){
    chomp; 
    my $foobar = $_;
    push (@dbnames, $foobar);
}
close TXT;
 
#####################################################
# FORMAT DBs

my $dbid;    # elements in array @dbnames

if ($runblast eq "N") {print "Skipping database formating\n";}	

elsif ($runblast eq "Y"){
	my $foo = 1;
	foreach $dbid (@dbnames){
	
    	print "Formating database $foo\t\($dbid\)\n";
    
    	if ($dbtype eq 'P'){
		system "$formatcommand -pT -i $dbid";
    	}
    
    	elsif ($dbtype eq 'N'){
		system "$formatcommand -pF -i $dbid";
    	}
    	else {print "enter database format type N or P";
      	}
    
    	$foo++;
	}
}
print "\n";

####################################################
# RUN BLAST
if ($runblast eq "N"){
print "Skipping BLAST searches\n";
system "cp blast_results/* .";
}


elsif ($runblast eq "Y"){
my $fa_query = "$queryfasta\.IN\.fa";
my $bar = 1;
foreach $dbid (@dbnames){

	print "Blasting database $bar\t\($dbid\)\n"; 
	system "$blastcommand -p $blastflava -d $dbid -i $fa_query -o $dbid\.blastout";
	$bar++;
	}
}
print "\n";


##########################################
# PARSE BLAST FOR USERS RECORDS

my $ken = 1;

foreach $dbid (@dbnames){
    print  "Parsing blast result $ken\t\($dbid\)\n"; 
    my $tmpid = "$dbid\.blastout";
    my $pout = BP($tmpid);
    $ken++;
}
print "\n";


##################################################
#INITIATE MATRIX OF RESULTS

open MATRIX, ">matrix\.$cutprint\.tmp";            # open file for writing results

foreach my $id (@queryids){
    print MATRIX "\t", $id;
}
print MATRIX "\n";

##################################################

##################################################
#PARSE PARSED BLAST RESULTS FOR MATRIX
print "Parsing Results into matrix\n";

foreach (@dbnames){
    
    my $in = $_;
    
    print MATRIX $in, "\t";
    open QHITSIN, "$in\.blastout\.$cutprint\.parsed";
    
    my $holding = "";
    my @holding = (); # will hold all queries
 
    while (<QHITSIN>){
	chomp;
	if ($_ =~ /(\w+)\t.*/){
	    push (@holding, $1);
	}

	else {print MATRIX "\n";}
    }	
	

    $holding = join "\n", @holding;
    my $out = PARSE4MATRIX($holding);
    print MATRIX "\n";
}


close MATRIX; 



###################################################
# FLIP MATRIX
#
my $flippedin = "matrix\.$cutprint\.tmp";
my $flipped = FLIPMATRIX($flippedin);


##################################################
#CLEAN UP HARD DRIVE
if ($runblast eq "Y"){
system "mkdir -p blast_results";
system "mkdir -p parsed_results";
}
system "rm -rf *.nsq *.nin *.nhr";
system "rm -rf *.psq *.pin *.phr";
system "mv *.blastout blast_results";
system "mv *.parsed parsed_results";
system "rm -rf $queryfasta\.IN\.fa";
print "Run complete\n";
exit;

###############################################
# SUBROUTINES #################################


 sub BP{

    use Bio::SearchIO;
    use Bio::AlignIO;
    
    
    
    my $in = shift or die;
    
    my $searchio = new Bio::SearchIO( -format => 'blast',
				      -file   => "$in" );
    open BLASTOUT, ">>$in\.$cutprint\.parsed";
    print BLASTOUT "Queryname\tPercent_identity\tEval\tBitScore\tQuerySt\tQueryEnd\tHitname\tHitSt\tHitEnd\n";
    my $countplus=0; #keeps track of the num of + strand HSPs
    my $countminus=0;
    my $qname="";
    my $querylen=0;
    my $hitname="";
    my $dbname="";
    my $hitdb="";
    my $hspID=0;
    
    while ( my $result = $searchio->next_result() ) {
	
	$qname = $result->query_name();
	$querylen = $result->query_length();
	$dbname = $result->database_name();
	
	while( my $hit = $result->next_hit ) {
	    $hitname = $hit->name();

	    while( my $hsp = $hit->next_hsp ) {
		
		my $t_al = $hsp->hit_string;
		my $q_al = $hsp->query_string;
		my $hom_str = $hsp->homology_string;
		my $hsp_len = $hsp->length('total');
		my $bitsc = $hsp->bits;
		my $pid = sprintf("%2.2d",$hsp->percent_identity);
		my $eval = $hsp->evalue;
		$eval =~ s/^e/1e/;
		my $frame = $hsp->frame;
		my $strand = $hsp->strand('sbjct');
		my ($qst,$qend) = $hsp->range('query');
		my ($tst,$tend) = $hsp->range('target');
		#my $testtwo = $qend-$qst;
		my $percentlength = ($hsp_len/$querylen)*100;
		
		if ($safe eq "Y")
			{
			if ($blastflava eq "blastn") 
				{
				if ($cutoff >= $eval && $percentlength >= $lengthcut && $pid >= $pidcut)
					{
		    			print BLASTOUT "$qname\t$pid\t$eval\t$bitsc\t$qst\t$qend\t$hitname\t$tst\t$tend\n";
					}
				}
			elsif ($blastflava ne "blastn") 
				{
				if ($cutoff >= $eval && $percentlength >= $lengthcut && $pid >= $pidcut)
					{
		    			print BLASTOUT "$qname\t$pid\t$eval\t$bitsc\t$qst\t$qend\t$hitname\t$tst\t$tend\n";
					}
				}
			}	
		
		elsif ($safe eq "N")
			{
			if ($cutoff >= $eval){
		    	print BLASTOUT "$qname\t$pid\t$eval\t$bitsc\t$qst\t$qend\t$hitname\t$tst\t$tend\n";
			}
		}
		
	    } #closes while hsp
	} #closes while hit
    } #closes while result
close BLASTOUT;    
}

##########################################################

sub PARSE4MATRIX{
    my $infile = shift or die "sub routine no input sub PARSE4MATRIX";
    my @blast = split "\n", $infile;
    my %dbid;
    my $blasthit;
    
   
    my $num = 0;
    foreach (@blast){
	$blasthit = $_;
	$dbid{$blasthit} = $num;
	$num++;

    }
    
    foreach (@queryids){
	
	my $tmp = $_;
	if (exists $dbid{$tmp}){
	    print MATRIX "1\t";
	}
	else {print MATRIX "0\t"; 
	  }
    }
    
}

sub FLIPMATRIX{
my $matrix = shift or die;

my $columns = "$querycount"+2;
my $stop = 1;


open FLIP, ">>$queryfasta\.$blastflava\.$cutprint\.L$lengthcut\.PID$pidcut\.matrix";

until ($stop == $columns)
{
system "cut -f$stop $matrix > $stop\.tmp" ;
$stop++
}

my @holding = ();
open INCOL, "<1\.tmp";

while (<INCOL>){
chomp;
push (@holding, $_);
}
close INCOL;


foreach (@holding){
my $foo = $_;
print FLIP "$foo\t";
}
print FLIP "\n";

my $num = 2;
until ($num == $columns)
{
my @hold = ();
open INCOL, "<$num\.tmp";

while (<INCOL>){
chomp;
push (@hold, $_);
}
close INCOL;


foreach (@hold){
my $foo = $_;
print FLIP "$foo\t";
}
print FLIP "\n";
$num++
}



close FLIP;
system "rm -rf *.tmp";
}


sub EXPCONVERT
{
my $eval = shift or die;

my $integer = 0;
my $posneg = 0;
my $expon = 0;
my $out = 0;


if ($eval =~ /(\d+)e([-]*)(\d+)/)
	{$integer = $1;
	$posneg = $2;
	$expon = $3;
	}

if ($posneg eq "-"){
	my $zeros = ($expon - 1); 
	my $zeroprint = "0";
	$out = "0.".$zeroprint x $zeros.$integer
	} 
	
else 	{my $zeros = ($expon); 
	my $zeroprint = "0";
	$out = $integer.$zeroprint x $zeros
	} 
return $out;
}