use strict; use warnings;
use Getopt::Long;

my $line;
my @array;
my $flag;
my %score;
my %turn;
my @total;
my $i;
my $help;
my @file;
my $check=0;
my $coral_path;
my $sym_path;
my $microbe_path;
my $output;
my $threshold_coral;
my $threshold_sym;

$help="USAGE:
  perl $0 --coral <Blastn alignments>  --sym <Blastn alignments>  --microbe <Megablast alignments>  --out <output> [options]*\n

ARGUMENTS:
  --coral <Blastn alignments>
    the transcripts alignment files with coral database
    files are separated by comma
  --sym <Blastn alignments>
    the transcripts alignment files with symbiodinium database
    files are separated by comma
  --microbe <Megablast alignments>
    the transcripts alignment files with microbe database
    files are separated by comma
  --out <output>
    output file
  --score_coral <score threshold>
    the threshold of alignment score for classifying coral
    the bigger of the threshold, the lower of sensitivity and the higher of specificity
    default: 45
  --score_sym <score threshold>
    the threshold of alignment score for classifying symbiodinium
    the bigger of the threshold, the lower of sensitivity and the higher of specificity
    default: 51
  --help|--h
    print help information\n";

$check=GetOptions("coral=s"=>\$coral_path,"sym=s"=>\$sym_path,"microbe=s"=>\$microbe_path,"out=s"=>\$output,"score_coral=i"=>\$threshold_coral,"score_sym=i"=>\$threshold_sym,"help|h"=>\$flag);
if($check==0||$flag)
{
	print $help;
	exit;
}
if(!$coral_path)
{
	print "Error! Please input the alignments of coral with --coral\n";
	exit;
}
if(!$sym_path)
{
	print "Error! Please input the alignments of symbiodinium with --sym\n";
	exit;
}
if(!$microbe_path)
{
	print "Error! Please input the alignments of microbe with --microbe\n";
	exit;
}
if(!$output)
{
	print "Error! Please input the name of output file with --out\n";
	exit;
}
if(!$threshold_coral)
{
	$threshold_coral=45;
}
if(!$threshold_sym)
{
	$threshold_sym=51;
}

###symbiodinium
@file=split(",",$sym_path);
for($i=0;$i<@file;$i++)
{
	if($file[$i]=~/gz$/)
	{
        	open(IN,"gzip -dc $file[$i] |") or die "Can't open $file[$i] file!\n";
	}
	else
	{
	        open(IN,"<$file[$i]") or die "Can't open $file[$i] file!\n";
	}
	while(<IN>)
	{
		chomp;
	        $line=$_;
	        if($line=~/^\#/)
	        {
	                $flag=0;
	                next;
	        }
	        if($flag==1)
	        {
	                next;
	        }
	        $flag=1;
	        @array=split("\t",$line);
		if(exists $score{$array[0]})
		{
			if($score{$array[0]}<$array[11])
			{
				$score{$array[0]}=$array[11];
			}
		}
		else
		{
			$turn{$array[0]}=2;
			$score{$array[0]}=$array[11];
		}
	}
	close IN;
}

##coral
@file=split(",",$coral_path);
for($i=0;$i<@file;$i++)
{
	if($file[$i]=~/gz$/)
	{
		open(IN,"gzip -dc $file[$i] |") or die "Can't open $file[$i] file!\n";
	}
	else
	{
		open(IN,"<$file[$i]") or die "Can't open $file[$i] file!\n";
	}
	while(<IN>)
	{
		chomp;
	        $line=$_;
	        if($line=~/^\#/)
	        {
	                $flag=0;
	                next;
	        }
	        if($flag==1)
	        {
	                next;
	        }
	        $flag=1;
	        @array=split("\t",$line);
		if(exists $score{$array[0]})
		{
			if($array[11]>$score{$array[0]})
			{
				$score{$array[0]}=$array[11];
				$turn{$array[0]}=1;
			}
		}
		else
		{
			$score{$array[0]}=$array[11];
			$turn{$array[0]}=1;
		}
	}
	close IN;
}

##microbe
@file=split(",",$microbe_path);
for($i=0;$i<@file;$i++)
{
	if($file[$i]=~/gz$/)
	{
		open(IN,"gzip -dc $file[$i] |") or die "Can't open $file[$i] file!\n";
	}
	else
	{
		open(IN,"<$file[$i]") or die "Can't open $file[$i] file!\n";
	}
	while(<IN>)
	{
		chomp;
		$line=$_;
		if($line=~/^\#/)
		{
			$flag=0;
			next;
		}
		if($flag==1)
		{
			next;
		}
		$flag=1;
		@array=split("\t",$line);
		if(exists $score{$array[0]})
		{
			if($array[11]>=$score{$array[0]})
			{
				delete($score{$array[0]});
			}
		}
	}
	close IN;
}

open(OUT,">$output") or die "Can't create $output file!\n";
print OUT "Transcripts\tSource\n";
$total[0]=0;
$total[1]=0;
foreach $line (keys %score)
{
	if($turn{$line}==1)
	{
		if($score{$line}>=$threshold_coral)
		{
			print OUT "$line\tCoral\n";
			$total[0]++;
			next;
		}
		next;
	}
	if($score{$line}>=$threshold_sym)
	{
		print OUT "$line\tSymbiodinium\n";
		$total[1]++;
	}
}
close OUT;
print "There are $total[0] transcripts of coral, $total[1] transcripts of symbiodinium.\n";
