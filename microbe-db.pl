use strict; use warnings;
use Getopt::Long;
use taxbuild;

my $dir;
my $line;
my @array;
my $taxDB;
my $i;
my $name;
my %tax;
my $flag;
my @value;
my $help;
my $nt;
my $output;
my $tax_path;
my %hash;

$help="USAGE:
 perl $0 --nt <Nt fasta> --out <output> --tax <accession2taxid file> [--dir <directory of *.dmp>]\n
ARGUMENTS:
  --nt <Nt fasta>
    the Nt fasta file downloaded from ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/, it can be in gzipped format
  --out <output>
    the output file
  --tax <accession2taxid file>
    the nucl_gb.accession2taxid file downloaded from ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/
  --dir <directory of *.dmp>
    the directory of names.dmp and nodes.dmp, dmp files are in taxdmp file downloaded from ftp://ftp.ncbi.nih.gov/pub/taxonomy
    default: current directory
  --help|--h
    the usage information\n";

GetOptions("nt=s"=>\$nt,"out=s"=>\$output,"tax=s"=>\$tax_path,"dir=s"=>\$dir,"help|h"=>\$flag);
if($flag)
{
	print $help;
	exit;
}
if(!$nt)
{
	print "Error! Please input the Nt fasta file with --nt\n";
	exit;
}
if(!$output)
{
	print "Error! Please input the output file with --out\n";
	exit;
}
if(!$tax_path)
{
	print "Error! Please input the accession2taxid file with --tax\n";
	exit;
}
if(!$dir)
{
	$dir=$ENV{'PWD'};
}
if($dir!~/\/$/)
{
	$dir=$dir."/";
}
$taxDB = taxbuild->new(nodes=>$dir."nodes.dmp",names=>$dir."names.dmp",save_mem=>0);
if($tax_path=~/gz$/)
{
	open(IN,"gzip -dc $tax_path |") or die "Can't open $tax_path file!\n";
}
else
{
	open(IN,"<$tax_path") or die "Can't open $tax_path file!\n";
}
while(<IN>)
{
	$line=$_;
	if($line=~/^accession/)
	{
		next;
	}
	@array=split("\t",$line);
	if(exists $tax{$array[2]})
	{
		$hash{$array[0]}=1;
		next;
	}
	$flag=$taxDB->get_term_at_level($array[2],"superkingdom");
	if($flag eq "Bacteria")
	{
		$tax{$array[2]}=1;
		$hash{$array[0]}=1;
		next;
	}
	if($flag eq "Archaea")
	{
		$tax{$array[2]}=1;
		$hash{$array[0]}=1;
		next;
	}
	if($flag eq "Viruses")
	{
		$tax{$array[2]}=1;
		$hash{$array[0]}=1;
		next;
	}
	$flag=$taxDB->get_term_at_level($array[0],"kingdom");
	if($flag eq "Fungi")
	{
		$tax{$array[2]}=1;
		$hash{$array[0]}=1;
	}
}
close IN;

if($nt=~/gz$/)
{
	open(IN,"gzip -dc $nt |") or die "Can't open $nt file!\n";
}
else
{
	open(IN,"<$nt") or die "Can't open $nt file!\n";
}
open(OUT,">$output") or die "Can't create $output file!\n";
while(<IN>)
{
	$line=$_;
	if($line=~/^\>/)
	{
		($line)=$line=~/^\>(.+)$/;
		$flag=0;
		@array=split("\x01",$line);
		for($i=0;$i<@array;$i++)
		{
			@value=split(" ",$array[$i]);
			if($value[0]!~/\.\d+$/)
			{
				next;
			}
			($name)=$value[0]=~/^(.+)\./;
			if(exists $hash{$name})
			{
				$flag++;
				last;
			}
		}
		if($flag==1)
		{
			print OUT "\>$name\n";
		}
		next;
	}
	if($flag==1)
	{
		print OUT $line;
	}
}
close IN;
close OUT;
