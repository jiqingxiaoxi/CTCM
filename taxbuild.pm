#!/usr/bin/perl

package taxbuild;

=head1 NAME

taxbuild - Get the taxonomy of a gene, protein or taxId

=head1 SYNOPSIS

  use taxbuild;

  my $taxDB = taxbuild->new(
                            nodes => $nodesFile,
                            names => $namesFile,
                            dict  => $dictFile,
                            save_mem => 0
                           );

  # Get the taxonomy given a GI identifier
  my @tax = $taxDB->get_taxonomy_from_gi("35961124");

  # Get the taxonomy term of a GI identifier at a given level
  my $term_at_level = $taxDB->get_term_at_level_from_gi("35961124","family");

  # Get the taxid of a GI identifier
  my $taxid = $taxDB->get_taxid("35961124");

  # Get the taxonomy given a taxid
  my @tax = $taxDB->get_taxonomy($taxid);

  # Get the taxonomy at a given level given a taxid
  my $taxid_at_level = $taxDB->get_term_at_level($taxid,"genus");

  # Get the level of a given taxonomical name
  my $level = $taxDB->get_level_from_name("Proteobacteria");

=head1 DESCRIPTION

This module has been designated to easily retrieve the taxonomy of a gene fast and with low memory usage.
It parses the taxonomy db that can be downloaded from NCBI at the following address: L<ftp://ftp.ncbi.nih.gov/pub/taxonomy/>
You need to pass to the constructor (the "new" method) the following:

=over 4

=item * The I<nodes> file of the taxonomy database. Both the name of the file or its filehandle if the file is already open are allowed. This file is located at I<taxdump.tar.gz>

=item * The I<names> file of the taxonomy database. Both the name of the file or its filehandle if the file is already open are allowed. This file is located at I<taxdump.tar.gz>

=item * (optional) The I<dict> file containing the correspondences between B<gi> and B<taxIds>. This file can be downloaded from the taxonomy database too (it is called I<gi_taxid_prot.dmp.gz>), but you can't use it directly with this module, instead you need to convert it to binary format. This conversion improves dramatically the speed of the script and saves a lot of memory. There should be an accompanion script with this module called I<tax2bin.pl> that can do this task.

=item * (optional) I<save_mem>. If you want to save some memory, the correspondencies between B<gi> and B<taxIds> will not be loaded into memory. This will save aprox 800Mb of system memory, but looking up for a taxId will be ~20% slower. Note that this option is I<off> by default and only has sense if you are using the I<dict> option too.

=back

=head1 METHODS

Currently this module provides the following methods:

=head2 C<get_taxid>

Given a I<gi> number, returns its taxId

  my $taxid = $taxDB->get_taxid($gi);

=head2 C<get_taxonomy_from_gi>

Given a I<gi> number, returns an array with its taxonomy

  my @taxonomy = $taxDB->get_taxonmy_from_gi($gi);
  my $taxonomy = join ";",@taxonomy;

=head2 C<get_taxonomy>

Given a taxId number, returns an array with its taxonomy

  my @taxonomy = $taxDB->get_taxonomy($taxid);

=head2 C<get_term_at_level_from_gi>

Returns the taxonomy term of the gene at a given level

  my $family = $taxDB->get_term_at_level($gi,"family");
  $families{$family++};

=head2 C<get_term_at_level>

Returns the parent term of a given taxonomy given the taxId and the level

If the term is not present at the taxonomy database an error is thrown.
At the moment of this writing, These are the taxonomy levels present at the taxonomy database:

C<class>
C<family>
C<forma>
C<genus>
C<infraclass>
C<infraorder>
C<kingdom>
C<order>
C<parvorder>
C<phylum>
C<species>
C<species group>
C<species subgroup>
C<subclass>
C<subfamily>
C<subgenus>
C<subkingdom>
C<suborder>
C<subphylum>
C<subspecies>
C<subtribe>
C<superclass>
C<superfamily>
C<superkingdom>
C<superorder>
C<superphylum>
C<tribe>
C<varietas>

But note that not all of them are present in all taxonomies.

=head2 C<get_level_from_name>

Given a taxonomical name, get its level

  my $level = $taxDB->get_level_from_name($name);

=head2 C<get_taxid_from_name>

Given a taxonomical name, get its taxId

  my $id = $taxDB->get_taxid_from_name($name);


=head1 AUTHOR

Miguel Pignatelli Moreno

Any comments should be addressed to: miguel.pignatelli@uv.es

=head1 LICENSE

Copyright 2008 Miguel Pignatelli, all rights reserved.

This library is free software; you may redistribute it and/or modify it under the same terms as Perl itself.

=cut

use strict;
use warnings;
use Carp qw/croak/;
use Data::Dumper;

our $VERSION = 0.04;

use constant FS => '\t\|\t';
use constant RS => '\t\|\n';

our %allowed_levels;

sub _check_level
  {
    my ($self, $level) = @_;
    croak "Level not defined" unless defined $level;;
    return $allowed_levels{$level};
  }
sub _print_levels
  {
    my ($self) = @_;
    print STDERR "$_\n" for sort keys %allowed_levels;
  }

sub new
  {
    my ($class, %args) = @_;
    my %opts;

    $args{'nodes'} or croak "Need the nodes.dmp file";
    $args{'names'} or croak "Need the names.dmp file";

    @opts{qw /nodesFile namesFile/} = @args{qw/nodes names/};
    my $save_mem = $args{'save_mem'} || 0;

    my $dictFile;
    if ($args{'dict'}) {
      $dictFile = $args{'dict'};
      croak "\nERROR\n$dictFile: File not found\n" unless -e $dictFile;
      croak "\nERROR\nThe file containing the gi <-> taxid correspondences must be converted to binary format.\nThis will increase dramatically the speed of this script.\nTo convert that file to binary use the tax2bin.pl script like:\n\nperl tax2bin.pl $dictFile > $dictFile.bin\n\nand use the resulting file instead\n" unless (-B $dictFile);
    }
    $opts{dict} = $dictFile;
    $opts{save_mem} = $save_mem;
    my $self = bless \%opts;

    $self -> _build_taxonomy();
    $self -> _name_nodes();
    $self -> _build_dict() if (defined $dictFile);
    return $self;
  }

sub get_taxonomy_from_gi
  {
    my ($self, $gi) = @_;
    my $taxid = $self->get_taxid ($gi);
    my @tax = $self->get_taxonomy ($taxid);
    return @tax;
  }

sub get_taxid
  {
    my ($self, $gi) = @_;
#    return $self->_binary_lookup ($gi);
    return $self->_direct_lookup ($gi);
  }

sub get_term_at_level_from_gi
  {
    my ($self, $gi, $level) = @_;
    do {
      print STDERR "Level $level not recognized\nAllowed levels:\n";
      $self->_print_levels;
      croak;
    } if (! defined $self->_check_level($level));
    my $taxid = $self->get_taxid($gi);
    return $self->get_term_at_level($taxid,$level);
  }

sub _build_dict
  {
    my ($self) = @_;
    my $dictFile = $self->{dict};
    my $data;
    open my $gi_FH, '<:raw', $dictFile or croak "$!:$dictFile";
    if ($self->{save_mem}){
      $self->{fh} = $gi_FH;
    } else {
      sysread( $gi_FH, $data, -s( $dictFile ) ) or croak $!;
      close $gi_FH;
      $self->{dict} = $data;
    }
  }

sub _build_taxonomy
  {
    my ($self) = @_;
    my $nodesFile = $self->{nodesFile};
    my $tax;
    if ((UNIVERSAL::isa($nodesFile, 'GLOB')) or (ref \$nodesFile eq 'GLOB')) {
      $tax = $nodesFile;
    } else {
      open $tax, "<", $nodesFile or croak $!;
    }
    while (<$tax>){
      chomp;
      _create_node (_parse_tax_rec($_));
    }
    close $tax;
  }

{
  my %nodes;

  sub _create_node
    {
      my ($node,$parent,$level) = @_;
      $allowed_levels{$level} = 1 if (! defined $allowed_levels{$level});
      @{$nodes{$node}}{qw/parent level/} = ($parent,$level);
    }

  sub _name_nodes
    {
      my ($self) = @_;
      my $namesFile = $self->{namesFile};
      my $nodesNames;
      if ((UNIVERSAL::isa($namesFile, 'GLOB')) or (ref \$namesFile eq 'GLOB')) {
	$nodesNames = $namesFile;
      } else {
	open $nodesNames, "<", $namesFile or croak $!;
      }
      while (<$nodesNames>){
	chomp;
	my ($taxId,$taxName,$comment) = _process_tax_name ($_);
	if ($comment eq "scientific name"){
	  ${$nodes{$taxId}}{name} = $taxName;
	}
      }
      close $nodesNames;
    }

  sub get_term_at_level
    {
      my ($self,$taxid,$level) = @_;
      do {
	print STDERR "Level $level not recognized\nAllowed levels:\n";
	$self->_print_levels;
	croak;
      } if (! defined $self->_check_level($level));
      return "" unless defined ${$nodes{$taxid}}{name};
      while (${$nodes{$taxid}}{name} ne "root"){
	return ${$nodes{$taxid}}{name} if (${$nodes{$taxid}}{level} eq $level);
	$taxid = ${$nodes{$taxid}}{parent};
      }
      return "undef";
    }

  sub get_taxonomy
    {
      my ($self, $taxid) = @_;
      return "" unless defined ${$nodes{$taxid}}{name};
      my @taxonomy;
      while (${$nodes{$taxid}}{name} ne "root"){
	push @taxonomy, ${$nodes{$taxid}}{name};
	$taxid = ${$nodes{$taxid}}{parent};
      }
      return reverse do{pop @taxonomy;@taxonomy};
    }

  sub get_taxonomy_with_levels
    {
      my ($self,$taxid) = @_;
      return "" unless defined ${$nodes{$taxid}}{name};
      my @taxonomy;
      while (${$nodes{$taxid}}{name} ne "root"){
	push @taxonomy, [${$nodes{$taxid}}{name},${$nodes{$taxid}}{level}];
	$taxid = ${$nodes{$taxid}}{parent};
      }
      return reverse do{pop @taxonomy;@taxonomy};
    }

  sub get_level_from_name
    {
      my ($self,$name) = @_;
      for (keys %nodes) {
	return ${$nodes{$_}}{level} if (${$nodes{$_}}{name} eq $name);
      }
      return undef;
    }

  sub get_taxid_from_name
    {
      my ($self,$name) = @_;
      for (keys %nodes){
	return $_ if (${$nodes{$_}}{name} eq $name)
      }
      return undef;
    }

  sub get_taxonomy_from_name
    {
      my ($self,$name) = @_;
      my $taxid = $self->get_taxid_from_name($name);
      return $self->get_taxonomy($taxid);
    }

}

sub _parse_tax_rec
{
    my $line = shift @_;
    return (split FS,$line)[0,1,2];
}


sub _process_tax_name
  {
    my $line = shift @_;
    my @fields = split FS, $line;
    $fields[3] =~ s/\t\|$//;
    return ($fields[0],$fields[1],$fields[3]);
  }

sub _binary_lookup {
  my ($self, $gi) = @_;
  my $target = pack 'N', $gi;
  my( $left, $right ) = ( 0, ( length( $self->{dict} ) ) / 8 );
  while( $left < $right ) {
    my $mid = int( ( $left + $right ) / 2 );
    my $key = substr $self->{dict}, $mid * 8, 4;
    if( $key lt $target ) {
      $left = $mid +1;
    }
    elsif( $key gt $target ) {
      $right = $mid;
    }
    elsif( $key eq $target ) {
      my( $key, $val ) = unpack 'NN', substr $self->{dict}, $mid * 8, 8;
      return $val;
    }
  }
}

sub _direct_lookup {
  my ($self,$gi) = @_;
  if ($self->{save_mem}){
    my $taxid;
    sysseek ($self->{fh},$gi*4,0);
    sysread($self->{fh},$taxid,4,);
    return (unpack "N",$taxid);
  } else {
    return (unpack "N",substr($self->{dict},$gi*4,4));
  }
}

1;
