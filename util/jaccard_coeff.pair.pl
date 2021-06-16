#!/usr/bin/perl -w

=head2 jaccard_coefficient.pair.pl
  
 this works with our current fingerprints, taking in a pair with identities
 output is a json structure with jaccard score + details on covarege
 
 ./jaccard_coefficient.pair.pl --id1 [identity1] --fin1 [.fin file 1] --id2 [identity1] --fin2 [.fin file 1] --use-uncovered [Optional flag to include locations with N flag]

=cut

use strict;
use Getopt::Long;
use constant DEBUG=>0;
use JSON::PP;

my($id1,$id2,$fin1,$fin2,$outfile,$use_uncovered,$short);
my $USAGE = "./jaccard_coefficient.pair.pl --id1 [identity1] --fin1 [.fin file 1] --id2 [identity1] --fin2 [.fin file 1] --use-uncovered [Optional flag to include locations with N flag]";
my $results = GetOptions("id1=s"      => \$id1,
                         "id2=s"      => \$id2,
                         "fin1=s"      => \$fin1,
                         "fin2=s"      => \$fin2,
                         "out=s"      => \$outfile,						 
                         "use-uncovered" => \$use_uncovered);

if ((!$id1 || !$id2 || !$fin1 || !$fin2)) {die $USAGE;}

my %finfiles=($fin1=>$id1,$fin2=>$id2);


my %snps   = ();
my %counts = ();

=head2
 check if the file exists
 and then use the values to populate snps hash
=cut

foreach my $file(keys %finfiles) {
 next if $file=~/^#/;
 chomp($file);
 if ( !-e $file || -S $file) {
   print STDERR "File[$file] is absent or empty or does not exist\n";
   die;
 }
 my @lines = `tail -n +2 $file | cut -f 3,5`;
 my $id = $finfiles{$file};

 ### this counts only sites with a variant, not all sites with coverage
 map{chomp;my @temp=split("\t");$snps{$id}->{$temp[0]} = $temp[1];$counts{$id}->{variant}++ if $temp[1]=~/[ACTG]/} @lines;
 ## this counts ll sites with coverage
 map{chomp;my @temp=split("\t");$snps{$id}->{$temp[0]} = $temp[1];$counts{$id}->{covered}++ if $temp[1]!~/N/} @lines;
 
 $counts{$id}->{total}=scalar @lines;

 
}

my ($jaccard_score,$match,$sites)=jaccard_coeff($snps{$id1},$snps{$id2});

$counts{jaccard_score}=$jaccard_score;
$counts{pair}{covered}=$sites;
$counts{pair}{match}=$match;

my $json=encode_json(\%counts);
print $json;

=head2 Coefficient function
 Depending on use_uncovered flag we calculate Jaccard coefficient
 as Intersect/Union
=cut

sub jaccard_coeff {
 my ($sample_one, $sample_two) = @_;
 my ($intersect,$union) = (0,0);

 ### go through the keys of sample_one
 foreach my $rs (keys %{$sample_one}) {
   print STDERR $rs."\n" if DEBUG;
   
   ### don't evaluate this site if it had not been interrogated in the second sample
   if(! $sample_two->{$rs}){
	   print STDERR "$rs not in sample_two" if DEBUG;
	   next;
   }
   
   
   
   if (($sample_one->{$rs} eq 'N' || $sample_two->{$rs} eq 'N') && !$use_uncovered) {
       print STDERR "Flag set to ignore uncovered and we have ".$sample_one->{$rs}." and ".$sample_two->{$rs}."\n" if DEBUG;
       next;
   }
   next if $sample_one->{$rs} eq 'N' && $sample_two->{$rs} eq 'N'; # Just skip these, we don't have any info here;

   

   ### union the total number of sites with coverage (or at least one covered, if use_uncovered is set)
   $union++;
   
   ### intersect : the number of sites with the same variant
   if ($sample_one->{$rs} eq $sample_two->{$rs}) {
     print STDERR "Have a match\n" if DEBUG;
     $intersect++;
   }
 }
 $union ||=1; # Avoid devision by zero
 
 my $score=sprintf "%.4f", $intersect/$union;
 return($score,$intersect,$union);
}
