#!usr/bin/perl
# This script takes an output file from nbinom.R and performs BH-FDR adjustment.
# usage: perl adjust_p_value.pl --in=[input file name] --out=[output file name] --fdr=[false discovery rate]

use strict;
use warnings;
use Getopt::Long;

my $in = "";
my $out = "";
my $fdr = 0.05; #default
GetOptions ("in=s" => \$in, "out=s" => \$out, "fdr=s" => \$fdr);

open IN, "<$in" or die "$!";
my @input = <IN>;
close IN;
chomp @input;

open OUT, ">$out" or die "$!";
print OUT "transcript_name_coordinate\tp-value\tadjusted_p-value\n";

my %hash = ();

for (my $i = 0; $i < @input; $i++)
{
	my @line = split(/\t/, $input[$i]);
	my $label = $line[0];
	my $p_value = $line[1];
	$hash{$label} = $p_value;
}

my @key = keys %hash;
# Sort the hash (transcriptName_coordinate => probability) by probability
# in ascending order.
my @sorted_key = sort { $hash{$a} <=> $hash{$b} } @key;
my $total = scalar @sorted_key;
my $threshold = 0;

# Go over the sorted hash from the one with the largest probability.
# Compute adjusted p-values.
# Store the index in the $threshould variable when the p-value gets below 0.05 or the specified FDR.
LINE1: for (my $i = $total; $i > 0; $i--)
{
	my $adjusted = $hash{$sorted_key[$i-1]} * ($total/($i));
	if ($adjusted < $fdr)
	{
		$threshold = $i-1;
		last LINE1;
	}
}

# Record the transcriptName_coordinate, probablity and adjusted p-values
# in a tab-delimited file.
for (my $i = 0; $i <= $threshold; $i++)
{
	my $adjusted = $hash{$sorted_key[$i]} * ($total/($i+1));
	print OUT "$sorted_key[$i]\t$hash{$sorted_key[$i]}\t$adjusted\n";
}

close OUT;
	

