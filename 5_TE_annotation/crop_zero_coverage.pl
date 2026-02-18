#!/usr/bin/perl -w
#
use strict;

chomp (my $length = `seqkit fx2tab --length $ARGV[1] | cut -f4`);
chomp (my $seq = `seqkit fx2tab --length $ARGV[1] | cut -f2`);

open IN,"$ARGV[0]" or die "Cannot open infile: $ARGV[0]\n";

my @coverage = (0)x$length;

while (my $row = <IN>) {
	chomp $row;
	my @line = split "\t",$row;
	foreach my $i ($line[6]..$line[7]) {
		$coverage[$i]++;
	}
}

close IN;

my $min = 9999999999;
my $max = 0;
foreach my $i (0..@coverage-1) {
	$max = $coverage[$i] if $coverage[$i] > $max;
	$min = $coverage[$i] if $coverage[$i] < $min;
}

printf "@coverage\n";
printf "$seq\n";
printf "$length\n";
printf "$min - $max\n";

my $threshold = $max*0.1;

printf "$threshold\n";
my @chars = split "",$seq;

# Shift from the left
foreach my $i (0..$length-1) {
	if ($coverage[$i] < $threshold) {
		shift @chars;
	} else {
		last;
	}
}

# Pop from the right
foreach my $j (0..$length-1) {
	my $i = $length-$j-1;
	if ($coverage[$i] < $threshold) {
		pop @chars;
	} else {
		last;
	}
}

printf join "",@chars;
printf "\n";

my $name = $ARGV[1];
$name =~ s/.maf.fa_consensus.fasta//;

open OUT,">$name.cons.fa";
printf OUT ">$name\n".(join "",@chars)."\n";
close OUT;
