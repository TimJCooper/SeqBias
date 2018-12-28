#!/usr/bin/env perl
#Package Version: 1.3

######################################################################################################
# Author(s): T.J.Cooper
# Updated: 15/2/2017
# Extracts genomic sequence surrounding user-specific peaks/features for sequence composition analysis
######################################################################################################

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use File::Basename qw(basename);
use List::Util qw(all);
my $scriptname = basename($0);	#Obtain script-name
my ($peaks, $fasta, $width, $mode, $outfile);
my $usage = "Usage: $scriptname -i <HistogramFile> -r <ReferenceFASTA> -w <Width> -m <Mode: Pos/Freq> -o <Output>";		#Error/usage message
GetOptions('i=s' => \$peaks,
			  'r=s' => \$fasta,
			  'w=i' => \$width,
			  'm=s' => \$mode,
			  'o=s' => \$outfile) or die("\n$usage\n");
die("\nError: Arguments or -flags are missing and/or incorrectly specified.\n\n$usage\n\n") unless all {defined} $peaks, $fasta, $width, $mode, $outfile;



open my $IN2, '<', $peaks or die "$!";
open my $OUT, '>', $outfile or die "$!";	#Create mapping-coordinate output file
my (@flankseq, $seqEX, %sequences);
my $seqio = Bio::SeqIO->new(-file => $fasta);
while(my $seqobj = $seqio->next_seq) {
	my $id  = $seqobj->display_id;
	my $seq = $seqobj->seq;
	$sequences{$id} = $seq;
}
sub seqstore {
	my @rcd = @_;
	if ($mode eq "Pos") {
		push (@flankseq, $rcd[0])
	} elsif ($mode eq "Freq") {
		push (@flankseq, $rcd[0]) for (1..$rcd[1]);
	}
	return;
}
<$IN2> for (1..1);
while (<$IN2>) {
	chomp $_;
	my(@F) = split("\t", $_);
	next if $F[1]-$width < 1;
	if ($F[3] > 0) {
		$seqEX = substr($sequences{$F[0]},$F[1]-1-$width,$width) #Extracts Xbp upstream (not incl. 5' end)
		.substr($sequences{$F[0]},$F[1]-1,$width+1);
		$seqEX =~ tr/GATC/CTAG/;
		$seqEX = reverse($seqEX);
		seqstore($seqEX,$F[3]);
	}
	if ($F[2] > 0) {
		$seqEX = substr($sequences{$F[0]},$F[1]-1-$width,$width) #Extracts Xbp upstream (not incl. 5' end)
		.substr($sequences{$F[0]},$F[1]-1,$width+1);
		seqstore($seqEX,$F[2]);
	}
}
my %freq;
my $inc = 1 / @flankseq;
for my $FS (@flankseq) {
	$freq{substr $FS, $_, 1}[$_] += $inc for 0 .. length($FS)-1;
}
printf $OUT " " . "%5s  " x (keys %freq) . "\n", sort keys %freq;
my $label = -$width-1;
for my $pos ( 1..length $flankseq[0] ) {
	$label++;
	printf $OUT "%1d" . "%7.4f" x (keys %freq) . "\n",
	$label, map { $freq{$_}[$pos - 1] // 0  } sort keys %freq;
}
my $run_time = time() - $^T;
print "\n-------------------------------------";
print "\nAnalysis Complete\n";
print "Execution Time: $run_time Seconds\n";
print "-------------------------------------\n\n";
