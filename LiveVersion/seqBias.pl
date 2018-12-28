#!/usr/bin/env perl
#Package Version: 3.1

######################################################################################################
# Author(s): T.J.Cooper
# Updated: 20/12/2017
# Extracts genomic sequence surrounding user-specific peaks/features for sequence bias analysis
######################################################################################################

use strict;
use warnings;
use Cwd;
use Bio::SeqIO;
use Getopt::Long;
use File::Basename;
use List::Util qw(all);
my ($name, $path, $suffix) = fileparse($0);	#Obtain script-name and directory
my $outext = '.txt';	#Output .file-extension
my ($peaks, $fasta, $width, $threshold, $mode, $exclusion, %list, $chr, $min, $max, %sequences);
my $usage = "Usage: $name -i <InputFolder> -r <ReferenceFASTA> -w <Width> -m <Mode: Pos/Freq>\nOptional: -t <FilterThreshold> -e <ExclusionList>";		#Error/usage message
GetOptions('i=s' => \$peaks,
			  'r=s' => \$fasta,
			  'e=s' => \$exclusion,
			  'w=i' => \$width,
			  't=i' => \$threshold,
			  'm=s' => \$mode) or die("\n$usage\n");
die("\nError: Arguments or -flags are missing and/or incorrectly specified.\n\n$usage\n\n") unless all {defined} $peaks, $fasta, $width, $mode;
my @files = glob($peaks."/*.txt");
my $chk = scalar(@files);
print "\nFailed to detect any valid files within the current directory.\n\n" if $chk == 0;
exit if $chk == 0;
chdir($peaks) or die "$!";
my $sub = cwd()."/Data";
my $sub2 = cwd()."/Plots";
mkdir("$sub");
mkdir("$sub2");
if (defined $exclusion) {
	open my $IN, '<', $exclusion or die "$!";
	while (<$IN>) {
		chomp $_;
		($chr, $min, $max) = split('\s+', $_);
		push @{$list{$chr}}, [$min,$max];
	}
	close $IN;
}
print "\n-------------------------------------\n";
print "Loading Genome Reference....\n";
print "-------------------------------------";
my $seqio = Bio::SeqIO->new(-file => $fasta);
while(my $seqobj = $seqio->next_seq) {
	my $id  = $seqobj->display_id;
	my $seq = $seqobj->seq;
	$sequences{$id} = $seq;
}
sub seqstore {
	my ($flankseq, @rcd) = @_;
	if ($mode eq "Pos") {
	  push (@$flankseq, $rcd[0]);
	} elsif ($mode eq "Freq") {
	  push (@$flankseq, $rcd[0]) for (1..$rcd[1]);
	}
	return;
}
print "\nCalculating Sequence Bias....\n";
print "-------------------------------------\n";
print "Currently processing:\n";
for my $file (@files) {
	open my $IN2, '<', $file or die "$!";
	(my $sample = basename($file)) =~ s/.txt//;
	$sample =~ s/FullMap.//;
	print "$sample\n";
	my (@flankseqW, @flankseqC, $seqEX);
	<$IN2> for (1..1);
	while (<$IN2>) {
		chomp $_;
		my(@F) = split("\t", $_);
		next if grep {$F[1] >= $_->[0] && $F[1] <= $_->[1]} @{ $list{$F[0]}};
		next if $F[1]-$width < 1;
		next if(defined ($threshold) && $F[2] < $threshold || defined ($threshold) && $F[3] < $threshold);
		if ($F[3] > 0) {
			$seqEX = substr($sequences{$F[0]},$F[1]-1-$width,$width) #Extracts Xbp upstream (not incl. 5' end)
			.substr($sequences{$F[0]},$F[1]-1,$width+1);
			$seqEX = uc($seqEX);
			$seqEX =~ tr/GATC/CTAG/;
			$seqEX = reverse($seqEX);
			seqstore(\@flankseqC,$seqEX,$F[3]);
		}
		if ($F[2] > 0) {
			$seqEX = substr($sequences{$F[0]},$F[1]-1-$width,$width) #Extracts Xbp upstream (not incl. 5' end)
			.substr($sequences{$F[0]},$F[1]-1,$width+1);
			$seqEX = uc($seqEX);
			seqstore(\@flankseqW,$seqEX,$F[2]);
		}
	}
	close $IN2;
	my (%freq, %freqW, %freqC, @charts);
	my @x = (-$width..$width);
	my $outfile = cwd()."/Data/"."SeqBias.Combined.".$sample."_"."$mode"."_".$threshold."_"."$width".$outext;
	my $outfile2 = cwd()."/Data/"."SeqBias.Watson.".$sample."_"."$mode"."_".$threshold."_"."$width".$outext;
	my $outfile3 = cwd()."/Data/"."SeqBias.Crick.".$sample."_"."$mode"."_".$threshold."_"."$width".$outext;
	my ($OUT, $OUT2, $OUT3);
	if (@flankseqW && @flankseqC) {
		my $incW = 1 / @flankseqW; my $incC = 1 / @flankseqC;
		open $OUT, '>', "$outfile" or die "$!";
		open $OUT2, '>', "$outfile2" or die "$!";
		open $OUT3, '>', "$outfile3" or die "$!";
		for my $FSW (@flankseqW) {
			$freqW{substr $FSW, $_, 1}[$_] += $incW for 0..length($FSW)-1;
			$freq{substr $FSW, $_, 1}[$_] += ($incW/2) for 0..length($FSW)-1;
		}
		for my $FSC (@flankseqC) {
			$freqC{substr $FSC, $_, 1}[$_] += $incC for 0..length($FSC)-1;
			$freq{substr $FSC, $_, 1}[$_] += ($incC/2) for 0..length($FSC)-1;
		}
		for my $ot ($OUT, $OUT2, $OUT3) {printf($ot "%s\t"."%s\t" x ((keys %freq)-1)."%s"."\n", "Pos",sort keys %freq)};
		for my $pos (1..length $flankseqW[0]) {
			printf($OUT "%d\t"."%7.4f\t" x ((keys %freq)-1)."%7.4f"."\n", $x[$pos-1], map{$freq{$_}[$pos-1] // 0} sort keys %freq);
			printf($OUT2 "%d\t"."%7.4f\t" x ((keys %freqW)-1)."%7.4f"."\n", $x[$pos-1], map{$freqW{$_}[$pos-1] // 0} sort keys %freqW);
			printf($OUT3 "%d\t"."%7.4f\t" x ((keys %freqC)-1)."%7.4f"."\n", $x[$pos-1], map{$freqC{$_}[$pos-1] // 0} sort keys %freqC);
		}
		close $OUT;
		close $OUT2;
		close $OUT3;
	} elsif (@flankseqW) {
		my $incW = 1 / @flankseqW;
		open $OUT2, '>', "$outfile2" or die "$!";
		for my $FSW (@flankseqW) {
			$freqW{substr $FSW, $_, 1}[$_] += $incW for 0..length($FSW)-1;
		}
		printf($OUT2 "%s\t"."%s\t" x ((keys %freqW)-1)."%s"."\n", "Pos",sort keys %freqW);
		for my $pos (1..length $flankseqW[0]) {
			printf($OUT2 "%d\t"."%7.4f\t" x ((keys %freqW)-1)."%7.4f"."\n", $x[$pos-1], map{$freqW{$_}[$pos-1] // 0} sort keys %freqW);
		}
		close $OUT2;
	} elsif (@flankseqC) {
		my $incC = 1 / @flankseqC;
		open $OUT3, '>', "$outfile3" or die "$!";
		for my $FSC (@flankseqC) {
			$freqC{substr $FSC, $_, 1}[$_] += $incC for 0..length($FSC)-1;
		}
		printf($OUT3 "%s\t"."%s\t" x ((keys %freqC)-1)."%s"."\n", "Pos",sort keys %freqC);
		for my $pos (1..length $flankseqC[0]) {
			printf($OUT3 "%d\t"."%7.4f\t" x ((keys %freqC)-1)."%7.4f"."\n", $x[$pos-1], map{$freqC{$_}[$pos-1] // 0} sort keys %freqC);
		}
		close $OUT3;
	}
	my $Routfile = $sample."_".$mode."_".$threshold."_".$width;
	my $RScript = $path."biasPlots.R";
	my $command = "Rscript '$RScript' $outfile2 $outfile3 $outfile $Routfile";
	system($command);
}
my $run_time = time() - $^T;
print "-------------------------------------";
print "\nAnalysis Complete\n";
print "Execution Time: $run_time Seconds\n";
print "-------------------------------------\n\n";
