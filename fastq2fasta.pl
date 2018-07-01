#!/usr/bin/perl

my $fastqFile = $ARGV[0];
my $outFile = $ARGV[1];

# variables
my $line = "";
my $index = 1;
my $header = "";
my @data = ();

# open output files
open OUT, ">$outFile" || die "could not open $outFile!";

# open input file
open IN, "<$fastqFile" || die "could not open $fastqFile!";
while ( $line = <IN> ) {

	chomp($line);

	# @NS500125:306:HJ2W7BGXY:1:11101:10394:1056 1:N:0:ACAGTG
	# ANAGATAACGCAGGTGTCCTAAGATGAGCTCAACGAGAACA
	# +
	# A#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
	# @NS500125:317:HG22YAFXX:1:11303:23520:1038 2:N:0:CGATGT
	# CTCACTAAATCATNC
	# +
	# EEE/EEEEEEEEE#E

	if ( 1 == $index && $line =~ /^@/ ) {
		$header = ">" . $line;
		# increment index
		$index++;
	}
	elsif ( 2 == $index ) {
		$seq = $line;
		# increment counter
		$index++;
	}
	elsif ( 3 == $index ) {
		# increment index
		$index++;
	}
	elsif ( 4 == $index ) {
		# write data
		print OUT $header . "\n";
		print OUT $seq . "\n";

		# reset index and variables
		$header = "";
		$seq = "";
		$index = 1;
	}
}

close IN;

close OUT;


