#!/usr/bin/perl

my $readsLengthFile = $ARGV[0];
my $blatOutputFile = $ARGV[1];
my $bestBlatHitsOutput = $ARGV[2];

# variables
my $line = "";
my %data;
my %length;
my %outHash;

open IN, "<$readsLengthFile" || die "could not open $readsLengthFile!";
while ( $line = <IN> ) {

	chomp($line);

	# @NS500125:306:HJ2W7BGXY:1:11101:10394:1056      41

	my @tmp = split(/\s+/,$line);
	my $readID = $tmp[0];
	my $read_length = $tmp[1];

	$length{$readID} = $read_length;
}

close IN;

my $id = "NaN";

open OUT, ">$bestBlatHitsOutput" || die "could not open $bestBlatHitsOutput!";

open IN, "<$blatOutputFile" || die "could not open $blatOutputFile!";
while ( $line = <IN> ) {

	chomp($line);

	#   0     						 1               2       3       4       5        6       7        8       9        10     11
	# Query id						 Subject id      %id aln.len mismatches gaps  q.strt   q.end   s.start  s.end    e-value score
	# @NS500125:306:HJ2W7BGXY:1:11101:16803:1062      4       100.00  34      0       0       1       34      5288493 5288460 9.9e-12 68.0
	# @NS500125:306:HJ2W7BGXY:1:11101:16803:1062      10      100.00  34      0       0       1       34      12839459        12839426        9.9e-12 68.0
	# @NS500125:306:HJ2W7BGXY:1:11101:16803:1062      10      100.00  34      0       0       1       34      9853457 9853424 9.9e-12 68.0
	# @NS500125:306:HJ2W7BGXY:1:11101:16803:1062      8       97.06   34      1       0       1       34      8385467 8385500 5.6e-11 66.0
	# @NS500125:306:HJ2W7BGXY:1:11101:16803:1062      8       97.06   34      1       0       1       34      6120310 6120277 5.6e-11 66.0
	# @NS500125:306:HJ2W7BGXY:1:11101:16803:1062      8       97.06   34      1       0       1       34      8384666 8384699 2.1e-10 64.0
	# @NS500125:306:HJ2W7BGXY:1:11101:16803:1062      9       97.06   34      1       0       1       34      4449729 4449696 2.1e-10 64.0

	my @tmp = split(/\s+/,$line);
	my $readID = $tmp[0];
	my $alignment_length = $tmp[3];
	my $score = $tmp[11];
	my $evalue = $tmp[10];
	
	my $coverage = ($alignment_length/$length{$readID});
	
	if ( exists $data{$readID}{"coverage"} ) {
		if ( $coverage > $data{$readID}{"coverage"} ) {
			$data{$readID}{"coverage"} = $coverage;
			$data{$readID}{"evalue"} = $evalue;
			$data{$readID}{"score"} = $score;
			$data{$readID}{"data"} = $line;
		}
		elsif ( $coverage == $data{$readID}{"coverage"} ) {
			if ( $evalue < $data{$readID}{"evalue"} ) {
				$data{$readID}{"coverage"} = $coverage;
				$data{$readID}{"evalue"} = $evalue;
				$data{$readID}{"score"} = $score;
				$data{$readID}{"data"} = $line;
			}
			elsif ( $evalue == $data{$readID}{"evalue"} ) {
				if ( $score > $data{$readID}{"score"} ) {
					$data{$readID}{"coverage"} = $coverage;
					$data{$readID}{"evalue"} = $evalue;
					$data{$readID}{"score"} = $score;
					$data{$readID}{"data"} = $line;
				}
				elsif ( $score == $data{$readID}{"score"} ) {
					$data{$readID}{"data"} = $data{$readID}{"data"} . "XXXXX" . $line;
				}
			}
		}
	}
	else {
		if ( $id !~ /NaN/ ) {
			if ( $id !~ /$readID/ ) {
				# flush the existing data to the output file
				my @outArr = split(/XXXXX/,$data{$id}{"data"});

				for ( my $i = 0; $i <= $#outArr; $i++ ) {
					print OUT $outArr[$i] . "\n";
				}

				# store the new readID data
				$data{$readID}{"coverage"} = $coverage;
				$data{$readID}{"evalue"} = $evalue;
				$data{$readID}{"score"} = $score;
				$data{$readID}{"data"} = $line;

				# update $id to be $readID
				$id = $readID;
			}
		}
		else {
			$data{$readID}{"coverage"} = $coverage;
			$data{$readID}{"evalue"} = $evalue;
			$data{$readID}{"score"} = $score;
			$data{$readID}{"data"} = $line;
			$id = $readID;
		}
	}
}

close IN;

close OUT;

