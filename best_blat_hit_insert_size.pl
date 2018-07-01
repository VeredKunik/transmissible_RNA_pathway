#!/usr/bin/perl

my $blatR1 = $ARGV[0];
my $blatR2 = $ARGV[1];
my $outputFile = $ARGV[2];

# variables
my $line = "";
my %data;
my %outHash;
my $id = "";


# open input file
open IN, "<$blatR1" || die "could not open $blatR1!";
while ( $line = <IN> ) {

	chomp($line);

	# @NS500125:317:HG22YAFXX:1:11101:11498:1611      gi|384444089|ref|NC_017248.1|   83.49   109     18      0       1       109     1801583 1801475 2.1e-37 153.0
	# @NS500125:317:HG22YAFXX:1:11101:11498:1611      gi|384444089|ref|NC_017248.1|   83.49   109     18      0       1       109     1602125 1602017 2.1e-37 153.0
	my @tmp = split(/\s+/,$line);
	my $readID = $tmp[0];
	my $genomeID = $tmp[1];
	my $genomicStart = (($tmp[8]<$tmp[9])?$tmp[8]:$tmp[9]);
	my $genomicEnd = (($tmp[8]<$tmp[9])?$tmp[9]:$tmp[8]);
	my $match_length = $tmp[3];
	my $homology = $tmp[2];
	my $evalue = $tmp[10];
	my $score = $tmp[11];
	my $strand = (($tmp[8]>$tmp[9])?-1:1);

	push(@{$data{$readID}{$genomeID}{"R1"}{"start"}},$genomicStart);
	push(@{$data{$readID}{$genomeID}{"R1"}{"end"}},$genomicEnd);
	push(@{$data{$readID}{$genomeID}{"R1"}{"match"}},$match_length);
	push(@{$data{$readID}{$genomeID}{"R1"}{"homology"}},$homology);
	push(@{$data{$readID}{$genomeID}{"R1"}{"evalue"}},$evalue);
	push(@{$data{$readID}{$genomeID}{"R1"}{"score"}},$score);
	push(@{$data{$readID}{$genomeID}{"R1"}{"strand"}},$strand);

}

close IN;


# open input file
open IN, "<$blatR2" || die "could not open $blatR2!";
while ( $line = <IN> ) {

	chomp($line);
	# R2 data
	# @NS500125:317:HG22YAFXX:1:11101:11498:1611      gi|384444089|ref|NC_017248.1|   73.39   124     33      0       1       124     1801350 1801473 2.2e-28 123.0
	# R1 data
	# @NS500125:317:HG22YAFXX:1:11101:11498:1611      gi|384444089|ref|NC_017248.1|   83.49   109     18      0       1       109     1801583 1801475 2.1e-37 153.0
	# @NS500125:317:HG22YAFXX:1:11101:11498:1611      gi|384444089|ref|NC_017248.1|   83.49   109     18      0       1       109     1602125 1602017 2.1e-37 153.0

	my @tmp = split(/\s+/,$line);
	my $readID = $tmp[0];
	my $genomeID = $tmp[1];
	my $genomicStart = (($tmp[8]<$tmp[9])?$tmp[8]:$tmp[9]);
	my $genomicEnd = (($tmp[8]<$tmp[9])?$tmp[9]:$tmp[8]);
	my $match_length = $tmp[3];
	my $homology = $tmp[2];
	my $evalue = $tmp[10];
	my $score = $tmp[11];
	my $strand = (($tmp[8]<$tmp[9])?-1:1);

	push(@{$data{$readID}{$genomeID}{"R2"}{"start"}},$genomicStart);
	push(@{$data{$readID}{$genomeID}{"R2"}{"end"}},$genomicEnd);
	push(@{$data{$readID}{$genomeID}{"R2"}{"match"}},$match_length);
	push(@{$data{$readID}{$genomeID}{"R2"}{"homology"}},$homology);
	push(@{$data{$readID}{$genomeID}{"R2"}{"evalue"}},$evalue);
	push(@{$data{$readID}{$genomeID}{"R2"}{"score"}},$score);

}

close IN;

# traverse the data and calculate the length of the insert
foreach my $readID ( keys %data ) {
	foreach my $genomeID ( keys %{$data{$readID}} ) {
		my $gap = 0;
		my $insert_size = 0;
		my $evalue = 1;

		if ( (exists $data{$readID}{$genomeID}{"R1"}{"start"}) && (exists $data{$readID}{$genomeID}{"R2"}{"start"}) ) {
			if ( $data{$readID}{$genomeID}{"R2"}{"strand"} == $data{$readID}{$genomeID}{"R1"}{"strand"} ) {
				for ( my $i = 0; $i <= $#{$data{$readID}{$genomeID}{"R1"}{"start"}}; $i++ ) {
					# reset variables for the current iteration
					$insert_size = 0;
					$evalue = 0;
					$gap = 0;

					my $R1_end = $data{$readID}{$genomeID}{"R1"}{"end"}->[$i];
					my $R1_start = $data{$readID}{$genomeID}{"R1"}{"start"}->[$i];
					my $R2_end = $data{$readID}{$genomeID}{"R2"}{"end"}->[$i];
					my $R2_start = $data{$readID}{$genomeID}{"R2"}{"start"}->[$i];

					if ( $R1_end < $genomicStart ) { # case 5: e1 < s2
						$gap = $genomicStart - $R1_end + 1;
						$match_length = $data{$readID}{$genomeID}{"R2"}{"match"}->[$i] + $data{$readID}{$genomeID}{"R1"}{"match"}->[$i];
					}
					elsif ( $R2_end < $R1_start ) { # case 6  e2 < s1
						$gap = $R1_start - $R2_end + 1;
						$match_length = $data{$readID}{$genomeID}{"R2"}{"match"}->[$i] + $data{$readID}{$genomeID}{"R1"}{"match"}->[$i];
					}
					elsif ( $R1_start < $R2_start ) { # cases 1 or 4 --> s1 < s2
						if ( $R2_end < $R1_end ) { # case 1 --> e2 < e1
							$gap = 0;
							$match_length = $data{$readID}{$genomeID}{"R1"}{"match"}->[$i];
						}
						elsif ( ( $R2_start < $R1_end) && ($R1_end < $R2_end)  ) { # case 4 --> s2 < e1 && e1 < e2
							$gap = 0;
							$match_length = $R2_end - $R1_start + 1;
						}
					}
					elsif ( $R2_start < $R1_start) { # cases 2 or 3 --> s2 < s1
						if ( $R1_end < $R2_end ) { # case 2 --> e1 < e2
							$gap = 0;
							$match_length = $data{$readID}{$genomeID}{"R2"}{"match"}->[$i];
						}
						elsif ( ($R1_start < $R2_end) && ($R2_end < $R1_end) ) { # case 3 --> s1 < e2 && e2 < e1
							$gap = 0;
							$match_length = $R1_end - $R2_start + 1;
						}
					}
	
					$insert_size = (($gap + $match_length)<300?($gap + $match_length):300);
					$evalue = ($data{$readID}{$genomeID}{"R1"}{"evalue"}+$data{$readID}{$genomeID}{"R2"}{"evalue"})/2;

					if ( exists $outHash{$readID}{"insert_size"} ) {
						if ( $insert_size > $outHash{$readID}{"insert_size"} ) {
							$outHash{$readID}{"insert_size"} = $insert_size;
							$outHash{$readID}{"evalue"} = $evalue;
						}
						elsif ( $insert_size == $outHash{$readID}{"insert_size"} ) {
							if ( $evalue < $outHash{$readID}{"evalue"} ) {
								$outHash{$readID}{"insert_size"} = $insert_size;
								$outHash{$readID}{"evalue"} = $evalue;
							}
							else { # reset variables for the next iteration
								$insert_size = 0;
								$evalue = 1;
								$gap = 0;
							}
						}
						else { # reset variables for the next iteration
							$insert_size = 0;
							$evalue = 1;
							$gap = 0;
						}
					}
					else {
						$outHash{$readID}{"insert_size"} = $insert_size;
						$outHash{$readID}{"evalue"} = $evalue;
					}
				}
			}		
		}
		elsif ( (exists $data{$readID}{$genomeID}{"R1"}{"start"}) && (!exists $data{$readID}{$genomeID}{"R2"}{"start"}) ) { # only R1 data exists
			for ( my $i = 0; $i <= $#{$data{$readID}{$genomeID}{"R1"}{"start"}}; $i++ ) {
				# reset variables for the current iteration
				$insert_size = $data{$readID}{$genomeID}{"R1"}{"match"}->[$i];
				$evalue = $data{$readID}{$genomeID}{"R1"}{"evalue"}->[$i];

				if ( $insert_size > $outHash{$readID}{"insert_size"} ) {
					$outHash{$readID}{"insert_size"} = $insert_size;
					$outHash{$readID}{"evalue"} = $evalue;
				}
				elsif ( $insert_size == $outHash{$readID}{"insert_size"} ) {
					if ( $evalue < $outHash{$readID}{"evalue"} ) {
						$outHash{$readID}{"insert_size"} = $insert_size;
						$outHash{$readID}{"evalue"} = $evalue;
					}
					else { # reset variables for the next iteration
						$insert_size = 0;
						$evalue = 1;
					}
				}
				else { # reset variables for the next iteration
					$insert_size = 0;
					$evalue = 1;
				}
			}
		}
		elsif ( (exists $data{$readID}{$genomeID}{"R2"}{"start"}) && (!exists $data{$readID}{$genomeID}{"R1"}{"start"}) ) { # only R2 data exists
			for ( my $i = 0; $i <= $#{$data{$readID}{$genomeID}{"R2"}{"start"}}; $i++ ) {
				# reset variables for the current iteration
				$insert_size = $data{$readID}{$genomeID}{"R2"}{"match"}->[$i];
				$evalue = $data{$readID}{$genomeID}{"R2"}{"evalue"}->[$i];

				if ( $insert_size > $outHash{$readID}{"insert_size"} ) {
					$outHash{$readID}{"insert_size"} = $insert_size;
					$outHash{$readID}{"evalue"} = $evalue;
				}
				elsif ( $insert_size == $outHash{$readID}{"insert_size"} ) {
					if ( $evalue < $outHash{$readID}{"evalue"} ) {
						$outHash{$readID}{"insert_size"} = $insert_size;
						$outHash{$readID}{"evalue"} = $evalue;
					}
					else { # reset variables for the next iteration
						$insert_size = 0;
						$evalue = 1;
					}
				}
				else { # reset variables for the next iteration
					$insert_size = 0;
					$evalue = 1;
				}
			}
		}
	}
}


open OUT, ">$outputFile" || die "could not open $outputFile!";

foreach my $readID ( keys %outHash ) {
	
	print OUT $readID . "\t";
	
	if ( 300 ==$outHash{$readID}{"insert_size"} ) {
		print OUT ">" . $outHash{$readID}{"insert_size"};
	}
	else {
		print OUT $outHash{$readID}{"insert_size"};
	}

	print OUT "\t" . $outHash{$readID}{"evalue"} . "\n";
}

close OUT;
