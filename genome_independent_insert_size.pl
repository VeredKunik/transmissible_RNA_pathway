#!/usr/bin/perl

my $fastaR1File = $ARGV[0];
my $fastaR2File = $ARGV[1];
my $sampleID = $ARGV[2];
my $outputFile = $ARGV[3];

# variables
my $line = "";
my $line1 = "";
my $line2 = "";

my %data;
my $id = "";
my $id1 = "";
my $id2 = "";
my %outHash;

open IN1, "<$fastaR1File" || die "could not open $fastaR1File!";
open IN2, "<$fastaR2File" || die "could not open $fastaR2File!";
open OUTPUT, ">$outputFile" || die "could not open $outputFile!";

while ( ($line1 = <IN1>) && ($line2 = <IN2>) ) {

	chomp($line1);
	chomp($line2);

	# >@NS500125:317:HG22YAFXX:1:11101:20922:1042 1:N:0:CGATGT
	# GAGTTTAAGCATATCAATAAGCGGAGGAAAAGAAA

	if ( ($line1 =~ /^>/) && ($line2 =~ /^>/) ) {
		my @tmp = split(/\s+/,$line1);
		$id1 = substr($tmp[0],1);

		my @tmp = split(/\s+/,$line2);
		$id2 = substr($tmp[0],1);
	}
	elsif ( ($line1 !~ /^>/) && ($line2 !~ /^>/) ) {
		if ( $id1 =~ /$id2/ && ($line1 !~ /N/) && ($line2 !~ /N/) ) {
			$id = $id1;
			$id =~ s/@//;

			$R1 = $line1;
			$R2 = &reverse_complement($line2);

			&create_fasta_file("$sampleID.R1.fasta",$id,$R1);
			&create_fasta_file("$sampleID.R2.fasta",$id,$R2);

			my $cmd = "blast-2.2.26/bin/bl2seq -i $sampleID.R2.fasta -j $sampleID.R1.fasta -p blastn -o $sampleID.$id.out";
			system($cmd);
			my $size = &parse_blast_output("$sampleID.$id.out",length($R1),length($R2));

			print OUTPUT $size . "\t" . $id . "\n";

			system("rm -f $sampleID.R1.fasta");
			system("rm -f $sampleID.R2.fasta");
			system("rm -f $sampleID.$id.out");
		}

		$id = "";
		$id1 = "";
		$id2 = "";
	}
}

close IN1;

close IN2;

close OUTPUT;

##################### subroutines #######################

sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}

sub parse_blast_output {
	my $blastFile = shift;
	my $R1_length = shift;
	my $R2_length = shift;

	my $match = 0;
	my $insert_size = 0;

	open IN, "<$blastFile" || die "could not open $blastFile!";
	while ( $line = <IN> ) {
		chomp($line);

		if ( $line =~ /Identities = / ) {
			#  Identities = 20/20 (100%)
			$line =~ s/ //g; # Identities=20/20(100%);

			my @tmp = split(/\//,$line);
			$match = $tmp[0]; # Identities=20
			$match =~ s/Identities=//; # 20 


		}
	}
	
	close IN;
	
	$insert_size = ($R1_length - $match) + $match + ($R2_length - $match);

	if ( 0 == $match ) {
		$insert_size = ">" . $insert_size;
	}

	return $insert_size;
}

sub create_fasta_file {

	my $file_name = shift;
	my $readID = shift;
	my $seq = shift;

	open OUT, ">$file_name" || die "could not open $file_name!";
	print OUT ">" . $readID . "\n";
	print OUT $seq . "\n";
	close OUT;
}