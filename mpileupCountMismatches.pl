#!/usr/bin/env perl
#  MORPhEUS. RNA polymerase errors cause splicing defects and can be regulated by differential expression of RNA polymerase subunits
#
# Lucas B. Carey.   eLife 2015
# lucas.carey@upf.edu
#
#
# Given the output of samtools mpileup
#   print the # ref & # errors at each position in the genome
#   optionally allows for filtering to remove positions with too many errors
#       ie: clonal/subclonal mutations
#    
use warnings;
use strict;
use Getopt::Long ;
use POSIX qw(floor);
use Switch ;

my $region_file = '';
my $debug = 0;
my %H ;
my %RB;
my %RC;
my %RefBaseCounts ;
my $max_mismatches_per_position = 100 ; 
my $max_percent_mismatches_per_position = 10 ; 
my $print_base_comp_flag = 1 ;
my $bedgraph_output_file = '';
my $mpileup_gzipped_file = '';
GetOptions(
	"region_file=s" => \$region_file,
	"mpileup_gzipped_file=s" => \$mpileup_gzipped_file,
	"bedgraph_output_file=s" => \$bedgraph_output_file,
	"print_base_comp" => \$print_base_comp_flag,
	"max_mismatches_per_position=i" => \$max_mismatches_per_position,
	"max_percent_mismatches_per_position=i" => \$max_percent_mismatches_per_position,
	"debug" => \$debug 
);

open(IN,"gunzip -c $mpileup_gzipped_file |") || die "$! \n $mpileup_gzipped_file \n";

# chrI	10	A	32	.............C....,............^S.	B[I4GGGAFG=A@1JGGF3DFDF4B;FC+@;@	39,20,19,41,18,14,14,15,10,22,27,41,19,20,15,20,29,26,8,8,8,7,7,7,7,10,5,5,4,4,3,2,2,2,1,1
# chrII	663002	C	43	.NNNNNNNNN.......NNNN......................	dooooooooodddddddooooddddddddddddddRRRRRRRR	47,47,47,39,33,32,32,32,32,32,31,30,30,29,29,29,28,24,24,24,21,21,18,18,15,15,14,13,12,12,12,12,12,12,12,12,12,10,9,9,9,9,9,9,9,9,5,5,4

# region_file is chrX:YYYY or chrX\tYYYY  (one position, one chr per line, locations list file for GATK)
if( $region_file){
	open(O,$region_file) || die "Cannot open $region_file $!\n";
	while(<O>){
		chomp;
		my @l = split(/[:\t-]/);
		if ($_ =~ m/\d[-\t]\d/){
			for(my $I=$l[1];$I<=$l[2];$I++){
				$H{$l[0] . $I}++;
			} 
		}else{
		$H{$l[0] . $l[1]}++;
	}
	}
	close(O);
}

$RefBaseCounts{'A'}=0; $RefBaseCounts{'C'}=0;$RefBaseCounts{'G'}=0;$RefBaseCounts{'T'}=0;

my $nref = 0;
my $nerr = 0;
#if($bedgraph_output_file){
#	print STDERR "Creating $bedgraph_output_file\n";
#	open(BGOFE,">$bedgraph_output_file.errors.bedgraph");
#	open(BGOFR,">$bedgraph_output_file.reads.bedgraph");
#}
while(<IN>){
	chomp;
	my @l = split(/\t/);
	if( ($region_file && $H{$l[0] . $l[1]}) || !$region_file){
		my $ref_base = uc($l[2]); 
		my $nreads = $l[3];
		$RB{$ref_base}++;
		$RC{$ref_base}+= $nreads;
		next if ($nreads == 0);
		next if ($l[4] =~ m/s[\+-]/); # no indels
		my $nref_fwd = () = $l[4] =~ /[\.]/gi;
	my $nref_rev = () = $l[4] =~ /[,]/gi;
	next unless( ($nref_fwd + $nref_rev) > 0 ) ;
	my $nNs = () = $l[4] =~ /[nN]/gi;
	my $A   = () = $l[4] =~ /[aA]/gi;
	my $C   = () = $l[4] =~ /[cC]/gi;
	my $T   = () = $l[4] =~ /[tT]/gi;
	my $G   = () = $l[4] =~ /[gG]/gi;
#	my $ne = $nreads - $nNs - $nr ;
	my $ne = $A + $C + $T + $G ;
	if ($ne>0 && $debug){
		print STDERR  "$region_file\t" if ($region_file) ;
		print STDERR  "$A,$C,$T,$G,$ne\t$_\n" ;
	}
	my $mpm = floor( ($nref_fwd + $nref_rev) * ( $max_percent_mismatches_per_position  / 100) );
	if( ($A>$max_mismatches_per_position || $C>$max_mismatches_per_position || $T>$max_mismatches_per_position || $G>$max_mismatches_per_position) || ($A>$mpm || $C>$mpm || $T>$mpm || $G>$mpm) ){
		print STDERR "$nref_fwd\t$nref_rev\t$nNs\t$ne\t$_\n" if ($debug);
	} else{
		$nref += ($nref_fwd + $nref_rev);
		$nerr += $ne ;
		$RefBaseCounts{$ref_base} += ($nref_fwd + $nref_rev);
		print "$l[0]\t" . ($l[1] -1) . "\t" . $l[1] . "\t$ref_base\t" . ($nref_fwd + $nref_rev) . "\t" . $ne . "\n"  ;
#		switch($ref_base){
#			case 'A' { $RefBaseCounts{'T'}+=$nref_rev ;}
#			case 'T' { $RefBaseCounts{'A'}+=$nref_rev ;}
#			case 'C' { $RefBaseCounts{'G'}+=$nref_rev ;}
#			case 'G' { $RefBaseCounts{'C'}+=$nref_rev ;}
#			}
	}
}
}
if($nref>0){
	print "$region_file\t" if ($region_file) ;
	print STDERR "$mpileup_gzipped_file\t$RefBaseCounts{'A'}\t$RefBaseCounts{'C'}\t$RefBaseCounts{'T'}\t$RefBaseCounts{'G'}\t" if($print_base_comp_flag);
	print STDERR "$mpileup_gzipped_file\t$nref\t$nerr\t" . $nerr / $nref . "\n"; 
}

foreach my $b (sort(keys %RB)){
	print STDERR "$mpileup_gzipped_file\t$region_file\t" if ($region_file) ;
	print STDERR "$mpileup_gzipped_file\t$b\t$RB{$b}\t$RC{$b}\n";
}
