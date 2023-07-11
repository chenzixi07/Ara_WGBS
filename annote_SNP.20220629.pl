#!/usr/bin/perl
use List::Util qw(min);
=pod
1       17111141        N       T,C     10
1       22451084        T       C       10
=cut
#usage: perl annote_SNP.pl <SNP.txt> <SNP.annote.all.xls> <SNP.annote.priority.xls>

my %priority=(
"CDS"=>1,
"intron"=>2,
"utr5"=>3,
"utr3"=>4,
"promoter"=>5,
"TEgene"=>6,
"TE"=>7,
"repeat"=>8,
"ncRNA"=>9,
"pseudogene"=>10,
"downstream"=>11,
"intergenic"=>12
);

my $bed="/mnt/d/SZU/Projects/031.LZG_WGBS/SNP/02.Filter_Merge/whole_genome.bed";
my $genenamefile="/mnt/d/SZU/Projects/031.LZG_WGBS/SNP/02.Filter_Merge/genenamefile.all.txt";
#my $bg="/mnt/d/SZU/Projects/031.LZG_WGBS/SNP/02.Filter_Merge/bg.pileup.allele_SNV.het.xls";
my $bg="/mnt/d/SZU/Projects/031.LZG_WGBS/SNP/02.Filter_Merge/bg.SNP.xls";

my %des;
open (GENE,"$genenamefile");
while(<GENE>){
	chomp;
	s/~//g;
	my @line=split /\t/;
	$des{$line[0]}=$line[1]."\t".$line[2];
}
close GENE;

my %bg;
open (BG,"$bg") || die $!;
while(<BG>){
    chomp;
    my @line = split /\t/;
    my $id = $line[0]."_".$line[1];
    $bg{$id} = $id;
}

open (TXT,"$ARGV[0]") || die $!;
open (TEMP,"> $ARGV[0].temp") || die $!;
while(<TXT>){
    chomp;
    my @line = split /\t/;
    my $start = $line[1]-1;
    my $id = $line[0]."_".$line[1];
    $line[1] = $start."\t".$line[1];
    $line[-1] = "+";
    print TEMP join("\t",@line)."\n" unless exists $bg{$id};
}
close TXT;
close TEMP;

my %result;
open (ALL,"> $ARGV[1]") || die $!;
open (PRI,"> $ARGV[2]") || die $!;
open (IN,"bedtools intersect -a $ARGV[0].temp -b $bed -wa -wb -nobuf|") || die $!;
while(<IN>){
	chomp;
	next if /Error/;
	$_=~s/\|/\t/g;
	my @line=split /\t/;
#	my $region_type=$line[14];
#	my $region_id=$line[13];
	my $id=$line[0]."_".$line[2];
	my $descrip;
	if (exists $des{$line[10]}){
		$descrip=$des{$line[10]};
	}else{
		$descrip="--\t--";
	}
	my @content=($line[0],$line[2],$line[3],$line[4],$line[6],$line[7],$line[8],$line[9],$line[10],$line[11],$line[13],$descrip);
	push @{$result{$id}{$priority{$line[11]}}},join("~",@content);
	print ALL join("\t",@content);
	print ALL "\n";
}
close IN;

foreach my $id (keys %result){
	my $type=min(keys %{$result{$id}});
	foreach my $content (@{$result{$id}{$type}}){
		$content =~ s/~/\t/g;
		print PRI "$content\n";
	}
}
close ALL;
close PRI;
`sort -nk1,2 $ARGV[1]|uniq > $ARGV[1].temp`;
`sort -nk1,2 $ARGV[2]|uniq > $ARGV[2].temp`;

`mv $ARGV[1].temp $ARGV[1]`;
`mv $ARGV[2].temp $ARGV[2]`;

`sed -i '1ichr\tpos\tref\talt\tregion_chr\tregion_start\tregion_end\tregion_id_alt\tregion_id\tregion_type\tregion_strand\tsymbol\tdescription' $ARGV[1]`;
`sed -i '1ichr\tpos\tref\talt\tregion_chr\tregion_start\tregion_end\tregion_id_alt\tregion_id\tregion_type\tregion_strand\tsymbol\tdescription' $ARGV[2]`;

`rm $ARGV[0].temp`