#!/usr/bin/perl
use List::Util qw(min);
=pod
1       21      22      +       0.3     0       0.3     CHH     Hyper   1       0       107 DR000001|DR000001|repeat .       -
1       22      23      +       0.272727272727273       0       0.272727272727273       CHH Hyper    1       0       107     DR000001|DR000001|repeat        .       -
=cut
#usage: perl annote_DML.pl <DML.bed> <DML.annote.all.xls> <DML.annote.priority.xls>

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

my $bed="/data4/LZG_WGBS/20210111_gexinghua_reannotation/00.script/annotation/final_bed/whole_genome.bed";
my $genenamefile="/data4/LZG_WGBS/20210111_gexinghua_reannotation/00.script/annotation/genenamefile.all.txt";

my %des;
open (GENE,"$genenamefile");
while(<GENE>){
	chomp;
	s/~//g;
	my @line=split /\t/;
	$des{$line[0]}=$line[1]."\t".$line[2];
}
close GENE;

my %result;
open (ALL,"> $ARGV[1]") || die $!;
open (PRI,"> $ARGV[2]") || die $!;
`cut -f1-8  $ARGV[0]|sed '1d'|sort -u|perl -alne '\$F[1]=(\$F[1]-1)."\\t".\$F[1];print join("\\t",\@F)' > $ARGV[0].bed`;
open (IN,"bedtools intersect -a $ARGV[0].bed -b $bed -wa -wb -nobuf|") || die $!;
while(<IN>){
	chomp;
	next if /Error/;
	$_=~s/\|/\t/g;
	my @line=split /\t/;
#	my $region_type=$line[14];
#	my $region_id=$line[13];
	my $id=$line[0]."_".$line[2];
	my $descrip;
	if (exists $des{$line[13]}){
		$descrip=$des{$line[13]};
	}else{
		$descrip="-\t-";
	}
	my @content=($line[0],$line[2],$line[3],$line[4],$line[5],$line[6],$line[7],$line[8],$line[9],$line[10],$line[11],$line[12],$line[13],$line[14],$line[16],$descrip);
	push @{$result{$id}{$priority{$line[14]}}},join("~",@content);
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

`sed -i '1ichr\tpos\tstrand\tML_case\tML_control\tML_diff\tContext\tDM_type\tregion_chr\tregion_start\tregion_end\tregion_id_alt\tregion_id\tregion_type\tregion_strand\tsymbol\tdescription' $ARGV[1]`;
`sed -i '1ichr\tpos\tstrand\tML_case\tML_control\tML_diff\tContext\tDM_type\tregion_chr\tregion_start\tregion_end\tregion_id_alt\tregion_id\tregion_type\tregion_strand\tsymbol\tdescription' $ARGV[2]`;
