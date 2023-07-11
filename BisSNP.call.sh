for i in {1..35}
do
  sample=`head -n $i all.bam.txt|tail -n1|cut -f1`
  bam=`head -n $i all.bam.txt|tail -n1|cut -f2`
  mkdir -p $sample
  cat > $sample/$sample.sh <<eof
fa=/data4/LZG_WGBS/20220518_gexinghua/genome.fa
fai=/data4/LZG_WGBS/20220518_gexinghua/genome.fa.fai
sample=$sample
#input=
sortbam=$bam
groupbam=\$sample.sorted.group.bam
methvct=\$sample.meth.vcf
snpvcf=\$sample.snp.vcf
snpsortvcf=\$sample.snp.sorted.vcf
snpfiltervcf=\$sample.snp.sorted.filtered.vcf
snpfilterlog=\$sample.snp.sorted.filter.summary.txt

#sort bam file and index
#samtools sort \$input \$sortbam
#samtools index \$sortbam

echo add_group
#add a new read group to the library
java -Xmx64g -jar /data1/wjx_software/picard/picard.jar AddOrReplaceReadGroups \
 I=\$sortbam \
 O=\$groupbam \
 RGID=4 \
 RGLB=lib1 \
 RGPL=illumina \
 RGPU=unit1 \
 RGSM=20

samtools index \$groupbam

echo call_snp
#call snp (attention parameters: nonDirectional for non directional library protocol)
java -Xmx64g -jar /data4/LZG_WGBS/20220518_gexinghua/BisSNP-1.0.0.jar \
 -R \$fa \
 -I \$groupbam \
 -T BisulfiteGenotyper \
 --trim_5_end_bp 0 \
 --trim_3_end_bp 0 \
 -vfn1 \$methvct -vfn2 \$snpvcf \
 -mbq 12 -minConv 1 -toCoverage 1000 -mmq 20 -nt 16
# --dbsnp /dbsnp/dbsnp147.GRCH37/All_20160601.vcf.gz  # -L chr11:7000000-7100000

#sort vcf and filter fake SNPs (By default, it filter out SNPs with quality score less than 20, reads coverage more than 250, strand bias more than -0.02, quality score by depth less than 1.0, mapping quality zero reads fraction more than 0.1 and 2 SNPs within the same 10 bp window.)
perl /data1/wjx_software/BisSNP/Utils/sortByRefAndCor.pl \
 --k 1 --c 2  \
 \$snpvcf \
 \$fai |perl -pe 's/^\\n//g' > \$snpsortvcf
java -Xmx64g -jar /data4/LZG_WGBS/20220518_gexinghua/BisSNP-1.0.0.jar \
 -R \$fa \
 -T VCFpostprocess \
 -oldVcf \$snpsortvcf \
 -newVcf \$snpfiltervcf \
 -snpVcf \$snpsortvcf \
 -o \$snpfilterlog
eof

done