# Count SNP number
for i in `ls *.snp.sorted.filtered.vcf`;do echo -n "$i ";grep -v "^#" $i|wc -l;done

# Count homozygous SNP number
for i in `ls *.snp.sorted.filtered.vcf`;do echo -n "$i ";bcftools filter -i 'GT="1/1" && DP>8' $i|grep -v "^#" |wc -l;done

# bggzip and index
for i in `ls *.snp.sorted.filtered.vcf`;do bgzip -i $i;bcftools index $i.gz;done

# rename and output homo.vcf
for i in `ls *.snp.sorted.filtered.vcf.gz`
do
sample=`echo $i|sed s/.snp.sorted.filtered.vcf.gz//`
cat > name.txt << eof
20 $sample
eof

bcftools filter -i 'GT="1/1" && DP>8' $i|bcftools reheader -s name.txt -o $sample.homo.vcf -
bgzip $sample.homo.vcf
bcftools index $sample.homo.vcf.gz

rm name.txt
done

# merge G0 samples to G0 group
bcftools isec G0_1.homo.vcf.gz G0_2.homo.vcf.gz G0_3.homo.vcf.gz -n +2 -o G0.homo.s2.txt
bcftools isec G0_1.homo.vcf.gz G0_2.homo.vcf.gz G0_3.homo.vcf.gz -n =3 -o G0.homo.s3.txt

bcftools merge --force-samples  -o G0_all.homo.vcf.gz -O z G0_1.homo.vcf.gz G0_2.homo.vcf.gz G0_3.homo.vcf.gz
bcftools view -T G0.homo.s2.txt G0_all.homo.vcf.gz -O z -o G0.homo.s2.vcf.gz
bcftools view -T G0.homo.s3.txt G0_all.homo.vcf.gz -O z -o G0.homo.s3.vcf.gz

bcftools index G0.homo.s2.vcf.gz
bcftools index G0.homo.s3.vcf.gz

wc -l G0.homo.s*.txt
# 261 G0.homo.s2.txt
#  97 G0.homo.s3.txt
zcat G0.homo.s2.vcf.gz|grep -v ^# |wc -l
#261
zcat G0.homo.s3.vcf|grep -v ^# |wc -l
#97

# use Venn to find SNPs in samples from experimental group but not in G0 samples
for i in `ls *.homo.vcf.gz|grep -v G0`;
do
  sample=`echo $i|cut -d'.' -f1`
  line=`bcftools isec -C $i G0.homo.s2.vcf.gz |wc -l`
  bcftools isec -C $i G0.homo.s2.vcf.gz -p ./
  mv 0000.vcf $sample.homo.s2.vcf
  echo $sample $line
done

for i in `ls *.homo.vcf.gz|grep -v G0`;
do
  sample=`echo $i|cut -d'.' -f1`
  line=`bcftools isec -C $i G0.homo.s3.vcf.gz |wc -l`
  bcftools isec -C $i G0.homo.s3.vcf.gz -p ./
  mv 0000.vcf $sample.homo.s3.vcf
  echo $sample $line
done

for i in `ls *.vcf|grep s2`;do grep -v ^# $i|wc -l;done
for i in `ls *.vcf|grep s3`;do grep -v ^# $i|wc -l;done

for i in `ls *.homo.s*.vcf`;do bgzip $i;bcftools index $i.gz;done

##### remove SNPs from GB SNP
bcftools concat GB.snp.vcf GB.indel.vcf > GB.all.vcf
bgzip GB.all.vcf
bcftools index GB.all.vcf.gz

for i in `ls *.homo.s*.vcf.gz|grep -v G0`
do
  sample=`echo $i|cut -d'.' -f1,2,3`
  bcftools isec -C $i GB.all.vcf.gz -p ./
  mv 0000.vcf $sample.GB.vcf
  bcftools isec -C $i GB.all.vcf.gz > $sample.GB.txt
done

wc -l *.homo.s2.GB.txt
wc -l *.homo.s3.GB.txt

##### annotation
#usage: perl annote_SNP.pl <SNP.txt> <SNP.annote.all.xls> <SNP.annote.priority.xls>
#for i in A21L1.homo.s2.GB.txt
for i in `ls *.GB.txt`
do
name=`echo $i|sed 's/.txt//'`
perl annote_SNP.20220627.pl $i $name.annote.all.xls $name.annote.priority.xls
done


### filter all G0 SNPs, not only homo but also hetero
bcftools filter -i 'GT="0/1" && DP>2' G0_1.snp.sorted.filtered.vcf.gz -O z -o G0_1.snp.sorted.filtered.hetero.vcf.gz
bcftools filter -i 'GT="0/1" && DP>2' G0_2.snp.sorted.filtered.vcf.gz -O z -o G0_2.snp.sorted.filtered.hetero.vcf.gz
bcftools filter -i 'GT="0/1" && DP>2' G0_3.snp.sorted.filtered.vcf.gz -O z -o G0_3.snp.sorted.filtered.hetero.vcf.gz
bcftools index G0_1.snp.sorted.filtered.hetero.vcf.gz
bcftools index G0_2.snp.sorted.filtered.hetero.vcf.gz
bcftools index G0_3.snp.sorted.filtered.hetero.vcf.gz

bcftools merge --force-samples  G0_1.snp.sorted.filtered.hetero.vcf.gz G0_2.snp.sorted.filtered.hetero.vcf.gz G0_3.snp.sorted.filtered.hetero.vcf.gz -O z -o G0_all.hetero.vcf.gz
bcftools index G0_all.hetero.vcf.gz

for i in `ls *.homo.*.GB.vcf`
do
  sample=`echo $i|cut -d'.' -f1-4`
  bgzip $i
  bcftools index $i.gz
  bcftools isec -C $i.gz G0_all.hetero.vcf.gz -p ./
  mv 0000.vcf $sample.G0_het.vcf
  bcftools isec -C $i.gz G0_all.hetero.vcf.gz > $sample.G0_het.txt
done

wc -l *.homo.s2.GB.G0_het.txt
wc -l *.homo.s3.GB.G0_het.txt

### filter het SNPs in GB, using bg.pileup.allele_SNV.het.xls
for i in `ls *.GB.G0_het.txt`
do
name=`echo $i|sed 's/.txt//'`
perl annote_SNP.20220629.pl $i $name.rmGB_SNP.annote.all.xls $name.rmGB_SNP.annote.priority.xls
done

for i in `ls *.s2.*.rmGB_SNP.annote.priority.xls`;do echo -n $i"  ";cut -f2 $i|sed '1d'|uniq|wc -l ;done
for i in `ls *.s3.*.rmGB_SNP.annote.priority.xls`;do echo -n $i"  ";cut -f2 $i|sed '1d'|uniq|wc -l ;done

### grep flower
for i in `ls *.GB.G0_het.txt`
do
name=`echo $i|sed 's/.txt//'`
grep -w -f flower.txt $name.rmGB_SNP.annote.all.xls > $name.rmGB_SNP.annote.all.flower.xls
done

wc -l *.rmGB_SNP.annote.all.flower.xls