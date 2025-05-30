#!/bin/bash

# lcwes version 4.0.0

# Set the directory where the FASTQ files are located
FASTQ_DIR="."

# Set paths to required tools and reference files
# Alignment-Variant Calling
REF_GENOME="/home/administrator/lifecode/genomes/databases/bwa_hg19/hg19.fa"
REF_GENOME_dict="/home/administrator/lifecode/genomes/databases/bwa_hg19/hg19.fa"
TARGETS="/home/administrator/lifecode/genomes/bed_files/WES_HG19/S33266340_Covered.adj.bed"
# Annotation
SNPEFF_JAR="/home/administrator/snpeff/snpEff/snpEff.jar"
CLINVAR_VCF="/home/administrator/lifecode/genomes/databases/clnvar_hg19/clinvar.chr.vcf.gz"
Mills_100G="/home/administrator/lifecode/genomes/databases/GATK_resources/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
# Downstream analysis
INTERVARDB="/home/administrator/lifecode/genomes/databases/intervar_humandb/hg19/intervar"
HUMANDB="/home/administrator/lifecode/genomes/databases/intervar_humandb/hg19/humandb"
FREEBAYES_REGIONS="/home/administrator/lifecode/genomes/databases/freebayes_regions_hg19/hg19_regions.txt"
# Computation
THREADS=28

# Function to process a single sample
process_sample() {
	local sample=$1
	local fastq1=$2
	local fastq2=$3

#-------------------------- Filtering & Alignment ---------------------------#

# Trimming (without UMI)
#java -jar /home/administrator/AGeNT/agent/lib/trimmer-3.1.2.jar \
#	-fq1 ${sample}_1.fq.gz \
#	-fq2 ${sample}_2.fq.gz \
#	-adaptor MGI -mbc null -IDEE_FIXE -out trimmed/${sample}

# Alignment
#bwa mem -R "@RG\tID:${sample}\tLB:exome_lib\tPL:MGISEQ\tPU:unit1\tSM:${sample}" -t $THREADS \
#	$REF_GENOME \
#	trimmed/${sample}_R1.fastq.gz \
#	trimmed/${sample}_R2.fastq.gz | samtools view -@ $THREADS -bS | samtools sort -@ $THREADS -o ${sample}_aligned_rg.bam

# Mark Duplicates
#gatk MarkDuplicates \
#	-I ${sample}_aligned_rg.bam \
#	-O ${sample}_aligned_marked.bam \
#	-M ${sample}_output.metrics.txt \
#	--ASSUME_SORT_ORDER coordinate \
#	--CREATE_INDEX true \
#	--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500

#samtools index ${sample}_aligned_marked.bam

#rm ${sample}_aligned_rg.bam*

# Base Quality Score Recalibration
#gatk BaseRecalibrator --java-options "-Xmx48G -XX:+UseParallelGC -XX:ParallelGCThreads=32" \
#	-R $REF_GENOME_dict \
#	-I ${sample}_aligned_marked.bam \
#	--known-sites $DBSNP_VCF \
#	--known-sites $Mills_100G \
#	-L $TARGETS \
#	-O ${sample}_recal_data.table

# Apply BQSR
#gatk --java-options "-Xmx48G -XX:+UseParallelGC -XX:ParallelGCThreads=$THREADS" ApplyBQSR \
#	-R $REF_GENOME_dict \
#	-I ${sample}_aligned_marked.bam \
#	--bqsr-recal-file ${sample}_recal_data.table \
#	-O ${sample}_aligned_marked_bqsr.bam

#rm ${sample}_aligned_marked.bam*

# Prepare annovar scripts
#ln -s /home/administrator/Annovar/annovar/*.pl .

#-------------------------- GATK Variant Calling ---------------------------#

# Variant calling GATK
#gatk HaplotypeCaller \
#	-R $REF_GENOME \
#	-I ${sample}_aligned_marked_bqsr.bam \
#	-O ${sample}_variants.vcf.gz \
#	--native-pair-hmm-threads $THREADS

# Extract SNPs
#gatk SelectVariants \
#	-V ${sample}_variants.vcf.gz \
#	-select-type SNP \
#	-O ${sample}_snps.vcf.gz

# Extract INDELs
#gatk SelectVariants \
#	-V ${sample}_variants.vcf.gz \
#	-select-type INDEL \
#	-O ${sample}_indels.vcf.gz

# variant filtering SNPS
#gatk VariantFiltration \
#	-V ${sample}_snps.vcf.gz \
#	-O ${sample}_snps.filtered.vcf.gz \
#	--filter-name "QD_filter" --filter-expression "QD < 2.0" \
#	--filter-name "FS_filter" --filter-expression "FS > 60.0" \
#	--filter-name "MQ_filter" --filter-expression "MQ < 40.0" \
#	--filter-name "MQRankSum_filter" --filter-expression "MQRankSum < -12.5" \
#	--filter-name "ReadPosRankSum_filter" --filter-expression "ReadPosRankSum < -8.0" \
#	--filter-name "QUAL_filter" --filter-expression "QUAL < 100.0" \
#	--filter-name "DP_filter" --filter-expression "DP < 10.0" \
#	--filter-name "SOR_filter" --filter-expression "SOR > 3.0"

#rm ${sample}_snps.vcf.gz*

# variant filtering INDELS
#gatk VariantFiltration \
#	-V ${sample}_indels.vcf.gz \
#	-O ${sample}_indels.filtered.vcf.gz \
#	--filter-name "QD_filter" --filter-expression "QD < 2.0" \
#	--filter-name "FS_filter" --filter-expression "FS > 200.0" \
#	--filter-name "ReadPosRankSum_filter" --filter-expression "ReadPosRankSum < -20.0" \
#	--filter-name "SOR_filter" --filter-expression "SOR > 10.0"

#rm ${sample}_indels.vcf.gz*

# Merge filtered SNPS and INDELS
#gatk MergeVcfs \
#	-I ${sample}_snps.filtered.vcf.gz \
#	-I ${sample}_indels.filtered.vcf.gz \
#	-O ${sample}_GATK.filtered.vcf.gz

#rm ${sample}_snps.filtered.vcf.gz* ${sample}_indels.filtered.vcf.gz*

#mkdir tmp
#mv ${sample}_variants.vcf.gz* tmp/

# Spilit mulitallelic sites / Normalize
#bcftools norm --threads $THREADS -m "-any" ${sample}_GATK.filtered.vcf.gz | \
#vt normalize - -n -r $REF_GENOME | \
#bgzip -@ $THREADS -c > ${sample}_GATK.filtered.norm.vcf.gz  && \
#tabix ${sample}_GATK.filtered.norm.vcf.gz

#mv ${sample}_GATK.filtered.vcf.gz* tmp/

# Filter only pass
#bcftools view -f PASS ${sample}_GATK.filtered.norm.vcf.gz -o ${sample}_GATK.filtered.norm.pass.vcf

#mv ${sample}_GATK.filtered.norm.vcf.gz* tmp/

#-------------------------- GATK Variant Annoation ---------------------------#

# Annotate with SnpEff
#java -jar $SNPEFF_JAR ann -v hg19 -canon ${sample}_GATK.filtered.norm.pass.vcf | \
#bcftools view --threads $THREADS -Oz -o ${sample}_GATK.filtered.norm.snpeff.vcf

#bgzip ${sample}_GATK.filtered.norm.snpeff.vcf
#bcftools index ${sample}_GATK.filtered.norm.snpeff.vcf.gz

#mv ${sample}_GATK.filtered.norm.pass.vcf tmp/

# Annotate with Clinvar
bcftools annotate --threads $THREADS -a $CLINVAR_VCF \
	-c CLNHGVS,CLNSIGCONF,ALLELEID,RS \
	-o ${sample}_GATK.filtered.norm.snpeff.clnvar.vcf ${sample}_GATK.filtered.norm.snpeff.vcf.gz

mv ${sample}_GATK.filtered.norm.snpeff.vcf.gz* tmp/
exit 0
# Annotate with gnomAD quality control fields
bcftools annotate --threads $THREADS \
	-a $GNOMAD_VCF \
	-c FILTER,INFO/segdup,INFO/lcr,INFO/rf_tp_probability,INFO/InbreedingCoeff \
	-o ${sample}_GATK_gnomad_qc.vcf \
	${sample}_GATK.filtered.norm.snpeff.vcf


# Convert to avinput
#perl convert2annovar.pl --format vcf4 \
#	--includeinfo \
#	--allsample \
#	--withfreq \
#	${sample}_GATK.vcf > ${sample}_GATK.avinput

#mv ${sample}_GATK.vcf tmp/

ln -s /home/administrator/Annovar/annovar/*.pl .

# Intervar/Annovar annotation
LCVar.py -b hg19 \
	-i ${sample}_GATK.avinput --input_type=AVinput \
	-o ${sample}_GATK.intervar \
	-t $INTERVARDB \
	-d $HUMANDB

#-------------------------- GATK Variants Processing ---------------------------#

# Convert 1->chr1 in intervar output file
python LCVarConv.py ${sample}_GATK.intervar.hg19_multianno.txt.intervar ${sample}_GATK.intervar.hg19_multianno.txt.chr.intervar

mv ${sample}_GATK.intervar.hg19_multianno.txt.intervar tmp/
mv ${sample}_GATK.intervar.hg19_multianno.txt.grl_p tmp/
mv ${sample}_GATK.avinput tmp/
mv ${sample}_GATK.intervar.hg19_multianno.txt tmp/

# Extract vcf headers
bcftools view -h tmp/${sample}_GATK.vcf > ${sample}_header.tmp

# Add INFO fields to vcf headers
awk '
/^##INFO/ && !added_info {
	print;
	print "##INFO=<ID=AVINPUTCHR,Number=1,Type=String,Description=\"Original ANNOVAR input chromosome\">";
	print "##INFO=<ID=AVINPUTSTART,Number=1,Type=Integer,Description=\"Original ANNOVAR input start position\">";
	print "##INFO=<ID=AVINPUTEND,Number=1,Type=Integer,Description=\"Original ANNOVAR input end position\">";
	print "##INFO=<ID=AVINPUTREF,Number=1,Type=String,Description=\"Original ANNOVAR input reference allele\">";
	print "##INFO=<ID=AVINPUTALT,Number=1,Type=String,Description=\"Original ANNOVAR input alternate allele\">";
	added_info=1;
	next;
}
{ print }
' ${sample}_header.tmp > ${sample}_header.txt
rm ${sample}_header.tmp

# Convert avinput to vcf
awk '{print $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $7 "\t" $15 "\t" $16 ";" "AVINPUTCHR=" $1 ";" "AVINPUTSTART=" $2 ";" "AVINPUTEND=" $3 ";" "AVINPUTREF=" $4 ";" "AVINPUTALT=" $5 ";" "\t" $17 "\t" $18}' tmp/${sample}_GATK.avinput > ${sample}_GATK.avinput.tmp

# Merge vcf headers to avinput (converted2vcf)
cat ${sample}_header.txt ${sample}_GATK.avinput.tmp > ${sample}_GATK.avinput.vcf
mv ${sample}_GATK.avinput.tmp tmp/
mv ${sample}_header.txt tmp/

# Snisift Info Extraction
conda run -n SNPSIFT SnpSift extractFields ${sample}_GATK.avinput.vcf CHROM POS REF ALT "ANN[0].GENE" "ANN[0].FEATUREID" "ANN[0].HGVS_P" "ANN[0].HGVS_C" "ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].RANK" DP AF "GEN[0].AD" CLNHGVS CLNSIGCONF ALLELEID FILTER RS AVINPUTSTART AVINPUTEND AVINPUTREF AVINPUTALT >  ${sample}_GATK.snpsift.tmp

# Re-order with avinput coordinates
awk -F'\t' '{print $1 "\t" $20 "\t" $21 "\t" $22 "\t" $23 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $19}' ${sample}_GATK.snpsift.tmp > ${sample}_GATK.snpsift.tsv

mv ${sample}_GATK.snpsift.tmp tmp/
mv ${sample}_GATK.avinput.vcf tmp/

# Merge Snpsift & Intervar
python LCVarMrg.py ${sample}_GATK.intervar.hg19_multianno.txt.chr.intervar ${sample}_GATK.snpsift.tsv ${sample}_GATK_merged.tsv ${sample}_GATK.unmatched.tsv

# Split Intervar Column to ACMG & ACMG Rules
python LCVarSplit.py ${sample}_GATK_merged.tsv ${sample}_GATK_merged_split.tsv

mv ${sample}_GATK.unmatched.tsv tmp/
mv ${sample}_GATK.intervar.hg19_multianno.txt.chr.intervar tmp/
mv ${sample}_GATK_merged.tsv tmp/
mv ${sample}_GATK.snpsift.tsv tmp/
${sample}_GATK.filtered.norm.snpeff.vcf.gz
#-------------------------- Variant Prioritization 1 ---------------------------#

# Prioritize Variants
python LCVarPrio.py -i ${sample}_GATK_merged_split.tsv -o ${sample}_variants.prioritized.tsv -n 25000 -f tsv

#-------------------------- MAGI ACMG implementation ---------------------------#

# Implement MAGI Vus sub classification
python LCVarMagi.py ${sample}_variants.prioritized.tsv ${sample}_variants.prioritized.magi.tsv

#-------------------------- Variant Prioritization 2 ---------------------------#
exit 0














#-------------------------- GAKT & Freebayes Merge ---------------------------#

# Merge GATK & Freebayes
#python lcwesmerge2.py ${sample}_GATK_merged_split.tsv ${sample}_freebayes_merged_split.tsv ${sample}_variants.tsv
exit 0
#mv ${sample}_GATK_merged_split.tsv tmp/
#mv ${sample}_freebayes_merged_split.tsv tmp/
#mv snpEff* tmp/
exit 0
#-------------------------- Variant Prioritization ---------------------------#


exit 0

}

# Main script
for fastq1 in $FASTQ_DIR/*_1.fq.gz; do
	fastq2="${fastq1/_1.fq.gz/_2.fq.gz}"
	if [ -f "$fastq2" ]; then
		sample=$(basename "$fastq1" | sed 's/_1.fq.gz//')
		process_sample "$sample" "$fastq1" "$fastq2"
	else
		echo "Warning: No matching read 2 file found for $fastq1"
	fi
done

echo "ALL DONE!"
