path="/mnt/g/arabidopsis/exp_dna/"
ids=($(find $path -name "*_1.fastq" -type f))
ref_file="/mnt/g/arabidopsis/ref/genomic.fna"
#samtools faidx $ref_file
#java -jar picard.jar CreateSequenceDictionary -R $ref_file -O /mnt/g/arabidopsis/ref/genomic.dict
for id in "${ids[@]}"
do
	bn=$(basename $id)
        pr="${bn%_1.fastq}"
	fastp -i $path$pr"_1.fastq" -I $path$pr"_2.fastq" -o $path$pr"_1.qc.fastq" -O $path$pr"_2.qc.fastq"
	bwa mem -t 4 -M -R "@RG\tID:$pr\tLB:RANDOM\tPL:ILLUMINA\tPU:unit1\tSM:run1" $ref_file $path$pr"_1.qc.fastq" $path$pr"_2.qc.fastq" > $path"out/"$pr".sam"
	samtools view -S -b $path"out/"$pr".sam" > $path"out/"$pr".bam" 
	samtools sort -l4 -o $path"out/"$pr".sort.bam" -T $path"out/"$pr -@4 $path"out/"$pr".bam"
	MarkDuplicates -I $path"out/"$pr".sort.bam" -O $path"out/"$pr".sort.markdup.bam" -M $path"out/"$pr".metrics.txt" -REMOVE_DUPLICATES true
	GenomeAnalysisTK -T RealignerTargetCreator -R $ref_file -I $path"out/"$pr".sort.markdup.bam" -o $path"out/"$pr".indelrealign.intervals"
	GenomeAnalysisTK -T IndelRealigner -R $ref_file -I $path"out/"$pr".sort.markdup.bam" -o $path"out/"$pr".realign.bam" -targetIntervals $path"out/"$pr".indelrealign.intervals"
	samtools index $path"out/"$pr".realign.bam"
	freebayes -f $ref_file $path"out/"$pr".realign.bam" >$path"out/"$pr".vcf"
	vcftools --vcf $path"out/"$pr".vcf" --min-alleles 2 --max-alleles 2 --max-missing 0.9 --mac 10 --min-meanDP 3 --recode --recode-INFO-all --minQ 20 --out $path"out/"$pr
	vcftools --vcf $path"out/"$pr".recode.vcf" --remove-indels --recode --recode-INFO-all --out $path"out/"$pr".noindel"
done
