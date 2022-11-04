
rule evaluate_mappings:
    input:
        touch("data/mapped_reads/{sample}/mapping_microbial_reads.done"),
         config["Rlibpath"],
        "data/kraken2_classification/{sample}.report",
        read_lengths = "data/read_lengths/{sample}/read_lengths.tsv"
    output:
        "data/tax_filtered_report/{sample}.report"
    script:
        "code/evaluate_mappings.R"


gunzip -r "/user_data/men/sepseq/databases/2022_10_18_kraken2_EUPATH_database/genomes/"

grep -r ">" *.fna | sed -r 's/[:>]+/\t/g' > /user_data/men/sepseq/clinical_studies/init_tests/kraken2_filter/data/mapped_reads/accession2genome.tsv


rule extract_read_lengths:
    input:
        reads = config["reads"]+"{sample}.fastq",
    output:
        read_lengths = "data/read_lengths/{sample}/read_lengths.tsv"
    params:
        seqkit = config["Seqkit"]
    shell:
        """    
        module load {params.seqkit}
        seqkit fx2tab -nl {input.reads} > {output.read_lengths}
        module purge
        """


##Workaround duplicate sam headers
WD=/user_data/men/sepseq/clinical_studies/init_tests/kraken2_filter
#Reads_50
WD=/user_data/men/sepseq/clinical_studies/init_tests/kraken2_filter
cd $WD
module load SAMtools/1.14-GCC-10.2.0
mkdir -p $WD/data/modified_sam/reads_50
mkdir -p $WD/data/sorted_sam/reads_50
mkdir -p $WD/data/sam_coverage/reads_50

for file in $WD/data/mapped_reads_header/reads_50/*.sam
do
BASENAME=$(basename $file)
awk 'BEGIN { i = 0; } /^@/ { if (/NZ_LR590469.1/) { if (i++ < 1) { print; } } else { print } } /^[^@]/ { print }' $file | awk 'BEGIN { i = 0; } /^@/ { if (/NZ_LR134116.1/) { if (i++ < 1) { print; } } else { print } } /^[^@]/ { print }' | awk 'BEGIN { i = 0; } /^@/ { if (/NZ_CP043727.1/) { if (i++ < 1) { print; } } else { print } } /^[^@]/ { print }' | awk 'BEGIN { i = 0; } /^@/ { if (/NZ_CP054043.1/) { if (i++ < 1) { print; } } else { print } } /^[^@]/ { print }' | awk 'BEGIN { i = 0; } /^@/ { if (/NZ_CP011975.1/) { if (i++ < 1) { print; } } else { print } } /^[^@]/ { print }' | awk 'BEGIN { i = 0; } /^@/ { if (/NZ_CP032487.1/) { if (i++ < 1) { print; } } else { print } } /^[^@]/ { print }' | awk 'BEGIN { i = 0; } /^@/ { if (/NZ_CP032488.1/) { if (i++ < 1) { print; } } else { print } } /^[^@]/ { print }' | awk 'BEGIN { i = 0; } /^@/ { if (/NC_003551.1/) { if (i++ < 1) { print; } } else { print } } /^[^@]/ { print }' > $WD/data/modified_sam/reads_50/${BASENAME}
done

for file in $WD/data/modified_sam/reads_50/*.sam
do
BASENAME=$(basename $file)
samtools sort -O sam -T sample50.sort -o $WD/data/sorted_sam/reads_50/${BASENAME} $file
done

for file in $WD/data/sorted_sam/reads_50/*.sam
do
BASENAME=$(basename $file)
samtools coverage $file > $WD/data/sam_coverage/reads_50/${BASENAME%.sam}\.tsv
done

for file in $WD/data/sam_coverage/reads_50/*.tsv
do
BASENAME=$(basename $file)
awk '$4 != 0' $file > $WD/data/sam_coverage/reads_50/${BASENAME%.tsv}\min.tsv
done


#Reads_70
WD=/user_data/men/sepseq/clinical_studies/init_tests/kraken2_filter/
cd $WD
module load SAMtools/1.14-GCC-10.2.0
mkdir -p $WD/data/modified_sam/reads_70/
mkdir -p $WD/data/sorted_sam/reads_70/
mkdir -p $WD/data/sam_coverage/reads_70/
for file in $WD/data/mapped_reads_header/reads_70/*.sam
do
BASENAME=$(basename $file)
awk 'BEGIN { i = 0; } /^@/ { if (/NZ_LR590469.1/) { if (i++ < 1) { print; } } else { print } } /^[^@]/ { print }' $file | awk 'BEGIN { i = 0; } /^@/ { if (/NZ_LR134116.1/) { if (i++ < 1) { print; } } else { print } } /^[^@]/ { print }' | awk 'BEGIN { i = 0; } /^@/ { if (/NZ_CP043727.1/) { if (i++ < 1) { print; } } else { print } } /^[^@]/ { print }' | awk 'BEGIN { i = 0; } /^@/ { if (/NZ_CP054043.1/) { if (i++ < 1) { print; } } else { print } } /^[^@]/ { print }' | awk 'BEGIN { i = 0; } /^@/ { if (/NZ_CP011975.1/) { if (i++ < 1) { print; } } else { print } } /^[^@]/ { print }' | awk 'BEGIN { i = 0; } /^@/ { if (/NZ_CP032487.1/) { if (i++ < 1) { print; } } else { print } } /^[^@]/ { print }' | awk 'BEGIN { i = 0; } /^@/ { if (/NZ_CP032488.1/) { if (i++ < 1) { print; } } else { print } } /^[^@]/ { print }' | awk 'BEGIN { i = 0; } /^@/ { if (/NC_003551.1/) { if (i++ < 1) { print; } } else { print } } /^[^@]/ { print }' > $WD/data/modified_sam/reads_70/${BASENAME}
done

for file in $WD/data/modified_sam/reads_70/*.sam
do
BASENAME=$(basename $file)
samtools sort -O sam -T sample70.sort -o $WD/data/sorted_sam/reads_70/${BASENAME} $file
done

for file in $WD/data/sorted_sam/reads_70/*.sam
do
BASENAME=$(basename $file)
samtools coverage $file > $WD/data/sam_coverage/reads_70/${BASENAME%.sam}\.tsv
done


for file in $WD/data/sam_coverage/reads_70/*.tsv
do
BASENAME=$(basename $file)
awk '$4 != 0' $file > $WD/data/sam_coverage/reads_70/${BASENAME%.tsv}\min.tsv
done

#Reads_90
WD=/user_data/men/sepseq/clinical_studies/init_tests/kraken2_filter/
cd $WD
module load SAMtools/1.14-GCC-10.2.0
mkdir -p $WD/data/modified_sam/reads_90/
mkdir -p $WD/data/sorted_sam/reads_90/
mkdir -p $WD/data/sam_coverage/reads_90/
for file in $WD/data/mapped_reads_header/reads_90/*.sam
do
BASENAME=$(basename $file)
awk 'BEGIN { i = 0; } /^@/ { if (/NZ_LR590469.1/) { if (i++ < 1) { print; } } else { print } } /^[^@]/ { print }' $file | awk 'BEGIN { i = 0; } /^@/ { if (/NZ_LR134116.1/) { if (i++ < 1) { print; } } else { print } } /^[^@]/ { print }' | awk 'BEGIN { i = 0; } /^@/ { if (/NZ_CP043727.1/) { if (i++ < 1) { print; } } else { print } } /^[^@]/ { print }' | awk 'BEGIN { i = 0; } /^@/ { if (/NZ_CP054043.1/) { if (i++ < 1) { print; } } else { print } } /^[^@]/ { print }' | awk 'BEGIN { i = 0; } /^@/ { if (/NZ_CP011975.1/) { if (i++ < 1) { print; } } else { print } } /^[^@]/ { print }' | awk 'BEGIN { i = 0; } /^@/ { if (/NZ_CP032487.1/) { if (i++ < 1) { print; } } else { print } } /^[^@]/ { print }' | awk 'BEGIN { i = 0; } /^@/ { if (/NZ_CP032488.1/) { if (i++ < 1) { print; } } else { print } } /^[^@]/ { print }' | awk 'BEGIN { i = 0; } /^@/ { if (/NC_003551.1/) { if (i++ < 1) { print; } } else { print } } /^[^@]/ { print }' > $WD/data/modified_sam/reads_90/${BASENAME}
done

for file in $WD/data/modified_sam/reads_90/*.sam
do
BASENAME=$(basename $file)
samtools sort -O sam -T sample.sort -o $WD/data/sorted_sam/reads_90/${BASENAME} $file
done

for file in $WD/data/sorted_sam/reads_90/*.sam
do
BASENAME=$(basename $file)
samtools coverage $file > $WD/data/sam_coverage/reads_90/${BASENAME%.sam}\.tsv
done

for file in $WD/data/sam_coverage/reads_90/*.tsv
do
BASENAME=$(basename $file)
awk '$4 != 0' $file > $WD/data/sam_coverage/reads_90/${BASENAME%.tsv}\min.tsv
done
