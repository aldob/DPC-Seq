#!/bin/bash

#SBATCH --cpus-per-task=9
#SBATCH --mem=100G
#SBATCH --time=0-10:00:00
#SBATCH -p general

rawdir=$PWD # rawdir is the location the script was run from and should only contain the bcl directory
readarray -t sample_ID < $rawdir/samples_$1.txt # instead of arguments, use a table of samples





################################################################# Initial setup and fastq unpacking ####################################################################

echo -e Setting up directory structure"\n"

mkdir $rawdir/log
mkdir $rawdir/processed
mkdir $rawdir/rawdata
mkdir $rawdir/rawdata/fastq
mkdir $rawdir/rawdata/fastqc
mkdir $rawdir/rawdata/pe_int
mkdir $rawdir/processed/covs
mkdir $rawdir/processed/wigs





echo -e "\n"Input bcl directory:
read bcldir # user inputs the bcl folder
echo -e Generating fastqs"\n"
bcl2fastq -R $bcldir -o $rawdir/rawdata/fastq -r 2 -w 2 -p 10 --no-lane-splitting --ignore-missing-bcls --barcode-mismatches 1
rm $rawdir/rawdata/fastq/Undetermined_* # remove the undetermined fastq files
echo -e Unzipping fastqs"\n"
for file in $rawdir/rawdata/fastq/*.gz # simultaneously unzip all fastq fiiles
do
  gzip -d ${file} & done
wait





############################################################### QC and trimming ####################################################################

echo -e Removing low quality reads and QC checking"\n"
for sample in ${sample_ID[@]}
do
  echo ${sample}
  fastp -g -w 4 -q 15 -u 50 -l 5 -A -h $rawdir/rawdata/fastqc/${sample}_fastp.html -i $rawdir/rawdata/fastq/${sample}_S_R1_001.fastq -I $rawdir/rawdata/fastq/${sample}_S_R2_001.fastq --out1 $rawdir/rawdata/fastq/${sample}_qc1.fq --out2 $rawdir/rawdata/fastq/${sample}_qc2.fq
done
wait





######################################################### Alignment and sorting ####################################################################

echo -e Alignment results"\n\n" >> $rawdir/log/alignments.log
echo -e "\n\n\n"Alignment results"\n"
for sample in ${sample_ID[@]} # align fastqs and pipe straight to samtools sort, outputting the bowtie alignment results to a log file
do
   (bowtie2 -t -p 8 --fr --maxins 2000 --dovetail -x /mnt/scratcha/sjlab/abader/genomes/hg38/hg38 \
   -1 $rawdir/rawdata/fastq/${sample}_qc1.fq -2 $rawdir/rawdata/fastq/${sample}_qc2.fq \
   | samtools sort -o $rawdir/rawdata/${sample}nsort.bam -n -O bam -l 0 -@ 8 -m 10G) \
   2> $rawdir/log/${sample}_align.log 

   echo -e ${sample} aligned successfully
   echo -e "\n"${sample} alignment: >> $rawdir/log/alignments.log
   tail -n 19 $rawdir/log/${sample}_align.log >> $rawdir/log/alignments.log
done
wait  



echo -e "\n"Indexing
for sample in ${sample_ID[@]} # simultaneously index all alignments 
do
  samtools index $rawdir/rawdata/${sample}nsort.bam & done
wait





######################################################### Coverage ####################################################################

echo -e "\n""\n"Running bamCoverage"\n"
for sample in ${sample_ID[@]}
do
    bamCoverage --bam $rawdir/rawdata/${sample}sort.bam -o $rawdir/processed/wigs/${sample}.bw -bs 20 --effectiveGenomeSize 2913022398 -p 1 --normalizeUsing CPM --extendReads 300 --maxFragmentLength 2000 & done
wait



dirsamps=(./chipseq/dpcseq1/processed/wigs/wt_0_merge.bw ./chipseq/dpcseq1/processed/wigs/wt_6_merge.bw ./chipseq/dpcseq1/processed/wigs/csb_0_merge.bw ./chipseq/dpcseq1/processed/wigs/csb_6_merge.bw)

computeMatrix scale-regions -p 8 -S ${dirsamps[*]} -R ./chipseq/dpcseq1/havana_transcripts2_mbt_csbd.bed -o ./chipseq/dpcseq1/processed/wigs/all_merge_csbd.tab -b 2000 -a 2000 -m 5000
plotProfile -m ./chipseq/dpcseq1/processed/wigs/all_merge_csbd.tab -o ./chipseq/dpcseq1/processed/wigs/new/all_merge_csbd.png --perGroup --plotWidth 9 --plotHeight 6 --dpi 1000 --colors '#00aaff' '#005ef6' '#ffad7b' '#ff5e00'



echo -e converting to bedpe"\n"
for sample in ${sample_ID[@]}
do
   bedtools bamtobed -bedpe -i $rawdir/rawdata/pe_int/${sample}nsort.bam > $rawdir/rawdata/pe_int/${sample}.bedpe & done
wait
echo -e cutting bedpe to bed"\n"
for sample in ${sample_ID[@]}
do
    $rawdir/bedpeTobed.py $rawdir/rawdata/pe_int/${sample}.bedpe $rawdir/rawdata/pe_int/${sample}cut.bed & done
wait
echo -e "\n""\n"Calculating bedcoverage"\n"
for sample in ${sample_ID[@]}
do
    bedtools coverage -a $rawdir/havana_transcripts2_mbt_csbd.bed -b $rawdir/rawdata/pe_int/${sample}_cut.bed > $rawdir/processed/covs/${sample}.cov & done
wait



echo -e calculating readcounts"\n"
>$rawdir/log/readcounts_total.log
for sample in ${sample_ID[@]}
do
   readcount=$(samtools view -c $rawdir/rawdata/pe_int/${sample}nsort.bam)
   echo -e ${sample},$readcount>>$rawdir/log/readcounts_total.log
   echo -e $readcount>>$rawdir/log/readcounts_$1.log
done
readarray -t readcounts < $rawdir/log/readcounts_$1.log

echo -e "\n""\n""\n"normalising covs"\n"
for ((i=0;i<${#sample_ID[@]};i++))
do
    awk -v x=${readcounts[$i]} '{print $1, $2, (($3+1)/(x/1000000))}' $rawdir/processed/covs/${sample}.cov > $rawdir/processed/covs/${sample}_norm.cov & done
wait







