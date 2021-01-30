#Build contigs from sample fastq files.
mkdir metaSPAdesContigs;
metaspades.py -o SampleName --pe1-1 SampleName-R1.fastq --pe1-2 SampleName-R2.fastq;
awk '/^>/{print ">SampleName_Contig_1" ++i; next}{print}' < SampleName/contigs.fasta > SampleName/SampleName_RenamedContigs.fasta;
cp SampleName/SampleName_RenamedContigs.fasta metaSPAdesContigs/;

#Concatenate contig fasta files. Deduplicate contigs using CD-HIT-EST.
cd-hit-est -i contigs.fasta -o deduplicated-contigs.fasta -c 0.95 -n 10 -G 0 -aS 0.95 -g 1 -r 1 -M 20000 -d 0

#Merge overlapping contigs using minimus2.
toAmos -s deduplicated-contigs.fasta -o deduplicated-contigs.afg
minimus2 deduplicated-contigs -D MINID=95

#Concatenate minimus2 output fasta file (merged contigs) and singletons.seq file (unmerged contigs).

#Filter contigs by length.
bbduk.sh in=merged-deduplicated-contigs.fasta out=filtered-contigs.fasta minlen=800

#Use blastx and taxonomizr (see separate taxonomizr R script) to identify candidate viral contigs.
blastx -db /path/to/ViralDatabase -query filtered-contigs.fasta -evalue 1e-3 -num_threads 60 -max_target_seqs 1 -outfmt "6 qseqid sseqid staxids evalue bitscore pident nident qcovs length mismatch qlen slen" -out blastx-viral-contigs.out

#Run megablast on candidate viral contigs called by blastx.
blastn -task megablast -db /path/to/ntDatabase -query blastx-viral-contigs.fasta -evalue 1e-10 -num_threads 60 -max_target_seqs 1 -outfmt "6 qseqid sseqid staxids evalue bitscore pident nident qcovs length mismatch qlen slen" -out megablast-viral-contigs.out

#Identify candidate bacteriophage contigs from full contig set using VirSorter.
./wrapper_phage_contigs_sorter_iPlant.pl -f filtered-contigs.fasta --virome --wdir output_folder --ncpu 60 --data-dir virsorter-data/

#Use blastx and taxonomizr to obtain taxonomy for candidate bacteriophage contigs identified by VirSorter. Use megablast and taxonomizr to identify VirSorter contigs with high percent identity to human sequences. Remove human contigs from VirSorter contig set.

#Combine contigs called viral by megablast (line 24) and contigs called by VirSorter into one fasta file.

#Index viral contigs.
bwa index ViralContigs.fasta

#Map QC reads from samples to viral contigs.
mkdir /Mapping;
mkdir /Counts;
bwa mem -M -L 97,97 ViralContigs.fasta /path/to/SampleName.fasta >/path/to/Mapping/SampleName-readsMapped.sam;
samtools view -h -F 0x900 /path/to/Mapping/SampleName-readsMapped.sam > /path/to/Mapping/SampleName-secondaryRemoved.sam;
samtools view -h -F 0x4 /path/to/Mapping/SampleName-secondaryRemoved.sam > /path/to/Mapping/SampleName-secondaryUnMappedRemoved.sam;
samtools view -S -b /path/to/Mapping/SampleName-secondaryUnMappedRemoved.sam > /path/to/Mapping/SampleName-secondaryUnMappedRemoved.bam;
samtools sort /path/to/Mapping/SampleName-secondaryUnMappedRemoved.bam > /path/to/Mapping/SampleName-secondaryUnMappedRemoved_sorted.bam;
samtools index /path/to/Mapping/SampleName-secondaryUnMappedRemoved_sorted.bam;
samtools idxstats /path/to/Mapping/SampleName-secondaryUnMappedRemoved_sorted.bam > /path/to/Counts/SampleName-counts.txt;

#Concatenate contig counts from each sample into one table to use for downstream analysis. 
#Normalize contig counts (RPK): (reads mapping to a contig in sample)(79,000/total QC reads of sample)/(contig length in kb). 
#Mask normalized counts <0.5.
#Convert normalized counts >0.5 to presence/absence for calculating Bray-Curtis dissimilarity.
