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

#Query length-filtered contigs against Gut Phage Database to identify candidate viral contigs.
tblastx -db /path/to/GPD -query filtered-contigs.fasta -evalue 1e-3 -num_threads 28 -max_target_seqs 1 -max_hsps 1 -outfmt "6 qseqid sseqid evalue bitscore pident nident qcovs length mismatch qlen slen" -out filtered-contigs-GPD.out

#Query contigs without a hit to GPD against Gut Virome Database.
tblastx -db /path/to/GVD -query filtered-contigs-notGPD.fasta -evalue 1e-3 -num_threads 28 -max_target_seqs 1 -max_hsps 1 -outfmt "6 qseqid sseqid evalue bitscore pident nident qcovs length mismatch qlen slen" -out filtered-contigs-GVD.out

#Combine GPD hits and GVD hits into one fasta file: candidate-viral-contigs.fasta.

#Run megablast on candidate viral contigs called by blastx.
blastn -task megablast -db /path/to/ntDatabase -query candidate-viral-contigs.fasta -evalue 1e-10 -num_threads 60 -max_target_seqs 1 -outfmt "6 qseqid sseqid staxids evalue bitscore pident nident qcovs length mismatch qlen slen" -out megablast-viral-contigs.out

#Index candidate viral contigs.
bwa index candidate-viral-contigs.fasta

#Map QC reads from samples to candidate viral contigs.
mkdir /path/to/Mapping;
mkdir /path/to/Counts;
bwa mem -M -L 97,97 candidate-viral-contigs.fasta /path/to/SampleName.fasta >/path/to/Mapping/SampleName-readsMapped.sam;
samtools view -h -F 0x900 /path/to/Mapping/SampleName-readsMapped.sam > /path/to/Mapping/SampleName-secondaryRemoved.sam;
samtools view -h -F 0x4 /path/to/Mapping/SampleName-secondaryRemoved.sam > /path/to/Mapping/SampleName-secondaryUnMappedRemoved.sam;
samtools view -S -b /path/to/Mapping/SampleName-secondaryUnMappedRemoved.sam > /path/to/Mapping/SampleName-secondaryUnMappedRemoved.bam;
samtools sort /path/to/Mapping/SampleName-secondaryUnMappedRemoved.bam > /path/to/Mapping/SampleName-secondaryUnMappedRemoved_sorted.bam;
samtools index /path/to/Mapping/SampleName-secondaryUnMappedRemoved_sorted.bam;
samtools idxstats /path/to/Mapping/SampleName-secondaryUnMappedRemoved_sorted.bam > /path/to/Counts/SampleName-counts.txt;

#Concatenate contig counts from each sample into one matrix. 
#Divide raw contig counts by contig length in kb. Use decontam (see decontam R script) to identify sequencing contaminants. Remove contaminants (default threshold) from contig count matrix.
#Use taxonomizr (see taxonomizr R script) to obtain taxonomy for megablast output. Remove contigs with hits to human sequences if percent identity and query coverage are high (both pident and qcovs >=90%, and at least one of pident or qcovs >=95%.
#Remove false positive virus families from contig count matrix.

#Normalize contig counts (RPK): (raw counts of a contig in sample)(79,000/total QC reads of sample)/(contig length in kb). 
#Mask normalized counts <0.5.
#Convert normalized counts >0.5 to presence/absence for unweighted beta diversity analyses.
#Log transform normalized counts >0.5 for Hellinger distance. Transformation = log10(x+1), where x = RPK count.
