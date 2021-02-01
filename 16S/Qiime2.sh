#Import sample fastq files to Qiime2
qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path ctl/ --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path seqs.qza

#Inspect quality. 
qiime demux summarize --i-data seqs.qza --o-visualization seqs.qzv

#Trim reads in fastq files using bbduk.
bbduk.sh in=SampleName.fastq out=trimmed_SampleName.fastq ftl=12 qtrim=r trimq=20 minavgquality=20 minlen=350 overwrite=t 1> SampleName.log.txt 2>&1;

#Import trimmed fastq files to Qiime2. Denoise using dada2.
qiime dada2 denoise-pyro --i-demultiplexed-seqs trimmed_seqs.qza --p-trunc-len 0 --o-table feature-table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza

#Create phylogenetic tree.
 qiime phylogeny align-to-tree-mafft-fasttree --i-sequences trimmed_seqs.qza --o-alignment alignment.qza  --o-masked-alignment masked-alignment.qza --o-tree tree.qza --o-rooted-tree rooted-tree.qza

 #Obtain taxonomy.
 qiime feature-classifier classify-sklearn --i-classifier classifier.qza --i-reads trimmed_seqs.qza --o-classification taxonomy.qza

 #Obtain core diversity metrics.
 qiime diversity core-metrics-phylogenetic --i-phylogeny  rooted-tree.qza --i-table  feature-table.qza --p-sampling-depth 2500 --m-metadata-file metadata.txt --output-dir diversity-metrics

 #Create relative abundance plots.
 qiime taxa barplot --i-table feature-table.qza --i-taxonomy taxonomy.qza --m-metadata-file metadata.txt --o-visualization barplot.qzv