#Decontam: identifies contaminants. Input: contig count matrix and metadata file. Output: csv file with contigs labeled as contaminants ("TRUE") or not contaminants ("FALSE").
library(decontam)
set.seed(100)
data <- read.delim("RawContigCounts_overKb.txt", row.names=1)
metadata<-read.delim("MetadataPreDecontam.txt", row.names = 1)
data<-as.matrix(t(data))
contam<-isContaminant(data, method = 'prevalence', neg =metadata$isNeg, threshold=0.1)
write.csv(contam, "Decontam.csv")

#Remove contaminants, false positives, and contigs with high percent identity to human sequences from contig count matrix.
#Normalize contig counts (RPK): (raw counts of a contig in sample)(79,000/total QC reads of sample)/(contig length in kb).
#Mask normalized counts <0.5. Convert to presence/absence for unweighted beta diversity analyses.

#Shannon diversity. Input: normalized contig counts table. Output: txt file with Shannon diversity for each sample.
library(vegan)
data<-read.delim("ContigCountsRPK.txt", row.names = 1)
dataTransposed<-t(data)
metadata<-read.delim("MetadataPostDecontam.txt")
data_vegan<-diversity(dataTransposed, index = 'shannon')
write.table(data_vegan,"ShannonDiversity.txt", sep = '\t')

#Linear mixed modeling. Input: txt file with Shannon diversity or richness, infant ID, and PMA for each sample.
library(nlme)
data <- read.delim("ShannonDiversity_ControlInfants.txt")
nlme.test <- lme(Shannon ~ PMA, random =~1|InfantID, data = data)
summary(nlme.test)

#Loess. Input: txt file with Shannon diversity or richness and PMA for each sample. Output: loess plot showing Shannon diversity or richness over time.
library(ggplot2)
data <- read.delim("ShannonDiversity_CaseInfants.txt")
loessPlot <- ggplot(data, aes(x = PMA, y = Alpha_div)) +
  geom_point() +
  stat_smooth(method="loess", se=TRUE, span=0.25, level=0.95) +
  scale_x_continuous(breaks=c(25,27,29,31,33,35)) +
  theme_bw()
loessPlot

#Sorensen dissimilarity. Input: contig presence/absence table. Output: Sorensen dissimilarity matrix.
library(vegan)
data<-read.delim("ContigPresAbs.txt", row.names = 1)
dataTransposed<-t(data)
dis <- vegdist(dataTransposed, method = "bray")
dis2<-as.matrix(dis)
write.table(dis2,"Sorensen.txt", sep = '\t')

#Bray-Curtis dissimilarity. Input: contig counts table. Output: Bray-Curtis dissimilarity matrix.
library(vegan)
data1<-read.delim("ContigsCountsRPK.txt", row.names = 1)
dataTransposed1<-t(data1)
dis <- vegdist(dataTransposed1, method = "bray")
dis2<-as.matrix(dis)
write.table(dis2,"WeightedBrayCurtis.txt", sep = '\t')

#Hellinger distance. Input: log-transformed contig count table (transformation = log10(x+1), where x = RPK count). Output: Hellinger distance matrix.
library(adespatial)
data5 <-read.delim("ContigCountsRPK_logTransf.txt", row.names=1)
dataTransposed5 <-t(data5)
dist.hel <-dist.ldc(dataTransposed5, method = "hellinger")
dist2 <-as.matrix(dist.hel)
write.table(dist2, "Hellinger.txt", sep = '\t')

#PCoA. Input: Contig count table and metadata file. Output: PCoA plot.
library(phyloseq)
library(ggplot2)
feature_table<-read.delim("ContigCountsRPK.txt", row.names = 1)
feature_table <- as.data.frame(feature_table)
dim(feature_table)
OTU=otu_table(feature_table,taxa_are_rows =TRUE)
sample_names(OTU)
metadata<-import_qiime_sample_data("MetadataPostDecontam.txt")
sample_names(metadata)
physeq =             phyloseq(OTU,metadata)
ord <- ordinate(physeq,
                method = "PCoA",
                distance = "bray")
PCoA_Plot  <- plot_ordination(physeq = physeq,
                                   ordination = ord,
                                   shape = "PMA_binned", # metadata variable
                                   color = "PMA_continuous", # metadata variable
                                   axes = c(1,2),
                                   title='Bray Curtis Axis1 vs Axis2') +
  theme_bw() +
  theme(text=element_text(size=20))+
  geom_point(size = 8)+
  scale_color_viridis_c(direction = -1)+
  scale_shape_manual(values = c( 19, 19, 19, 19, 19, 19, 19, 19, 19, 19))
PCoA_Plot

#PERMANOVA. Input: Contig count table and metadata file. PMA is a continuous variable. Case/control status and subject ID are categorical.
library(vegan)
set.seed(1)
data1<-read.delim("ContigsCountsRPK_controls.txt", row.names = 1)
dataTransposed1<-t(data1)
dist.1 <- vegdist(dataTransposed1, method = "bray")
metadata <- read.delim("Metadata_controls.txt")
adonis.test1.cont.subjID <- adonis(dist.1 ~ PMA, data = metadata, permutations = 999, strata = metadata$SubjectID)
adonis.test1.cont.subjID ##Gives p values for PMA, see figures 1 and 2.

data2<-read.delim("ContigsCountsRPK.txt", row.names = 1)
dataTransposed2<-t(data2)
dist.2 <- vegdist(dataTransposed2, method = "bray")
metadata2 <- read.delim("MetadataPostDecontam.txt")
adonis.test2 <- adonis(dist.2 ~ Case_or_control, data = metadata2, permutations = 999, strata = metadata2$SubjectID)
adonis.test2 ##Gives p values for case/control status, see figure 2.
