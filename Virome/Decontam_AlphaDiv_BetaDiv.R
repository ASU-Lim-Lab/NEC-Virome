#Decontam: identifies contaminants. Input: species output table from MEGAN and metadata file. Output: csv file with species labeled as contaminants ("TRUE") or not contaminants ("FALSE").
library(vegan)
set.seed(100)
data <- read.delim("SpeciesTable.txt", row.names=1)
data<-as.matrix(t(data))
metadata<-read.table("Metadata.txt", sep = "\t", header = T, row.names = 1, stringsAsFactors = TRUE, comment.char = "")
contam<-isContaminant(data, method = 'prevalence', neg =metadata$Neg, threshold=0.1)
write.csv(contam, "Decontam.csv")

#Remove contaminants and false positives from species table.
#Normalize species counts: (79,000/total QC reads of sample)(counts of virus species in sample).

#Shannon diversity. Input: normalized species table. Output: txt file with Shannon diversity for each sample.
library(vegan)
data<-read.delim("SpeciesPostDecontamNormalized.txt", row.names = 1)
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

#Bray-Curtis dissimilarity. Input: species presence/absence table. Output: Bray-Curtis dissimilarity matrix.
library(vegan)
data<-read.delim("SpeciesPresAbs.txt", row.names = 1)
dataTransposed<-t(data)
dis <- vegdist(dataTransposed, method = "bray")
dis2<-as.matrix(dis)
write.table(dis2,"BrayCurtis.txt", sep = '\t')

#PCoA. Input: Species presence/absence table and metadata file. Output: PCoA plot.
library(phyloseq)
library(ggplot2)
feature_table<-read.delim("SpeciesPresAbs.txt", row.names = 1)
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
                                   shape = "PMA", # metadata variable
                                   color = "PMA", # metadata variable
                                   axes = c(1,2),
                                   title='Bray Curtis Axis1 vs Axis2') +
  theme_bw() +
  theme(text=element_text(size=20))+
  geom_point(size = 8)+
  scale_color_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a"))+
  scale_shape_manual(values = c( 19, 19, 19, 19, 19, 19, 19, 19, 19, 19))
PCoA_Plot