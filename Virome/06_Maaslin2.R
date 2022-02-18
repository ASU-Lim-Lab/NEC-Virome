library(Maaslin2)
library(gplots)
library(RColorBrewer)
##Input: data = contig count table and metadata = bacterial count table.
##Filtered to only include contigs and bacterial genera present in at least 10% of samples. Infant ID set as random effect.

#1: NEC-associated contig interactions in case infants (Figure 4C, left heatmap).
data.cases.lefse<-read.delim("Contig_Counts_0.5_cases_lefse.txt", row.names = 1)
metadata.cases<-read.delim("16S-cases-filtered.txt", row.names = 1)
fitData.cases.lefse= Maaslin2(
  input_data = data.cases.lefse, 
  input_metadata = metadata.cases,
  output = "cases_lefse",
  normalization = "none",
  transform = "none",
  fixed_effects = c("f__Enterobacteriaceae.__", "g__Staphylococcus", "g__Klebsiella", "g__Escherichia", "f__Enterococcaceae.__", "g__Streptococcus", "g__Propionibacterium", "g__Clostridium", "g__Citrobacter", "g__Corynebacterium", "g__Dialister", "g__Proteus", "g__Acinetobacter", "g__Bifidobacterium", "g__Lactobacillus", "g__Haemophilus"),
  random_effects = c("InfantID"))

#2: Control time-associated contig interactions in control infants (Extended Data Figure 4D, left heatmap).
data.ctrls.lefse<-read.delim("Contig_Counts_0.5_ctrls_lefse.txt", row.names = 1)
metadata.ctrls<-read.delim("16S-level-6-ctrls-filtered-0.1.txt", row.names = 1)
fitData.ctrls.lefse= Maaslin2(
  input_data = data.ctrls.lefse, 
  input_metadata = metadata.ctrls,
  output = "level6_ctrls_lefse",
  normalization = "none",
  transform = "none",
  fixed_effects = c("f__Enterobacteriaceae.__", "f__Enterococcaceae.__", "g__Enterococcus", "g__Staphylococcus", "g__Veillonella", "g__Escherichia", "g__Clostridium", "g__Streptococcus", "f__Pseudomonadaceae.__", "g__Actinomyces", "g__Corynebacterium", "g__Propionibacterium", "g__Lactobacillus", "f__Clostridiaceae.__", "f__Peptostreptococcaceae.g__.Clostridium.", "g__Klebsiella", "g__Haemophilus"),
  random_effects = c("InfantID"))

#3: NEC-associated contig interactions in control infants (Figure 4C, right heatmap).
data.ctrls.lefse2<-read.delim("case_lefse_contig_counts_in_ctrls.txt", row.names = 1)
metadata.ctrls<-read.delim("16S-level-6-ctrls-filtered-0.1.txt", row.names = 1)
fitData.ctrls.lefse= Maaslin2(
  input_data = data.ctrls.lefse2, 
  input_metadata = metadata.ctrls,
  output = "level6_case_lefse_in_ctrls",
  normalization = "none",
  transform = "none",
  fixed_effects = c("f__Enterobacteriaceae.__", "f__Enterococcaceae.__", "g__Enterococcus", "g__Staphylococcus", "g__Veillonella", "g__Escherichia", "g__Clostridium", "g__Streptococcus", "f__Pseudomonadaceae.__", "g__Actinomyces", "g__Corynebacterium", "g__Propionibacterium", "g__Lactobacillus", "f__Clostridiaceae.__", "f__Peptostreptococcaceae.g__.Clostridium.", "g__Klebsiella", "g__Haemophilus"),
  random_effects = c("InfantID"))

#4: Control time-associated contig interactions in case infants (Extended Data Figure 4D, right heatmap).
data.cases.lefse2<-read.delim("ctrl_lefse_in_cases.txt", row.names = 1)
metadata.cases<-read.delim("16S-level-6-cases-filtered-0.1.txt", row.names = 1)
fitData.cases.lefse= Maaslin2(
  input_data = data.cases.lefse2, 
  input_metadata = metadata.cases,
  output = "level6_ctrl_lefse_in_cases",
  normalization = "none",
  transform = "none",
  fixed_effects = c("f__Enterobacteriaceae.__", "g__Staphylococcus", "g__Klebsiella", "g__Escherichia", "f__Enterococcaceae.__", "g__Streptococcus", "g__Propionibacterium", "g__Clostridium", "g__Citrobacter", "g__Corynebacterium", "g__Dialister", "g__Proteus", "g__Acinetobacter", "g__Bifidobacterium", "g__Lactobacillus", "g__Haemophilus"),
  random_effects = c("InfantID"))

#Plot heatmaps. Repeat for each Maaslin2 analysis. Input: table with coefficient values of significant Maaslin2 comparisons. Columns = bacterial genera and rows = viral contigs.
data <- as.matrix(read.delim("cases-Maaslin2-coefficients.txt", header=TRUE, row.names = 1))
heatmap <- heatmap.2(data, Rowv = TRUE, Colv = TRUE, trace = "none", dendrogram = "both", col = c("#053061","#2166ac","#4393c3","#92c5de","#d1e5f0","#ffffff", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f"), breaks = c(-1100,-400,-20,-10,-5,-0.01,0.01,5,10,20,400,1100), key = FALSE)

