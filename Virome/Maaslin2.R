library(Maaslin2)
##Input: data = contig count table and metadata = bacterial count table.
##Filtered to only include contigs and bacterial genera present in at least 10% of samples. Infant ID set as random effect.

#1: Input = all contigs present in at least 10% of samples.
data3<-read.delim("ContigCounts_Cases_Filtered.txt", row.names = 1)
metadata10<-read.delim("16S-cases-filtered.txt", row.names = 1)
fitData10= Maaslin2(
  input_data = data3, 
  input_metadata = metadata10,
  output = "cases_output",
  normalization = "none",
  transform = "none",
  fixed_effects = c("f__Enterobacteriaceae.__", "g__Staphylococcus", "g__Klebsiella", "g__Escherichia", "f__Enterococcaceae.__", "g__Streptococcus", "g__Propionibacterium", "g__Clostridium", "g__Citrobacter", "g__Corynebacterium", "g__Dialister", "g__Proteus", "g__Acinetobacter", "g__Bifidobacterium", "g__Lactobacillus", "g__Haemophilus"),
  random_effects = c("InfantID"))

data4<-read.delim("ContigCounts_Controls_filtered.txt", row.names = 1)
metadata11<-read.delim("16S-controls-filtered.txt", row.names = 1)
fitData11= Maaslin2(
  input_data = data4, 
  input_metadata = metadata11,
  output = "controls_output",
  normalization = "none",
  transform = "none",
  fixed_effects = c("f__Enterobacteriaceae.__", "f__Enterococcaceae.__", "g__Enterococcus", "g__Staphylococcus", "g__Veillonella", "g__Escherichia", "g__Clostridium", "g__Streptococcus", "f__Pseudomonadaceae.__", "g__Actinomyces", "g__Corynebacterium", "g__Propionibacterium", "g__Lactobacillus", "f__Clostridiaceae.__", "f__Peptostreptococcaceae.g__.Clostridium.", "g__Klebsiella", "g__Haemophilus"),
  random_effects = c("InfantID"))

#2: Input = LEfSe discriminant contigs
data.cases.lefse<-read.delim("DiscriminantContigCounts_cases.txt", row.names = 1)
metadata.cases<-read.delim("16S-cases-filtered.txt", row.names = 1)
fitData.cases.lefse= Maaslin2(
  input_data = data.cases.lefse, 
  input_metadata = metadata.cases,
  output = "cases_lefse",
  normalization = "none",
  transform = "none",
  fixed_effects = c("f__Enterobacteriaceae.__", "g__Staphylococcus", "g__Klebsiella", "g__Escherichia", "f__Enterococcaceae.__", "g__Streptococcus", "g__Propionibacterium", "g__Clostridium", "g__Citrobacter", "g__Corynebacterium", "g__Dialister", "g__Proteus", "g__Acinetobacter", "g__Bifidobacterium", "g__Lactobacillus", "g__Haemophilus"),
  random_effects = c("InfantID"))

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



