## Virome analysis workflow
1. Quality-filter reads. See 01_QC.sh
2. Build contigs, map reads from samples to contigs. See 02A_Contig-Building-Mapping.sh and 02B_taxonomizr.R
3. Run decontam. Remove contaminants and false positives, normalize contig counts. Perform alpha and beta diversity analyses (manuscript Figures 1, 2 and 3). See 03_Decontam-AlphaDiv-BetaDiv.R
4. LEfSe analysis (manuscript Figure 3). See 04_LEfSe.sh
5. Predict phage lifestyles (supplementary figure 3). See 05_Prodigal-PHACTS.sh
6. Transkingdom analysis (manuscript figure 4). See 06_Maaslin2.sh
