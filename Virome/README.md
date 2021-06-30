## Virome analysis workflow
1. Quality-filter reads. See QC.sh
2. Build contigs, map reads from samples to contigs, normalize contig counts. See Contig-Building-Mapping.sh
3. Run decontam and remove contaminants and false positives. Normalize contig counts. See Decontam-AlphaDiv-BetaDiv.R
4. Alpha and beta diversity analyses (manuscript Figures 1, 2 and 3). See Decontam-AlphaDiv-BetaDiv.R
5. LEfSe analysis (manuscript Figure 3). See LEfSe.sh
6. Predict phage lifestyles (supplementary figure 3). See Prodigal-PHACTS.sh
7. Transkingdom analysis (manuscript figure 4).