## Virome analysis workflow
1. QC \[add code\]
2. Query quality-filtered reads against viral database using BLASTx \[add code\]
3. Parse BLASTx output files in MEGAN. LCA parameters: min score: 100; min support percent: 0; min support: 1.
3. Run decontam and remove contaminants and false positives. Normalize species counts.
4. Alpha and beta diversity analyses (manuscript Figures 1 & 2). See Decontam_AlphaDiv_BetaDiv.R
5. Build contigs, map reads from samples to contigs, normalize contig counts. See Contig_Building_Mapping.sh
6. Beta diversity analysis (manuscript Figure 3).
7. LEfSe analysis (manuscript Figure 3), using normalized contig counts and browser-based LEfSe module available at: http://huttenhower.sph.harvard.edu/galaxy
8. Prevalence and abundance analysis (manuscript Figure 3 and Supplementary Figure 3). See Prevalence_Abundance.R