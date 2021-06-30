#Use LEfSe to identify discriminant contigs in case and control samples +/- 10 days before NEC onset.
#Input: table with contig relative abundance per sample. Filtered to contigs present in at least 10% of samples.
format_input.py Contig_RelativeAbundance_cases.txt case_lefse_input.in -c 1 -u 2 -o 1000000
run_lefse.py case_lefse_input.in case_lefse_output.res