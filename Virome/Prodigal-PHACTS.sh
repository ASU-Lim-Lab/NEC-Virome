#Use prodigal to predict open reading frames and get amino acid sequences for each contig.
#Input: individual nucleotide fasta file for each contig. Output: individual faa file for each contig.
prodigal -i ContigName.fasta -a ContigName.faa -p meta;

#Use PHACTS to predict lifestyles (temperate/lytic) for contigs with at least 5 predicted open reading frames.
perl phacts.pl --file ContigName.faa --classes classes_lifestyle;

#The script phacts.pl was used as published (https://github.com/deprekate/PHACTS), with lines 309-316 modified to the following:
}else{
	open(FH, '>', "$file_path-predictions.txt");
	foreach (sort {$ave_predictions{$b} <=> $ave_predictions{$a} } keys %ave_predictions){
		print FH "$file_path","\t",$_,"\t",$ave_predictions{$_},"\t",$std_predictions{$_},"\n";
	}
	close(FH);
}