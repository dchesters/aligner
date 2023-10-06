# aligner
#
# Citation:
# 	Chesters, Douglas (2019) aligner.pl, high-throughput DNA barcode sequence alignment. Computer software available from github.com/dchesters/aligner.
#  
# DEPENDENCIES:
#	NCBI Blast+ (command line tools).
# 	Emboss (command 'needle' is used). Tested on EMBOSS version 6.6.0.0
# 
#  QUICK START
# 	Aligner is run in 2 stages, the first you provide the profile alignment and the queries, and the script will process these:
# 	perl aligner.pl -process_input -references $reference_alignment -subject $unaligned -delete_short 475
# 
# 	The output is 2 files named processed_profile_alignment, processed_queries.
# 	Next you need to do a Blast+ search, this matches each query to its most similar profile:
# 	blastn-2.12.0 -task blastn -query processed_queries -subject processed_profile_alignment -out BlastOUT -word_size 10 -perc_identity 60 -max_target_seqs 5 -evalue 1e-8 -dust no -strand plus -outfmt '6 qseqid sseqid evalue pident length'
# 
# 	The single output file is tabular format Blast result. 
# 	Finally, using the Blast result and the file of queries, pairwise NW alignment is conducted:
# 	perl ~/usr_scripts/aligner.pl -blast_results BlastOUT -references processed_profile_alignment -subject processed_queries -outfile alignerOUT -delete_short 475 -gapopen 20.0 -gapextend 0.5 -endweight Y
# 
