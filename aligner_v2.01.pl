
###################################################################################################################################
#
#
#	aligner.pl
#
#    	Copyright (C) 2019-2023 Douglas Chesters
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#	contact address dc0357548934@live.co.uk
#
#
#
#
###################################################################################################################################
# 
# 
# 
# 
# 
# Citation:
# 	Chesters, Douglas (2019) aligner.pl, high-throughput DNA barcode sequence alignment. Computer software available from github.com/dchesters/aligner.
# 
# 
# 
# 
# 
# DEPENDENCIES:
#	NCBI Blast+ (command line tools).
# 	Emboss (command 'needle' is used). Tested on EMBOSS version 6.6.0.0
# 
# 
# 
# 
# 
# 
# QUICK START
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
# 
# 
# 
# 
# 
# MORE DETAILS:
# 	Basic familiarity with use of Linux terminal is required to use this software.
# 	This program requires a curated reference alignment, and a set of query unaligned sequences, each in fasta format.
#	 First, a Blast search is done of queries against the references, 
# 	 this is to find the most similar reference, for each sequence.
# 	Then for each pair (query and its similar reference), pairwise Needleman-Wunsch alignment is conducted.
# 	Default pair alignment settings that seem to work well were determined during 
#	a study using Chinese bee COI barcodes: gapopen=20.0, gapextend=0.5, endweight=Y.
# 	These can be changed from within the script (below) if neccesary.
# 	Next, the query sequence is processed as neccessary, to maintain the integrity of the reference alignment.
# 	Aligned queries are then written to output file.
# 	I have tried several other software which do this same job, 
#	 however, each of them either crashed or output major alignment errors, under certain conditions.
# 	Program was developed for Linux platform (so will probably work also with UNIX or OSX).
# 
# 
# 
# 
# ######################################################################################################
# 
# 
# CHANGE LOG:
# 
# 2016-12-12: 	Development started while working on large insect DNA barcode data, 
# 		due to various issues with existing barcode alignment tools.
# 		basically a wrapper for Blast to references followed by pairwise NW alignment
# 2017-02-17: 	Error messages given when references provided are not aligned, or pairwise alignment fails
# 2017-04-18: 	Changed default gap paremeters, previous ones gave unnecessary gaps in bee barcode study.
# 2017-04-28: 	Bug fix.
# 2019-04-14: 	Wrote in GNU GPL, uploaded to source forge. 
#		has been used in 3 studies currently being submitted:		
# 		Chesters. The phylogeny of insects in the data-driven era.
# 		Chesters et al. Climatic and Vegetational Drivers of Insect Beta Diversity at the Continental Scale.
# 		Wang et al. Multiple components of plant diversity loss determine herbivore phylogenetic diversity in a subtropical forest experiment.
# 2019-05-30: 	Checks that queries in the blast file are included in the sequence file, if not, prints error message and quits. 
# 2019-11-01: 	Now fully invoked by command line.
# 2021-03-09: 	95% into aligning a million barcodes and crashed due to seperate job taking all memory,
# 		so if needle doesnt produce output, script now continues
# 		Noticed many sequences end up very short after alignment, so added option to remove short sequences
# 2021-03-11: 	Short alignments apparently mostly just short seqs in the first place. count these and report.
# 2022-02-23: 	Noticed crashes if run on profile alignment with dos lines breaks. error message now given to this effect.
# 		Prints log file file pairwise alignment details. (currently confused why taxon specific reference gives exactly same as generalized one)
# 		Some temp files deleted after completion.
# 2022-02-24:	Permits self hits, since subject may be of same species as something in the profile alignment.
# 2023-01-21:	Added a few comments to code and printed some more details.
# 		Finally figured out how to run the PW aligner without its 1 line screenoutput (which was printed to screen for each sequence).
# 		Prints log file (low_quality_alignments.txt) with details on queries that didnt hit any references, 
#		and those which produced a poor quality pair alignment.
# 2023-06-09	prints some useful info if user inputs profile alignment with gaps.
# 2023-07-09	Minor change to log file.
# 2023-07-12	Split running this script into 2 so that a more thorough check of sequence quality can be done.
# 		I noticed on a large Lep COI dataset that a lot of the failed alignemnts are due to blocks of ambiguous sites.
# 		Thus procedure is changed, first run a check of references to remove seqs with gaps or ambiguous blocks,
# 		and check queries for ambiguous blocks, short sequences, and presence of gaps.
# 2023-07-13	Improved reference selection for where Blast is run with max_target_seqs > 1.
# 2023-10-06	Updated instructions in script (now run 2-part). Uploaded to github (v2.01).
# 
# 
# 
# 
# 
# 
# 
# 
#
# 
# 
# 
# 
# 
# 
######################################################################################################


$arg_string  = join ' ', @ARGV;
$processed_profile_filename = "processed_profile_alignment";
$processed_queries_filename = "processed_queries";

#####################################
read_command_arguments($arg_string);#
#####################################



# $blast_file 		= $ARGV[0];
# $refs			= $ARGV[1];
# $qs			= $ARGV[2];
# $outfile_name 		= $ARGV[3];



	# settings work well for COI of Chinese bees: $gapopen 20.0, $gapextend 0.5, $endweight = Y
	# other Needle PW alignment options see http://emboss.sourceforge.net/apps/cvs/emboss/apps/needle.html


						# default 10.0
# $gapopen 			= "20.0";# (Floating point number from 1.0 to 100.0)
						# CN bees COI, 10.0 created unneccessary gaps


						# default 0.5
# $gapextend 			= "0.5"; # (Floating point number from 0.0 to 10.0)


						# 
# $endweight 			= "Y"; # boolean    [N] Apply end gap penalties.



unless($blast_file =~ /\w/ && $refs =~ /\w/ && $qs =~ /\w/ && $outfile_name =~ /\w/)
	{die "\ncommand error. could not read 4 file names in your command. quitting.\n"};



##############################################################################################################################


# first count hits for each query

my %count_hits_for_each_query = ();

open(IN, $blast_file) || die "\nerror 162. cant open file Blast output ($blast_file)\n";
print "\n opened file:$blast_file\n";
while (my $line = <IN>)
	{
	#Daphnia_longispina_JN903664	Daphnia_longispina_JX134334	0.0	97.35	528
	if($line =~ /^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
		{
		my $query = $1;my $hit = $2; my $pident= $4; my $alignment_length = $5;
		# print "query:$query hit:$hit pident:$pident alignment_length:$alignment_length\n";
		if($pident >= 50)
			{
			#		i dont remember the rationale for excluding self hits.
			#		but for current purpose i need this so must deactivate
			#unless($query eq $hit)

				if($alignment_length =~ /\d/)
					{
					if($hit_of_longest_alignment{$query} =~ /\d/)
						{
						$autodetect_max_target_seqs = 1;
						if($alignment_length > $hit_of_longest_alignment{$query})
							{
							$count_hits_for_each_query{$query}++;
							$align_query_to_these_hits{$query} = $hit;
							$hit_of_longest_alignment{$query} = $alignment_length;
							$better_alignment_found++;
							}else{
							$prev_alignment_preferred++;
							};
						$multiple_hits_for_query++;
						}else{
						$count_hits_for_each_query{$query}++;
						$align_query_to_these_hits{$query} = $hit;
						$hit_of_longest_alignment{$query} = $alignment_length;
						$first_encounter_query++;
						};

					}else{
					die "\nerror 200. no alignment length parsed from Blast output:$line\n";
					};

			}else{
			print "\nerror why is similarity so low:$line\n";
			}
		}
	#print $line;

	}

close IN;

my @all_queries = keys %count_hits_for_each_query;#@all_queries = sort {$b cmp $a} @all_queries;
@all_queries = sort @all_queries;
print "\nread file Blast result file ($blast_file), counted $#all_queries queries.
\n";
if($autodetect_max_target_seqs == 1)
	{
	print "multiple hits found per query (Blast run with max_target_seqs > 1)
	first_encounter_query:$first_encounter_query
	multiple_hits_for_query:$multiple_hits_for_query
	better_alignment_found:$better_alignment_found
	prev_alignment_preferred:$prev_alignment_preferred

";
	};



################################################################################################

my %seqs = ();
$convert=0;

##################################################
$reference_alignemnt_length = read_refs();#
##################################################


################################################################################################

################################################################################################

# read subject file specified by user and put IDs into %q_seqs
# if($delete_short_seqs == 1) && too short, then ID stored in %queries_shorter_than_cutoff

#######################
read_qs();#
#######################

my @user_queries = keys %q_seqs;
my $count_user_queries = scalar @user_queries;
print "count_user_queries:$count_user_queries\n";

#######################################################################################

################################################################################################


open(LOG2 ,">low_quality_alignments.txt");


foreach my $key(@all_queries) # these queries stored while reading blast results.
	{
	if(exists($q_seqs{$key})) # whereas these stored from the fasta file.
		{
	#	print "\tfound in fasta:$key\n"
		}else{
	#	print "$key\n";
		$error_instances++;
		print "\terror, query from blast file is not in fasta file:$key\n"
		};
	};
#####################
# also do inverse, find which entries in the query find had no blast results

foreach my $key(@user_queries) # all in query fasta file
	{
	if(exists($count_hits_for_each_query{$key})) # queries in blast results file
		{
		if($queries_shorter_than_cutoff{$key} ==1) # short to begin with, ignore.
			{
			}else{
			# query hit something in the references,
			};
		}else{
		$user_queries_no_blast_hit++;
		if($queries_shorter_than_cutoff{$key} ==1) # short to begin with, ignore.
			{
			}else{
			# length not an issue, still no blast hit
			print LOG2 "$key\n";
			};

		};
	};
print "user_queries_no_blast_hit:$user_queries_no_blast_hit\n";

#######################################################################################

if($error_instances >= 1)
	{
	die "\nthere are sequences missing from your fasta file, probably wrong sequence file specified. quitting.\n\n";
	};




################################################################################################


################################################################################################








open(OUTPUT , ">$outfile_name") || die "\nerror 165\n";
open(LOG ,">aligner_LOG.txt");


foreach my $key(@all_queries) # these queries stored while reading blast results.
	{
	$acounter++;
	my $hits_for_this_query = $align_query_to_these_hits{$key};
	$hits_for_this_query =~ s/\n$//;
#	my @hits = split /\n/ , $hits_for_this_query;

	if($acounter =~ /00$/)
		{
		print "$acounter of $#all_queries key:$key , aligning to $hits_for_this_query\n";
		};

		$q_seq = $q_seqs{$key};$ref_seq = $ref_seqs{$hits_for_this_query};
		unless($ref_seq =~ /[actg]/i && $q_seq =~ /[actg]/i)
			{

			print "ref_seq $hits_for_this_query:($ref_seq)\nq_seq $key:($q_seq)\n";
			print "\n$acounter of $#all_queries key:$key , aligning to $hits_for_this_query\n";

			unless($ref_seq =~ /[actg]/i)
				{
				print "error, no reference seqeunce while processing $key\n " , " this error can be caused by dos line breaks in profile alignment.\n"
				};
			unless($q_seq =~ /[actg]/i){print "error, no query seqeunce while processing $key\n"};
			die "\nerror 163\n";

			};
		
		system("rm alignerQUERY alignerREFERENCE needleOUT");

		open(OUT, ">alignerQUERY") || die "\nerror 27.\n";
		print OUT ">$key\n$q_seq\n";
		close(OUT);

		open(OUT2, ">alignerREFERENCE") || die "\nerror 27.\n";
		print OUT2 ">$hits_for_this_query\n$ref_seq\n";
		close(OUT2);

		# emboss option -auto runs the software without usual printout of the function being run 
		my $command = "needle -asequence alignerQUERY -bsequence alignerREFERENCE -outfile needleOUT -gapopen $gapopen -gapextend $gapextend -endweight $endweight -awidth3 10000 -aformat3 fasta -auto";
		system("$command");

		# not needed now anyways
		# this still gives the needle screen output:
		# $unwanted_screenoutput = `$command`;
		# so does this:
		# my $unwanted_screenoutput = qx{$command};


		#######################
		parse_needle_output( $hits_for_this_query );#
		#######################

#	if($key =~ /Acyrthosiphon_spKU374293/){die ""};



	}


close OUTPUT;
close LOG;
close LOG2;

system("rm alignerQUERY alignerREFERENCE needleOUT");

# ref_seqs

# 
my $count_aligned_queries  = scalar @all_queries;
print "
count_user_queries:$count_user_queries
	user_queries_no_blast_hit:$user_queries_no_blast_hit
count_aligned_queries:$count_aligned_queries
PW_alignment_failed:$PW_alignment_failed
";

if($delete_short_seqs == 1)
	{
	print "user specified to remove aligned subjects if < $delete_short.\n";
	print "input file has $subjects_shorter_than_limit_denovo sequence that are shorter than this in the first place\n";
	print "short_seqs_deleted after alignment:$short_seqs_deleted\n";
	}

print "\n\nFIN.\n";
exit;



############################################################################################################


# SUBS follow


############################################################################################################


sub parse_needle_output
{

my $currentrefID = shift;

open(NEEDLE_OUT , "needleOUT") || print "\nerror 107:no needle output\n";
my $filestring = "";
my $line_no = 0;
while (my $input = <NEEDLE_OUT>)
	{
	$line_no++;
#	print "line_no:$line_no input:$input\n";
	$filestring .= $input;
	}
close(NEEDLE_OUT);

if($line_no >= 1)
	{

my @split = split />/ , $filestring;

my $query_ID = "NA";
if($split[1] =~ s/^(.+)//){$query_ID = $1};
$split[2] =~ s/^.+//;
$split[1] =~ s/[\n\r\s]+//g;
$split[2] =~ s/[\n\r\s]+//g;

my $seqlength = length($split[1]);
my $trimmed_query_seq = "";

for my $i(0 .. ($seqlength-1))
	{
	my $base_query = substr $split[1] , $i , 1;my $base_reference = substr $split[2] , $i , 1;
	unless($base_reference =~ "-"){$trimmed_query_seq .= $base_query};
	};

# print "0:$split[0]\n";
# print "1:($split[1])\n"; # query 
# print "2:($split[2])\n"; # reference
# print "trimmed_query_seq:$trimmed_query_seq\n";

my $length_remaining; my $short_deleted=0;my $short_deleted_print="Length_OK";
if($delete_short_seqs == 1)
	{
	my $trimmed_query_seq_test = $trimmed_query_seq;$trimmed_query_seq_test =~ s/[nN\-\?]//g;$length_remaining = length($trimmed_query_seq_test);
	# print "$length_remaining,$delete_short\n";
	if($length_remaining >= $delete_short)
		{
		print OUTPUT ">$query_ID\n$trimmed_query_seq\n";
		}elsif($length_remaining  <= 150)
		{
		# extremly short
		$short_deleted=1;$short_deleted_print = "Very_short";

		if($queries_shorter_than_cutoff{$query_ID} == 1) # short to begin with
			{
			print "subject $query_ID extremly short (<150bp), even before attempting alighment\n";
			}else{
			print "warning, subject $query_ID has extremly short ($length_remaining) sequence after alignment to top hit refereence ($currentrefID).\n";
			print LOG2 "Short after alignment\t$query_ID\n";
			};
		# die "";
		$short_seqs_deleted++;
		}else{

		if($queries_shorter_than_cutoff{$query_ID} == 1) # short to begin with
			{
			# print "subject $query_ID extremly short (<150bp), even before attempting alighment\n";
			}else{
			# print "warning, subject $query_ID has extremly short ($length_remaining) sequence after alignment to top hit refereence ($currentrefID).\n";
			print LOG2 "after alignment 150-$delete_short\t$query_ID\n";
			};


		$short_seqs_deleted++;
		$short_deleted=1;$short_deleted_print = "Too_short";
		};
	}else{
	print OUTPUT ">$query_ID\n$trimmed_query_seq\n";
	};



my $new_aligned_length = length($trimmed_query_seq);

 print LOG "$query_ID\t$currentrefID\t$length_remaining\t$short_deleted_print\n";

unless($new_aligned_length == $reference_alignemnt_length)
	{
	print "\nalignment error .... 
	reference alignment is length $reference_alignemnt_length, however length of $query_ID is $new_aligned_length
";
die "";
	};


	}else{
	$PW_alignment_failed++;
	};


}; # sub parse_needle_output



#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################





sub read_refs
{


my $file_as_string = "";
open(IN_FILTER , "$refs") || die "\nerror. cant open file ($refs)\n";
while (my $line = <IN_FILTER>)	{$file_as_string .= $line};
close(IN_FILTER);
my @all_lines = split />/, $file_as_string;
print "read reference alignment, " , scalar @all_lines , " seqs in file\n";
my $length_first_entry;
for my $each_line(1 .. $#all_lines)
	{
	my $line = $all_lines[$each_line];
	if($line =~ /^(.+)/)
		{
		my $speciesid = $1;	#print "Species ID:$speciesid\n";
		$line =~ s/^.+\n//;
		$line =~ s/\012\015?|\015\012?//g;$line =~ s/\n//g;$line =~ s/\r//g;
		$line =~ s/[\s\t]//g;
# print "Species ID 2:$speciesid\n";
		$ref_seqs{$speciesid}  = $line;  # print "stored_refseq:$speciesid, entry:$line\n";
		if($line =~ /\-/)
			{
			print "\n\nerror. current implementation requires reference alignment with no gaps.\n" , 
				"perhaps run following on your profile alignment: remove_fasta_entries_with_gaps.pl\n\n";
			die "";
			};


		if($each_line == 1)
			{
			$length_first_entry = length($line);
			}else{
			unless($length_first_entry == length($line))
				{
				print "\nerror. reference alignment you provided seems not aligned\n";
				print "length of entry $speciesid is " , length($line);
				print ", which differs from length ($length_first_entry) of first entry in file\n";
				die  "";
				};
			};
	#	print " speciesid:$speciesid";

		}
	}

print "
all in reference alignment file are length $length_first_entry
";
return($length_first_entry);




}; # sub read_refs



#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################




sub read_command_arguments
{

my $arguments = shift;
$arguments =~ s/\n//;$arguments =~ s/\r//;


# $blast_file 		= $ARGV[0];
# $refs			= $ARGV[1];
# $qs			= $ARGV[2];
# $outfile_name 		= $ARGV[3];


	# settings work well for COI of Chinese bees: $gapopen 20.0, $gapextend 0.5, $endweight = Y
	# other Needle PW alignment options see http://emboss.sourceforge.net/apps/cvs/emboss/apps/needle.html

						# default 10.0
# $gapopen 			= "20.0";# (Floating point number from 1.0 to 100.0)
						# CN bees COI, 10.0 created unneccessary gaps

						# default 0.5
# $gapextend 			= "0.5"; # (Floating point number from 0.0 to 10.0)

						# 
# $endweight 			= "Y"; # boolean    [N] Apply end gap penalties.








if($arguments =~ /-references\s+(\S+)/)
	{
	$refs = $1;
	}else{
	print "\nerror reading command arguments\n";die""
	}
if($arguments =~ /-subject\s+(\S+)/)
	{
	$qs = $1;
	}else{
	print "\nerror reading command arguments\n";die""
	}
if($arguments =~ /-delete_short\s+(\d+)/)
	{
	$delete_short = $1;$delete_short_seqs = 1;
	}else{
	};
# if just processing input files, dont need to check settings below
if($arguments =~ /\-process_input/)
	{
	process_inputs();
	}else{
	};



if($arguments =~ /-blast_results\s+(\S+)/)
	{
	$blast_file = $1;
	}else{
	print "\nerror reading command arguments\n";die""
	}

if($arguments =~ /-outfile\s+(\S+)/)
	{
	$outfile_name = $1;
	}else{
	print "\nerror reading command arguments\n";die""
	}
if($arguments =~ /-gapopen\s+(\S+)/)
	{
	$gapopen = $1;
	}else{
	print "\nerror reading command arguments\n";die""
	}
if($arguments =~ /-gapextend\s+(\S+)/)
	{
	$gapextend = $1;
	}else{
	print "\nerror reading command arguments\n";die""
	}
if($arguments =~ /-endweight\s+(\S+)/)
	{
	$endweight = $1;
	}else{
	print "\nerror reading command arguments\n";die""
	}



};


#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################



# Acyrthosiphon_spKU374293 is too short

# >Acsala_anomala_KJ374863
# ------TGAGCTGGAATAATAGGAACCTCATTAAGATTATTAATTCGAGCTGAATTAGGTAACCCTGGATCATTAATTGGAGATGATCAAATTTATAATACTATTGTAACTGCTCATGCTTTTATTATAATTTTTTTTATAGTTATACCTATTATAATTGGTGGATTTGGTAATTGATTAGTACCTTTAATATTAGGAGCTCCCGATATAGCCTTTCCTCGAATAAATAATATAAGTTTTTGACTTCTTCCCCCCTCATTAACTTTATTAATTTCGAGAAGAATTGTAGAAAATGGAGCAGGAACAGGATGAACAGTGTACCCCCCACTTTCATCCAATATTGCCCACAGAGGTAGTTCAGTAGATTTAGCTATTTTTTCCCTACATTTAGCAGGTATTTCATCTATTCTCGGAGCTATTAATTTTATTACAACAATTATTAATATACGATTAAATAGATTATCCTTTGATCAAATACCTTTATTTGTTTGAGCTGTAGGTATTACAGCTTTCTTATTCTT-CTTTCATTACCTGTCCTAGCAGGAGCTATTACTATACTATTAACTGATCGAAATTTAAATACATCCTTTTTTGACCCAGCAGGAGGTGGAGATCCTATTCT-
# >Acyrthosiphon_spKU374293
# GATCATCTCTTAGAATTTTAATTCGTTTAGAATTAAGAAAT------TCTATTATTAATAATAATCAATTATATAATGTAATTGTTACAATTCATGCTTTTATTATAATTTTTTTTATAACAATACCAATTGTAATTGGTGGTTTTGGAAACTGATTAATTCCTATAATAATAGGATGTCCTGATATATCATTTCCACGTTTAAATAATATTAGATTTTGATTATTACCACCATCATTAATAATAATAATTTGTAGACTAATT----AAAATGGAACAGGAACAGGATGAACTATTTATCCACCTTTATCAAATAATATTGCACA-----CAATATCAGTTGATTTAACTATCTTCTCTTTACATTTAGCAGGAATTTCATCAATTTTAGGAGCAATTAACTTTATTTGTACAATTCTTAATATA--ATCAAAA--TATAAAATTAAATCAAATTCCCCTTTTCCCTTGATCAATTTTAATTACAGCATT-TTATTAATTTTATCTTTACCAGTTTTAGCTGGTGCTATTACAATATTATTAACTGATCGTAATTTAAATACATCATTTTTTGATCCAGCAGGAGGAGGAGATCCTATTTTA


#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################


sub read_qs
{
	my $file_as_string = "";


			# this file given by user for switch -subject
	open(IN_FILTER , "$qs") || die "\nerror. cant open file($qs)\n"; 

	while (my $line = <IN_FILTER>)	{$file_as_string .= $line};close(IN_FILTER);
	unless($file_as_string =~ />/){die "\nerror 138, expecting fasta format. quitting.\n"};
	my @all_lines = split />/, $file_as_string;
	print "read queries from file $qs, there are " , scalar @all_lines , " seqs.\n";

	for my $each_line(1 .. $#all_lines)
		{
		my $line = $all_lines[$each_line];
		if($line =~ /^(.+)/)
			{
			my $speciesid = $1;	$speciesid =~ s/\n//;$speciesid =~ s/\r//;


		#	if($speciesid =~ /BJ2012_1495_711/)
		#		{
		#		print "test2\n";
		#		print "$speciesid\n";
		#		};


			$line =~ s/^.+\n//;
			$line =~ s/\012\015?|\015\012?//g;$line =~ s/\n//g;$line =~ s/\r//g;
			$line =~ s/[\s\t]//g;


			if($delete_short_seqs == 1)
				{
				my $trimmed_query_seq_test = $line;$trimmed_query_seq_test =~ s/[nN\-\?]//g;
				my $length_remaining = length($trimmed_query_seq_test);
				# print "$length_remaining,$delete_short\n";
				if($length_remaining <= $delete_short)
					{
					$subjects_shorter_than_limit_denovo++;
					$queries_shorter_than_cutoff{$speciesid} = 1;
					}
				}



# BJ2012_1495_711
#			if($speciesid =~ /BJ2012_1495_711/)
#				{
#				print "storeing:$speciesid $line\n";
#				};


			$q_seqs{$speciesid}  = $line;

			}
		}

if($delete_short_seqs == 1)
	{
	print "user specified to remove aligned subjects if < $delete_short.\n";
	print "input file has $subjects_shorter_than_limit_denovo sequence that are shorter than this in the first place\n";
	};


} # sub read_qs



#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################


sub process_inputs
{
print "

processing inputs .... 
	profile:$refs
	queries:$qs
	length cutoff:$delete_short

";

$database_file = $refs;



###################################################################

open(FASTA_IN, $database_file) || die "Cant open $database_file.\n";
open(OUT11, ">$processed_profile_filename") || die"";
my $fasta_entry = "";
while(my $fasta_line= <FASTA_IN>)
	{
	if($fasta_line =~ /^>.+/)
		{
		unless(length($fasta_entry)<=2)
			{
			$fasta_entry =~ s/^>//;
			########################################
			process_entry($fasta_entry);#
			########################################
			}
		$fasta_entry = $fasta_line;
		}else{
		$fasta_entry .= $fasta_line;
		}
	};
close(FASTA_IN);

# do the last entry!:
unless(length($fasta_entry)<=2)
	{
	$fasta_entry =~ s/^>//;
	########################################
	process_entry($fasta_entry);#
	########################################
	}



	sub process_entry
		{
		my $line = shift;
		my $current_id = "";
		if ($line =~ /^(.+)\n/ )
			{
			$current_id = $1;
			$line =~ s/^.+\n//;#$line =~ s/\n//;$line =~ s/\r//;my $current_seq_length = length($line);
			}else{
			die "\nerror 5631.\n"
			}
		$references_checked++;

		my $print_reference = 1;
		if($line =~ /[ACTG]N{10,}[ACTG]/i)
			{
			$REF_entries_with_ambiguous_blocks++;print "Ambigous block in REF $current_id\n";
			$print_reference = 0;
			};

		if($line =~ /\-/)
			{
			$entry_with_gaps_removed++;
			$print_reference = 0;
			};

		if($print_reference == 1)
			{
			print OUT11 ">$current_id\n$line";$REF_printed++;
			}else{
			$REF_not_printed++;
			};
		}

print "
references_checked:$references_checked
	entry_with_gaps_removed:$entry_with_gaps_removed
	REF_entries_with_ambiguous_blocks:$REF_entries_with_ambiguous_blocks
	REF_printed:$REF_printed
	REF_not_printed:$REF_not_printed
";
close OUT11;
print "new processed profile alignment printed to file name $processed_profile_filename\nnow checking queries $qs .....\n";
###################################################################

$database_file = $qs;

open(FASTA_IN, $database_file) || die "Cant open $database_file.\n";
open(OUT12, ">$processed_queries_filename") || die"";
my $fasta_entry = "";
while(my $fasta_line= <FASTA_IN>)
	{
	if($fasta_line =~ /^>.+/)
		{
		unless(length($fasta_entry)<=2)
			{
			$fasta_entry =~ s/^>//;
			########################################
			process_Qentry($fasta_entry);#
			########################################
			}
		$fasta_entry = $fasta_line;
		}else{
		$fasta_entry .= $fasta_line;
		}
	};
close(FASTA_IN);
unless(length($fasta_entry)<=2)
	{
	$fasta_entry =~ s/^>//;
	########################################
	process_Qentry($fasta_entry);#
	########################################
	}



	sub process_Qentry
		{
		my $line = shift;
		my $current_id = "";
		if ($line =~ /^(.+)\n/ )
			{
			$current_id = $1;
			$line =~ s/^.+\n//;#$line =~ s/\n//;$line =~ s/\r//;my $current_seq_length = length($line);
			}else{
			die "\nerror 5631.\n"
			}

		if($line =~ /\-/){die "\nerror 935. Query sequences contains gap chararcters, please remove prior to alignment.\n"};
		$queries_checked++;
		my $print_query = 1;
		if($line =~ /[ACTG]N{10,}[ACTG]/i)
			{
			# print "ambiguous block:$line\n";
			$QUERY_entries_with_ambiguous_blocks++;#print "Ambigous block in $current_id\n";
			$print_query = 0;
			};

	#	if(length($line) < $delete_short)
	#		{
	#		$Query_too_short++; # print "Short query $current_id\n";
	#		$print_query = 0;
	#		};
			if($delete_short_seqs == 1)
				{
				my $trimmed_query_seq_test = $line;$trimmed_query_seq_test =~ s/[nN\-\?]//g;
				my $length_remaining = length($trimmed_query_seq_test);	# print "$length_remaining,$delete_short\n";
				if($length_remaining <= $delete_short)
					{
					$Query_too_short++;$print_query = 0;
					};
				}




		if($print_query == 1)
			{
			print OUT12 ">$current_id\n$line";$Q_printed++;
			}else{
			$Q_not_printed++;
			};
		}


close OUT12;

print "
queries_checked:$queries_checked
	Query_entries_with_ambiguous_blocks:$QUERY_entries_with_ambiguous_blocks
	Query_too_short:$Query_too_short
	Query_printed:$Q_printed
	Query_not_printed:$Q_not_printed

Wrote 2 files ($processed_queries_filename, $processed_profile_filename) that will be used in next step, Blast.\n\nFIN.\n
";
###################################################################

exit;
};

#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################

