=head
This script is provided in support of Bush, et al. Genomic diversity affects the accuracy of bacterial SNP calling pipelines. https://www.biorxiv.org/content/10.1101/653774v1
It is provided for reproducibility purposes only. No attempt has been made to optimise the code.

Prerequisites:
1. Download the data archive from: https://ora.ox.ac.uk/objects/uuid:8f902497-955e-4b84-9b85-693ee0e4433e. Its contents are referred to here as $in_dir, $fq_dir, $rep_dir, and $locs_dir.
2. This archive also contains information regarding third-party tools, in the file "Installation_instructions.txt". This script will test for the presence of some of these tools.

Fuller details are available in the readme provided with the archive.
=cut

use strict;
use warnings;
use Acme::Tools qw(avg);

# REQUIREMENTS
my $root     = '/home/ndm.local/steveb';
my $path     = "$root/varcall_evaluation";
my $progs    = "$root/programs"; # all third-party tools in "Installation_instructions.txt" are installed here
my $in_dir   = "$path/filtered_vcfs"; # created by 3.filter_VCFs.pl
my $fq_dir   = "$path/illumina_reads_for_nanopore_illumina_hybrid_assemblies";
my $rep_dir  = "$path/representative_genomes_for_nanopore_illumina_hybrid_assemblies";
my $locs_dir = "$path/nanopore_illumina_hybrid_assemblies";

# CHECK THAT ALL DEPENDENCIES ARE AVAILABLE
my $fatal = 0; # ALL THE FOLLOWING PATHS ARE HARD-CODED PREREQUISITES. WE WILL ABORT IF THEY CANNOT BE FOUND.
if (!(-d($path)))     { $fatal++; print "ERROR: cannot find $path\n";     }
if (!(-d($progs)))    { $fatal++; print "ERROR: cannot find $progs\n";    }
if (!(-d($in_dir)))   { $fatal++; print "ERROR: cannot find $in_dir\n";   }
if (!(-d($fq_dir)))   { $fatal++; print "ERROR: cannot find $fq_dir\n";   }
if (!(-d($rep_dir)))  { $fatal++; print "ERROR: cannot find $rep_dir\n";  }
if (!(-d($locs_dir))) { $fatal++; print "ERROR: cannot find $locs_dir\n"; }
exit 1 if ($fatal > 0);

# PARAMETERS
my @aligners = (qw/bbmap bowtie2 bwa-mem bwa-sw bwa+stampy cushaw3 gassst gem hisat2 minimap2 mosaik ngm novoalign smalt snap snippy spandx speedseq stampy yara/);
my @callers	 = (qw/16GT deepvariant freebayes gatk lofreq mpileup octopus pilon platypus snippy snver snvsniffer solsnp spandx speedseq strelka varscan/);
my %aligners = map {$_ => 1} @aligners;
my %callers  = map {$_ => 1} @callers;
my $exclude_repetitive_regions = 1; # repetitive regions are identified by self-self BLASTn of each representative genome. The BLASTn reports are available in $rep_dir. The criteria by which we consider a region as repetitive are hard-coded on line 70 (≥ 95% identity over length ≥ 100bp, with no more than 1 gap, and an E-value < 0.05, not including the match of the entire genome against itself).

# OUTPUT
my $out_dir = '';
if ($exclude_repetitive_regions == 0)
	{ $out_dir = "$path/evaluationOutput"; }
elsif ($exclude_repetitive_regions == 1)
	{ $out_dir = "$path/evaluationOutput-repeatmasked"; }
if (!(-d($out_dir))) { mkdir $out_dir or die $!; }
my $out_file1 = "$out_dir/summary_statistics.tsv";
my $out_file2 = "$out_dir/pipeline_rankings.tsv";
open(OUT1,'>',$out_file1) or die $!; open(OUT2,'>',$out_file2) or die $!;
print OUT1 "Genome\tNumber of cleaned (that is, quality- and adapter-trimmed) reads (millions)\tAligner\tCaller\tAligner/caller\tTotal number of true SNPs\tNo. of true positives\tNo. of true positives as a % of no. of true SNPs\tNo. of false positives\tNo. of false positives (that are also heterozygous calls)\tNo. of false positives as a % of no. of true SNPs\tNo. of false negatives\tNo. of false negatives as a % of no. of true SNPs\tTotal number of errors (=FP+FN) per million sequenced bases\tTotal number of FPs per million sequenced bases\tTotal number of FNs per million sequenced bases\tPrecision (positive predicitive value)\tRecall (sensitivity)\tF-score (=2 . precision. recall/(precision + recall))\tMiss rate\n";
print OUT2 "Rank\tAligner/caller\tSum of individual ranks (F-score, precision, recall, total number of true positives, false positives and false negatives, total number of errors per million sequenced bases)\tRanks and scores for each metric\n";

# RECAP OF WHAT WE'VE ALREADY DONE SO FAR: WE'VE ALIGNED RAW READS (FASTQs IN $fq_dir) TO A REPRESENTATIVE GENOME OF THE SPECIES SEQUENCED (FASTAs IN $rep_dir), AND THEN CALLED SNPs (VCFs IN $in_dir).
# WE'VE ALSO ALREADY CREATED CLOSED ASSEMBLIES FROM THESE RAW READS, FROM WHICH WE IDENTIFIED THE LOCATION OF SNPs RELATIVE TO THE REPRESENTATIVE (CORE) GENOME USING TWO WHOLE GENOME ALIGNERS, NUCMER AND PARSNP (THESE LOCATIONS ARE IN $locs_dir).
# AS SUCH, WE KNOW WHAT SNPs TO EXPECT IN EACH VCF.
# TO EVALUATE THESE VCFs, WE'LL FIRST STORE AS OUR "TRUTH SET" THOSE CALLS MADE BY BOTH NUCMER AND PARSNP.
my %true_snp_locs = (); my %repetitive_regions = ();
opendir(DIR,$locs_dir) or die $!;
my @genomes = readdir(DIR);
closedir(DIR) or die $!;
foreach my $genome (@genomes)
	{ next if (($genome eq '.') or ($genome eq '..'));
	  
	  my $vcf1 = "$locs_dir/$genome/$genome.varcall_relative_to_representative.nucmer.vcf";
	  my $vcf2 = "$locs_dir/$genome/$genome.varcall_relative_to_representative.parsnp.vcf";
	  next if ( (!(-e($vcf1))) or (!(-e($vcf2))) );
	  
	  print "obtaining SNP positions from $genome...\n";
	  
	  # store positions in internally repetitive regions, identified by self-self BLASTn (the output format is described at http://www.metagenomics.wiki/tools/blast/blastn-output-format-6)
	  my $blastn_out = "$locs_dir/$genome/$genome.blastn.outfmt6";
	  if (-e($blastn_out))
		{ open(IN,$blastn_out) or die $!;
		  while(<IN>)
			{ next if ($. <= 2); # we want to ignore the first (and best) blast hit, which is of the entire sequence against itself. We're concerned only with internally repetitive regions.
			  my $line = $_; chomp($line);
		      my @line = split(/\t/,$line);
			  my $q_chr = $line[0]; my $t_chr = $line[1]; my $pc_identity = $line[2]; my $alignment_length = $line[3]; my $gap_opens = $line[5]; my $q_start = $line[6]; my $q_end = $line[7]; my $t_start = $line[8]; my $t_end = $line[9]; my $e_value = $line[10]; # q = query, t = target
			  if (($alignment_length >= 100) and ($pc_identity >= 95) and ($e_value <= 0.05) and ($gap_opens <= 1))
				{ for(my $pos=$q_start;$pos<=$q_end;$pos++)	{ $repetitive_regions{$genome}{$pos}++; }
				  for(my $pos=$t_start;$pos<=$t_end;$pos++)	{ $repetitive_regions{$genome}{$pos}++; }
				}
			}
		  close(IN) or die $!;
		}
	  my $num_of_masked_positions = scalar keys %{$repetitive_regions{$genome}};
	  
	  # store nucmer calls
	  open(IN,$vcf1) or die $!;
	  while(<IN>)
		{ next if ($. == 1);
		  my $line = $_; chomp($line);
		  next if ($line =~ /^\#/);
		  my @line = split(/\t/,$line);
		  my $chr = $line[0]; my $pos = $line[1]; my $ref = $line[3]; my $alt = $line[4];
		  if ( ((length($ref))==1) && ((length($alt))==1) )
			{ next if ( ($exclude_repetitive_regions == 1) and (exists($repetitive_regions{$genome}{$pos})) );
			  $true_snp_locs{$genome}{"$chr:$pos"}{ref}{nucmer} = $ref;
			  $true_snp_locs{$genome}{"$chr:$pos"}{alt}{nucmer} = $alt;
			}
		}
	  close(IN) or die $!;
	  
	  # store ParSnp calls
	  open(IN,$vcf2) or die $!;
	  while(<IN>)
		{ next if ($. == 1);
		  my $line = $_; chomp($line);
		  next if ($line =~ /^\#/);
		  my @line = split(/\t/,$line);
		  my $chr = $line[0]; my $pos = $line[1]; my $ref = $line[3]; my $alt = $line[4];
		  if ( ((length($ref))==1) && ((length($alt))==1) )
			{ next if ( ($exclude_repetitive_regions == 1) and (exists($repetitive_regions{$genome}{$pos})) );
			  $true_snp_locs{$genome}{"$chr:$pos"}{ref}{parsnp} = $ref;
			  $true_snp_locs{$genome}{"$chr:$pos"}{alt}{parsnp} = $alt;
			}
		}
	  close(IN) or die $!;
	}

# PARSE VCFs
opendir(DIR,$in_dir) or die $!;
my @genome_ids = readdir(DIR);
closedir(DIR) or die $!;
my @sorted_genome_ids = sort {$a cmp $b} @genome_ids;
my $genome_ids_seen = 0; my $genome_ids_total = @genome_ids; $genome_ids_total = $genome_ids_total-2;
my %scoring_metrics = ();
foreach my $genome_id (@sorted_genome_ids)
	{ next if (($genome_id eq '.') or ($genome_id eq '..'));
	  next if (!(-d("$in_dir/$genome_id")));
	  next if (!(exists($true_snp_locs{$genome_id})));
	  $genome_ids_seen++;
	  opendir(DIR,"$in_dir/$genome_id") or die $!;
	  my @vcfs = readdir(DIR);
	  closedir(DIR) or die $!;
	  my %snps_detected = ();
	  my $files_seen = 0; my $files_total = @vcfs; $files_total = $files_total-2;
	  foreach my $vcf (@vcfs)
		{ next if (($vcf eq '.') or ($vcf eq '..'));
		  $files_seen++;
		  next if (-d($vcf));
		  next if ($vcf !~ /^.*?\.vcf$/);
		  if ($vcf =~ /^(.*?)\.(.*?)\.(.*?)\.vcf$/)
			{ my $genome = $1; my $aligner = $2; my $caller = $3;
			  next if ( (!(exists($aligners{$aligner}))) || (!(exists($callers{$caller}))) );
			  print "$genome_ids_seen of $genome_ids_total --> $files_seen of $files_total\n";
			  open(IN,"$in_dir/$genome_id/$vcf") or die $!;
			  while(<IN>)
				{ my $line = $_; chomp($line);
				  next if ($line =~ /^\#/);
				  my @line = split(/\t/,$line);
				  my $chr = $line[0]; my $pos = $line[1]; my $ref = $line[3]; my $alt = $line[4]; my $qual = $line[5]; my $filter = $line[6]; my $info = $line[7];
				  next if ( ($exclude_repetitive_regions == 1) and (exists($repetitive_regions{$genome}{"$chr:$pos"})) );
				  my $format = ''; my $type = '';
				  if (defined($line[8])) { $format = $line[8]; }
				  if (defined($line[9])) { $type   = $line[9]; }
				  my @format = split(/\:/,$format);
				  my @type   = split(/\:/,$type);
				  my $snp_status = '';
				  for(my $x=0;$x<@format;$x++)
					{ my $format_value = $format[$x];
					  my $type_value   = $type[$x];
					  if ($format_value eq 'GT')
						{ $snp_status = $type_value; }
					}
				  next if (($snp_status eq '0/0') or ($snp_status eq '0|0')); # CHECKPOINT 1: exclude 0/0 or 0|0 calls
				  next if ((length($ref))!=(length($alt))); # CHECKPOINT 2: this cannot be a SNP; rather, it's likely an indel (or MNP)
				  next if (($ref !~ /^[ATCG]{1,}$/) || ($alt !~ /^[ATCG]{1,}$/)); # CHECKPOINT 3: this cannot be a biallelic SNP or MNP
				  next if (($ref eq 'N') or ($alt eq 'N') or ($ref eq '.') or ($alt eq '.')); # CHECKPOINT 4: exclude N calls (SolSNP makes these) and no-calls (Pilon makes these)
				  my $snp_type = '';
				  if 	(($snp_status eq '0/1') or ($snp_status eq '0|1')) { $snp_type = 'heterozygote';	   }
				  elsif (($snp_status eq '1/1') or ($snp_status eq '1|1')) { $snp_type = 'homozygous variant'; }
				  else						   							   { $snp_type = 'other';			   } # there shouldn't be any of these...
				  next if ($filter ne 'PASS'); # CHECKPOINT 5: this does not meet various exclusion criteria either specifically applied re: minimum mapping quality and depth of coverage (see 3.filter_VCFs.pl), or the internal criteria of that pipeline (which is particularly notable in the case of, e.g., Strelka)
				  if ( ((length($ref))==(length($alt))) && ((length($ref))>1) ) # MNP
					{ my @ref = split(//,$ref);
					  my @alt = split(//,$alt);
					  for(my $x=0;$x<@ref;$x++)
						{ my $ref = $ref[$x]; my $alt = $alt[$x];
						  if ($x > 0) { $pos++; }
						  if ($ref ne $alt)
							{ $snps_detected{$genome}{$aligner}{$caller}{"$chr:$pos"}{ref}      = $ref;
							  $snps_detected{$genome}{$aligner}{$caller}{"$chr:$pos"}{alt}      = $alt;
							  $snps_detected{$genome}{$aligner}{$caller}{"$chr:$pos"}{snp_type} = $snp_type;
							}
						}
					}
				  else # SNP
					{ $snps_detected{$genome}{$aligner}{$caller}{"$chr:$pos"}{ref}      = $ref;
					  $snps_detected{$genome}{$aligner}{$caller}{"$chr:$pos"}{alt}      = $alt;
					  $snps_detected{$genome}{$aligner}{$caller}{"$chr:$pos"}{snp_type} = $snp_type;
					}
				}
			  close(IN) or die $!;
			}
		}
	
	  # OUTPUT PIPELINE PERFORMANCE SCORES
	  foreach my $aligner (@aligners)
		{ foreach my $caller (@callers)
			{ next if (($aligner eq 'snippy')   and ($caller ne 'snippy'));
			  next if (($aligner ne 'snippy')   and ($caller eq 'snippy'));
			  next if (($aligner eq 'spandx')   and ($caller ne 'spandx'));
			  next if (($aligner ne 'spandx')   and ($caller eq 'spandx'));
			  next if (($aligner eq 'speedseq') and ($caller ne 'speedseq'));
			  next if (($aligner ne 'speedseq') and ($caller eq 'speedseq'));
			  next if (!(exists($snps_detected{$genome_id}{$aligner}{$caller})));
			  
			  my $num_of_true_snps = 0;
			  my $true_positives = 0; my $false_positives = 0; my $false_negatives = 0;
			  my $false_positives_that_are_heterozygotes = 0;
			  print "calculating F-score for: $genome_id, $aligner, $caller...\n";
			  while((my $snp_loc,my $irrel)=each(%{$true_snp_locs{$genome_id}}))
				{ my $true_ref_nucmer = ''; my $true_alt_nucmer = '';
				  if (exists($true_snp_locs{$genome_id}{$snp_loc}{ref}{nucmer})) { $true_ref_nucmer = $true_snp_locs{$genome_id}{$snp_loc}{ref}{nucmer}; }
				  if (exists($true_snp_locs{$genome_id}{$snp_loc}{alt}{nucmer})) { $true_alt_nucmer = $true_snp_locs{$genome_id}{$snp_loc}{alt}{nucmer}; }
				  my $true_ref_parsnp = ''; my $true_alt_parsnp = '';
				  if (exists($true_snp_locs{$genome_id}{$snp_loc}{ref}{parsnp})) { $true_ref_parsnp = $true_snp_locs{$genome_id}{$snp_loc}{ref}{parsnp}; }
				  if (exists($true_snp_locs{$genome_id}{$snp_loc}{alt}{parsnp})) { $true_alt_parsnp = $true_snp_locs{$genome_id}{$snp_loc}{alt}{parsnp}; }
				  
				  # the purpose of the following checkpoints is to ensure we're looking only at SNP calls made at unambiguous sites, i.e. those where nucmer & ParSnp are concordant
				  next if (($true_ref_nucmer eq '') or ($true_ref_parsnp eq '') or ($true_alt_nucmer eq '') or ($true_alt_parsnp eq ''));
				  next if ($true_ref_nucmer ne $true_ref_parsnp);
				  next if ($true_alt_nucmer ne $true_alt_parsnp);

				  $num_of_true_snps++;
				  my $true_ref = $true_ref_nucmer; my $true_alt = $true_alt_nucmer;
				  if (exists($snps_detected{$genome_id}{$aligner}{$caller}{$snp_loc}))
					{ my $ref = $snps_detected{$genome_id}{$aligner}{$caller}{$snp_loc}{ref};
					  my $alt = $snps_detected{$genome_id}{$aligner}{$caller}{$snp_loc}{alt};
					  if (($ref eq $true_ref) and ($alt eq $true_alt))
						{ $true_positives++; } # the position exists in the true set, and we've called the SNP correctly: this is a true positive
					  else
						{ $false_negatives++; } # the position exists in the true set, but we've not found the SNP: this is a false negative
					}
				  else
					{ $false_negatives++; } # the position exists in the true set, but we've not found the SNP: this is a false negative
				}
			  while((my $snp_loc,my $irrel)=each(%{$snps_detected{$genome_id}{$aligner}{$caller}}))
				{ next if (exists($true_snp_locs{$genome_id}{$snp_loc}));
				  # by skipping any position in %true_snp_locs, we are ALSO skipping the ambiguous positions, i.e. those where $true_snp_locs{$genome_id}{$snp_loc}{ref}{nucmer} ne $true_snp_locs{$genome_id}{$snp_loc}{ref}{parsnp} (or if either of these values don't exist). As such, we won't be counting as FPs those calls made at sites where we can't impute the reference position.
				  # to elaborate, the %true_snp_locs hash contains both ambiguous AND non-ambiguous positions - the former are those where nucmer and ParSnp calls are not concordant, and so are skipped; see the checkpoints in the above loop. As such, by the time we get to HERE, every call at a position we're confident of has been classified as either TP or FN (and we won't have counted calls made at ambiguous positions as being anything). The purpose of this "next" command is to ensure we won't be counting as FPs those calls made at positions we're either sure OR unsure of. It would otherwise be possible for us to call as a false positive a SNP made at a position where nucmer and ParSnp don't agree - which isn't a fair judgement to make as we aren't ourselves sure what the true reference position is.
				  
				  $false_positives++; # we've called something in a position that isn't in the true set: this is a false positive
				  my $type = $snps_detected{$genome_id}{$aligner}{$caller}{$snp_loc}{snp_type};
				  if ($type eq 'heterozygote') { $false_positives_that_are_heterozygotes++; }
				}
			
		      # calculate precision (positive predictive value), recall (sensitivity) and miss rate, as in https://www.nature.com/articles/srep17875
			  my $precision    = 0; if (($true_positives+$false_positives) > 0) { $precision = sprintf("%.4f",($true_positives/($true_positives+$false_positives)));  }
			  my $recall	   = 0; if (($true_positives+$false_negatives) > 0) { $recall    = sprintf("%.4f",($true_positives/($true_positives+$false_negatives)));  }
			  my $f1_score     = 0; if (($precision+$recall) 			   > 0) { $f1_score  = sprintf("%.4f",(2*(($recall*$precision)/($precision+$recall))));		  }
			  my $miss_rate    = 1; if (($true_positives+$false_negatives) > 0) { $miss_rate = sprintf("%.4f",($false_negatives/($true_positives+$false_negatives))); }
			  my $total_errors = $false_positives+$false_negatives;
			  my $num_of_million_bases_sequenced = 'NA';
			  if (-e("$fq_dir/$genome_id.count"))
				{ open(IN,"$fq_dir/$genome_id.count") or die $!;
				  while(<IN>)
					{ my $line = $_; chomp($line);
					  $num_of_million_bases_sequenced = $line;
					}
				  close(IN) or die $!;
				  $num_of_million_bases_sequenced = sprintf("%.2f",($num_of_million_bases_sequenced/1000000));
				}
			  my $total_fp_per_million_bases = sprintf("%.4f",($false_positives/$num_of_million_bases_sequenced));
			  my $total_fn_per_million_bases = sprintf("%.4f",($false_negatives/$num_of_million_bases_sequenced));
			  my $total_errors_per_million_bases = 'NA';
			  if ($num_of_million_bases_sequenced !~ /^NA$/)
				{ $total_errors_per_million_bases = sprintf("%.2f",($total_errors/$num_of_million_bases_sequenced)); }
			  $total_errors = $total_errors_per_million_bases;
			  
			  # express true positives, false positives and false negatives as a % of the number of true SNPs
			  my $true_positives_as_pc_of_true_snps  = sprintf("%.4f",(($true_positives/$num_of_true_snps)*100));
			  my $false_positives_as_pc_of_true_snps = sprintf("%.4f",(($false_positives/$num_of_true_snps)*100));
			  my $false_negatives_as_pc_of_true_snps = sprintf("%.4f",(($false_negatives/$num_of_true_snps)*100));
			  
			  print OUT1 "$genome_id\t$num_of_million_bases_sequenced\t$aligner\t$caller\t$aligner/$caller\t$num_of_true_snps\t$true_positives\t$true_positives_as_pc_of_true_snps\t$false_positives\t$false_positives_that_are_heterozygotes\t$false_positives_as_pc_of_true_snps\t$false_negatives\t$false_negatives_as_pc_of_true_snps\t$total_errors\t$total_fp_per_million_bases\t$total_fn_per_million_bases\t$precision\t$recall\t$f1_score\t$miss_rate\n";
			  push(@{$scoring_metrics{caller}{f1}{$caller}},$f1_score);
			  push(@{$scoring_metrics{caller}{recall}{$caller}},$recall);
			  push(@{$scoring_metrics{caller}{precision}{$caller}},$precision);
			  push(@{$scoring_metrics{caller}{total_errors}{$caller}},$total_errors);
			  push(@{$scoring_metrics{caller}{true_positives}{$caller}},$true_positives);
			  push(@{$scoring_metrics{caller}{false_positives}{$caller}},$false_positives);
			  push(@{$scoring_metrics{caller}{false_negatives}{$caller}},$false_negatives);
			  push(@{$scoring_metrics{caller}{total_fp_per_million_bases}{$caller}},$total_fp_per_million_bases);
			  push(@{$scoring_metrics{caller}{total_fn_per_million_bases}{$caller}},$total_fn_per_million_bases);			  
			  push(@{$scoring_metrics{aligner}{f1}{$aligner}},$f1_score);
			  push(@{$scoring_metrics{aligner}{recall}{$aligner}},$recall);
			  push(@{$scoring_metrics{aligner}{precision}{$aligner}},$precision);
			  push(@{$scoring_metrics{aligner}{total_errors}{$aligner}},$total_errors);
			  push(@{$scoring_metrics{aligner}{true_positives}{$aligner}},$true_positives);
			  push(@{$scoring_metrics{aligner}{false_positives}{$aligner}},$false_positives);
			  push(@{$scoring_metrics{aligner}{false_negatives}{$aligner}},$false_negatives);
			  push(@{$scoring_metrics{aligner}{total_fp_per_million_bases}{$aligner}},$total_fp_per_million_bases);
			  push(@{$scoring_metrics{aligner}{total_fn_per_million_bases}{$aligner}},$total_fn_per_million_bases);			  
			  push(@{$scoring_metrics{"aligner/caller"}{f1}{"$aligner/$caller"}},$f1_score);
			  push(@{$scoring_metrics{"aligner/caller"}{recall}{"$aligner/$caller"}},$recall);
			  push(@{$scoring_metrics{"aligner/caller"}{precision}{"$aligner/$caller"}},$precision);
			  push(@{$scoring_metrics{"aligner/caller"}{total_errors}{"$aligner/$caller"}},$total_errors);
			  push(@{$scoring_metrics{"aligner/caller"}{true_positives}{"$aligner/$caller"}},$true_positives);
			  push(@{$scoring_metrics{"aligner/caller"}{false_positives}{"$aligner/$caller"}},$false_positives);
			  push(@{$scoring_metrics{"aligner/caller"}{false_negatives}{"$aligner/$caller"}},$false_negatives);
			  push(@{$scoring_metrics{"aligner/caller"}{total_fp_per_million_bases}{"$aligner/$caller"}},$total_fp_per_million_bases);
			  push(@{$scoring_metrics{"aligner/caller"}{total_fn_per_million_bases}{"$aligner/$caller"}},$total_fn_per_million_bases);
			}
		}
	}
close(OUT1) or die $!;

if ($genome_ids_seen == 0)
	{ print "ERROR: no VCFs detected in $in_dir; unable to proceed\n";
	  exit 1;
	}
	
# OBTAIN PERFORMANCE RANKS FOR EACH PIPELINE
my %aligner_caller_ranks = ();
my @variables = ("aligner","caller","aligner/caller");
foreach my $variable (@variables)
	{ my @values = ();
	  if 	($variable eq 'aligner') { @values = @aligners; }
	  elsif ($variable eq 'caller')  { @values = @callers;  }
	  elsif ($variable eq 'aligner/caller')
		{ foreach my $aligner (@aligners)
			{ foreach my $caller (@callers)
				{ next if (($aligner eq 'snippy')   and ($caller ne 'snippy'));
				  next if (($aligner ne 'snippy')   and ($caller eq 'snippy'));
				  next if (($aligner eq 'spandx')   and ($caller ne 'spandx'));
				  next if (($aligner ne 'spandx')   and ($caller eq 'spandx'));
			      next if (($aligner eq 'speedseq') and ($caller ne 'speedseq'));
				  next if (($aligner ne 'speedseq') and ($caller eq 'speedseq'));
				  push(@values,"$aligner/$caller");
				}
			}
		}
	  if (($variable eq 'aligner') or ($variable eq 'caller') or ($variable eq 'aligner/caller'))
		{ my $highest_f1 = ''; my $highest_precision = ''; my $highest_recall = ''; my $highest_true_positives = ''; my $lowest_false_positives = ''; my $lowest_false_negatives = ''; my $lowest_total_errors = '';
		  my $lowest_f1 = ''; my $lowest_precision = ''; my $lowest_recall = ''; my $lowest_true_positives = ''; my $highest_false_positives = ''; my $highest_false_negatives = ''; my $highest_total_errors = '';
		  my @rank_order_f1 = (); my @rank_order_precision = (); my @rank_order_recall = (); my @rank_order_true_positives = (); my @rank_order_false_positives = (); my @rank_order_false_negatives = (); my @rank_order_total_errors = ();
		  my @metrics = (qw/f1 precision recall true_positives false_positives false_negatives total_errors/);
		  my $total_num_of_observations = 0;
		  foreach my $metric (@metrics)
			{ my @avg = ();
			  next if (!(exists($scoring_metrics{$variable}{$metric})));
			  my $num_of_observations = 0;
			  while((my $value,my $irrel)=each(%{$scoring_metrics{$variable}{$metric}}))
				{ my $num = @{$scoring_metrics{$variable}{$metric}{$value}};
				  my $avg_for_this_value = avg(@{$scoring_metrics{$variable}{$metric}{$value}});
				  push(@avg,[$avg_for_this_value,$value]);
				  $num_of_observations += $num;
				}
			  my @sorted_avg = ();
			  if (($metric eq 'total_errors') or ($metric eq 'false_negatives') or ($metric eq 'false_positives'))
				{ @sorted_avg = map { $_->[0] } sort { $a->[1] <=> $b->[1] } map { [$_, $_->[0]] } @avg; }
			  else
				{ @sorted_avg = map { $_->[0] } sort { $b->[1] <=> $a->[1] } map { [$_, $_->[0]] } @avg; }
			  my $best_score  = $sorted_avg[0][0]; $best_score = sprintf("%.4f",$best_score);
			  my $best_value  = $sorted_avg[0][1];
			  my $worst_score = $sorted_avg[$#sorted_avg][0]; $worst_score = sprintf("%.4f",$worst_score);
			  my $worst_value = $sorted_avg[$#sorted_avg][1];
			  if 	($metric eq 'f1')		  	   { $highest_f1 			 = "$best_value ($best_score)"; $lowest_f1 			     = "$worst_value ($worst_score)"; $total_num_of_observations = $num_of_observations; @rank_order_f1 			 = @sorted_avg; }
			  elsif ($metric eq 'recall')		   { $highest_recall 		 = "$best_value ($best_score)"; $lowest_recall           = "$worst_value ($worst_score)"; $total_num_of_observations = $num_of_observations; @rank_order_recall 		 = @sorted_avg; }
			  elsif ($metric eq 'precision')	   { $highest_precision 	 = "$best_value ($best_score)"; $lowest_precision        = "$worst_value ($worst_score)"; $total_num_of_observations = $num_of_observations; @rank_order_precision 	 	 = @sorted_avg; }
			  elsif ($metric eq 'total_errors')	   { $lowest_total_errors	 = "$best_value ($best_score)"; $highest_total_errors    = "$worst_value ($worst_score)"; $total_num_of_observations = $num_of_observations; @rank_order_total_errors    = @sorted_avg; }
			  elsif ($metric eq 'true_positives')  { $highest_true_positives = "$best_value ($best_score)"; $lowest_true_positives   = "$worst_value ($worst_score)"; $total_num_of_observations = $num_of_observations; @rank_order_true_positives  = @sorted_avg; }
			  elsif ($metric eq 'false_positives') { $lowest_false_positives = "$best_value ($best_score)"; $highest_false_positives = "$worst_value ($worst_score)"; $total_num_of_observations = $num_of_observations; @rank_order_false_positives = @sorted_avg; }
			  elsif ($metric eq 'false_negatives') { $lowest_false_negatives = "$best_value ($best_score)"; $highest_false_negatives = "$worst_value ($worst_score)"; $total_num_of_observations = $num_of_observations; @rank_order_false_negatives = @sorted_avg; }
			}
		  my $num = 0;
		  for(my $x=0;$x<@rank_order_f1;$x++)
			{ $num = $x+1;
			  my $rank_f1 			   = $rank_order_f1[$x][0]; 			 my $rank_f1_score 			    = $rank_order_f1[$x][1];
			  my $rank_recall 		   = $rank_order_recall[$x][0]; 		 my $rank_recall_score 		    = $rank_order_recall[$x][1];
			  my $rank_precision 	   = $rank_order_precision[$x][0]; 	     my $rank_precision_score 	    = $rank_order_precision[$x][1];
			  my $rank_total_errors    = $rank_order_total_errors[$x][0];    my $rank_total_errors_score    = $rank_order_total_errors[$x][1];
			  my $rank_true_positives  = $rank_order_true_positives[$x][0];  my $rank_true_positives_score  = $rank_order_true_positives[$x][1];
			  my $rank_false_positives = $rank_order_false_positives[$x][0]; my $rank_false_positives_score = $rank_order_false_positives[$x][1];
			  my $rank_false_negatives = $rank_order_false_negatives[$x][0]; my $rank_false_negatives_score = $rank_order_false_negatives[$x][1];
			  $rank_f1 = sprintf("%.4f",$rank_f1);
			  $rank_recall = sprintf("%.4f",$rank_recall);
			  $rank_precision = sprintf("%.4f",$rank_precision);
			  $rank_total_errors = sprintf("%.4f",$rank_total_errors);
			  $rank_true_positives = sprintf("%.4f",$rank_true_positives);
			  $rank_false_positives = sprintf("%.4f",$rank_false_positives);
			  $rank_false_negatives = sprintf("%.4f",$rank_false_negatives);
			  if ($variable eq 'aligner/caller')
				{ push(@{$aligner_caller_ranks{'recall'}},[$rank_recall,$rank_recall_score]);
				  push(@{$aligner_caller_ranks{'f-score'}},[$rank_f1,$rank_f1_score]);
				  push(@{$aligner_caller_ranks{'precision'}},[$rank_precision,$rank_precision_score]);
				  push(@{$aligner_caller_ranks{'total errors'}},[$rank_total_errors,$rank_total_errors_score]);
				  push(@{$aligner_caller_ranks{'true positives'}},[$rank_true_positives,$rank_true_positives_score]);
				  push(@{$aligner_caller_ranks{'false positives'}},[$rank_false_positives,$rank_false_positives_score]);
				  push(@{$aligner_caller_ranks{'false negatives'}},[$rank_false_negatives,$rank_false_negatives_score]);
				}
			}
		}
	  else
		{ foreach my $value (@values)
			{ my $highest_f1 = ''; my $highest_precision = ''; my $highest_recall = ''; my $highest_true_positives = ''; my $lowest_false_positives = ''; my $lowest_false_negatives = ''; my $lowest_total_errors = '';
			  my $lowest_f1 = ''; my $lowest_precision = ''; my $lowest_recall = ''; my $lowest_true_positives = ''; my $highest_false_positives = ''; my $highest_false_negatives = ''; my $highest_total_errors = '';
			  my @rank_order_f1 = (); my @rank_order_precision = (); my @rank_order_recall = (); my @rank_order_true_positives = (); my @rank_order_false_positives = (); my @rank_order_false_negatives = (); my @rank_order_total_errors = ();
			  my @metrics = (qw/f1 precision recall true_positives false_positives false_negatives total_errors/);
			  my $total_num_of_observations = 0;
			  foreach my $metric (@metrics)
				{ my @avg = ();
				  next if (!(exists($scoring_metrics{$variable}{$value}{$metric})));
				  my $num_of_observations = 0;
				  while((my $combo,my $irrel)=each(%{$scoring_metrics{$variable}{$value}{$metric}}))
					{ my $num = @{$scoring_metrics{$variable}{$value}{$metric}{$combo}};
					  my $avg_for_this_combo = avg(@{$scoring_metrics{$variable}{$value}{$metric}{$combo}});
					  push(@avg,[$avg_for_this_combo,$combo]);
					  $num_of_observations += $num;
					}
				  my @sorted_avg = ();
				  if (($metric eq 'total_errors') or ($metric eq 'false_negatives') or ($metric eq 'false_positives'))
					{ @sorted_avg = map { $_->[0] } sort { $a->[1] <=> $b->[1] } map { [$_, $_->[0]] } @avg; }
				  else
					{ @sorted_avg = map { $_->[0] } sort { $b->[1] <=> $a->[1] } map { [$_, $_->[0]] } @avg; }
				  my $best_score  = $sorted_avg[0][0]; $best_score = sprintf("%.4f",$best_score);
				  my $best_value  = $sorted_avg[0][1];
				  my $worst_score = $sorted_avg[$#sorted_avg][0]; $worst_score = sprintf("%.4f",$worst_score);
				  my $worst_value = $sorted_avg[$#sorted_avg][1];
				  if 	($metric eq 'f1')		  	   { $highest_f1 			 = "$best_value ($best_score)"; $lowest_f1 			     = "$worst_value ($worst_score)"; $total_num_of_observations = $num_of_observations; @rank_order_f1 			 = @sorted_avg; }
				  elsif ($metric eq 'recall')		   { $highest_recall 		 = "$best_value ($best_score)"; $lowest_recall           = "$worst_value ($worst_score)"; $total_num_of_observations = $num_of_observations; @rank_order_recall 		 = @sorted_avg; }
				  elsif ($metric eq 'precision')	   { $highest_precision 	 = "$best_value ($best_score)"; $lowest_precision        = "$worst_value ($worst_score)"; $total_num_of_observations = $num_of_observations; @rank_order_precision 	 	 = @sorted_avg; }
				  elsif ($metric eq 'total_errors')	   { $lowest_total_errors	 = "$best_value ($best_score)"; $highest_total_errors    = "$worst_value ($worst_score)"; $total_num_of_observations = $num_of_observations; @rank_order_total_errors    = @sorted_avg; }
				  elsif ($metric eq 'true_positives')  { $highest_true_positives = "$best_value ($best_score)"; $lowest_true_positives   = "$worst_value ($worst_score)"; $total_num_of_observations = $num_of_observations; @rank_order_true_positives  = @sorted_avg; }
				  elsif ($metric eq 'false_positives') { $lowest_false_positives = "$best_value ($best_score)"; $highest_false_positives = "$worst_value ($worst_score)"; $total_num_of_observations = $num_of_observations; @rank_order_false_positives = @sorted_avg; }
				  elsif ($metric eq 'false_negatives') { $lowest_false_negatives = "$best_value ($best_score)"; $highest_false_negatives = "$worst_value ($worst_score)"; $total_num_of_observations = $num_of_observations; @rank_order_false_negatives = @sorted_avg; }
				}
			  my $num = 0;
			  for(my $x=0;$x<@rank_order_f1;$x++)
				{ $num = $x+1;
				  my $rank_f1 			   = $rank_order_f1[$x][0]; 			 my $rank_f1_score 			    = $rank_order_f1[$x][1];
				  my $rank_recall 		   = $rank_order_recall[$x][0]; 		 my $rank_recall_score 		    = $rank_order_recall[$x][1];
				  my $rank_precision 	   = $rank_order_precision[$x][0]; 	     my $rank_precision_score 	    = $rank_order_precision[$x][1];
				  my $rank_total_errors    = $rank_order_total_errors[$x][0];    my $rank_total_errors_score    = $rank_order_total_errors[$x][1];
				  my $rank_true_positives  = $rank_order_true_positives[$x][0];  my $rank_true_positives_score  = $rank_order_true_positives[$x][1];
				  my $rank_false_positives = $rank_order_false_positives[$x][0]; my $rank_false_positives_score = $rank_order_false_positives[$x][1];
				  my $rank_false_negatives = $rank_order_false_negatives[$x][0]; my $rank_false_negatives_score = $rank_order_false_negatives[$x][1];
				  $rank_f1 = sprintf("%.4f",$rank_f1);
				  $rank_recall = sprintf("%.4f",$rank_recall);
				  $rank_precision = sprintf("%.4f",$rank_precision);
				  $rank_total_errors = sprintf("%.4f",$rank_total_errors);
				  $rank_true_positives = sprintf("%.4f",$rank_true_positives);
				  $rank_false_positives = sprintf("%.4f",$rank_false_positives);
				  $rank_false_negatives = sprintf("%.4f",$rank_false_negatives);
				}
			}
		}
	}

# WHAT IS THE ALL-ROUND BEST PERFORMING PIPELINE? IT WILL HAVE THE LOWEST SUM OF ALL RANKS.
my %summed_aligner_and_caller_ranks = ();
my @metrics = ("f-score","precision","recall","true positives","false positives","false negatives","total errors");
foreach my $metric (@metrics)
	{ my @arr = @{$aligner_caller_ranks{$metric}};
	  my @sorted_arr = ();
	  if (($metric eq 'total errors') or ($metric eq 'false negatives') or ($metric eq 'false positives'))
		{ @sorted_arr = map { $_->[0] } sort { $a->[1] <=> $b->[1] } map { [$_, $_->[0]] } @arr; }
	  else
		{ @sorted_arr = map { $_->[0] } sort { $b->[1] <=> $a->[1] } map { [$_, $_->[0]] } @arr; }
	  my $rank = 0;
	  for(my $x=0;$x<@sorted_arr;$x++)
		{ my $score = $sorted_arr[$x][0]; my $aligner_and_caller = $sorted_arr[$x][1];
		  my $prev_score;
		  if (defined($sorted_arr[$x-1]))
			{ $prev_score = $sorted_arr[$x-1][0]; }
		  $rank++ unless ((defined($prev_score)) and ($prev_score == $score));
		  $summed_aligner_and_caller_ranks{$aligner_and_caller}{overall_rank} += $rank;
		  push(@{$summed_aligner_and_caller_ranks{$aligner_and_caller}{rank_per_metric}},"$rank ($metric: $score)");
		}
	}
my @arr = ();
while((my $aligner_and_caller,my $irrel)=each(%summed_aligner_and_caller_ranks))
	{ my $rank 			  =   $summed_aligner_and_caller_ranks{$aligner_and_caller}{overall_rank};
	  my @rank_per_metric = @{$summed_aligner_and_caller_ranks{$aligner_and_caller}{rank_per_metric}};
	  my $rank_per_metric = join(", ",@rank_per_metric);
	  push(@arr,[$rank,$aligner_and_caller,$rank_per_metric]);
	}
my @sorted_arr = map { $_->[0] } sort { $a->[1] <=> $b->[1] } map { [$_, $_->[0]] } @arr;
for(my $x=0;$x<@sorted_arr;$x++)
	{ my $rank = $x+1;
	  my $score = $sorted_arr[$x][0]; my $aligner_and_caller = $sorted_arr[$x][1]; my $rank_per_metric = $sorted_arr[$x][2];
	  print OUT2 "$rank\t$aligner_and_caller\t$score\t$rank_per_metric\n";
	}
close(OUT2) or die $!;	
exit 1;