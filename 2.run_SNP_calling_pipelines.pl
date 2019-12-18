=head
This script is provided in support of Bush, et al. Genomic diversity affects the accuracy of bacterial SNP calling pipelines. https://www.biorxiv.org/content/10.1101/653774v1
It is provided for reproducibility purposes only. No attempt has been made to optimise the code.

Prerequisites:
1. Download the data archive from: https://ora.ox.ac.uk/objects/uuid:8f902497-955e-4b84-9b85-693ee0e4433e. Its contents are referred to here as $fq_dir and $rep_dir.
2. This archive also contains information regarding third-party tools, in the file "Installation_instructions.txt". This script will test for the presence of these tools. Note that Novoalign is proprietary software and, if available, needs to be accessible via $PATH.

Fuller details are available in the readme provided with the above archive.
=cut

use strict;
use warnings;
use POSIX qw(ceil);

# REQUIREMENTS
my $root    = '/home/ndm.local/steveb';
my $path    = "$root/varcall_evaluation"; # in this directory we will create (a) a bash script containing the commands necessary to populate (b) a sub-directory to hold the VCFs
my $progs   = "$root/programs"; # all third-party tools in "Installation instructions.txt" are installed here
my $fq_dir  = "$path/illumina_reads_for_nanopore_illumina_hybrid_assemblies";
my $rep_dir = "$path/representative_genomes_for_nanopore_illumina_hybrid_assemblies"; # a directory of representative genomes for each of the sequenced species, with choices informed by Kaiju, Kraken and MALDI TOF predictions

# CHECK THAT ALL DEPENDENCIES ARE AVAILABLE
my $fatal = 0; # ALL THE FOLLOWING PATHS ARE HARD-CODED PREREQUISITES. WE WILL ABORT IF THEY CANNOT BE FOUND.
if (!(-d($path)))    { $fatal++; print "error: cannot find $path\n";    }
if (!(-d($progs)))   { $fatal++; print "error: cannot find $progs\n";   }
if (!(-d($fq_dir)))  { $fatal++; print "error: cannot find $fq_dir\n";  }
if (!(-d($rep_dir))) { $fatal++; print "error: cannot find $rep_dir\n"; }
## PATHS TO INDEXERS
my $BGSBuild_path		= "$progs/4.2-amazon/BGS-Build";
my $gem_build_path      = "$progs/gem3-mapper/bin/gem-indexer";
my $hisat2build_path    = "$progs/hisat2-2.1.0/hisat2-build";
my $bowtie2build_path   = "$progs/bowtie2-2.3.4.1-linux-x86_64/bowtie2-build";
my $yara_indexer_path   = "$progs/yara-build/bin/yara_indexer";
my $soap3dpbuilder_path = "$progs/4.2-amazon/soap3-dp-builder";
if (!(-e($BGSBuild_path)))       { $fatal++; print "error: cannot find $BGSBuild_path\n";       }
if (!(-e($gem_build_path)))      { $fatal++; print "error: cannot find $gem_build_path\n";      }
if (!(-e($hisat2build_path)))    { $fatal++; print "error: cannot find $hisat2build_path\n";    }
if (!(-e($bowtie2build_path)))   { $fatal++; print "error: cannot find $bowtie2build_path\n";   }
if (!(-e($yara_indexer_path)))   { $fatal++; print "error: cannot find $yara_indexer_path\n";   }
if (!(-e($soap3dpbuilder_path))) { $fatal++; print "error: cannot find $soap3dpbuilder_path\n"; }
## PATHS TO ALIGNERS
my $bwa_dir 	 	   = "$progs/bwa-0.7.17";
my $bwa_path 	 	   = "$progs/bwa-0.7.17/bwa";
my $gem_path		   = "$progs/gem3-mapper/bin/gem-mapper";
my $ngm_path		   = "$progs/Cibiv-NextGenMap-33e92fb/bin/ngm-0.5.5/ngm";
my $snap_path		   = "$progs/snap-aligner";
my $yara_path		   = "$progs/yara-build/bin/yara_mapper";
my $smalt_path		   = "$progs/bin/smalt";
my $gassst_path 	   = "$progs/Gassst_v1.28/Gassst";
my $hisat2_path		   = "$progs/hisat2-2.1.0/hisat2";
my $stampy_path  	   = "$progs/stampy-1.0.32/stampy.py";
my $bowtie2_dir 	   = "$progs/bowtie2-2.3.4.1-linux-x86_64";
my $bowtie2_path 	   = "$progs/bowtie2-2.3.4.1-linux-x86_64/bowtie2";
my $cushaw3_path	   = "$progs/cushaw3";
my $minimap2_path	   = "$progs/minimap2/minimap2";
my $mosaik_pe_ann	   = "$progs/MOSAIK/src/networkFile/2.1.78.pe.ann";
my $mosaik_se_ann	   = "$progs/MOSAIK/src/networkFile/2.1.78.se.ann";
my $MosaikBuild_path   = "$progs/MOSAIK/bin/MosaikBuild";
my $MosaikAligner_path = "$progs/MOSAIK/bin/MosaikAligner";
if (!(-e($bwa_path))) 			{ $fatal++; print "error: cannot find $bwa_path\n"; 		  }
if (!(-e($gem_path))) 			{ $fatal++; print "error: cannot find $gem_path\n"; 		  }
if (!(-e($ngm_path))) 			{ $fatal++; print "error: cannot find $ngm_path\n"; 		  }
if (!(-e($snap_path))) 			{ $fatal++; print "error: cannot find $snap_path\n"; 		  }
if (!(-e($yara_path))) 			{ $fatal++; print "error: cannot find $yara_path\n"; 		  }
if (!(-e($smalt_path))) 		{ $fatal++; print "error: cannot find $smalt_path\n"; 		  }
if (!(-e($gassst_path))) 		{ $fatal++; print "error: cannot find $gassst_path\n"; 		  }
if (!(-e($hisat2_path))) 		{ $fatal++; print "error: cannot find $hisat2_path\n"; 		  }
if (!(-e($stampy_path))) 		{ $fatal++; print "error: cannot find $stampy_path\n"; 	 	  }
if (!(-e($bowtie2_path))) 		{ $fatal++; print "error: cannot find $bowtie2_path\n"; 	  }
if (!(-e($cushaw3_path))) 		{ $fatal++; print "error: cannot find $cushaw3_path\n"; 	  }
if (!(-e($minimap2_path))) 	    { $fatal++; print "error: cannot find $minimap2_path\n"; 	  }
if (!(-e($mosaik_pe_ann)))		{ $fatal++; print "error: cannot find $mosaik_pe_ann\n";	  }
if (!(-e($mosaik_se_ann)))		{ $fatal++; print "error: cannot find $mosaik_se_ann\n";	  }
if (!(-e($MosaikBuild_path)))	{ $fatal++; print "error: cannot find $MosaikBuild_path\n";	  }
if (!(-e($MosaikAligner_path)))	{ $fatal++; print "error: cannot find $MosaikAligner_path\n"; }
## PATHS TO CALLERS
my $I6gt_path      	   	   = "$progs/16GT";
my $gatk_path      		   = "$progs/gatk-4.0.1.2/gatk";
my $snver_path	   		   = "$progs/SNVerIndividual.jar";
my $pilon_path			   = "$progs/pilon-1.23.jar";
my $prokka_dir			   = "$progs/prokka/bin"; # is required to be on the $PATH in order for Snippy to work
my $snippy_path			   = "$progs/snippy/bin/snippy";
my $spandx_path			   = "$progs/SPANDx/SPANDx.sh";
my $lofreq_path	   		   = "$progs/lofreq_star-2.1.2/bin/lofreq";
my $solsnp_path			   = "$progs/SolSNP/SolSNP.jar";
my $strelka_dir   		   = "$progs/strelka-2.9.3.centos6_x86_64/bin";
my $varscan_path   		   = "$progs/VarScan.v2.3.9.jar";
my $octopus_path  		   = "$progs/octopus/bin/octopus";
my $bcftools_dir  		   = "$progs/bcftools";
my $bcftools_path  		   = "$progs/bcftools/bcftools";
my $platypus_path  		   = "$progs/Platypus_0.8.1/Platypus.py";
my $speedseq_path		   = "$progs/speedseq/bin/speedseq";
my $freebayes_path 		   = "$progs/freebayes/bin/freebayes";
my $SNVSniffer_path		   = "$progs/SNVSniffer";
my $bam2snapshot_path 	   = "$I6gt_path/bam2snapshot";
my $snapshotSnpcaller_path = "$I6gt_path/snapshotSnpcaller";
if (!(-e($I6gt_path)))      	    { $fatal++; print "error: cannot find $I6gt_path\n";	  		  }
if (!(-e($gatk_path)))      	    { $fatal++; print "error: cannot find $gatk_path\n";	  		  }
if (!(-e($snver_path))) 			{ $fatal++; print "error: cannot find $snver_path\n";  	  		  }
if (!(-e($pilon_path))) 			{ $fatal++; print "error: cannot find $pilon_path\n";  	  		  }
if (!(-d($prokka_dir))) 			{ $fatal++; print "error: cannot find $prokka_dir\n";  	  		  }
if (!(-e($snippy_path))) 			{ $fatal++; print "error: cannot find $snippy_path\n";  	  	  }
if (!(-e($spandx_path))) 			{ $fatal++; print "error: cannot find $spandx_path\n";  	  	  }
if (!(-e($lofreq_path))) 			{ $fatal++; print "error: cannot find $lofreq_path\n";    		  }
if (!(-e($solsnp_path))) 			{ $fatal++; print "error: cannot find $solsnp_path\n";    		  }
if (!(-d($strelka_dir))) 			{ $fatal++; print "error: cannot find $strelka_dir\n";   		  }
if (!(-e($varscan_path))) 			{ $fatal++; print "error: cannot find $varscan_path\n";   		  }
if (!(-e($octopus_path))) 			{ $fatal++; print "error: cannot find $octopus_path\n";   		  }
if (!(-e($bcftools_path))) 			{ $fatal++; print "error: cannot find $bcftools_path\n";  		  }
if (!(-e($platypus_path)))  		{ $fatal++; print "error: cannot find $platypus_path\n";  		  }
if (!(-e($speedseq_path)))  		{ $fatal++; print "error: cannot find $speedseq_path\n";  		  }
if (!(-e($freebayes_path))) 	    { $fatal++; print "error: cannot find $freebayes_path\n"; 		  }
if (!(-e($SNVSniffer_path)))		{ $fatal++; print "error: cannot find $SNVSniffer_path\n";		  }
if (!(-e($bam2snapshot_path))) 		{ $fatal++; print "error: cannot find $bam2snapshot_path\n"; 	  } # a constituent of 16GT
if (!(-e($snapshotSnpcaller_path))) { $fatal++; print "error: cannot find $snapshotSnpcaller_path\n"; } # a constituent of 16GT
## PATHS TO ANCILLARY SOFTWARE
my $vt_dir    	  		  	  = "$progs/vt";
my $seqtk_dir    	  		  = "$progs/seqtk";
my $seqtk_path    	  		  = "$progs/seqtk/seqtk";
my $wgsim_path	  			  = "$progs/samtools-1.7/bin/wgsim";
my $picard_path   			  = "$progs/picard.jar";
my $fqtools_path			  = "$progs/fqtools/bin/fqtools";
my $samtools_dir			  = "$progs/samtools-1.7/bin/";
my $samtools_path 			  = "$progs/samtools-1.7/bin/samtools";
my $sam_1pt99_path			  = "$progs/sam-1.99.jar"; # required to post-process MOSAIK BAMs, to fix an error in their indexing bins
my $gassst_to_sam_path		  = "$progs/Gassst_v1.28/gassst_to_sam";
my $vcfallelicprimitives_path = "$progs/vcflib/bin/vcfallelicprimitives";
if (!(-d($vt_dir)))    			   	   { $fatal++; print "error: cannot find $vt_dir\n";    			    }
if (!(-d($seqtk_dir)))    			   { $fatal++; print "error: cannot find $seqtk_dir\n";    			    }
if (!(-e($seqtk_path)))    			   { $fatal++; print "error: cannot find $seqtk_path\n";    			}
if (!(-e($wgsim_path)))    			   { $fatal++; print "error: cannot find $wgsim_path\n";    			}
if (!(-e($picard_path)))   			   { $fatal++; print "error: cannot find $picard_path\n";   			}
if (!(-e($fqtools_path)))   		   { $fatal++; print "error: cannot find $fqtools_path\n";   			}
if (!(-d($samtools_dir))) 			   { $fatal++; print "error: cannot find $samtools_dir\n"; 			    }
if (!(-e($samtools_path))) 			   { $fatal++; print "error: cannot find $samtools_path\n"; 			}
if (!(-e($sam_1pt99_path))) 		   { $fatal++; print "error: cannot find $sam_1pt99_path\n"; 			}
if (!(-e($gassst_to_sam_path))) 	   { $fatal++; print "error: cannot find $gassst_to_sam_path\n"; 		}
if (!(-e($vcfallelicprimitives_path))) { $fatal++; print "error: cannot find $vcfallelicprimitives_path\n"; }
exit if ($fatal > 0);

# OUTPUT
my $out_dir = "$path/vcfs";
if (!(-d($out_dir))) { mkdir $out_dir or die $!; }
my $sh_file = "$path/run_SNP_calling_pipelines.sh";
open(SH,'>',$sh_file) or die $!;
print SH "#!/bin/bash\n";
print SH "export C_INCLUDE_PATH=$progs/include\n"; # required for mpileup, platypus and fqtools
print SH "export LD_LIBRARY_PATH=$progs/lib\n"; # required for mpileup, platypus and fqtools
print SH "export LIBRARY_PATH=$progs/lib\n"; # required for mpileup, platypus and fqtools
print SH "export PATH=:$progs/bin:\$PATH\n"; # required for Snippy (specifically, to obtain GNU Parallel)
print SH "export PATH=:$bwa_dir:\$PATH\n"; # required for Snippy
print SH "export PATH=:$bcftools_dir:\$PATH\n"; # required for Snippy
print SH "export PATH=:$samtools_dir:\$PATH\n"; # required for Snippy. This prepend to PATH is because on the cluster used for this study, /usr/bin already contained a lesser version of samtools
print SH "PATH=\$PATH:$prokka_dir\n"; # required for Snippy
print SH "PATH=\$PATH:$seqtk_dir\n"; # required for Snippy
print SH "PATH=\$PATH:$vt_dir\n"; # required for Snippy
print SH "PATH=\$PATH:$progs/bin\n"; # required for samtools/wgsim
print SH "PATH=\$PATH:$progs/lib\n";
print SH "PATH=\$PATH:$progs/include\n";

# PARAMETERS
my $num_procs = 10;
## the arrays @aligners and @callers list all possible combinations evaluated in this study
my @aligners  = (qw/bbmap bowtie2 bwa-mem bwa-sw bwa+stampy cushaw3 gassst gem hisat2 minimap2 mosaik ngm novoalign smalt snap snippy spandx speedseq stampy yara/);
my @callers   = (qw/16GT deepvariant freebayes gatk lofreq mpileup octopus pilon platypus snippy snver snvsniffer solsnp spandx speedseq strelka varscan/);
=cut
## for the example $out_file provided in this archive, commands were output only for the following combinations, and only for $genome_id eq 'cft073'. These commands were run to populate the $out_dir provided in this archive.
my @aligners  = (qw/bwa-mem snippy speedseq/);
my @callers   = (qw/freebayes gatk lofreq mpileup snippy speedseq/);
=cut
my %aligners  = map {$_ => 1} @aligners;
my %callers   = map {$_ => 1} @callers;

opendir(DIR,$fq_dir) or die $!;
my @files = readdir(DIR);
closedir(DIR) or die $!;
my %genome_ids = ();
foreach my $file (@files)
	{ next if (($file eq '.') or ($file eq '..'));
	  if ($file =~ /^(.*?)\.\d+\.fq.gz$/)
		{ $genome_ids{$1}++; }
	}
my @genome_ids = ();
while((my $genome_id,my $irrel)=each(%genome_ids))
	{ push(@genome_ids,$genome_id); }
my @sorted_genome_ids = sort {$a cmp $b} @genome_ids;
foreach my $genome_id (@sorted_genome_ids)
	{ my $dict = ''; my $ref = ''; my $idx = '';
	  if (($genome_id eq 'rbhstw00029') or ($genome_id eq 'rbhstw00053') or ($genome_id eq 'rbhstw00127') or ($genome_id eq 'rbhstw00350'))
		{ $dict = "$rep_dir/Citrobacter_freundii.dict";
		  $ref  = "$rep_dir/Citrobacter_freundii.fa";
		  $idx  = "$rep_dir/Citrobacter_freundii";
		}
	  elsif (($genome_id eq 'rbhstw00059') or ($genome_id eq 'rbhstw00131') or ($genome_id eq 'rbhstw00340'))	  
		{ $dict = "$rep_dir/Enterobacter_cloacae.dict";
		  $ref  = "$rep_dir/Enterobacter_cloacae.fa";
		  $idx  = "$rep_dir/Enterobacter_cloacae";
		}
	  elsif ($genome_id eq 'rbhstw00309')
		{ $dict = "$rep_dir/Enterobacter_kobei.dict";
		  $ref  = "$rep_dir/Enterobacter_kobei.fa";
		  $idx  = "$rep_dir/Enterobacter_kobei";
		}
	  elsif (($genome_id eq 'cft073') or ($genome_id eq 'rbhstw00122') or ($genome_id eq 'rbhstw00277') or ($genome_id eq 'rhb10c07') or ($genome_id eq 'rhb11c04'))
		{ $dict = "$rep_dir/Escherichia_coli.dict";
		  $ref  = "$rep_dir/Escherichia_coli.fa";
		  $idx  = "$rep_dir/Escherichia_coli";
		}
	  elsif (($genome_id eq 'rbhstw00167') or ($genome_id eq 'rbhstw00189'))
		{ $dict = "$rep_dir/Klebsiella_oxytoca.dict";
		  $ref  = "$rep_dir/Klebsiella_oxytoca.fa";
		  $idx  = "$rep_dir/Klebsiella_oxytoca";
		}
	  elsif (($genome_id eq 'mgh78578') or ($genome_id eq 'rbhstw00128') or ($genome_id eq 'rhb14c01'))
		{ $dict = "$rep_dir/Klebsiella_pneumoniae.dict";
		  $ref  = "$rep_dir/Klebsiella_pneumoniae.fa";
		  $idx  = "$rep_dir/Klebsiella_pneumoniae";
		}
	  next if (!(-e($ref))); # we're not concerned at this point that $dict and $idx do not exist as commands to create them will be output as part of this script
	  
	  my $root = $genome_id;
  	  my $fq_1 = "$fq_dir/$genome_id.1.fq.gz";
	  my $fq_2 = "$fq_dir/$genome_id.2.fq.gz";
	  next if ( (!(-e($fq_1))) or (!(-e($fq_2))) ); # CHECKPOINT: there are no reads to align!

	  if (!(-d("$out_dir/$genome_id"))) { mkdir "$out_dir/$genome_id" or die $!; }
	  
	  # HOW MANY READS ARE IN EACH FQ?
	  print SH "$fqtools_path count $fq_1 > $fq_dir/$genome_id.count\n" unless (-e("$fq_dir/$genome_id.count"));
	  
	  # IS THERE ANYTHING STILL TO DO FOR THIS $root?
	  my %to_run = ();
	  foreach my $aligner (@aligners)
		{ foreach my $caller (@callers)
			{ next if (-e("$out_dir/$root/$root.$aligner.$caller.vcf")); # CHECKPOINT: we've produced output from this particular combination of aligner and caller before
			  
			  # Snippy, SpeedSeq and SPANDx are all-in-one self-contained pipelines, so cannot be paired with any other $aligner or $caller
			  next if (($aligner eq 'snippy') 	and ($caller  ne 'snippy'));
			  next if (($caller  eq 'snippy') 	and ($aligner ne 'snippy'));
			  next if (($aligner eq 'speedseq') and ($caller  ne 'speedseq'));
			  next if (($caller  eq 'speedseq') and ($aligner ne 'speedseq'));
			  next if (($aligner eq 'spandx') 	and ($caller  ne 'spandx'));
			  next if (($caller  eq 'spandx') 	and ($aligner ne 'spandx'));
			  
			  # for various reasons, these combinations were not found to work on our system; see 'operating notes' in supplementary text 1
			  next if (($aligner eq 'bbmap')    and ($caller  eq '16GT'));
			  next if (($aligner eq 'bbmap')    and ($caller  eq 'pilon'));
			  next if (($aligner eq 'bbmap')    and ($caller  eq 'snver'));
			  next if (($aligner eq 'bbmap')    and ($caller  eq 'solsnp'));
			  next if (($aligner eq 'bbmap')    and ($caller  eq 'strelka'));
			  next if (($aligner eq 'gassst')   and ($caller  eq '16GT'));
			  next if (($aligner eq 'gassst')   and ($caller  eq 'deepvariant'));
			  next if (($aligner eq 'gassst')   and ($caller  eq 'freebayes'));
			  next if (($aligner eq 'gassst')   and ($caller  eq 'gatk'));
			  next if (($aligner eq 'gassst')   and ($caller  eq 'lofreq'));
			  next if (($aligner eq 'gassst')   and ($caller  eq 'octopus'));
			  next if (($aligner eq 'gassst')   and ($caller  eq 'platypus'));
			  next if (($aligner eq 'gassst')   and ($caller  eq 'snver'));
			  next if (($aligner eq 'gassst')   and ($caller  eq 'solsnp'));
			  next if (($aligner eq 'gassst')   and ($caller  eq 'strelka'));
			  next if (($aligner eq 'hisat2')   and ($caller  eq 'deepvariant'));
			  next if (($aligner eq 'hisat2')   and ($caller  eq 'octopus'));
			  next if (($aligner eq 'hisat2')   and ($caller  eq 'platypus'));
			  next if (($aligner eq 'mosaik')   and ($caller  eq 'lofreq'));
			  next if (($aligner eq 'mosaik')   and ($caller  eq 'mpileup'));
			  next if (($aligner eq 'mosaik')   and ($caller  eq 'pilon'));
			  next if (($aligner eq 'mosaik')   and ($caller  eq 'platypus'));
			  next if (($aligner eq 'mosaik')   and ($caller  eq 'snver'));
			  next if (($aligner eq 'mosaik')   and ($caller  eq 'snvsniffer'));
			  next if (($aligner eq 'mosaik')   and ($caller  eq 'strelka'));
			  next if (($aligner eq 'mosaik')   and ($caller  eq 'varscan'));
		      $to_run{$aligner}{$caller}++;
			}
		}
	  next if (scalar keys %to_run == 0); # CHECKPOINT: don't go on if there's nothing to do
	  
	  # CREATE INDEXES, IF NECESSARY
	  print SH "cd $rep_dir\n";
	  print SH "$bwa_path index $ref\n" 		   unless ( (-e("$ref.amb")) && (-e("$ref.ann")) && (-e("$ref.bwt")) && (-e("$ref.pac")) && (-e("$ref.sa")) );
	  print SH "$gem_build_path -i $ref $idx -t $num_procs\n" unless ( (-e("$idx.gem")) && ("$idx.info") );
	  print SH "$snap_path index $ref $idx -t$num_procs -bSpace\n" unless (-d($idx));
	  print SH "$yara_indexer_path $ref -o $ref.yara_idx\n" unless ( (-e("$ref.yara_idx.lf.drv")) && (-e("$ref.yara_idx.lf.drs")) && (-e("$ref.yara_idx.lf.drp")) && (-e("$ref.yara_idx.lf.pst")) && (-e("$ref.yara_idx.sa.ind")) && (-e("$ref.yara_idx.sa.val")) && (-e("$ref.yara_idx.sa.len")) && (-e("$ref.yara_idx.txt.limits")) && (-e("$ref.yara_idx.rid.limits")) && (-e("$ref.yara_idx.rid.concat")) && (-e("$ref.yara_idx.txt.size")) && (-e("$ref.yara_idx.txt.concat")) );
	  print SH "$hisat2build_path $ref $idx\n"     unless ( (-e("$idx.1.ht2")) && (-e("$idx.2.ht2")) && (-e("$idx.3.ht2")) && (-e("$idx.4.ht2")) && (-e("$idx.5.ht2")) && (-e("$idx.6.ht2")) && (-e("$idx.7.ht2")) && (-e("$idx.8.ht2")) );
	  print SH "$bowtie2build_path $ref $idx\n"    unless ( (-e("$idx.1.bt2")) && (-e("$idx.2.bt2")) && (-e("$idx.3.bt2")) && (-e("$idx.4.bt2")) && (-e("$idx.rev.1.bt2")) && (-e("$idx.rev.2.bt2")) );
	  print SH "$soap3dpbuilder_path $ref\n"	   unless ( (-e("$ref.index.amb")) && (-e("$ref.index.amb")) && (-e("$ref.index.ann")) &&(-e("$ref.index.bwt")) && (-e("$ref.index.fmv")) && (-e("$ref.index.lkt")) && (-e("$ref.index.pac")) && (-e("$ref.index.rev.bwt")) && (-e("$ref.index.rev.fmv")) && (-e("$ref.index.rev.lkt")) && (-e("$ref.index.rev.pac")) && (-e("$ref.index.sa")) && (-e("$ref.index.tra")) );
	  print SH "$BGSBuild_path $ref.index\n"	   unless ( (-e("$ref.index.fmv.gpu")) && (-e("$ref.index.fmv.mmap")) && (-e("$ref.index.pac.mmap")) && (-e("$ref.index.rev.fmv.gpu")) && (-e("$ref.index.rev.fmv.mmap")) );
	  print SH "$cushaw3_path index $ref\n"		   unless ( (-e("$ref.ann")) && (-e("$ref.map")) && (-e("$ref.pac")) && (-e("$ref.rbwt")) && (-e("$ref.rsa")) );
	  print SH "$smalt_path index $idx $ref\n"	   unless ( (-e("$idx.sma")) && (-e("$idx.smi")) );
	  print SH "$samtools_path faidx $ref\n" 	   unless   (-e("$ref.fai"));
	  print SH "$minimap2_path -d $idx.mmi $ref\n" unless   (-e("$idx.mmi"));
	  print SH "$stampy_path -G $idx $ref\n" 	   unless   (-e("$idx.stidx"));
	  print SH "$stampy_path -g $idx -H $idx\n"    unless   (-e("$idx.sthash"));
	  print SH "novoindex -t $num_procs $idx.ndx $ref\n" unless (-e("$idx.ndx")); # IMPORTANT: Novoalign is proprietary software and, if available, needs to be accessible via $PATH.
	  print SH "java -jar $picard_path CreateSequenceDictionary REFERENCE=$ref OUTPUT=$dict\n" unless (-e($dict));
	  print SH "cd $out_dir\n";
	  
	  # to use DeepVariant (via Docker), make Docker-mountable input and output directories, then copy into the former $ref and $ref.fai
	  if (exists($callers{'deepvariant'}))
		{ if (!(-d("$out_dir/deepvariant_input")))  { print SH "mkdir $out_dir/deepvariant_input\n";  }
		  if (!(-d("$out_dir/deepvariant_output"))) { print SH "mkdir $out_dir/deepvariant_output\n"; }
		  print SH "cp $ref $out_dir/deepvariant_input\n";
		  print SH "cp $ref.fai $out_dir/deepvariant_input\n";
		}
	  
	  # ALIGN READS
	  foreach my $aligner (@aligners)
		{ next if (!(exists($to_run{$aligner})));
		  if ($aligner eq 'bbmap')
			{ print SH "bbmap.sh ref=$ref in=$fq_1 in2=$fq_2 out=$out_dir/$root.$aligner.unsorted.no_read_groups.sam nodisk\n";
			  print SH "$samtools_path view -bS $out_dir/$root.$aligner.unsorted.no_read_groups.sam > $out_dir/$root.$aligner.unsorted.no_read_groups.bam\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.no_read_groups.sam\n";
			  print SH "java -jar $picard_path CleanSam INPUT=$out_dir/$root.$aligner.unsorted.no_read_groups.bam OUTPUT=$out_dir/$root.$aligner.unsorted.no_read_groups.cleaned.bam TMP_DIR=$out_dir\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.no_read_groups.bam\n";
			  print SH "java -jar $picard_path FixMateInformation INPUT=$out_dir/$root.$aligner.unsorted.no_read_groups.cleaned.bam OUTPUT=$out_dir/$root.$aligner.unsorted.no_read_groups.cleaned.fixmated.bam TMP_DIR=$out_dir\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.no_read_groups.cleaned.bam\n";
			  print SH "java -jar $picard_path AddOrReplaceReadGroups INPUT=$out_dir/$root.$aligner.unsorted.no_read_groups.cleaned.fixmated.bam OUTPUT=$out_dir/$root.$aligner.unsorted.bam RGID=$root RGLB=lib RGPL=Illumina RGPU=unit RGSM=sample TMP_DIR=$out_dir\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.no_read_groups.cleaned.fixmated.bam\n";
			}
		  elsif ($aligner eq 'bowtie2')
			{ print SH "$bowtie2_path --rg-id $root --rg SM:sample --rg LIB:lib --rg PL:Illumina -p $num_procs -x $idx -1 $fq_1 -2 $fq_2 -S $out_dir/$root.$aligner.unsorted.sam\n";
			  print SH "$samtools_path view -bS $out_dir/$root.$aligner.unsorted.sam > $out_dir/$root.$aligner.unsorted.bam\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.sam\n";
			  print SH "java -jar $picard_path CleanSam INPUT=$out_dir/$root.$aligner.unsorted.bam OUTPUT=$out_dir/$root.$aligner.unsorted.cleaned.bam TMP_DIR=$out_dir\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.bam\n";
			  print SH "mv $out_dir/$root.$aligner.unsorted.cleaned.bam $out_dir/$root.$aligner.unsorted.bam\n";
			}
		  elsif ($aligner eq 'bwa-mem')
			{ print SH "$bwa_path mem -R '\@RG\\tID:group_$root\\tSM:sample_$root\\tPL:Illumina\\tLIB:lib_$root\\tPU:unit_$root' -t $num_procs -M $ref $fq_1 $fq_2 | $samtools_path view -Shb - > $out_dir/$root.$aligner.unsorted.bam\n";
			  print SH "java -jar $picard_path CleanSam INPUT=$out_dir/$root.$aligner.unsorted.bam OUTPUT=$out_dir/$root.$aligner.unsorted.cleaned.bam TMP_DIR=$out_dir\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.bam\n";
			  print SH "mv $out_dir/$root.$aligner.unsorted.cleaned.bam $out_dir/$root.$aligner.unsorted.bam\n";
			}
		  elsif ($aligner eq 'bwa-sw')
			{ my $concat_fq = "$fq_dir/$genome_id.fq.gz";
			  print SH "cat $fq_1 $fq_2 > $concat_fq\n";
			  print SH "$bwa_path bwasw -t $num_procs -f $out_dir/$root.$aligner.unsorted.no_read_groups.sam $ref $concat_fq\n";
			  print SH "rm $concat_fq\n";
			  print SH "$samtools_path view -bS $out_dir/$root.$aligner.unsorted.no_read_groups.sam > $out_dir/$root.$aligner.unsorted.no_read_groups.bam\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.no_read_groups.sam\n";
			  print SH "java -jar $picard_path CleanSam INPUT=$out_dir/$root.$aligner.unsorted.no_read_groups.bam OUTPUT=$out_dir/$root.$aligner.unsorted.no_read_groups.cleaned.bam TMP_DIR=$out_dir\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.no_read_groups.bam\n";
			  print SH "java -jar $picard_path AddOrReplaceReadGroups INPUT=$out_dir/$root.$aligner.unsorted.no_read_groups.cleaned.bam OUTPUT=$out_dir/$root.$aligner.unsorted.bam RGID=$root RGLB=lib RGPL=Illumina RGPU=unit RGSM=sample TMP_DIR=$out_dir\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.no_read_groups.cleaned.bam\n";
			}
		  elsif ($aligner eq 'bwa+stampy') # using command lines as given in http://www.well.ox.ac.uk/~gerton/README.txt
			{ print SH "$bwa_path aln -t $num_procs $ref $fq_1 > $out_dir/$root.$aligner.1.sai\n";
			  print SH "$bwa_path aln -t $num_procs $ref $fq_2 > $out_dir/$root.$aligner.2.sai\n";
			  print SH "$bwa_path sampe $ref $out_dir/$root.$aligner.1.sai $out_dir/$root.$aligner.2.sai $fq_1 $fq_2 | $samtools_path view -Sb - > $out_dir/$root.$aligner.unsorted.pre_stampy.bam\n";
			  print SH "$stampy_path -g $idx -h $idx -t $num_procs --bamkeepgoodreads -M $out_dir/$root.$aligner.unsorted.pre_stampy.bam | $samtools_path view -Sb - > $out_dir/$root.$aligner.unsorted.no_read_groups.bam\n";
			  print SH "rm $out_dir/$root.$aligner.1.sai $out_dir/$root.$aligner.2.sai $out_dir/$root.$aligner.unsorted.pre_stampy.bam\n";
			  print SH "java -jar $picard_path CleanSam INPUT=$out_dir/$root.$aligner.unsorted.no_read_groups.bam OUTPUT=$out_dir/$root.$aligner.unsorted.no_read_groups.cleaned.bam TMP_DIR=$out_dir\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.no_read_groups.bam\n";
			  print SH "java -jar $picard_path AddOrReplaceReadGroups INPUT=$out_dir/$root.$aligner.unsorted.no_read_groups.cleaned.bam OUTPUT=$out_dir/$root.$aligner.unsorted.bam RGID=$root RGLB=lib RGPL=Illumina RGPU=unit RGSM=sample TMP_DIR=$out_dir\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.no_read_groups.cleaned.bam\n";
			}
		  elsif ($aligner eq 'cushaw3')
			{ print SH "$cushaw3_path align -rgid $root -rgsm sample -rglb lib -rgpl Illumina -t $num_procs -r $ref -q $fq_1 $fq_2 -o $out_dir/$root.$aligner.unsorted.sam\n";
			  print SH "$samtools_path view -bS $out_dir/$root.$aligner.unsorted.sam > $out_dir/$root.$aligner.unsorted.bam\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.sam\n";
			  print SH "java -jar $picard_path CleanSam INPUT=$out_dir/$root.$aligner.unsorted.bam OUTPUT=$out_dir/$root.$aligner.unsorted.cleaned.bam TMP_DIR=$out_dir\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.bam\n";
			  print SH "mv $out_dir/$root.$aligner.unsorted.cleaned.bam $out_dir/$root.$aligner.unsorted.bam\n";
			}
		  elsif ($aligner eq 'gassst')
			{ print SH "$seqtk_path seq -A $fq_1 > $out_dir/$genome_id.1.fa\n";
			  print SH "$seqtk_path seq -A $fq_2 > $out_dir/$genome_id.2.fa\n";
			  print SH "cat $out_dir/$genome_id.1.fa $out_dir/$genome_id.2.fa > $out_dir/$genome_id.fa\n";
			  print SH "rm $out_dir/$genome_id.1.fa $out_dir/$genome_id.2.fa\n";
			  print SH "$gassst_path -d $ref -i $out_dir/$genome_id.fa -p 95 -n $num_procs -o $out_dir/$root.$aligner.unsorted.no_read_groups.gassst\n";
			  print SH "$gassst_to_sam_path $out_dir/$root.$aligner.unsorted.no_read_groups.gassst $out_dir/$root.$aligner.unsorted.no_read_groups.sam\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.no_read_groups.gassst\n";
			  print SH "$samtools_path view -bS $out_dir/$root.$aligner.unsorted.no_read_groups.sam > $out_dir/$root.$aligner.unsorted.no_read_groups.bam\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.no_read_groups.sam\n";
			  print SH "java -jar $picard_path CleanSam INPUT=$out_dir/$root.$aligner.unsorted.no_read_groups.bam OUTPUT=$out_dir/$root.$aligner.unsorted.no_read_groups.cleaned.bam TMP_DIR=$out_dir\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.no_read_groups.bam\n";
			  print SH "java -jar $picard_path AddOrReplaceReadGroups INPUT=$out_dir/$root.$aligner.unsorted.no_read_groups.cleaned.bam OUTPUT=$out_dir/$root.$aligner.unsorted.bam RGID=$root RGLB=lib RGPL=Illumina RGPU=unit RGSM=sample TMP_DIR=$out_dir\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.no_read_groups.cleaned.bam\n";
			  print SH "rm $out_dir/$genome_id.fa\n";
			}
		  elsif ($aligner eq 'gem')
			{ print SH "$gem_path -I $idx.gem -1 $fq_1 -2 $fq_2 -r '\@RG\\tID:group_$root\\tSM:sample_$root\\tPL:Illumina\\tLIB:lib_$root\\tPU:unit_$root' -z -p -o $out_dir/$root.$aligner.unsorted.sam\n";
			  print SH "$samtools_path view -bS $out_dir/$root.$aligner.unsorted.sam > $out_dir/$root.$aligner.unsorted.bam\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.sam\n";
			  print SH "java -jar $picard_path CleanSam INPUT=$out_dir/$root.$aligner.unsorted.bam OUTPUT=$out_dir/$root.$aligner.unsorted.cleaned.bam TMP_DIR=$out_dir\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.bam\n";
			  print SH "mv $out_dir/$root.$aligner.unsorted.cleaned.bam $out_dir/$root.$aligner.unsorted.bam\n";
			}
		  elsif ($aligner eq 'hisat2')
			{ print SH "$hisat2_path --rg-id $root --rg SM:sample --rg LB:lib --rg PL:Illumina -p $num_procs -x $idx -1 $fq_1 -2 $fq_2 -S $out_dir/$root.$aligner.unsorted.sam\n";
			  print SH "$samtools_path view -bS $out_dir/$root.$aligner.unsorted.sam > $out_dir/$root.$aligner.unsorted.bam\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.sam\n";
			  print SH "java -jar $picard_path CleanSam INPUT=$out_dir/$root.$aligner.unsorted.bam OUTPUT=$out_dir/$root.$aligner.unsorted.cleaned.bam TMP_DIR=$out_dir\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.bam\n";
			  print SH "mv $out_dir/$root.$aligner.unsorted.cleaned.bam $out_dir/$root.$aligner.unsorted.bam\n";
			}
		  elsif ($aligner eq 'minimap2')
			{ print SH "$minimap2_path -ax sr $ref -R '\@RG\\tID:group_$root\\tSM:sample_$root\\tPL:Illumina\\tLIB:lib_$root\\tPU:unit_$root' -t $num_procs $fq_1 $fq_2 | $samtools_path view -Shb - > $out_dir/$root.$aligner.unsorted.bam\n";
			  print SH "java -jar $picard_path CleanSam INPUT=$out_dir/$root.$aligner.unsorted.bam OUTPUT=$out_dir/$root.$aligner.unsorted.cleaned.bam TMP_DIR=$out_dir\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.bam\n";
			  print SH "mv $out_dir/$root.$aligner.unsorted.cleaned.bam $out_dir/$root.$aligner.unsorted.bam\n";
			}
	      elsif ($aligner eq 'mosaik') # note that the MOSAIK aligner often produce BAMs with errors in their indexing bins (this is an optional field but should be correct): "bin field of BAM record does not equal value computed based on alignment start and end, and length of sequence to which read is aligned". This causes downstream problems with Picard MarkDuplicates and BAM indexing unless VALIDATION_STRINGENCY is set to "lenient" - or the file is corrected. See https://sourceforge.net/p/samtools/mailman/message/31853465/ and https://gatkforums.broadinstitute.org/gatk/discussion/4290/sam-bin-field-error-for-the-gatk-run
			{ print SH "$MosaikBuild_path -fr $ref -oa $out_dir/$root.dat\n";
			  print SH "$MosaikBuild_path -q $fq_1 -q2 $fq_2 -st illumina -out $out_dir/$root.$aligner.mkb\n";
			  print SH "$MosaikAligner_path -in $out_dir/$root.$aligner.mkb -ia $out_dir/$root.dat -out $out_dir/$root.$aligner.unsorted -p $num_procs -annpe $mosaik_pe_ann -annse $mosaik_se_ann\n";
			  print SH "rm $out_dir/$root.dat $out_dir/$root.$aligner.mkb $out_dir/$root.$aligner.unsorted.stat\n";
			  print SH "java -classpath $sam_1pt99_path net.sf.samtools.FixBAMFile $out_dir/$root.$aligner.unsorted.bam $out_dir/$root.$aligner.unsorted.fixed.bam\n";
			  print SH "mv $out_dir/$root.$aligner.unsorted.fixed.bam $out_dir/$root.$aligner.unsorted.bam\n";
			  print SH "java -jar $picard_path CleanSam INPUT=$out_dir/$root.$aligner.unsorted.bam OUTPUT=$out_dir/$root.$aligner.unsorted.cleaned.bam TMP_DIR=$out_dir\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.bam\n";
			  print SH "mv $out_dir/$root.$aligner.unsorted.cleaned.bam $out_dir/$root.$aligner.unsorted.bam\n";
			}
		  elsif ($aligner eq 'ngm')
			{ print SH "$ngm_path -t $num_procs -r $ref -1 $fq_1 -2 $fq_2 --rg-id $root --rg-sm sample --rg-lb lib --rg-pl Illumina -o $out_dir/$root.$aligner.unsorted.sam\n"; # why do we output SAM instead of taking the -b option to output BAM? Because of an error when validating index bins (https://gatkforums.broadinstitute.org/gatk/discussion/4290/sam-bin-field-error-for-the-gatk-run): "bin field of BAM record does not equal value computed based on alignment start and end, and length of sequence to which read is aligned". SAM does not have a field for indexing bin, so if we output in this format, there's no problem.
			  print SH "$samtools_path view -bS $out_dir/$root.$aligner.unsorted.sam > $out_dir/$root.$aligner.unsorted.bam\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.sam\n";
			  print SH "java -jar $picard_path CleanSam INPUT=$out_dir/$root.$aligner.unsorted.bam OUTPUT=$out_dir/$root.$aligner.unsorted.cleaned.bam TMP_DIR=$out_dir\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.bam\n";
			  print SH "mv $out_dir/$root.$aligner.unsorted.cleaned.bam $out_dir/$root.$aligner.unsorted.bam\n";
			}
		  elsif ($aligner eq 'novoalign') # IMPORTANT: Novoalign is proprietary software and, if available, needs to be accessible via $PATH.
			{ print SH "novoalign -c $num_procs -d $idx.ndx -F STDFQ -o SAM \$'\@RG\\tSM:sample\\tID:group\\tPL:Illumina\\tPU:unit\\tLB:lib' -f $fq_1 $fq_2 > $out_dir/$root.$aligner.unsorted.sam\n";
			  print SH "$samtools_path view -bS $out_dir/$root.$aligner.unsorted.sam > $out_dir/$root.$aligner.unsorted.bam\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.sam\n";
			  print SH "java -jar $picard_path CleanSam INPUT=$out_dir/$root.$aligner.unsorted.bam OUTPUT=$out_dir/$root.$aligner.unsorted.cleaned.bam TMP_DIR=$out_dir\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.bam\n";
			  print SH "mv $out_dir/$root.$aligner.unsorted.cleaned.bam $out_dir/$root.$aligner.unsorted.bam\n";
			}
		  elsif ($aligner eq 'smalt')
			{ print SH "$smalt_path map -n $num_procs -o $out_dir/$root.$aligner.unsorted.sam $idx $fq_1 $fq_2\n";
			  print SH "$samtools_path view -bS $out_dir/$root.$aligner.unsorted.sam > $out_dir/$root.$aligner.unsorted.no_read_groups.bam\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.sam\n";
			  print SH "java -jar $picard_path CleanSam INPUT=$out_dir/$root.$aligner.unsorted.no_read_groups.bam OUTPUT=$out_dir/$root.$aligner.unsorted.no_read_groups.cleaned.bam TMP_DIR=$out_dir\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.no_read_groups.bam\n";
			  print SH "java -jar $picard_path AddOrReplaceReadGroups INPUT=$out_dir/$root.$aligner.unsorted.no_read_groups.cleaned.bam OUTPUT=$out_dir/$root.$aligner.unsorted.bam RGID=$root RGLB=lib RGPL=Illumina RGPU=unit RGSM=sample TMP_DIR=$out_dir\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.no_read_groups.cleaned.bam\n";
			}
		  elsif ($aligner eq 'snap')
			{ print SH "$snap_path paired $idx $fq_1 $fq_2 -R '\@RG\\tID:group\\tSM:sample\\tPL:Illumina\\tPU:unit\\tLB:lib' -o $out_dir/$root.$aligner.unsorted.sam\n"; # -t $num_procs. Performance may be affected by multithreading: https://github.com/amplab/snap/issues/61
			  print SH "$samtools_path view -bS $out_dir/$root.$aligner.unsorted.sam > $out_dir/$root.$aligner.unsorted.bam\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.sam\n";
			  print SH "java -jar $picard_path CleanSam INPUT=$out_dir/$root.$aligner.unsorted.bam OUTPUT=$out_dir/$root.$aligner.unsorted.cleaned.bam TMP_DIR=$out_dir\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.bam\n";
			  print SH "mv $out_dir/$root.$aligner.unsorted.cleaned.bam $out_dir/$root.$aligner.unsorted.bam\n";
			}
		  elsif ($aligner eq 'speedseq')
			{ print SH "$speedseq_path align -R '\@RG\\tID:group\\tSM:sample\\tPL:Illumina\\tPU:unit\\tLB:lib' -t $num_procs -o $out_dir/$root.$aligner $ref $fq_1 $fq_2\n";
			}
		  elsif ($aligner eq 'stampy')
			{ print SH "$stampy_path -g $idx -h $idx -t $num_procs -M $fq_1 $fq_2 | $samtools_path view -Sb - > $out_dir/$root.$aligner.unsorted.no_read_groups.bam\n";
			  print SH "java -jar $picard_path CleanSam INPUT=$out_dir/$root.$aligner.unsorted.no_read_groups.bam OUTPUT=$out_dir/$root.$aligner.unsorted.no_read_groups.cleaned.bam TMP_DIR=$out_dir\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.no_read_groups.bam\n";
			  print SH "java -jar $picard_path AddOrReplaceReadGroups INPUT=$out_dir/$root.$aligner.unsorted.no_read_groups.cleaned.bam OUTPUT=$out_dir/$root.$aligner.unsorted.bam RGID=$root RGLB=lib RGPL=Illumina RGPU=unit RGSM=sample TMP_DIR=$out_dir\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.no_read_groups.cleaned.bam\n";
			}
		  elsif ($aligner eq 'yara')
			{ print SH "$yara_path $ref.yara_idx $fq_1 $fq_2 --threads $num_procs -rg '\@RG\\tID:group_$root\\tSM:sample_$root\\tPL:Illumina\\tLIB:lib_$root\\tPU:unit_$root' -output-file $out_dir/$root.$aligner.unsorted.uncleaned.bam\n";
			  print SH "java -jar $picard_path CleanSam INPUT=$out_dir/$root.$aligner.unsorted.uncleaned.bam OUTPUT=$out_dir/$root.$aligner.unsorted.bam TMP_DIR=$out_dir\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.uncleaned.bam\n";
			}
			
		  # POST-PROCESSING OF BAM FILES
		  if (($aligner ne 'snippy') and ($aligner ne 'spandx') and ($aligner ne 'speedseq'))
			{ print SH "java -jar $picard_path SortSam INPUT=$out_dir/$root.$aligner.unsorted.bam OUTPUT=$out_dir/$root.$aligner.sorted.bam SORT_ORDER=coordinate TMP_DIR=$out_dir\n"; # the lenient validation stringency is solely for the use of MOSAIK, and for the reason of a recently added validation of the index bin in Picard Tools v1.9: https://gatkforums.broadinstitute.org/gatk/discussion/4290/sam-bin-field-error-for-the-gatk-run
			  print SH "java -jar $picard_path MarkDuplicates INPUT=$out_dir/$root.$aligner.sorted.bam OUTPUT=$out_dir/$root.$aligner.bam METRICS_FILE=$out_dir/$root.$aligner.metrics ASSUME_SORTED=true\n";
			  print SH "java -jar $picard_path BuildBamIndex INPUT=$out_dir/$root.$aligner.bam\n";
			  print SH "rm $out_dir/$root.$aligner.unsorted.bam $out_dir/$root.$aligner.sorted.bam $out_dir/$root.$aligner.metrics\n";
			}
		  
		  # CALL VARIANTS
		  foreach my $caller (@callers)
			{ next if (-e("$out_dir/$root.$aligner.$caller.vcf")); # CHECKPOINT: we've seen this particular combination of aligner and caller before
			  next if (!(exists($to_run{$aligner}{$caller})));
			  
			  if ($caller eq '16GT')
				{ print SH "cd $I6gt_path\n";
				  print SH "$bam2snapshot_path -i $ref.index -b $out_dir/$root.$aligner.bam -o $out_dir/$root.$aligner.$caller\n";
				  print SH "$snapshotSnpcaller_path -i $ref.index -o $out_dir/$root.$aligner.$caller\n";
				  print SH "perl txt2vcf.pl $out_dir/$root.$aligner.$caller.txt $root.$aligner.$caller $ref > $out_dir/$root.$aligner.$caller.vcf\n";
				  print SH "rm $out_dir/$root.$aligner.$caller.alignmentQC.txt $out_dir/$root.$aligner.$caller.snapshot $out_dir/$root.$aligner.$caller.txt\n";
				  print SH "cd $out_dir\n";
				}
			  elsif ($caller eq 'deepvariant') # see https://github.com/google/deepvariant/blob/master/docs/deepvariant-quick-start.md
				{ print SH "cp $out_dir/$root.$aligner.bam $out_dir/deepvariant_input\n";
				  print SH "cp $out_dir/$root.$aligner.bai $out_dir/deepvariant_input\n";
				  my $ref_minus_root = $ref; $ref_minus_root =~ s/$rep_dir\///;
				  print SH "docker run -v $out_dir/deepvariant_input:/input -v $out_dir/deepvariant_output:/output gcr.io/deepvariant-docker/deepvariant:0.8.0 /opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=/input/$ref_minus_root --reads=/input/$root.$aligner.bam --output_vcf=/output/$root.$aligner.$caller.vcf --num_shards=$num_procs\n";
				  print SH "rm $out_dir/deepvariant_input/$root.$aligner.bam $out_dir/deepvariant_input/$root.$aligner.bai\n";
				  print SH "cp $out_dir/deepvariant_output/$root.$aligner.$caller.vcf $out_dir/$root.$aligner.$caller.vcf\n";
				}
			  elsif ($caller eq 'freebayes')
				{ print SH "$freebayes_path -f $ref --ploidy 1 $out_dir/$root.$aligner.bam > $out_dir/$root.$aligner.$caller.vcf\n"; # note that you must first install VCFlib (and make vcflib/bin available on $PATH) and then remove the hard-coded paths in freebayes-parallel! See https://github.com/ekg/freebayes/issues/376
				}
			  elsif ($caller eq 'gatk') # note: "to get into the granular details on what is going on for specific loci, remember you can use HaplotypeCaller's -bamout option +/- the --emitDroppedReads option" (see https://gatkforums.broadinstitute.org/gatk/discussion/8094/no-varaint-calling-with-wgsim-simulated-reads and https://software.broadinstitute.org/gatk/documentation/article.php?id=1235)
				{ print SH "$gatk_path HaplotypeCaller -R $ref -I $out_dir/$root.$aligner.bam -O $out_dir/$root.$aligner.$caller.vcf\n"; # the default output is only variant sites; to output all sites, append to this command line "--emit-ref-confidence BP_RESOLUTION"
				  print SH "rm $out_dir/$root.$aligner.$caller.vcf.idx\n";
				}
			  elsif ($caller eq 'lofreq') # see http://csb5.github.io/lofreq/commands/
				{ print SH "$lofreq_path call-parallel --pp-threads $num_procs -f $ref -o $out_dir/$root.$aligner.$caller.vcf $out_dir/$root.$aligner.bam\n";
				}
			  elsif ($caller eq 'mpileup') # see https://samtools.github.io/bcftools/howtos/variant-calling.html
				{ print SH "$bcftools_path mpileup -Ou -f $ref $out_dir/$root.$aligner.bam | $bcftools_path call --threads $num_procs --ploidy 1 -mv -Ov -o $out_dir/$root.$aligner.$caller.vcf\n"; # the "v" in "-mv" specifies that the output file only contain variant sites; omit this to output all sites
				}
			  elsif ($caller eq 'octopus')
				{ print SH "$octopus_path --reference $ref --reads $out_dir/$root.$aligner.bam --legacy --threads $num_procs > $out_dir/$root.$aligner.$caller.vcf\n";
				}
			  elsif ($caller eq 'pilon')
				{ print SH "java -jar $pilon_path --genome $ref --bam $out_dir/$root.$aligner.bam --threads $num_procs --outdir $out_dir/$root.$aligner.$caller --output $root.$aligner.$caller --vcf\n";
				  print SH "mv $out_dir/$root.$aligner.$caller/$root.$aligner.$caller.vcf $out_dir/$root.$aligner.$caller.vcf\n";
				  print SH "rm -r $out_dir/$root.$aligner.$caller\n";
				}
			  elsif ($caller eq 'platypus') # note: an error is possible with Platypus resulting in a segmentation fault: "Exception OverflowError: 'value too large to convert to short' in 'htslibWrapper.ReadIterator.get' ignored". This is due to CIGAR strings with H operations longer than the maximum value of short. See https://groups.google.com/forum/#!topic/platypus-users/xHp_ZwhyxuM
				{ print SH "python $platypus_path callVariants --bamFiles=$out_dir/$root.$aligner.bam --logFileName=$out_dir/$root.$aligner.$caller.log --refFile=$ref --output=$out_dir/$root.$aligner.$caller.vcf\n";
			      print SH "rm $out_dir/$root.$aligner.$caller.log\n";
				}
			  elsif (($aligner eq 'snippy') and ($caller eq 'snippy'))
				{ print SH "$snippy_path --cpus $num_procs --outdir $out_dir/$root-snippy --prefix $root --cleanup --ref $ref --R1 $fq_1 --R2 $fq_2\n";
				  print SH "mv $out_dir/$root-snippy/$root.filt.vcf $out_dir/$root.$aligner.$caller.vcf\n";
				  print SH "rm -r $out_dir/$root-snippy\n";
				}
			  elsif ($caller eq 'snver')
				{ print SH "java -jar $snver_path -i $out_dir/$root.$aligner.bam -r $ref -o $out_dir/$root.$aligner.$caller\n"; # if the BAM contains no valid alignments, SNVer will only produce an empty VCF, $out_dir/$root.$aligner.$caller.raw.vcf, not the filtered one: $out_dir/$root.$aligner.$caller.filter.vcf
				  print SH "if [ -f $out_dir/$root.$aligner.$caller.filter.vcf ]; then mv $out_dir/$root.$aligner.$caller.filter.vcf $out_dir/$root.$aligner.$caller.vcf; else mv $out_dir/$root.$aligner.$caller.raw.vcf $out_dir/$root.$aligner.$caller.vcf; fi\n";
				  print SH "rm $out_dir/$root.$aligner.$caller.failed.log $out_dir/$root.$aligner.$caller.indel.filter.vcf $out_dir/$root.$aligner.$caller.indel.raw.vcf $out_dir/$root.$aligner.$caller.filter.vcf $out_dir/$root.$aligner.$caller.raw.vcf\n";
				}
			  elsif ($caller eq 'snvsniffer')
				{ print SH "$samtools_path view -H $out_dir/$root.$aligner.bam > $out_dir/$root.$aligner.header.sam\n";
				  print SH "$SNVSniffer_path snp -f 2 -g $ref -o $out_dir/$root.$aligner.$caller.vcf $out_dir/$root.$aligner.header.sam $out_dir/$root.$aligner.bam\n";
				  print SH "rm $out_dir/$root.$aligner.header.sam\n";
				}
			  elsif ($caller eq 'solsnp')
				{ print SH "java -jar $solsnp_path INPUT=$out_dir/$root.$aligner.bam OUTPUT=$out_dir/$root.$aligner.$caller.vcf R=$ref OUTPUT_FORMAT=VCF\n";
				}
			  elsif (($aligner eq 'spandx') and ($caller eq 'spandx'))
				{ print SH "mkdir $out_dir/$root.$aligner.$caller\n";
				  print SH "cd $out_dir/$root.$aligner.$caller\n";
				  print SH "cp $fq_1 strain_1_sequence.fastq.gz\n";
				  print SH "cp $fq_2 strain_2_sequence.fastq.gz\n";
				  print SH "cp $ref $genome_id.fasta\n";
				  print SH "$spandx_path -r $genome_id -m yes -t Illumina -p PE -z yes\n";
				  print SH "cp $out_dir/$root.$aligner.$caller/Outputs/SNPs_indels_PASS/strain.snps.PASS.vcf $out_dir/$root.$aligner.$caller.vcf\n";
				  print SH "cd $out_dir\n";
				  print SH "rm -r $out_dir/$root.$aligner.$caller\n";
				}
			  elsif (($aligner eq 'speedseq') and ($caller eq 'speedseq'))
				{ print SH "$speedseq_path var -t $num_procs -o $out_dir/$root.$aligner.$caller $ref $out_dir/$root.$aligner.bam\n";
				  print SH "gunzip $out_dir/$root.$aligner.$caller.vcf.gz\n";
				  print SH "rm $out_dir/$root.$aligner.$caller.vcf.gz.tbi\n";
				}
			  elsif ($caller eq 'strelka')
				{ print SH "python $strelka_dir/configureStrelkaGermlineWorkflow.py --bam $out_dir/$root.$aligner.bam --referenceFasta $ref --runDir $out_dir/$root.$aligner.$caller\n";
				  print SH "python $out_dir/$root.$aligner.$caller/runWorkflow.py -m local -j $num_procs\n";
				  print SH "gunzip $out_dir/$root.$aligner.$caller/results/variants/genome.S1.vcf.gz\n";
				  print SH "mv $out_dir/$root.$aligner.$caller/results/variants/genome.S1.vcf $out_dir/$root.$aligner.$caller.vcf\n";
				  print SH "rm -r $out_dir/$root.$aligner.$caller\n";
				}
			  elsif ($caller eq 'varscan') # see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4278659/
				{ print SH "$samtools_path mpileup -B -q 1 -f $ref $out_dir/$root.$aligner.bam > $out_dir/$root.$aligner.$caller.mpileup\n";
				  print SH "java -jar $varscan_path mpileup2snp $out_dir/$root.$aligner.$caller.mpileup --output-vcf 1 > $out_dir/$root.$aligner.$caller.vcf\n";
				  print SH "rm $out_dir/$root.$aligner.$caller.mpileup\n";
				}
			
			  # REGULARIZE VCFs SO THAT DIFFERENT REPRESENTATIONS OF THE SAME INDEL OR COMPLEX VARIANT ARE NOT COUNTED AS DIFFERENT VARIANTS
			  print SH "$vcfallelicprimitives_path $out_dir/$root.$aligner.$caller.vcf > $out_dir/$root.$aligner.$caller.regularised.vcf\n";
			  print SH "rm $out_dir/$root.$aligner.$caller.vcf\n";
			  print SH "mv $out_dir/$root.$aligner.$caller.regularised.vcf $out_dir/$genome_id/$root.$aligner.$caller.vcf\n";
			}
		  if (($aligner ne 'snippy') and ($aligner ne 'spandx') and ($aligner ne 'speedseq'))
			{ print SH "rm $out_dir/$root.$aligner.bam $out_dir/$root.$aligner.bai\n"; }
  		  if ($aligner eq 'speedseq')
			{ print SH "rm $out_dir/$root.$aligner.bam $out_dir/$root.$aligner.bam.bai $out_dir/$root.$aligner.discordants.bam $out_dir/$root.$aligner.discordants.bam.bai $out_dir/$root.$aligner.splitters.bam $out_dir/$root.$aligner.splitters.bam.bai\n"; }
		}
	  if (exists($callers{'deepvariant'}))
		{ print SH "rm -r $out_dir/deepvariant_output\n";
		  print SH "rm -r $out_dir/deepvariant_input\n";
		}
	}
close(SH) or die $!;
exit 1;