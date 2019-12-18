=head
This script is provided in support of Bush, et al. Genomic diversity affects the accuracy of bacterial SNP calling pipelines. https://www.biorxiv.org/content/10.1101/653774v1
It is provided for reproducibility purposes only. No attempt has been made to optimise the code.

Prerequisites:
1. Download the data archive from: https://ora.ox.ac.uk/objects/uuid:8f902497-955e-4b84-9b85-693ee0e4433e. Its contents are referred to here as $rep_dir.
2. This archive also contains information regarding third-party tools, in the file "Installation_instructions.txt". This script will test for the presence of some of these tools.

Fuller details are available in the readme provided with the archive.
=cut

use strict;
use warnings;

# REQUIREMENTS
my $root    = '/home/ndm.local/steveb';
my $path    = "$root/varcall_evaluation";
my $progs   = "$root/programs"; # all third-party tools in "Installation_instructions.txt" are installed here
my $in_dir  = "$path/nanopore_illumina_hybrid_assemblies"; # a directory of assemblies (fasta files), one in each sub-directory. This script will populate this directory with additional files (VCFs).
my $rep_dir = "$path/representative_genomes_for_nanopore_illumina_hybrid_assemblies"; # a directory of representative genomes for each of the sequenced species, with choices informed by Kaiju, Kraken and MALDI TOF predictions

# CHECK THAT ALL DEPENDENCIES ARE AVAILABLE
my $fatal = 0; # ALL THE FOLLOWING PATHS ARE HARD-CODED PREREQUISITES. WE WILL ABORT IF THEY CANNOT BE FOUND.
if (!(-d($path)))    { $fatal++; print "error: cannot find $path\n";    }
if (!(-d($progs)))   { $fatal++; print "error: cannot find $progs\n";   }
if (!(-d($in_dir)))  { $fatal++; print "error: cannot find $in_dir\n";  }
if (!(-d($rep_dir))) { $fatal++; print "error: cannot find $rep_dir\n"; }
## PATHS TO WHOLE GENOME ALIGNERS
my $mash_path			= "$progs/mash-Linux64-v2.1/mash";
my $nucmer_dir 			= "$progs/mummer-4.0.0beta2";
my $parsnp_path			= "$progs/Parsnp-Linux64-v1.2/parsnp";
my $nucmer_path  		= "$progs/mummer-4.0.0beta2/nucmer";
my $dnadiff_path 		= "$progs/mummer-4.0.0beta2/dnadiff";
my $showsnps_path 		= "$progs/mummer-4.0.0beta2/show-snps";
my $harvesttools_path	= "$progs/harvesttools-Linux64-v1.2/harvesttools";
my $MUMmerSNPs2VCF_path = "$progs/MUMmerSNPs2VCF.py"; # from https://github.com/liangjiaoxue/PythonNGSTools
if (!(-e($mash_path)))  		 { $fatal++; print "error: cannot find $mash_path\n";  			}
if (!(-d($nucmer_dir)))  		 { $fatal++; print "error: cannot find $nucmer_dir\n";  		}
if (!(-e($nucmer_path)))  		 { $fatal++; print "error: cannot find $nucmer_path\n";  		}
if (!(-e($dnadiff_path))) 		 { $fatal++; print "error: cannot find $dnadiff_path\n"; 		}
if (!(-e($showsnps_path))) 		 { $fatal++; print "error: cannot find $showsnps_path\n"; 		}
if (!(-e($harvesttools_path))) 	 { $fatal++; print "error: cannot find $harvesttools_path\n";   }
if (!(-e($MUMmerSNPs2VCF_path))) { $fatal++; print "error: cannot find $MUMmerSNPs2VCF_path\n"; }
exit if ($fatal > 0);

# PARAMETERS
my $num_procs = 10;

# OUTPUT
my $sh_file = "$path/align_nanopore_assemblies_to_representative_genomes_and_call_variants.sh";
open(SH,'>',$sh_file) or die $!;
print SH "#!/bin/bash\n";
print SH "PATH=\$PATH:$nucmer_dir\n";

opendir(DIR,$in_dir) or die $!;
my @strain_ids = readdir(DIR);
closedir(DIR) or die $!;
my @sorted_strain_ids = sort {$a cmp $b} @strain_ids;
foreach my $strain (@sorted_strain_ids)
	{ next if (($strain eq '.') or ($strain eq '..'));
	  my $original_genome = "$in_dir/$strain/$strain.fasta";
	  
	  my $representative = ''; my $representative_genome = '';
	  if (($strain eq 'rbhstw00029') or ($strain eq 'rbhstw00053') or ($strain eq 'rbhstw00127') or ($strain eq 'rbhstw00350'))
		{ $representative = 'Citrobacter_freundii'; $representative_genome = "$rep_dir/Citrobacter_freundii.fa"; }
	  elsif (($strain eq 'rbhstw00059') or ($strain eq 'rbhstw00131') or ($strain eq 'rbhstw00340'))
		{ $representative = 'Enterobacter_cloacae'; $representative_genome = "$rep_dir/Enterobacter_cloacae.fa"; }
	  elsif ($strain eq 'rbhstw00309')
		{ $representative = 'Enterobacter_kobei'; $representative_genome = "$rep_dir/Enterobacter_kobei.fa"; }
	  elsif (($strain eq 'cft073') or ($strain eq 'rbhstw00122') or ($strain eq 'rbhstw00277') or ($strain eq 'rhb10c07') or ($strain eq 'rhb11c04'))
		{ $representative = 'Escherichia_coli'; $representative_genome = "$rep_dir/Escherichia_coli.fa"; }
	  elsif (($strain eq 'rbhstw00167') or ($strain eq 'rbhstw00189'))
		{ $representative = 'Klebsiella_oxytoca'; $representative_genome = "$rep_dir/Klebsiella_oxytoca.fa"; }
	  elsif (($strain eq 'mgh78578') or ($strain eq 'rbhstw00128') or ($strain eq 'rhb14c01'))
		{ $representative = 'Klebsiella_pneumoniae'; $representative_genome = "$rep_dir/Klebsiella_pneumoniae.fa"; }
	  
	  if (!(-e($original_genome))) 		 { print "WARNING: cannot find $original_genome\n"; 	  }
	  if (!(-e($representative_genome))) { print "WARNING: cannot find $representative_genome\n"; }
	  next if ( (!(-e($original_genome))) or (!(-e($representative_genome))) );
	  
	  # CREATE A VERSION OF THE REPRESENTATIVE GENOME CONTAINING ONLY THE CORE SEQUENCE, I.E. THE "COMPLETE GENOME", WITH NO PLASMIDS
	  ## we're doing this to avoid downstream complications with ParSnp (detailed at, e.g., https://github.com/marbl/parsnp/issues/16 and https://github.com/marbl/parsnp/issues/33)
	  my $header = ''; my %seqs = ();
	  open(IN,$representative_genome) or die $!;
	  while(<IN>)
		{ my $line = $_; chomp($line);
		  if ($line =~ /^\>(.*?)$/)
			{ $header = $1;
			}
		  next if ($line =~ /^\>.*?$/);
		  $seqs{$header} .= $line;
		}
	  close(IN) or die $!;
	  my $representative_core_genome = "$rep_dir/$representative.core.fa";
	  if (!(-e($representative_core_genome)))
		{ open(OUT,'>',$representative_core_genome) or die $!;
		  while((my $header,my $seq)=each(%seqs))
			{ if ($header =~ /complete genome$/)
				{ print OUT ">$header\n$seq\n";
				}
			}
		  close(OUT) or die $!;
		}
	  
	  # WHAT IS THE MASH DISTANCE BETWEEN THE ORIGINAL AND REPRESENTATIVE GENOME?
	  print SH "cd $rep_dir\n";
	  print SH "$mash_path sketch -p $num_procs -o $strain $original_genome\n" 				 unless (-e("$rep_dir/$strain.msh"));
	  print SH "$mash_path sketch -p $num_procs -o $representative $representative_genome\n" unless (-e("$rep_dir/$representative.msh"));
	  print SH "$mash_path dist -p $num_procs $rep_dir/$strain.msh $rep_dir/$representative.msh > $rep_dir/$strain.mash-dist_output.tsv\n" unless (-e("$rep_dir/$strain.mash-dist_output.tsv"));

	  # ALIGN THE ORIGINAL GENOME (NANOPORE/ILLUMINA HYBRID ASSEMBLY) TO THE REPRESENTATIVE GENOME OF THAT SPECIES, SO AS TO KNOW THE LOCATION OF ALL SNPs RELATIVE TO THIS REFERENCE. THESE ARE WHAT WE WILL LOOK FOR WHEN ALIGNING THE ORIGINAL ILLUMINA READS (USED TO MAKE EACH NANOPORE/ILLUMINA HYBRID ASSEMBLY) TO THE REPRESENTATIVE GENOME OF THAT SPECIES.
	  
	  print SH "cd $in_dir/$strain\n";

	  ### NUCMER ###
	  if (!(-e("$in_dir/$strain/$strain.varcall_relative_to_representative.nucmer.vcf")))
		{ print SH "$nucmer_path -p $strain $representative_genome $original_genome\n";
		  print SH "$dnadiff_path -p $strain -d $in_dir/$strain/$strain.delta\n";
		  print SH "$showsnps_path -Clr -x 1 -T $in_dir/$strain/$strain.delta > $in_dir/$strain/$strain.snps.filter\n";
		  print SH "rm $in_dir/$strain/$strain.1coords $in_dir/$strain/$strain.delta $in_dir/$strain/$strain.1delta $in_dir/$strain/$strain.mcoords $in_dir/$strain/$strain.mdelta $in_dir/$strain/$strain.qdiff $in_dir/$strain/$strain.rdiff\n";
		  print SH "rm $in_dir/$strain/$strain.snps $in_dir/$strain/$strain.report\n";
		  print SH "rm $in_dir/$strain/$strain.unref $in_dir/$strain/$strain.unqry\n"; # NOTE: these files do not always get created by $dnadiff_path, so this step may report a warning ("unable to delete non-existent file...")
		  print SH "python $MUMmerSNPs2VCF_path $in_dir/$strain/$strain.snps.filter $in_dir/$strain/$strain.varcall_relative_to_representative.nucmer.vcf\n";
		  print SH "rm $in_dir/$strain/$strain.snps.filter\n";
		}
		
	  ### PARSNP ###
	  if (!(-e("$in_dir/$strain/$strain.varcall_relative_to_representative.parsnp.vcf")))
		{ print SH "mkdir $in_dir/$strain/genome\n";
		  print SH "cp $original_genome $in_dir/$strain/genome/$strain.fasta\n";
		  print SH "$parsnp_path -r $representative_core_genome -d $in_dir/$strain/genome -o $in_dir/$strain/parsnp\n";
		  print SH "$harvesttools_path -i $in_dir/$strain/parsnp/parsnp.ggr -V $in_dir/$strain/$strain.varcall_relative_to_representative.parsnp.vcf\n";
		  print SH "rm -r $in_dir/$strain/genome\n";
		  print SH "rm -r $in_dir/$strain/parsnp\n";
		}
	}
close(SH) or die $!;
exit 1;