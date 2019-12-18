=head
This script is provided in support of Bush, et al. Genomic diversity affects the accuracy of bacterial SNP calling pipelines. https://www.biorxiv.org/content/10.1101/653774v1
It is provided for reproducibility purposes only. No attempt has been made to optimise the code.

Prerequisites:
1. Download the data archive from: https://ora.ox.ac.uk/objects/uuid:8f902497-955e-4b84-9b85-693ee0e4433e.
2. This archive also contains information regarding third-party tools, in the file "Installation_instructions.txt". This script will test for the presence of some of these tools.

Fuller details are available in the readme provided with the archive.
=cut

use strict;
use warnings;

# REQUIREMENTS
my $root   		  = '/home/ndm.local/steveb';
my $path   		  = "$root/varcall_evaluation";
my $progs  		  = "$root/programs";
my $in_dir 		  = "$path/vcfs"; # from 2.run_SNP_calling_pipelines.pl
my $bgzip_path    = "$progs/bin/bgzip";
my $tabix_path    = "$progs/bin/tabix"; # tabix and bgzip are components of the SAMtools package, commonly used to compress and index VCFs; see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3042176/
my $vcfsort_path  = "$progs/bin/vcf-sort"; # a component of VCFtools; see https://vcftools.github.io/index.html
my $bcftools_path = "$progs/bcftools/bcftools";

# CHECK THAT ALL DEPENDENCIES ARE AVAILABLE
my $fatal = 0; # ALL THE FOLLOWING PATHS ARE HARD-CODED PREREQUISITES. WE WILL ABORT IF THEY CANNOT BE FOUND.
if (!(-d($path)))     	   { $fatal++; print "ERROR: cannot find $path\n";    		}
if (!(-d($progs)))    	   { $fatal++; print "ERROR: cannot find $progs\n";   		}
if (!(-d($in_dir)))   	   { $fatal++; print "ERROR: cannot find $in_dir\n";  		}
if (!(-e($bgzip_path)))    { $fatal++; print "ERROR: cannot find $bgzip_path\n";    }
if (!(-e($tabix_path)))    { $fatal++; print "ERROR: cannot find $tabix_path\n";    }
if (!(-e($vcfsort_path)))  { $fatal++; print "ERROR: cannot find $vcfsort_path\n";  }
if (!(-e($bcftools_path))) { $fatal++; print "ERROR: cannot find $bcftools_path\n"; }
exit 1 if ($fatal > 0);

# PARAMETERS
my @aligners = (qw/bbmap bowtie2 bwa-mem bwa-sw bwa+stampy cushaw3 gassst gem hisat2 minimap2 mosaik ngm novoalign smalt snap snippy spandx speedseq stampy yara/);
my @callers  = (qw/16GT deepvariant freebayes gatk lofreq mpileup octopus pilon platypus snippy snver snvsniffer solsnp spandx speedseq strelka varscan/);
my %aligners = map {$_ => 1} @aligners;
my %callers  = map {$_ => 1} @callers;

# OUTPUT
my $out_dir = "$path/filtered_vcfs";
if (!(-d($out_dir))) { mkdir $out_dir or die $!; }
my $sh_file = "$path/filter_VCFs.sh";
open(SH,'>',$sh_file) or die $!;
print SH "#!/bin/bash\n";

opendir(DIR,$in_dir) or die $!;
my @strain_ids = readdir(DIR);
closedir(DIR) or die $!;
my $strain_ids_seen = 0; my $strain_ids_total = @strain_ids; $strain_ids_total = $strain_ids_total-2;
foreach my $strain_id (@strain_ids)
	{ next if (($strain_id eq '.') or ($strain_id eq '..'));
	  $strain_ids_seen++;
	  next if (!(-d("$in_dir/$strain_id")));
	  opendir(DIR,"$in_dir/$strain_id") or die $!;
	  my @files = readdir(DIR);
	  closedir(DIR) or die $!;
	  my $files_seen = 0; my $files_total = @files; $files_total = $files_total-2;
	  my @sorted_files = sort {$a cmp $b} @files;
	  foreach my $file (@sorted_files)
		{ next if (($file eq '.') or ($file eq '..'));
		  $files_seen++;
		  print "reading VCFs for $strain_id (strain $strain_ids_seen of $strain_ids_total): $files_seen of $files_total\n";
		  if ($file =~ /^(.*?)\.(.*?)\.(.*?)\.vcf$/)
			{ my $strain_id = $1; my $aligner = $2; my $caller = $3;
			  my $file_name = "$strain_id.$aligner.$caller";
			  next if ( (!(exists($aligners{$aligner}))) || (!(exists($callers{$caller}))) );
			  next if (-e("$out_dir/$strain_id/$file")); # CHECKPOINT: we have seen this already
			  if (!(-d("$out_dir/$strain_id"))) { mkdir "$out_dir/$strain_id" or die $!; }
			  if ($caller eq 'solsnp')
				{ print SH "cat $in_dir/$strain_id/$file | awk '{ sub(/ID\=GQ\,Number\=1\,Type\=Integer/,\"ID=GQ,Number=1,Type=Float\"); print }' > $out_dir/$strain_id/$file_name.original.vcf\n"; }
			  else
				{ print SH "cp $in_dir/$strain_id/$file $out_dir/$strain_id/$file_name.original.vcf\n"; }
			  print SH "$vcfsort_path $out_dir/$strain_id/$file_name.original.vcf > $out_dir/$strain_id/$file_name.original.sorted.vcf\n";
			  print SH "mv $out_dir/$strain_id/$file_name.original.sorted.vcf $out_dir/$strain_id/$file_name.original.vcf\n";
			  print SH "$bgzip_path $out_dir/$strain_id/$file_name.original.vcf\n";
			  print SH "$tabix_path -p vcf $out_dir/$strain_id/$file_name.original.vcf.gz\n";
			  if    ($caller eq '16GT')		   { print SH "$bcftools_path filter -S . -s Q20 -e '%QUAL<20' -Ov $out_dir/$strain_id/$file_name.original.vcf.gz | $bcftools_path filter -S . -s OneEachWay -e '(DP8[0] == 0 && DP8[2] == 0 && DP8[4] == 0 && DP8[6] == 0) || (DP8[1] == 0 && DP8[3] == 0 && DP8[5] == 0 && DP8[7] == 0)' -m+ -Ov | $bcftools_path filter -S . -s Consensus75 -e 'FRAC<0.75' -m+ -Ov | $bcftools_path filter -S . -s HQDepth5 -e 'DP<=5' -m+ -Ov -o $out_dir/$strain_id/$file\n"; }
			  elsif ($caller eq 'deepvariant') { print SH "$bcftools_path filter -S . -s Q20 -e '%QUAL<20' -Ov $out_dir/$strain_id/$file_name.original.vcf.gz | $bcftools_path filter -S . -s Consensus75 -e 'VAF<0.75' -m+ -Ov | $bcftools_path filter -S . -s HQDepth5 -e 'DP<=5' -m+ -Ov -o $out_dir/$strain_id/$file\n"; }
			  elsif ($caller eq 'freebayes')   { print SH "$bcftools_path filter -S . -s Q20 -e '%QUAL<20' -Ov $out_dir/$strain_id/$file_name.original.vcf.gz | $bcftools_path filter -S . -s OneEachWay -e 'SAF == 0 || SAR == 0' -m+ -Ov | $bcftools_path filter -S . -s Consensus75 -e '((SAF+SAR)/(SRF+SRR))<3' -m+ -Ov | $bcftools_path filter -S . -s HQDepth5 -e 'DP<=5' -m+ -Ov -o $out_dir/$strain_id/$file\n"; }
			  elsif ($caller eq 'gatk')		   { print SH "$bcftools_path filter -S . -s Q20 -e '%QUAL<20' -Ov $out_dir/$strain_id/$file_name.original.vcf.gz | $bcftools_path filter -S . -s HQDepth5 -e 'DP<=5' -m+ -Ov -o $out_dir/$strain_id/$file\n"; }
			  elsif ($caller eq 'lofreq')	   { print SH "$bcftools_path filter -S . -s Q20 -e '%QUAL<20' -Ov $out_dir/$strain_id/$file_name.original.vcf.gz | $bcftools_path filter -S . -s OneEachWay -e 'DP4[2] == 0 || DP4[3] == 0' -m+ -Ov | $bcftools_path filter -S . -s Consensus75 -e '((DP4[2]+DP4[3])/(DP4[0]+DP4[1]))<3' -m+ -Ov | $bcftools_path filter -S . -s HQDepth5 -e '(DP4[2]+DP4[3])<=5' -m+ -Ov -o $out_dir/$strain_id/$file\n"; }
			  elsif ($caller eq 'mpileup')     { print SH "$bcftools_path filter -S . -s Q20 -e '%QUAL<20' -Ov $out_dir/$strain_id/$file_name.original.vcf.gz | $bcftools_path filter -S . -s OneEachWay -e 'DP4[2] == 0 || DP4[3] == 0' -m+ -Ov | $bcftools_path filter -S . -s Consensus75 -e '((DP4[2]+DP4[3])/(DP4[0]+DP4[1]))<3' -m+ -Ov | $bcftools_path filter -S . -s HQDepth5 -e '(DP4[2]+DP4[3])<=5' -m+ -Ov -o $out_dir/$strain_id/$file\n"; }
			  elsif ($caller eq 'octopus')	   { print SH "$bcftools_path filter -S . -s Q20 -e '%QUAL<20' -Ov $out_dir/$strain_id/$file_name.original.vcf.gz | $bcftools_path filter -S . -s HQDepth5 -e 'DP<=5' -m+ -Ov -o $out_dir/$strain_id/$file\n"; }
			  elsif ($caller eq 'pilon')	   { print SH "$bcftools_path filter -S . -s Q20 -e '%QUAL<20' -Ov $out_dir/$strain_id/$file_name.original.vcf.gz | $bcftools_path filter -S . -s Consensus75 -e 'AF<0.75' -m+ -Ov | $bcftools_path filter -S . -s HQDepth5 -e 'DP<=5' -m+ -Ov -o $out_dir/$strain_id/$file\n"; }
			  elsif ($caller eq 'platypus')	   { print SH "$bcftools_path filter -S . -s Q20 -e '%QUAL<20' -Ov $out_dir/$strain_id/$file_name.original.vcf.gz | $bcftools_path filter -S . -s OneEachWay -e 'NF == 0 || NR == 0' -m+ -Ov | $bcftools_path filter -S . -s Consensus75 -e '(TR/(TC-TR))<3' -m+ -Ov | $bcftools_path filter -S . -s HQDepth5 -e 'TR<=5' -m+ -Ov -o $out_dir/$strain_id/$file\n"; }
			  elsif ($caller eq 'snippy')	   { print SH "$bcftools_path filter -S . -s Q20 -e '%QUAL<20' -Ov $out_dir/$strain_id/$file_name.original.vcf.gz | $bcftools_path filter -S . -s Consensus75 -e '(AO/RO)<3' -m+ -Ov | $bcftools_path filter -S . -s HQDepth5 -e 'DP<=5' -m+ -Ov -o $out_dir/$strain_id/$file\n"; }
			  elsif ($caller eq 'snver')	   { print SH "$bcftools_path filter -S . -s OneEachWay -e 'AC1 == 0 || AC2 == 0' -Ov $out_dir/$strain_id/$file_name.original.vcf.gz | $bcftools_path filter -S . -s Consensus75 -e '((AC1+AC2)/(RC1+RC2))<3' -m+ -Ov | $bcftools_path filter -S . -s HQDepth5 -e 'DP<=5' -m+ -Ov -o $out_dir/$strain_id/$file\n"; }
			  elsif ($caller eq 'solsnp')	   { print SH "$bcftools_path filter -S . -s Q20 -e '%QUAL<20' -Ov $out_dir/$strain_id/$file_name.original.vcf.gz -o $out_dir/$strain_id/$file\n"; } # note that there are errors in SolSNP VCFs such that vcfallelicprimitives reports: (1) "Invalid character '.' in 'GQ' FORMAT field" (the value is a float but the FORMAT field defines it as an integer), AND (2) "INFO 'AR' is not defined in the header, assuming Type=String". The latter triggers a warning, but the former an error. This would result in empty output unless corrected (as above).
			  elsif ($caller eq 'spandx')	   { print SH "$bcftools_path filter -S . -s Q20 -e '%QUAL<20' -Ov $out_dir/$strain_id/$file_name.original.vcf.gz | $bcftools_path filter -S . -s HQDepth5 -e 'DP<=5' -m+ -Ov -o $out_dir/$strain_id/$file\n"; }
			  elsif ($caller eq 'speedseq')	   { print SH "$bcftools_path filter -S . -s Q20 -e '%QUAL<20' -Ov $out_dir/$strain_id/$file_name.original.vcf.gz | $bcftools_path filter -S . -s OneEachWay -e 'SAF == 0 || SAR == 0' -m+ -Ov | $bcftools_path filter -S . -s Consensus75 -e '(AO/RO)<3' -m+ -Ov | $bcftools_path filter -S . -s HQDepth5 -e 'DP<=5' -m+ -Ov -o $out_dir/$strain_id/$file\n"; }
			  elsif ($caller eq 'strelka')	   { print SH "$bcftools_path filter -S . -s Q20 -e '%QUAL<20' -Ov $out_dir/$strain_id/$file_name.original.vcf.gz | $bcftools_path filter -S . -s OneEachWay -e 'ADF[*:1] == 0 || ADR[*:1] == 0' -m+ -Ov | $bcftools_path filter -S . -s Consensus75 -e '((ADF[*:1]+ADR[*:1])/(ADF[*:0]+ADR[*:0]))<3' -m+ -Ov | $bcftools_path filter -S . -s HQDepth5 -e 'DP<=5' -m+ -Ov -o $out_dir/$strain_id/$file\n"; }
			  elsif ($caller eq 'snvsniffer')  { print SH "$bcftools_path filter -S . -s HQDepth5 -e 'DP<=5' -Ov $out_dir/$strain_id/$file_name.original.vcf.gz -o $out_dir/$strain_id/$file\n"; }
			  elsif ($caller eq 'varscan')	   { print SH "$bcftools_path filter -S . -s OneEachWay -e 'ADF == 0 || ADR == 0' -Ov $out_dir/$strain_id/$file_name.original.vcf.gz | $bcftools_path filter -S . -s Consensus75 -e '((ADF+ADR)/(RDF+RDR))<3' -m+ -Ov | $bcftools_path filter -S . -s HQDepth5 -e 'DP<=5' -m+ -Ov -o $out_dir/$strain_id/$file\n"; }
			  print SH "rm $out_dir/$strain_id/$file_name.original.vcf.gz $out_dir/$strain_id/$file_name.original.vcf.gz.tbi\n";
			}
		}
	}

close(SH) or die $!;
exit 1;