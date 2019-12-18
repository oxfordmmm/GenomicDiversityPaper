# Genomic Diversity Paper Scripts
The Scripts in the repo is for "Genomic diversity affects the accuracy of bacterial SNP calling pipelines"

This archive contains the Perl scripts necessary to reproduce the results of the SNP calling evaluation available as:

Bush, et al. Genomic diversity affects the accuracy of bacterial SNP calling pipelines. https://www.biorxiv.org/content/10.1101/653774v1

Perl scripts are to be used in numerical order. The purpose of, and prerequisites for, each script are detailed in the header of each and briefly recapitulated here. See also in-code comments for further information.
All scripts are provided as is, for reproducibility purposes only. No attempt has been made to optimise the code.

All scripts rely on supporting data available in the archive: https://ora.ox.ac.uk/objects/uuid:8f902497-955e-4b84-9b85-693ee0e4433e. This archive also contains the fuller set of supporting documentation, although a brief functional description of each script follows:

Script 1 is used to align each Nanopore/Illumina hybrid assembly (in "nanopore_illumina_hybrid_assemblies") to its representative genome (in "representative_genomes_for_nanopore_illumina_hybrid_assemblies"), and then create the nucmer & ParSnp SNP calls.

Script 2 generates a .sh that contains the bash commands necessary to run every pairwise combination of read aligner and variant caller, alongside 'all-in-one' pipelines. These bash commands then populate an output directory ("vcfs").

Script 3 filters the contents of the "vcfs" directory, produced by script 2, in order to populate the directory "filtered_vcfs".

Script 4 parses the contents of the "filtered_vcfs" directory to generate two summary files, summary_statistics.tsv and pipeline_rankings.tsv, within the output directory "evaluationOutput" or "evaluationOutput-repeatmasked" (a parameter in this script determines which is created).

