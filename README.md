# SemHelpers
A repo of perl scripts used to improve genome assemblies (tuned for Illumina Synthetic Long Reads), born from the work on assembling allotetraploid species Trifolium repens (White Clover).

correct_small_gaps_v2.pl: This script will take a FASTA file (reference sequence) and an mpileup file, and will do indel correction on said reference according to the alignment (goes by consensus, it is tuned to work with Illumina Synthetic Long Reads).

merge_and_replace.pl: Alternatively, this script will take a MAF alignment (such as one done with last or lastZ) and will replace the reference genome with the aligned reads based on similarity (user-defined).

maf_masking_by_feature.pl: This script will take a MAF alignment and a GFF file (or more), and mask the alignment where those features are located. Useful for, for example, align whole genomes taking advantage of synteny, and then mask exons or genes to do mutation rate estimates.

sort_scaffolds_by_LD.pl
place_scaffolds_list_simplified.pl : These two scripts are merely aids to create chained scaffolds/pseudomolecules based on positions to a close relative reference sequence. The lists can be created via LD mapping, simple Megablast, or a combination or both.
