# poly-div_sfs.pl
Estimates polymorphism and divergence parameters from FASTA alignments

This Perl script can estimate polymorphism and divergence summary statistics from a fasta alignment. It will separate synonymous and nonsynonymous sites if the alignment is for coding sequences. You can control for reading frames with the -r argument, and provide one or multiple outgroups with the -o argument. In sites with more than two alleles, the two most frequent ones are considered, while the least frequent one is ignored. All estimates are based on site frequency spectra (SFS). For an example see Zeng et al. (2006) and Campos et al. (2012). It will also print SFS for all, synonymous, and nonsynonymous sites if the -S 1/s/n argument is given. If an outgroup is provided the SFS will be unfolded. The most recent version also outputs estimates for 4-fold sites, and calculates the codon adaptation index (CAI).

Start by downloading the script:

    
