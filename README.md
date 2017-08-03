# poly-div_sfs.pl
Estimates polymorphism and divergence parameters from FASTA alignments

This Perl script can estimate polymorphism and divergence summary statistics from a fasta alignment. It will separate synonymous and nonsynonymous sites if the alignment is for coding sequences. You can control for reading frames with the -r argument, and provide one or multiple outgroups with the -o argument. In sites with more than two alleles, the two most frequent ones are considered, while the least frequent one is ignored. All estimates are based on site frequency spectra (SFS). For an example see Ronen et al. (2013) and Campos et al. (2014). It will also print SFS for all, synonymous, and nonsynonymous sites if the -S 1/s/n argument is given. If an outgroup is provided the SFS will be unfolded. The most recent version also outputs estimates for 4-fold sites, and calculates the codon adaptation index (CAI).

# Installation

Start by downloading the script:

    wget https://raw.githubusercontent.com/santiagosnchez/poly-div_sfs.pl/master/poly%2Bdiv_sfs.pl
    
Then run it with the `-h` (help) flag:

    perl poly+div_sfs.pl -h

If you encounter this error:

    Can't locate Statistics/ChisqIndep.pm in @INC (you may need to install the Statistics::ChisqIndep module) (@INC contains: /usr/local/Cellar/perl/5.24.1/lib/perl5/site_perl/5.24.1/darwin-thread-multi-2level /usr/local/Cellar/perl/5.24.1/lib/perl5/site_perl/5.24.1 /usr/local/Cellar/perl/5.24.1/lib/perl5/5.24.1/darwin-thread-multi-2level /usr/local/Cellar/perl/5.24.1/lib/perl5/5.24.1 /usr/local/lib/perl5/site_perl/5.24.1 .) at poly+div_sfs.pl line 18.
    BEGIN failed--compilation aborted at poly+div_sfs.pl line 18.

You will need to install the `Statistics::ChisqIndep.pm` module. You can easily do this uing `cpanm`. On Mac OSX you can use `brew` to install `cpanm`:

    brew install cpanm
    
On Ubuntu/Linux use `apt-get`:

    sudo apt-get install cpanm

Once `cpanm` is installed install the module:

    sudo cpanm Statistics::ChisqIndep

Now rerun `poly+div_sfs.pl`:

You should see:

    perl poly+div_sfs.pl -h
    
    Try:
    perl poly+div.pl -f your_fasta_file  [required]   FASTA file alignment
                     -s sp1[:sp2:sp3]    [optional]   Header tags that can distinguish between populations/species
                     -e id1[:id2:id3]    [optional]   Exclude sequences based on header tags
                     -c                  [optional]	  Takes the alignment as inframe coding sequences
                                                      producing estimates for synonymous and nonsynonymous sites
                     -o spX[:spY:spZ]    [optional]   Invokes MKT and specifies outgroup(s) of each species in the -s array labels.
                                                      If multiple sequences are found it will take only one sequence; the one with less missing data.
                     -v                  [optional]   Verbose mode
                     -t 3[0.5]           [optional]   Sample size threshold for segregating sites. Could be specified through
                                                      an integer or a fraction of N. Will exclude sites with insufficient data.
                     -F 1[0.04]          [optional]   Will exclude sites with minor allele frequency specified with am integer or a 
                                                      fraction of N.
                     -r 1[2]             [optional]   Start reading frame in a position other than 0 (or 1st frame). The default is 0.
                     -S 0[1:s:n]         [optional]   Only print the site frequency spectrum. If -o is given the SFS will be unfolded.
                                                      If used in combination with -c, it will print SFS for synonymous(s) or replacement(n) sites.


# Running the program



