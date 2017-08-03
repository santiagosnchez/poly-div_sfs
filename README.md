# poly-div_sfs.pl
Estimates polymorphism and divergence parameters from FASTA alignments

This Perl script can estimate polymorphism and divergence summary statistics from a fasta alignment. It will separate synonymous and nonsynonymous sites if the alignment is for coding sequences. You can control for reading frames with the -r argument, and provide one or multiple outgroups with the -o argument. In sites with more than two alleles, the two most frequent ones are considered, while the least frequent one is ignored. All estimates are based on site frequency spectra (SFS). For an example see Ronen et al. (2013) and Campos et al. (2014). It will also print SFS for all, synonymous, and nonsynonymous sites if the -S 1/s/n argument is given. If an outgroup is provided the SFS will be unfolded. The most recent version also outputs estimates for 4-fold sites, and calculates the codon adaptation index (CAI).

## Installation

Start by downloading the repository:

    git clone https://github.com/santiagosnchez/poly-div_sfs
    
You can also just download the script with:

    wget -Nq https://raw.githubusercontent.com/santiagosnchez/poly-div_sfs.pl/master/poly%2Bdiv_sfs.pl
    
Go to the `poly-div_sfs` dirctory and run the script with the `-h` (help) flag:

    cd poly-div_sfs
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

If you wish to have the script accessible from any location, you can do:

    chmod +x poly+div_sfs.pl
    sudo cp poly+div_sfs.pl /usr/local/bin

## Running the program

The example file includes protein-coding sequences (e.g. no introns) for 7 *Amanita* species in the *Amanita jacksonii* complex (Sánchez-Ramírez et al. 2015). The sampling includes multiple individuals per species.

To get **diversity** estimates simply speficy the file and the name tag of a species that matches **all** records of the headers in that group.

    perl poly+div_sfs.pl -f example_Amanita_scaffold-0.g11.fasta -s jacksonii

Results in:

    example_Amanita_scaffold-0.g11.fasta,jacksonii,44,137,0.00974,0.00269,-2.63807,0.03774

As you can see, the results are in CSV format, and includes the name of the file, the species, the number of haplotyes, etc.

To see, what each field means, run with the `-v` verbose argument.

    perl poly+div_sfs.pl -f example_WholeGene_Amanita_scaffold-0.g11.fasta -s jacksonii -v
    # Reading frame is 0 + 1
    Gene,Pop,N,S,Theta,Pi,TajimasD,Psing
    example_WholeGene_Amanita_scaffold-0.g11.fasta,jacksonii,44,137,0.00974,0.00269,-2.63807,0.03774
    
Now you have more information and headers for each field.

### Multiple species/populations

If you want estimates for different populations or species, you simply need to provide a list of the species you want in the `-s` argument.

    perl poly+div_sfs.pl -f example_WholeGene_Amanita_scaffold-0.g11.fasta -s jacksonii:sp_F11:sp_T31:sp_jack2 -v

### Divergence

By providing an outgroup `-o` you can estimate **divergence**-based statistics. By default the Kimura-2-parameter distance will be calculated.

    perl poly+div_sfs.pl -f example_WholeGene_Amanita_scaffold-0.g11.fasta -s jacksonii:sp_jack2 -o sp_F11 -v

If multiple sequences match the outgroup tag, either the sequence with least missing data or one at random will be picked.

### Site-class (synonymous/nonsynonymous) separation

Simply add the `-c` argument.

    perl poly+div_sfs.pl -f example_Coding_Amanita_scaffold-0.g11.fasta -s jacksonii:sp_F11:sp_T31:sp_jack2 -v -c

### Divergence at silent and replacement sites

For **divergence** estimates by site-class separation, specify an outgroup to the `-o` argument together with the `-c` argument.

    perl poly+div_sfs.pl -f example_WholeGene_Amanita_scaffold-0.g11.fasta -s jacksonii:sp_jack2 -o sp_F11  -v -c

    

