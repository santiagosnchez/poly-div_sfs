#!/usr/bin/perl
# ©Santiago Sanchez-Ramirez, University of Toronto
# poly+div.pl
#
# If you get this error right away: "Can't locate Statistics/ChisqIndep.pm in @INC",
# try installing CPAN and installing the modules Statistics::Distributions and Statistics::ChisqIndep
# on Mac OS X you can do by simply typing on the command line:
# $ sudo cpan App::cpanminus
# Add you password and follow the instructions as best as you can. I would recommend going for the "local::lib" approach.
# After cpanm is installed, you can install any module easily. So do:
# $ sudo cpanm Statistics::Distributions
# $ sudo cpanm Statistics::ChisqIndep

#use warnings;
#use strict;
use List::MoreUtils qw(any all uniq);
use POSIX qw(ceil floor);
use Statistics::ChisqIndep;
my @lines=();
my $file;
my $gene;
my @spp=();
my @exclude=();
my @out=();
my $v=0;
my $c=0;
my $m=0;
my $SFS=0;
my $thresh=0;
my $frac=0;
my $start_frame=0;

if (grep { /^-he{0,1}l{0,1}p{0,1}$/ } @ARGV){
	die "
Try:
perl poly+div_sfs.pl -h

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
                 -F 1[0.04]          [optional]   Will exclude sites with minor allele frequency specified with an integer or a 
                                                  fraction of N.
                 -r 1[2]             [optional]   Start reading frame in a position other than 0 (or 1st frame). The default is 0.
                 -S 0[1:s:n]         [optional]   Only print the site frequency spectrum. If -o is given the SFS will be unfolded.
                                                  If used in combination with -c, it will print SFS for synonymous(s) or replacement(n) sites. 
\n";
}

if (grep { /^-v$/ } @ARGV){
	$v = 1;
}

if (my ($indSFS) = grep { $ARGV[$_] =~ /^-S$/ } 0 .. $#ARGV){
	$SFS = $ARGV[$indSFS+1];
	$v = 0;
}

if (grep { /^-c$/ } @ARGV){
	$c = 1;
	if ($v == 1){
		print "# Estimating values for synonymous and nonsynonymous sites\n";
	}
}

if (my ($indF) = grep { $ARGV[$_] =~ /^-f$/ } 0 .. $#ARGV){
	$file = $ARGV[$indF+1];
	$gene = $ARGV[$indF+1];
	$gene =~ s/.*\///g;
	open(F,"<",$ARGV[$indF+1]) or die "Cannot open $file.\n";
	while(<F>){
		next if (/^\s*$/);
		chomp($_);
		$_ =~ s/\s//g;
		push @lines, $_;
	}
} else {
	die "-f flag not found. Try with: -f your_file\n";
}

if (my ($indS) = grep { $ARGV[$_] =~ /^-s$/ } 0 .. $#ARGV){
	if ($ARGV[$indS+1] =~ m/:/){
		@spp = split(":", $ARGV[$indS+1]);
	} else {
		push @spp, $ARGV[$indS+1];
	}
	if (scalar(@spp) == 0){
		die "The -s flag was found, but with no pattern. Try -s spp1[:spp2:spp3] [optional]\n";
	}	
} else {
	if ($v == 1){
		print "# Values are over all individuals\n";
	}
}

if (my ($indE) = grep { $ARGV[$_] =~ /^-e$/ } 0 .. $#ARGV){
	if ($ARGV[$indE+1] =~ m/:/){
		@exclude = split(":", $ARGV[$indE+1]);
	} else {
		push @exclude, $ARGV[$indE+1];
	}
	if (scalar(@exclude) == 0){
		die "The -e flag was found, but with no pattern. Try -e ind1[:ind2:ind3] [optional]\n";
	}
	if ($v == 1){
		print "# excluding records that match: ";
		print "@exclude\n";
	}
}

if (my ($indR) = grep { $ARGV[$_] =~ /^-r$/ } 0 .. $#ARGV){
	$start_frame = $ARGV[$indR+1];
	if ($v == 1){
		print "# Reading frame is $start_frame + 1\n";
	}
} else {
	if ($v == 1){
		print "# Reading frame is $start_frame + 1\n";
	}
}

if (my ($indO) = grep { $ARGV[$_] =~ /^-o$/ } 0 .. $#ARGV){
	if ($ARGV[$indO+1] =~ m/:/){
		@out = split(":", $ARGV[$indO+1]);
	} else {
		push @out, $ARGV[$indO+1];
	}
	if ($v == 1){
		print "# Estimating divergence based on outgroup sequences\n";
	}
	if (scalar(@out) == 0){
		die "The -o flag was found, but with no pattern. Try -o sppX[:sppY:sppZ] [optional]\n";
	}	
}

if (my ($indT) = grep { $ARGV[$_] =~ /^-t$/ } 0 .. $#ARGV){
	$thresh = $ARGV[$indT+1];
	if ($v == 1){
		if ($thresh =~ m/\./){
			print "# Filtering SNPs wth sample sizes below (excludes missing data): $thresh*N indiv\n";
		} else {
			print "# Filtering SNPs with sample sizes below (excludes missing data): $thresh indiv\n";
		}
	}
}

if (my ($indf) = grep { $ARGV[$_] =~ /^-F$/ } 0 .. $#ARGV){
	$frac = $ARGV[$indf+1];
	if ($v == 1){
		if ($frac =~ m/\./){
			print "# Filtering SNPs with minor allele freq below: $frac*N indiv\n";
		} else {
			print "# Filtering SNPs with minor allele freq below: $frac indiv\n";
		}
	}
}


my @indH = grep { $lines[$_] =~ /^>/ } 0 .. $#lines;
if (scalar(@indH) == 0){
	die "No data or wrong format\n";
}
my %data=();
my $sl = length(join('', @lines[ $indH[0]+1 .. $indH[0+1]-1 ]));

for my $i (0 .. $#indH){
	if ($i != $#indH){
		if (length(join('', @lines[ $indH[$i]+1 .. $indH[$i+1]-1 ])) != $sl){
			die "Sequence \"@lines[$indH[$i]]\" has a different length in file $file\n";
		}
		$data{@lines[$indH[$i]]} = join('', @lines[ $indH[$i]+1 .. $indH[$i+1]-1 ]);
		$data{@lines[$indH[$i]]} =~ tr/a-z/A-Z/;
	} else {
		if (length(join('', @lines[ $indH[$i]+1 .. $#lines ])) != $sl){
			die "Sequence \"@lines[$indH[$i]]\" has a different length in file $file\n";
		}
		$data{@lines[$indH[$i]]} = join('', @lines[ $indH[$i]+1 .. $#lines ]);
		$data{@lines[$indH[$i]]} =~ tr/a-z/A-Z/;
	}
}

if (scalar(@exclude) > 0){
	map { $x = $_; @exh = grep { /$x/ } (keys %data); for $exh (@exh){ delete $data{$exh} } } @exclude;
}

my %matchOut=();
my @ingroup=();
my @labels = sort {$a cmp $b} keys %data;
my $outone;

if (scalar(@out) != 0){
	if (scalar(@out) > 1){
		if (scalar(@spp) == scalar(@out)){
			for my $i (0 .. $#out){
				my @tmpout = grep { /$out[0]/ } @labels;
				$matchOut{$spp[$i]} = $tmpout[getbest(@tmpout)];
			}
		} else {
			die "the legth of arrays in -s and -o is not equal\n";
		}
	} 
	elsif (scalar(@out) == 1){
		if (scalar(@spp) > 1){
			my @tmpout = grep { /$out[0]/ } @labels;
			for my $i (0 .. $#spp){
				$matchOut{$spp[$i]} = $tmpout[getbest(@tmpout)];
			}
		}
		elsif (scalar(@spp) == 1){
			my @tmpout = grep { /$out[0]/ } @labels;
			$matchOut{$spp[0]} = $tmpout[getbest(@tmpout)];
		}
		else {
			my @tmpout = grep { /$out[0]/ } @labels;
			if (scalar(@tmpout) == 1){
				$outone = $tmpout[0];
			} else {
				sort {$a cmp $b} @tmpout;
				$outone = $tmpout[getbest(@tmpout)];
			}
			@ingroup = grep { !/$out[0]/ } keys %data;
			if ($v == 1){
				print "# Because -s is absent, the ingroup will be defined as N minus the outgroup\n";
			}
		}
	}
}
elsif (scalar(@out) == 0){
	@ingroup = sort {$a cmp $b} keys %data;
}

if ($c == 0){
	if (scalar(@spp) > 0){
		for $spp (@spp){
			my @pop = grep { /$spp/ } sort {$a cmp $b} keys %data;
			my $Nind = scalar(@pop);
			if ($Nind == 0){
				die "# $spp not found on FASTA labels\n";
			}
			my ($segs,$sfs_f,$sfs_u,$co_f,$co_u,$ds,$do,$di);
			if (exists $matchOut{$spp[0]}){
				($segs,$sfs_f,$sfs_u,
				$co_f,$co_u,$ds,$do,$di) = AllSegSites(\$sl,\@pop,\$Nind,\%data,\$thresh,\$frac,\$data{$matchOut{$spp}});
				if (($SFS == 1) and (scalar(@spp) == 1)){
					print "@$sfs_u\n" and exit;
				}
				my @res = polymorphism(\$sl,\$Nind,\@$sfs_f,@$sfs_u,1);
				my $dxy = divergence(\$sl,\@$di,\@$do);
				if (($v == 1) and ($spp eq $spp[0])){
					print "# Dxy is K2P corrected\n";
					print "Gene,Pop,N,S,Theta,Pi,TajimasD,Fay&Wu_H,Psing,Dxy\n";
				}
				print "$gene,$spp,$Nind," . join(',',@res) . ",$dxy" . "\n";
			} else {
				($segs,$sfs_f,$co) = AllSegSites(\$sl,\@pop,\$Nind,\%data,\$thresh,\$frac);
				if (($SFS == 1) and (scalar(@spp) == 1)){
					print "@$sfs_f\n" and exit;
				}
				my @res = polymorphism(\$sl,\$Nind,\@$sfs_f,0);
				if (($v == 1) and ($spp eq $spp[0])){
					print "Gene,Pop,N,S,Theta,Pi,TajimasD,Psing\n";
				}
				print "$gene,$spp,$Nind," . join(',',@res) . "\n";
			}
		}
	} else {
		my @pop = @ingroup;
		my $Nind = scalar(@pop);
		my ($segs,$sfs_f,$sfs_u,$co,$ds,$do,$di);
		if (exists $data{$outone}){
			($segs,$sfs_f,$sfs_u,
			$co_f,$co_u,$ds,$do,$di) = AllSegSites(\$sl,\@pop,\$Nind,\%data,\$thresh,\$frac,\$data{$outone});			
			if ($SFS == 1){
				print "@$sfs_u\n" and exit;
			}
			my @res = polymorphism(\$sl,\$Nind,\@$sfs_f,\@$sfs_u,1);
			my $dxy = divergence(\$sl,\@$di,\@$do);
			if (($v == 1) and ($spp eq $spp[0])){
				print "# Dxy is K2P corrected\n";
				print "Gene,N,S,Theta,Pi,TajimasD,Fay&Wu_H,Psing,Dxy\n";
			}
			print "$gene,$Nind," . join(',',@res) . ",$dxy" . "\n";
		} else {
			($segs,$sfs_f,$co) = AllSegSites(\$sl,\@pop,\$Nind,\%data,\$thresh,\$frac);		
			if ($SFS == 1){
				print "@$sfs_f\n" and exit;
			}
			my @res = polymorphism(\$sl,\$Nind,\@$sfs_f,0);
			if (($v == 1) and ($spp eq $spp[0])){
				print "Gene,N,S,Theta,Pi,TajimasD,Psing\n";
			}
			print "$gene,$Nind," . join(',',@res) . "\n";
		}
	}
} else {
	if (scalar(@spp) > 0){
		for $spp (@spp){
			my @pop = grep { /$spp/ } sort {$a cmp $b} keys %data;
			my $Nind = scalar(@pop);
			if ($Nind == 0){
				die "# $spp not found on FASTA labels\n";
			}
			if (exists $matchOut{$spp[0]}){
				my ($segs,$sfs_f,$sfs_u,
				$co_f,$co_u,$ds,$do,$di) = AllSegSites(\$sl,\@pop,\$Nind,\%data,\$thresh,\$frac,\$data{$matchOut{$spp}});
				my ($ss4,$ss3,$ss2,$ss0,
				$scfp4,$scfp3,$scfp2,$ncfp,
				$ds4,$ds3,$ds2,$ds0,
				$scfd4,$scfd3,$scfd2,$ncfd,$cai,$stop) = SegSitesCoding(\$sl,\@pop,\%data,\@$segs,\$data{$matchOut{$spp}},\@$ds);
				my @synp = sort {$a <=> $b} (@$ss4,@$ss3,@$ss2);
				my @synd = sort {$a <=> $b} (@$ds4,@$ds3,@$ds2);
				my $sisp = $$scfp4+$$scfp3+$$scfp2;
				my $sisd = $$scfd4+$$scfd3+$$scfd2;
				my @repco_p = @$co_f[@$ss0];
				my @synco_p = @$co_f[@synp];
				my @synco_p4f = @$co_f[@$ss4];
				my @syncou_p = @$co_u[@synp];
				my @sfsu_s = sfs(\@syncou_p,\$Nind,1);
				my @repcou_p = @$co_u[@$ss0];
				my @sfsu_n = sfs(\@repcou_p,\$Nind,1);
				if (($SFS =~ /s/i) and (scalar(@spp) == 1)){
					unshift @sfsu_s, ($sisp-scalar(@synp));
					unshift @sfsu_s, $sisp;
					push @sfsu_s, scalar(@synd);
					print "@sfsu_s\n" and exit;
				}
				elsif (($SFS =~ /n/i) and (scalar(@spp) == 1)){
					unshift @sfsu_n, ($$ncfp-scalar(@$ss0));
					unshift @sfsu_n, $$ncfp;
					push @sfsu_n, scalar(@$ds0);
					print "@sfsu_n\n" and exit;
				}
				my @sfs_s = sfs(\@synco_p,\$Nind,0);
				my @sfs_n = sfs(\@repco_p,\$Nind,0);
				my @sfs_4f = sfs(\@synco_p4f,\$Nind,0);
				my @res_s = polymorphism(\$sisp,\$Nind,\@sfsu_s,1);
				my @res_n = polymorphism(\$$ncfp,\$Nind,\@sfsu_n,1);
				my @res_4f = polymorphism(\$$scfp4,\$Nind,\@sfs_4f,0);
				my @dorep = @$do[@$ds0];
				my @direp = @$di[@$ds0];
				my @dosyn = @$do[@synd];
				my @disyn = @$di[@synd];
				my @dosyn_4f = @$do[@$ds4];
				my @disyn_4f = @$di[@$ds4];
				my $dxy_s = divergence(\$sisp,\@disyn,\@dosyn);
				my $dxy_n = divergence(\$$ncfd,\@direp,\@dorep);
				my $dxy_s4f = divergence(\$$scfd4,\@disyn_4f,\@dosyn_4f);
				my @mkt = mkt(\@res_s,\@res_n,\scalar(@synd),\scalar(@$ds0),\$dxy_s,\$dxy_n);
				my @mkt_4f = mkt(\@res_4f,\@res_n,\scalar(@$ds4),\scalar(@$ds0),\$dxy_s4f,\$dxy_n);
				my @res = ($res_s[0],$res_4f[0],$res_n[0],
					   $res_s[1],$res_4f[1],$res_n[1],
					   $res_s[2],$res_4f[2],$res_n[2],
					   $res_s[3],$res_4f[3],$res_n[3],
					   $res_s[4],$res_n[4],
					   $res_s[5],$res_4f[4],$res_n[5],
					   $dxy_s,$dxy_s4f,$dxy_n,@mkt,@mkt_4f,$$cai,$$stop);
				if (($v == 1) and ($spp eq @spp[0])){
					print "# Using K2P-corrected divergence for MKT\n";
					print "Gene,Pop,N,S_s,S_s4f,S_n,Theta_s,Theta_s4f,Theta_n,".
					      "Pi_s,Pi_s4f,Pi_n,TajimasD_s,TajimasD_s4f,TajimasD_n,Fay&Wu_H_s, Fay&Wu_H_n,Psing_s,Psing_s4f,Psing_n,".
					      "Dxy_s,Dxy_s4f,Dxy_n,NI,alpha,p-val,NI_4f,alpha_4f,p-val_4f,CAI,StopCod\n";
				}
				print "$gene,$spp,$Nind," . join(',',@res) . "\n";
			} else {
				my ($segs,$sfs_f,$co) = AllSegSites(\$sl,\@pop,\$Nind,\%data,\$thresh,\$frac);
				my ($ss4,$ss3,$ss2,$ss0,
				$scfp4,$scfp3,$scfp2,$ncfp,$cai,$stop) = SegSitesCoding(\$sl,\@pop,\%data,\@$segs);
				my @syn = sort {$a <=> $b} (@$ss4,@$ss3,@$ss2);
				my $sis = $$scfp4+$$scfp3+$$scfp2;
				my @repco = @$co[@$ss0];
				my @synco = @$co[@syn];
				my @synco_p4f = @$co[@$ss4];
				my @sfs_s = sfs(\@synco,\$Nind,0);
				my @sfs_n = sfs(\@repco,\$Nind,0);
				my @sfs_4f = sfs(\@synco_p4f,\$Nind,0);
				if (($SFS =~ /s/i) and (scalar(@spp) == 1)){
					print "@sfs_s\n" and exit;
				}
				elsif (($SFS =~ /n/i) and (scalar(@spp) == 1)){
					print "@sfs_n\n" and exit;
				}
				my @res_s = polymorphism(\$sis,\$Nind,\@sfs_s,0);
				my @res_n = polymorphism(\$$ncfp,\$Nind,\@sfs_n,0);
				my @res_4f = polymorphism(\$$scfp4,\$Nind,\@sfs_4f,0);
				my @res = ($res_s[0],$res_4f[0],$res_n[0],
					   $res_s[1],$res_4f[1],$res_n[1],
					   $res_s[2],$res_4f[2],$res_n[2],
					   $res_s[3],$res_4f[3],$res_n[3],
					   $res_s[4],$res_4f[4],$res_n[4],$$cai,$$stop);
				if (($v == 1) and ($spp eq @spp[0])){
					print "# No outgroup data found, or label not found. Only reporting polymorphism data.\n";
					print "Gene,Pop,N,S_s,S_s4f,S_n,Theta_s,Theta_s4f,Theta_n,".
				      	  "Pi_s,Pi_s4f,Pi_n,TajimasD_s,TajimasD_s4f,TajimasD_n,Psing_s,Psing_s4f,Psing_n,CAI,StopCod\n";
				}
				print "$gene,$spp,$Nind," . join(',',@res) . "\n";
			}
		}
	} else {
		my @pop = @ingroup;
		my $Nind = scalar(@pop);
		if (exists $data{$outone}){
			my ($segs,$sfs_f,$sfs_u,
			$co_f,$co_u,$ds,$do,$di) = AllSegSites(\$sl,\@pop,\$Nind,\%data,\$thresh,\$frac,\$data{$outone});
			my ($ss4,$ss3,$ss2,$ss0,
			$scfp4,$scfp3,$scfp2,$ncfp,
			$ds4,$ds3,$ds2,$ds0,
			$scfd4,$scfd3,$scfd2,$ncfd,$cai,$stop) = SegSitesCoding(\$sl,\@pop,\%data,\@$segs,\$data{$outone},\@$ds);
			my @synp = sort {$a <=> $b} (@$ss4,@$ss3,@$ss2);
			my @synd = sort {$a <=> $b} (@$ds4,@$ds3,@$ds2);
			my $sisp = $$scfp4+$$scfp3+$$scfp2;
			my $sisd = $$scfd4+$$scfd3+$$scfd2;
			my @repco_p = @$co_f[@$ss0];
			my @synco_p = @$co_f[@synp];
			my @synco_p4f = @$co_f[@$ss4];
			my @syncou_p = @$co_u[@synp];
			my @sfsu_s = sfs(\@syncou_p,\$Nind,1);
			my @repcou_p = @$co_u[@$ss0];
			my @sfsu_n = sfs(\@repcou_p,\$Nind,1);
			if ($SFS =~ /s/i){
				unshift @sfsu_s, ($sisp-scalar(@synp));
				unshift @sfsu_s, $sisp;
				push @sfsu_s, scalar(@synd);
				print "@sfsu_s\n" and exit;
			}
			elsif ($SFS =~ /n/i){
				unshift @sfsu_n, ($$ncfp-scalar(@$ss0));
				unshift @sfsu_n, $$ncfp;
				push @sfsu_n, scalar(@$ds0);
				print "@sfsu_n\n" and exit;
			}
			my @sfs_s = sfs(\@synco_p,\$Nind,0);
			my @sfs_n = sfs(\@repco_p,\$Nind,0);
			my @sfs_4f = sfs(\@synco_p4f,\$Nind,0);
			my @res_s = polymorphism(\$sisp,\$Nind,\@sfsu_s,1);
			my @res_n = polymorphism(\$$ncfp,\$Nind,\@sfsu_n,1);
			my @res_4f = polymorphism(\$$scfp4,\$Nind,\@sfs_4f,0);
			my @dorep = @$do[@$ds0];
			my @direp = @$di[@$ds0];
			my @dosyn = @$do[@synd];
			my @disyn = @$di[@synd];
			my @dosyn_4f = @$do[@$ds4];
			my @disyn_4f = @$di[@$ds4];
			my $dxy_s = divergence(\$sisp,\@disyn,\@dosyn);
			my $dxy_n = divergence(\$$ncfd,\@direp,\@dorep);
			my $dxy_s4f = divergence(\$$scfd4,\@disyn_4f,\@dosyn_4f);
			my @mkt = mkt(\@res_s,\@res_n,\scalar(@synd),\scalar(@$ds0),\$dxy_s,\$dxy_n);
			my @mkt_4f = mkt(\@res_4f,\@res_n,\scalar(@$ds4),\scalar(@$ds0),\$dxy_s4f,\$dxy_n);
			my @res = ($res_s[0],$res_4f[0],$res_n[0],
				   $res_s[1],$res_4f[1],$res_n[1],
				   $res_s[2],$res_4f[2],$res_n[2],
				   $res_s[3],$res_4f[3],$res_n[3],
				   $res_s[4],$res_n[4],
				   $res_s[5],$res_4f[4],$res_n[5],
				   $dxy_s,$dxy_s4f,$dxy_n,@mkt,@mkt_4f,$$cai,$$stop);
			if ($v == 1){
				print "# Using K2P-corrected divergence for MKT\n";
				print "Gene,N,S_s,S_s4f,S_n,Theta_s,Theta_s4f,Theta_n,".
				      "Pi_s,Pi_s4f,Pi_n,TajimasD_s,TajimasD_s4f,TajimasD_n,Fay&Wu_H_s,Fay&Wu_H_n,Psing_s,Psing_s4f,Psing_n,".
				      "Dxy_s,Dxy_s4f,Dxy_n,NI,alpha,p-val,NI_4f,alpha_4f,p-val_4f,CAI,StopCod\n";
			}
			print "$gene,$Nind," . join(',',@res) . "\n";
		} else {
			my ($segs,$sfs_f,$co) = AllSegSites(\$sl,\@pop,\$Nind,\%data,\$thresh,\$frac);
			my ($ss4,$ss3,$ss2,$ss0,
			$scfp4,$scfp3,$scfp2,$ncfp,$cai,$stop) = SegSitesCoding(\$sl,\@pop,\%data,\@$segs);
			my @syn = sort {$a <=> $b} (@$ss4,@$ss3,@$ss2);
			my $sis = $$scfp4+$$scfp3+$$scfp2;
			my @repco = @$co[@$ss0];
			my @synco = @$co[@syn];
			my @synco_p4f = @$co[@$ss4];
			my @sfs_s = sfs(\@synco,\$Nind,0);
			my @sfs_n = sfs(\@repco,\$Nind,0);
			my @sfs_4f = sfs(\@synco_p4f,\$Nind,0);
			if ($SFS =~ /s/i){
				print "@sfs_s\n" and exit;
			}
			elsif ($SFS =~ /n/i){
				print "@sfs_n\n" and exit;
			}
			my @res_s = polymorphism(\$sis,\$Nind,\@sfs_s,0);
			my @res_n = polymorphism(\$$ncfp,\$Nind,\@sfs_n,0);
			my @res_4f = polymorphism(\$$scfp4,\$Nind,\@sfs_4f,0);
			my @res = ($res_s[0],$res_4f[0],$res_n[0],
				   $res_s[1],$res_4f[1],$res_n[1],
				   $res_s[2],$res_4f[2],$res_n[2],
				   $res_s[3],$res_4f[3],$res_n[3],
				   $res_s[4],$res_4f[4],$res_n[4],$$cai,$$stop);
			if ($v == 1){
				print "# No outgroup data found, or label not found. Only reporting polymorphism data.\n";
				print "Gene,N,S_s,S_s4f,S_n,Theta_s,Theta_s4f,Theta_n,".
			      	  "Pi_s,Pi_s4f,Pi_n,TajimasD_s,TajimasD_s4f,TajimasD_n,Psing_s,Psing_s4f,Psing_n,CAI,StopCod\n";
			}
			print "$gene,$Nind," . join(',',@res) . "\n";
		}
	}
}

sub iupac {
	my ($nuc) = @_;
	$nuc =~ tr/a-z/A-Z/;
	my %ambig = (
'R' => ['A','G'],
'W' => ['A','T'],
'M' => ['A','C'],
'S' => ['C','G'],
'Y' => ['C','T'],
'K' => ['T','G']);
	if ($ambig{$nuc}){
		return(@{$ambig{$nuc}});
	} else {
		return($nuc);
	}
}

sub getbest {
	my @x = @_;
	my @count=();
	for my $i (0 .. $#x){
		my $se = $data{$x[$i]};
		my $co = $se =~ tr/[ATGCatgc]//;
		push @count, $co;
	}
	my $max = max(@count);
	my ($best) = grep { $count[$_] =~ /$max/ } 0 .. $#count;
	return($best);
}

sub max {
	my @x = sort {$a <=> $b} @_;
	return(@x[0]);
}

sub round {
	my ($val) = @_;
	my $down = floor($val);
	my $up = ceil($val);
	my $mid = $down + 0.5;
	if ($mid > $val){
		return($up);
	} else {
		return($down);
	}
}

sub siteSortCount {
	my @c = @_;
	my @cu = uniq(@c);
	my @f = ();
	for $cod (@cu){
		push @f, scalar(grep { /$cod/ } @c);
	}
	my @idx = sort {$f[$a] <=> $f[$b]} 0 .. $#f;
	my @fs = @f[@idx];
	my @cus = @cu[@idx];
	return(\@cus,\@fs);
}

sub countDiffTsTv {
	my ($i,$j) = @_;
	my $p=0;
	my $q=0;
	for my $x (0 .. $#$i){
		next if (($$i[$x] !~ m/[ATGC]/i) or ($$j[$x] !~ m/[ATGC]/i));
		if ($$i[$x] ne $$j[$x]){
			if (( scalar(grep { /[AG]/i } ($$i[$x],$$j[$x])) == 2 ) or 
			   ( scalar(grep { /[TC]/i } ($$i[$x],$$j[$x])) == 2 )) {
				++$p;
			}
			elsif (( scalar(grep { /[AT]/i } ($$i[$x],$$j[$x])) == 2) or 
			      ( scalar(grep { /[AC]/i } ($$i[$x],$$j[$x])) == 2) or 
			      ( scalar(grep { /[TG]/i } ($$i[$x],$$j[$x])) == 2) or
			      ( scalar(grep { /[CG]/i } ($$i[$x],$$j[$x])) == 2)) {
				++$q;
			}
		}
	}
	return(\$p,\$q);
}

sub checkSynSites {
	my ($cod) = @_;
	map { $_ =~ tr/a-z/A-Z/ } @$cod;
	if ((scalar(grep { /CT[ACTG]/i } @$cod) >= 1) or 
	    (scalar(grep { /TC[ACTG]/i } @$cod) >= 1) or 
	    (scalar(grep { /CC[ACTG]/i } @$cod) >= 1) or 
	    (scalar(grep { /AC[ACTG]/i } @$cod) >= 1) or 
	    (scalar(grep { /GC[ACTG]/i } @$cod) >= 1) or 
	    (scalar(grep { /GT[ACTG]/i } @$cod) >= 1) or 
	    (scalar(grep { /CG[ACTG]/i } @$cod) >= 1) or 
	    (scalar(grep { /GG[ACTG]/i } @$cod) >= 1)){
		return((2,4));
	}
	elsif (scalar(grep { /AT[ACT]/i } @$cod) >= 1){
		return((2,3));
	}
	elsif ((scalar(grep { /TT[CT]/i } @$cod) >= 1) or 
	       (scalar(grep { /TT[AG]/i } @$cod) >= 1) or 
	       (scalar(grep { /TA[CT]/i } @$cod) >= 1) or 
	       (scalar(grep { /CA[CT]/i } @$cod) >= 1) or 
	       (scalar(grep { /CA[AG]/i } @$cod) >= 1) or 
	       (scalar(grep { /AA[CT]/i } @$cod) >= 1) or 
	       (scalar(grep { /AA[AG]/i } @$cod) >= 1) or 
	       (scalar(grep { /GA[CT]/i } @$cod) >= 1) or
	       (scalar(grep { /GA[AG]/i } @$cod) >= 1) or 
	       (scalar(grep { /TG[CT]/i } @$cod) >= 1) or
	       (scalar(grep { /AG[CT]/i } @$cod) >= 1) or
	       (scalar(grep { /AG[AG]/i } @$cod) >= 1) or
	       (scalar(grep { /GA[AG]/i } @$cod) >= 1)){
		return((2,2));
	}
	elsif ((scalar(grep { /[CT]TA/i } @$cod) == 2) or 
	       (scalar(grep { /[CT]TG/i } @$cod) == 2) or 
	       (scalar(grep { /[AC]GA/i } @$cod) == 2) or 
	       (scalar(grep { /[AC]GG/i } @$cod) == 2)){
		return((0,2));
	}
	elsif ((scalar(grep { /ATG/i } @$cod) >= 1) or 
	       (scalar(grep { /TGG/i } @$cod) >= 1)){
		return((2,0));
	}
	else {
		return((0,0));
	}
}

sub trans {
	my ($cod) = @_;
	$cod =~ tr/a-z/A-Z/;
	my %gc = (
'TTT' => 'F','TCT' => 'S','TAT' => 'Y','TGT' => 'C',
'TTC' => 'F','TCC' => 'S','TAC' => 'Y','TGC' => 'C',
'TTA' => 'L','TCA' => 'S','TAA' => '*','TGA' => '*',
'TTG' => 'L','TCG' => 'S','TAG' => '*','TGG' => 'W',

'CTT' => 'L','CCT' => 'P','CAT' => 'H','CGT' => 'R',
'CTC' => 'L','CCC' => 'P','CAC' => 'H','CGC' => 'R',
'CTA' => 'L','CCA' => 'P','CAA' => 'Q','CGA' => 'R',
'CTG' => 'L','CCG' => 'P','CAG' => 'Q','CGG' => 'R',

'ATT' => 'I','ACT' => 'T','AAT' => 'N','AGT' => 'S',
'ATC' => 'I','ACC' => 'T','AAC' => 'N','AGC' => 'S',
'ATA' => 'I','ACA' => 'T','AAA' => 'K','AGA' => 'R',
'ATG' => 'M','ACG' => 'T','AAG' => 'K','AGG' => 'R',

'GTT' => 'V','GCT' => 'A','GAT' => 'D','GGT' => 'G',
'GTC' => 'V','GCC' => 'A','GAC' => 'D','GGC' => 'G',
'GTA' => 'V','GCA' => 'A','GAA' => 'E','GGA' => 'G',
'GTG' => 'V','GCG' => 'A','GAG' => 'E','GGG' => 'G'
	);
	if ($gc{$cod}){
		return($gc{$cod});
	} else {
		return('?');
	}
}

sub trans2 {
	my ($cod) = @_;
	$cod =~ tr/a-z/A-Z/;
	my %gc = (
'TTT' => 'F','TCT' => 'S4','TAT' => 'Y','TGT' => 'C',
'TTC' => 'F','TCC' => 'S4','TAC' => 'Y','TGC' => 'C',
'TTA' => 'L2','TCA' => 'S4','TAA' => '*','TGA' => '*',
'TTG' => 'L2','TCG' => 'S4','TAG' => '*','TGG' => 'W',

'CTT' => 'L4','CCT' => 'P','CAT' => 'H','CGT' => 'R4',
'CTC' => 'L4','CCC' => 'P','CAC' => 'H','CGC' => 'R4',
'CTA' => 'L4','CCA' => 'P','CAA' => 'Q','CGA' => 'R4',
'CTG' => 'L4','CCG' => 'P','CAG' => 'Q','CGG' => 'R4',

'ATT' => 'I','ACT' => 'T','AAT' => 'N','AGT' => 'S2',
'ATC' => 'I','ACC' => 'T','AAC' => 'N','AGC' => 'S2',
'ATA' => 'I','ACA' => 'T','AAA' => 'K','AGA' => 'R2',
'ATG' => 'M','ACG' => 'T','AAG' => 'K','AGG' => 'R2',

'GTT' => 'V','GCT' => 'A','GAT' => 'D','GGT' => 'G',
'GTC' => 'V','GCC' => 'A','GAC' => 'D','GGC' => 'G',
'GTA' => 'V','GCA' => 'A','GAA' => 'E','GGA' => 'G',
'GTG' => 'V','GCG' => 'A','GAG' => 'E','GGG' => 'G'
	);
	if ($gc{$cod}){
		return($gc{$cod});
	} else {
		return('?');
	}
}

sub codfreq3pos {
	my ($cod) = @_;
	$cod =~ tr/a-z/A-Z/;
	my %gcf = (
'TTT' => 0.5,'TCT' => 0.75,'TAT' => 0.5,'TGT' => 0.5,
'TTC' => 0.5,'TCC' => 0.75,'TAC' => 0.5,'TGC' => 0.5,
'TTA' => 0.5,'TCA' => 0.75,'TAA' => 0,  'TGA' => 0,
'TTG' => 0.5,'TCG' => 0.75,'TAG' => 0,  'TGG' => 0,

'CTT' => 0.75,'CCT' => 0.75,'CAT' => 0.5,'CGT' => 0.75,
'CTC' => 0.75,'CCC' => 0.75,'CAC' => 0.5,'CGC' => 0.75,
'CTA' => 0.75,'CCA' => 0.75,'CAA' => 0.5,'CGA' => 0.75,
'CTG' => 0.75,'CCG' => 0.75,'CAG' => 0.5,'CGG' => 0.75,

'ATT' => 2/3, 'ACT' => 0.75,'AAT' => 0.5,'AGT' => 0.5,
'ATC' => 2/3, 'ACC' => 0.75,'AAC' => 0.5,'AGC' => 0.5,
'ATA' => 2/3, 'ACA' => 0.75,'AAA' => 0.5,'AGA' => 0.5,
'ATG' => 0,   'ACG' => 0.75,'AAG' => 0.5,'AGG' => 0.5,

'GTT' => 0.75,'GCT' => 0.75,'GAT' => 0.5,'GGT' => 0.75,
'GTC' => 0.75,'GCC' => 0.75,'GAC' => 0.5,'GGC' => 0.75,
'GTA' => 0.75,'GCA' => 0.75,'GAA' => 0.5,'GGA' => 0.75,
'GTG' => 0.75,'GCG' => 0.75,'GAG' => 0.5,'GGG' => 0.75
	);
	if ($gcf{$cod}){
		return($gcf{$cod});
	}
}

sub codfreq1pos {
	my ($cod) = @_;
	$cod =~ tr/a-z/A-Z/;
	my %gcf = (
'TTT' => 0,  'TCT' => 0,'TAT' => 0,'TGT' => 0,
'TTC' => 0,  'TCC' => 0,'TAC' => 0,'TGC' => 0,
'TTA' => 0.5,'TCA' => 0,'TAA' => 0,'TGA' => 0,
'TTG' => 0.5,'TCG' => 0,'TAG' => 0,'TGG' => 0,

'CTT' => 0,  'CCT' => 0,'CAT' => 0,'CGT' => 0,
'CTC' => 0,  'CCC' => 0,'CAC' => 0,'CGC' => 0,
'CTA' => 0.5,'CCA' => 0,'CAA' => 0,'CGA' => 0.5,
'CTG' => 0.5,'CCG' => 0,'CAG' => 0,'CGG' => 0.5,

'ATT' => 0, 'ACT' => 0,'AAT' => 0,'AGT' => 0,
'ATC' => 0, 'ACC' => 0,'AAC' => 0,'AGC' => 0,
'ATA' => 0, 'ACA' => 0,'AAA' => 0,'AGA' => 0.5,
'ATG' => 0, 'ACG' => 0,'AAG' => 0,'AGG' => 0.5,

'GTT' => 0,'GCT' => 0,'GAT' => 0,'GGT' => 0,
'GTC' => 0,'GCC' => 0,'GAC' => 0,'GGC' => 0,
'GTA' => 0,'GCA' => 0,'GAA' => 0,'GGA' => 0,
'GTG' => 0,'GCG' => 0,'GAG' => 0,'GGG' => 0
	);
	if ($gcf{$cod}){
		return($gcf{$cod});
	}
}

sub codonClass {
	my @codons = @_;
	my %codon_class = (
		'1' => ['ATG','TGG'],
		'2' => ['TTT','TTC','TAT','TAC','CAT','CAC','CAA','CAG']
	);
}

sub AllSegSites {
	my ($sl,$pop,$Nind,$data,$thresh,$frac,$out) = @_;
	my @ss = ();
	my @co_u = ();
	my @co_f = ();
	my $n2 = floor($$Nind/2);
	my @sfs_u = map 0, 0 .. ($n2+($n2-1));
	my @sfs_f = map 0, 0 .. $n2;
	my @ds = ();
	my @do = ();
	my @di = ();
	for my $i (0 .. ($$sl-1)){
		my @site=();
		foreach(@$pop){
			if (substr($$data{$_},$i,1) =~ m/[ACTG]/i){
				push @site, substr($$data{$_},$i,1);
			}
		}
		my ($alleles,$count) = siteSortCount(@site);
		if (scalar(@$alleles) == 2){
			if ((($$thresh < 1) and (scalar(@site) >= round($$Nind*$$thresh))) and 
			((($$frac < 1) and (($$count[0]/scalar(@site)) > $$frac)) or 
			(($$frac >= 1) and ($$count[0] > $$frac)))){
				push @ss, $i;
				if ((length($$out) != 0) and (substr($$out,$i,1) =~ m/[ATGC]/i)){
					if ($$alleles[0] ne substr($$out,$i,1)){
						$sfs_u[$$count[0]] += 1;
						push @co_u, $$count[0];
					} else {
						$sfs_u[$$Nind-$$count[0]] += 1;
						push @co_u, ($$Nind-$$count[0]);
					}
				} else {
					$sfs_u[$$count[0]] += 1;
					push @co_u, $$count[0];
				}
				$sfs_f[$$count[0]] += 1;
				push @co_f, $$count[0];
			}
			elsif ((($$thresh > 1) and (scalar(@site) > $$thresh)) and 
			((($$frac < 1) and (($$count[0]/scalar(@site)) > $$frac)) or 
			(($$frac >= 1) and ($$count[0] > $$frac)))){
				push @ss, $i;
				if ((length($$out) != 0) and (substr($$out,$i,1) =~ m/[ATGC]/i)){
					if ($$alleles[0] ne substr($$out,$i,1)){
						$sfs_u[$$count[0]] += 1;
						push @co_u, $$count[0];
					} else {
						$sfs_u[$$Nind-$$count[0]] += 1;
						push @co_u, ($$Nind-$$count[0]);
					}
				} else {
					$sfs_u[$$count[0]] += 1;
					push @co_u, $$count[0];
				}
				$sfs_f[$$count[0]] += 1;
				push @co_f, $$count[0];
			}
		}
		elsif (scalar(@$alleles) > 2){
			if ((($$thresh < 1) and (scalar(@site) >= round($$Nind*$$thresh))) and 
			(($$frac < 1) and (($$count[0]/scalar(@site)) > $$frac))){
				push @ss, $i;
				@allelesL = @$alleles[($#$alleles-1) .. $#$alleles];
				@countL = @$count[($#$count-1) .. $#$count];
				if ((length($$out) != 0) and (substr($$out,$i,1) =~ m/[ATGC]/i)){
					if ($allelesL[0] ne substr($$out,$i,1)){
						$sfs_u[$countL[0]] += 1;
						push @co_u, $countL[0];
					} else {
						$sfs_u[$$Nind-$countL[0]] += 1;
						push @co_u, ($$Nind-$countL[0]);
					}
				} else {
					$sfs_u[$countL[0]] += 1;
					push @co_u, $countL[0];
				}
				$sfs_f[$countL[0]] += 1;
				push @co_f, $countL[0];
			}
			elsif ((($$thresh >= 1) and (scalar(@site) >= $$thresh)) and 
			(($$frac >= 1) and ($$count[0] > $$frac))){
				push @ss, $i;
				@allelesL = @$alleles[($#$alleles-1) .. $#$alleles];
				@countL = @$count[($#$count-1) .. $#$count];
				if ((length($$out) != 0) and (substr($$out,$i,1) =~ m/[ATGC]/i)){
					if ($allelesL[0] ne substr($$out,$i,1)){
						$sfs_u[$countL[0]] += 1;
						push @co_u, $countL[0];
					} else {
						$sfs_u[$$Nind-$countL[0]] += 1;
						push @co_u, ($$Nind-$countL[0]);
					}
				} else {
					$sfs_u[$countL[0]] += 1;
					push @co_u, $countL[0];
				}
				$sfs_f[$countL[0]] += 1;
				push @co_f, $countL[0];
			}
		}
		elsif (scalar(@$alleles) == 1){
			if (length($$out) != 0){
				if ($$alleles[0] ne substr($$out,$i,1)){
					push @ds, $i;
					push @do, substr($$out,$i,1);
					push @di, $$alleles[0];
				}
			}
			
		}
	}
	if (length($$out) != 0){
		shift @sfs_f;
		shift @sfs_u;
		return(\@ss,\@sfs_f,\@sfs_u,\@co_f,\@co_u,\@ds,\@do,\@di);
	} else {
		shift @sfs_f;
		return(\@ss,\@sfs_f,\@co_f);
	}
}

sub SegSitesCoding {
	my ($sl,$pop,$data,$ss,$out,$ds) = @_;
	my @csp4=();
	my @csp3=();
	my @csp2=();
	my @csd4=();
	my @csd3=();
	my @csd2=();
	my $totalcod;
	my $stop=0;
	my @fourfold_p=();
	my @twofold_p=();
	my @threefold_p=();
	my @zerofold_p=();
	my @fourfold_d=();
	my @twofold_d=();
	my @threefold_d=();
	my @zerofold_d=();
	my %aa_d=();
	for (my $i = $start_frame; $i < $$sl; $i += 3){
		my @cod=();
		my $codOut;
		my $stop_pos;
		foreach(@$pop){
			if (substr($$data{$_},$i,3) =~ m/[ACTG]{3}/i){
				if (substr($$data{$_},$i,3) =~ m/TAA|TAG|TGA/i){
					$stop_pos = 1;
				} else {
					push @cod, substr($$data{$_},$i,3);
				}
			}
		}
		my ($codcheck,$codfreq) = siteSortCount(@cod);
		map { $p = trans2($_); push @{$aa_d{$p}}, $_ } @$codcheck;
		my @syn = checkSynSites(\@$codcheck);
		if ($syn[1] == 4){
			my $pos3 = $i+$syn[0];
			my $pos1 = $i;
			my $pos2 = $i+1;
			if (my ($idx) = grep { $$ss[$_] == $pos3 } 0 .. $#$ss){
				push @fourfold_p, $idx;
			}
			if (my ($idx) = grep { $$ss[$_] == $pos1 } 0 .. $#$ss){
				push @zerofold_p, $idx;
			}
			if (my ($idx) = grep { $$ss[$_] == $pos2 } 0 .. $#$ss){
				push @zerofold_p, $idx;
			}
			push @csp4, @$codcheck[0];
		}
		elsif ($syn[1] == 3){
			my $pos3 = $i+$syn[0];
			my $pos1 = $i;
			my $pos2 = $i+1;
			if (my ($idx) = grep { $$ss[$_] == $pos3 } 0 .. $#$ss){
				push @threefold_p, $idx;
			}
			if (my ($idx) = grep { $$ss[$_] == $pos1 } 0 .. $#$ss){
				push @zerofold_p, $idx;
			}
			if (my ($idx) = grep { $$ss[$_] == $pos2 } 0 .. $#$ss){
				push @zerofold_p, $idx;
			}
			push @csp3, @$codcheck[0];
		}
		elsif ($syn[1] == 2){
			my ($pos1,$pos2,$pos3);
			if ($syn[0] == 0){
				$pos1 = $i+$syn[0];
				$pos2 = $i+1;
				if (my ($idx) = grep { $$ss[$_] == $pos1 } 0 .. $#$ss){
					push @twofold_p, $idx;
				}
				if (my ($idx) = grep { $$ss[$_] == $pos2 } 0 .. $#$ss){
					push @zerofold_p, $idx;
				}
				push @csp2, @$codcheck[0];
			}
			elsif ($syn[0] == 2) {
				$pos3 = $i+$syn[0];
				$pos2 = $i+1;
				$pos1 = $i;
				if (my ($idx) = grep { $$ss[$_] == $pos3 } 0 .. $#$ss){
					push @twofold_p, $idx;
				}
				if (my ($idx) = grep { $$ss[$_] == $pos2 } 0 .. $#$ss){
					push @zerofold_p, $idx;
				}
				if (my ($idx) = grep { $$ss[$_] == $pos1 } 0 .. $#$ss){
					push @zerofold_p, $idx;
				}
				push @csp2, @$codcheck[0];
			}
		}
		elsif ($syn[1] == 0){
			$pos = $i+$syn[0];
			if (my ($idx) = grep { $$ss[$_] == $pos } 0 .. $#$ss){
				push @zerofold_p, $idx;
			}
			$pos2 = $i+1;
			if (my ($idx) = grep { $$ss[$_] == $pos2 } 0 .. $#$ss){
				push @zerofold_p, $idx;
			}
		}
		if (length($$out) != 0){
			if (substr($$out,$i,3) =~ m/[ACTG]{3}/i){
				if (substr($$out,$i,3) =~ m/TAA|TAG|TGA/i){
					if ($stop_pos != $i){
						$stop_pos = 1;
					}
				} else {
					$codOut = substr($$out,$i,3);
				}
			}
			if (length($codOut) == 3){
				if (scalar(@$codcheck) == 1){
					my @c = ($$codcheck[0],$codOut);
					my @syn = checkSynSites(\@c);
					if ($syn[1] == 4){
						my $pos3 = $i+$syn[0];
						my $pos1 = $i;
						my $pos2 = $i+1;
						if (my ($idx) = grep { $$ds[$_] == $pos3 } 0 .. $#$ds){
							push @fourfold_d, $idx;
						}
						if (my ($idx) = grep { $$ds[$_] == $pos1 } 0 .. $#$ds){
							push @zerofold_d, $idx;
						}
						if (my ($idx) = grep { $$ds[$_] == $pos2 } 0 .. $#$ds){
							push @zerofold_d, $idx;
						}
						push @csd4, @c[0];
					}
					elsif ($syn[1] == 3){
						my $pos3 = $i+$syn[0];
						my $pos1 = $i;
						my $pos2 = $i+1;
						if (my ($idx) = grep { $$ds[$_] == $pos3 } 0 .. $#$ds){
							push @threefold_d, $idx;
							push @csd3, @c[0];
						}
						if (my ($idx) = grep { $$ds[$_] == $pos1 } 0 .. $#$ds){
							push @zerofold_d, $idx;
						}
						if (my ($idx) = grep { $$ds[$_] == $pos2 } 0 .. $#$ds){
							push @zerofold_d, $idx;
						}
						push @csd3, @c[0];
					}
					elsif ($syn[1] == 2){
						my ($pos1,$pos2,$pos3);
						if ($syn[0] == 0){
							$pos1 = $i+$syn[0];
							$pos2 = $i+1;
							if (my ($idx) = grep { $$ds[$_] == $pos1 } 0 .. $#$ds){
								push @twofold_d, $idx;
							}
							if (my ($idx) = grep { $$ds[$_] == $pos2 } 0 .. $#$ds){
								push @zerofold_d, $idx;
							}
							push @csd2, @c;
						}
						elsif ($syn[0] == 2) {
							$pos3 = $i+$syn[0];
							$pos2 = $i+1;
							$pos1 = $i;
							if (my ($idx) = grep { $$ds[$_] == $pos3 } 0 .. $#$ds){
								push @twofold_d, $idx;
							}
							if (my ($idx) = grep { $$ds[$_] == $pos2 } 0 .. $#$ds){
								push @zerofold_d, $idx;
							}
							if (my ($idx) = grep { $$ds[$_] == $pos1 } 0 .. $#$ds){
								push @zerofold_d, $idx;
							}
							push @csd2, @c;
						}
					}
					elsif ($syn[1] == 0){
						$pos = $i+$syn[0];
						if (my ($idx) = grep { $$ds[$_] == $pos } 0 .. $#$ds){
							push @zerofold_d, $idx;
						}
						$pos2 = $i+1;
						if (my ($idx) = grep { $$ds[$_] == $pos2 } 0 .. $#$ds){
							push @zerofold_d, $idx;
						}
					}
				}
			}
		}
		$totalcod += 1;
		$stop += $stop_pos;
	}
    my $scfp4=0;
	foreach my $codi (@csp4){
		$scfp4 += codfreq3pos($codi);
	}
	my $scfp3=0;
	foreach my $codi (@csp3){
		$scfp3 += codfreq3pos($codi);
    	}
	my $scfp2=0;
	foreach my $codi (@csp2){
		$scfp2 += (codfreq3pos($codi) + codfreq1pos($codi));
    	}
	my $all_scfp=$scfp4+$scfp3+$scfp2;
	my $ncfp = ($totalcod*3)-$all_scfp;
	#my $y;
	#foreach(keys %aa_d){
	#	print $_ . " " . "@{$aa_d{$_}}\n";
	#	$x = F_ENC(\$_, \%aa_d);
	#	next if ($x == 0); 
	#	$y += 1/$x;
	#}
	#print $y . "\n";
	#Ne_ENC(%aa_d);
	#exit;
	my $CAI = sprintf("%.5f",cai(%aa_d));
	if (length($$out) != 0){
		my $scfd4=0;
		foreach my $codi (@csd4){
			$scfd4 += codfreq3pos($codi);
	    	}
		my $scfd3=0;
		foreach my $codi (@csd3){
			$scfd3 += codfreq3pos($codi);
	    	}
		my $scfd2=0;
		foreach my $codi (@csd2){
			$scfd2 += (codfreq3pos($codi) + codfreq1pos($codi));
	    	}
		my $all_scfd=$scfd4+$scfd3+$scfd2;
	    my $ncfd = ($totalcod*3)-$all_scfd;
		return(\@fourfold_p,\@threefold_p,\@twofold_p,\@zerofold_p,
			\$scfp4,\$scfp3,\$scfp2,\$ncfp,
			\@fourfold_d,\@threefold_d,\@twofold_d,\@zerofold_d,
			\$scfd4,\$scfd3,\$scfd2,\$ncfd,\$CAI,
			\$stop);
	} else {
		return(\@fourfold_p,\@threefold_p,\@twofold_p,\@zerofold_p,
			\$scfp4,\$scfp3,\$scfp2,\$ncfp,\$CAI,
			\$stop);
	}
}

sub F_ENC {
	my ($aa,$cod_aa) = @_;
	my %cods_in_aa = (
		'F'=>['TTT','TTC'],
		'L2'=>['TTA','TTG'],
		'L4'=>['CTT','CTC','CTA','CTG'],
		'I'=>['ATT','ATC','ATA'],
		'M'=>['ATG'],
		'V'=>['GTT','GTC','GTA','GTG'],
		'S4'=>['TCT','TCC','TCA','TCG'],
		'S2'=>['AGT','AGC'],
		'P'=>['CCT','CCC','CCA','CCG'],
		'T'=>['ACT','ACC','ACA','ACG'],
		'A'=>['GCT','GCC','GCA','GCG'],
		'Y'=>['TAT','TAC'],
		'H'=>['CAT','CAC'],
		'Q'=>['CAA','CAG'],
		'N'=>['AAT','AAC'],
		'K'=>['AAA','AAG'],
		'D'=>['GAT','GAC'],
		'E'=>['GAA','GAG'],
		'C'=>['TGT','TGC'],
		'W'=>['TGG'],
		'R2'=>['CGT','CGC'],
		'R4'=>['CGA','CGG','AGA','AGG'],
		'G'=>['GGT','GGC','GGA','GGG']);
	my @cod_count=();
	for $AA (@{$cods_in_aa{$$aa}}){
		push @cod_count, scalar(grep { /$AA/ } @{$$cod_aa{$$aa}});
	}
	my $m_enc = scalar(@{$cods_in_aa{$$aa}});
	my $n_enc = eval(join('+',@cod_count));
	my @p1_enc = map { $_ / $n_enc } @cod_count;
	my @p2_enc = map { $_**2 } @p1_enc;
	my @p3_enc = map { (($_+1)/($n_enc+$m_enc))**2 } @cod_count;
	my $sum_pi_enc = eval(join('+',@p3_enc));
	#my $F_enc;
	#eval { $F_enc = (($n_enc * $sum_pi_enc)-1)/($n_enc - 1) } or $F_enc = 0;
	#print $F_enc . ' ' . $sum_pi_enc . "\n";
	return($sum_pi_enc);
}

sub Ne_ENC {
	my %cods = @_;
	my %cod_class = (
		'M'=>1,'W'=>1,
		'F'=>2,'Y'=>2,'H'=>2,'Q'=>2,'N'=>2,'K'=>2,'D'=>2,'E'=>2,'C'=>2,'L2'=>5,'S2'=>5,'R2'=>5,
		'I'=>3,
		'V'=>4,'P'=>4,'T'=>4,'A'=>4,'G'=>4,'L4'=>5,'S4'=>5,'R4'=>5);
	my %class_Fs;
	foreach (keys %cods){
		if ($cod_class{$_} == 1){
			push @{$class_Fs{1}}, F_ENC(\$_,\%cods);
		}
		elsif ($cod_class{$_} == 2){
			push @{$class_Fs{2}}, F_ENC(\$_,\%cods);
		}
		elsif ($cod_class{$_} == 3){
			push @{$class_Fs{3}}, F_ENC(\$_,\%cods);
		}
		elsif ($cod_class{$_} == 4){
			push @{$class_Fs{4}}, F_ENC(\$_,\%cods);
		}
	}
	my $enc;
	foreach (1, 2, 3, 4, 5){
		my @Fs = @{$class_Fs{$_}};
		my $Fst = scalar(@Fs);
		my $Fs_average = eval(join('+',@Fs))/$Fst;
		my $Fs_ratio = $Fst/$Fs_average;
#		print $Fst . " " . $Fs_average . " " . $Fs_ratio . "\n";
		$enc += $Fs_ratio;
	}
	print $enc . "\n";
}

sub cai {
	my %amino = @_;
	my @AA = keys %amino;
	@AA = grep { $_ !~ /[\*|M]/ } @AA;
	my @CAIi=();
	for $aa (@AA){
		my %seen;
		$seen{$_}++ for @{$amino{$aa}};
		@tmp = sort { $seen{$b} <=> $seen{$a} } (keys %seen);
		if (scalar(@tmp) > 1){
			foreach(@tmp){
				$W = sprintf("%.2f",$seen{$_}/$seen{$tmp[0]});
				push @CAIi, $W;
			}
		}
	}
	my $CAI = eval(join("*",@CAIi))**(1/scalar(keys %amino));
	return($CAI);
}

sub sfs {
	my ($co,$Nind,$fol) = @_;
	my %seen;
	++$seen{$_} for @$co;
	my $n2 = floor($$Nind/2);
	my @sfs_u = map 0, 0 .. ($n2+($n2-1));
	my @sfs_f = map 0, 0 .. $n2;
	if ($fol == 0){
		for my $i (0 .. $n2){
			if ($seen{$i}){
				$sfs_f[$i] = $seen{$i};
			}
		}
		shift(@sfs_f);
		return(@sfs_f);
	} 
	elsif ($fol == 1){
		for my $i (0 .. ($n2+($n2-1))){
			if ($seen{$i}){
				$sfs_u[$i] = $seen{$i};
			}
		}
		shift(@sfs_u);
		return(@sfs_u);
	}
}

sub polymorphism {
	my ($sl,$Nind,$sfs,$f) = @_;
	my @sfs_scaled = map { $i = $_+1; $$sfs[$_]*$i } 0 .. $#$sfs;
	my @sfs_ps = map { $_/$$sl } @sfs_scaled;
	my $sums = eval(join("+",@$sfs));
	my @P = map { $i = $_+1; $sfs_ps[$_]*($$Nind-$i) } 0 .. $#sfs_ps;
	my @Pt = map { $i = $_+1; $sfs_scaled[$_]*($$Nind-$i) } 0 .. $#sfs_scaled;
	my $pis = (2/($$Nind*($$Nind-1)))*(eval(join("+",@P)));
	my $pit = (2/($$Nind*($$Nind-1)))*(eval(join("+",@Pt)));
	my ($a1,$Dvar,$Hvar) = Dvar(\$$Nind,\$sums,\$f);
	my $ths = ($sums/$$sl)/$$a1;
	my $tht = $sums/$$a1;
	my $D;
	eval { $D = ($pit-$tht)/$$Dvar } or $D = "NA";
	my @res;
	if ($f == 1){
		my @h = map { $i = $_+1; $sfs_scaled[$_]*((((2*$$Nind)-(4*$i))/($$Nind*($$Nind-1)))) } 0 .. $#sfs_scaled;
		my $H;
		eval { $H = eval(join("+",@h))/$$Hvar } or $H = "NA";
		@res = ($ths,$pis,$D,$H);
	} else {
		@res = ($ths,$pis,$D);
	}
	push @res, $sfs_ps[0];
	if ($sums != 0){
		map { $_ = sprintf("%.5f",$_) } @res;
	}
	unshift @res, $sums;
	return(@res);
}

sub Dvar {
	my ($n,$s,$f) = @_;
	my $a1=0;
	my $a2=0;
	foreach(1 .. ($$n-1)){
		$a1 += 1/$_;
		$a2 += 1/($_**2);
	}
	foreach(1 .. $$n){
		$bn1 += 1/($_**2);
	}
	my $b1 = ($$n+1)/(3*($$n-1));
	my $b2 = (2*(($$n**2)+$$n+3))/(9*$$n*($$n-1));
	my $c1 = $b1 - (1/$a1);
	my $c2 = $b2 - (($$n+2)/($a1*$$n)) + ($a2/($a1**2));
	my $e1 = $c1/$a1;
	my $e2 = $c2/(($a1**2)+$a2);
	my $Dvar = sqrt(($e1*$$s)+($e2*$$s*($$s-1)));
	if ($$f == 1){
		# See Zeng et al. 2006 Genetics 174: 1431–1439
		my $left = (($$n-2)/(6*($$n-1)))*($$s/$a1);
		my $right_up = (18*($$n**2)*((3*$$n)+2)*$bn1)-((88*($$n**3))+(9*($$n**2))-(13*$$n)+6);
		my $right_down = 9*$$n*(($$n-1)**2);
		my $Hvar = sqrt($left+(($right_up/$right_down)*($$s*($$s-1)/(($a1**2)+$a2))));
		return(\$a1,\$Dvar,\$Hvar);
	} else {
		return(\$a1,\$Dvar);
	}
}

sub divergence {
	my ($sl,$di,$do) = @_;
	my ($dp,$dq) = countDiffTsTv(\@$di,\@$do);
	my $dxy;
	eval { $dxy = K2Pd(\($$dp/$$sl),\($$dq/$$sl)) or $dxy = 0 };
	$dxy = sprintf("%.5f",$dxy); 
	return($dxy);
}

sub mkt {
	my ($res_s,$res_n,$div_s,$div_n,$dxy_s,$dxy_n) = @_;
	my ($S_s,$S_n) = ($$res_s[0],$$res_n[0]);
	my ($pi_s,$pi_n) = ($$res_s[2],$$res_n[2]);
	my $alpha;
	my $ni;
	my $pval;
	if (any { $_ == 0 } ($pi_s,$pi_n,$$div_s,$$div_n) ){
		$alpha = 'NA';
		$ni = 'NA';
		$pval = 'NA';
		return (($ni,$alpha,$pval));
	} else {
		my @obs = ([$S_s,$$div_s],[$S_n,$$div_n]);
		my $chi = new Statistics::ChisqIndep;
		$chi->load_data(\@obs);
		$pval = $chi->{p_value};
		eval { $alpha = 1 - (($$dxy_s*$pi_n)/($$dxy_n*$pi_s)) } or $alpha = 'NA';
		eval { $ni = ($pi_n/$pi_s)/($$dxy_n/$$dxy_s) } or $ni = 'NA';
		if (any { $_ ne 'NA' } ($ni,$alpha)){
			map { $_ = sprintf("%.5f",$_) } ($ni,$alpha,$pval);
			return(($ni,$alpha,$pval));
		} else {
			if ($pval != 0){
				$pval = sprintf("%.5f",$pval);
			}
			return (($ni,$alpha,$pval));
		}
	}
}

sub JCd {
	my ($p1) = @_;
	my $d1;
	eval { $d1 = -(3/4)*log(1-((4*$$p1)/3)) } or $d1 = 0;
	return($d1);
}

sub K2Pd {
	my ($p1,$q1) = @_;
	my $K1;
	eval { $K1 = -(1/2)*log((1-(2*$$p1)-$$q1)*sqrt(1-(2*$$q1))) } or $K1 = 0;
	return($K1);
}





