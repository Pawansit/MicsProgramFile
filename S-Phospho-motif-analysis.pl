#!/usr/bin/perl -w


use strict;
use warnings; 
#use diagnostics;

############## Reading data file #######################################


open(A,"Plasmodb_Summary-v1.csv");
system("cut -d',' -f11 ../Processed-data-schizont/Sp-All.csv > tmp.txt");
system("cut -d';' -f1 tmp.txt > Sp-Motif.txt");
system("rm tmp.txt");
open(B,"Sp-Motif.txt");


#################################### Open Motif data ####################

my@data_Nonp = <A>;
my@slice_data_Nonp = @data_Nonp[1..200];
my@data_Pho  = <B>;

############# Non Phospo-motif #########################################
my@NonPho_motif = Motif_finder(@slice_data_Nonp);
my @unique_Nonp_motif = do { my %seen; grep { !$seen{$_}++ } @NonPho_motif };
my@slice_data_Nonp_unique = UniqueNonPmotif(\@data_Pho,\@unique_Nonp_motif);
my@AA_Nonp_freq = AA_count(\@slice_data_Nonp_unique);

##################
my@AA_Pho_freq =Serine_Array_OUT_File(\@slice_data_Nonp_unique,\@data_Pho);
Serine_Matix_OUT_File(\@slice_data_Nonp_unique,\@data_Pho,\@AA_Nonp_freq,\@AA_Pho_freq);


########### Motif Correction ###########################################
print scalar@data_Pho,"\n";
Smotif13mer(@data_Pho);



###################### Sub routines ####################################
sub Smotif13mer
{
	my@data = @_;
	open(OUT,">Sp-Motif-final.txt");
	
	for(my$i=0;$i<scalar@data;$i++)
	{
		my$r = length$data[$i];
		if($r == 16)
		{
			my$e1 = uc(substr($data[$i],1,13));
			print OUT "$e1\n";
		}
		#print OUT "\n";
	}
}






sub Serine_Array_OUT_File
{
	my($slice_data_Nonp_unique,$data_Pho) = @_;
	open(A,">Phopho-S-Array.csv");
	open(B,">Phopho-S-Array_Relative.csv");


	###############Phospho-motif########################################
	my@AA_Pho_freq = AA_count($data_Pho);
	my@AA = ("A,","C,","D,","E,","F,","G,","H,","I,","K,","L,","M,","N,","P,","Q,","R,","S,","T,","W,","Y,","V,");
	print A @AA,"\n";
	for(my$i = 0;$i <scalar@AA_Pho_freq;$i++)
	{
		print A "$AA_Pho_freq[$i],";
	}
	################## Normalize value Overall #########################
	my@overall_Nor = logfraction(\@AA_Pho_freq,\@AA_Nonp_freq);
	print B @AA,"\n";
	for(my$j = 0;$j<scalar@overall_Nor;$j++)
	{
		print B "$overall_Nor[$j],";
	}
	close(A);
	close(B);
	return(@AA_Pho_freq);
}


sub Serine_Matix_OUT_File
{
	my($slice_data_Nonp_unique,$data_Pho,$AA_Nonp_freq,$AA_Pho_freq) = @_;
	open(A,">Phopho-S-InD.csv");
	open(B,">Phopho-S-InD_Relative.csv");
	my@AA = ("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","W","Y","V");
	print A "AA,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6\n";
	print B "AA,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6\n";
	for(my$i=0;$i<scalar@AA;$i++){
		my@AA_ind_Nonp = AA_Ind_count($slice_data_Nonp_unique,$AA[$i]);
		my@AA_ind_Pho  = AA_Ind_count($data_Pho,$AA[$i]);
		print A "$AA[$i],";
		print B "$AA[$i],";
		for(my$p =1;$p<scalar@AA_ind_Nonp-1;$p++)
		{
			my$perNonp = substr((($AA_ind_Nonp[$p]/$$AA_Nonp_freq[$i])*100),0,4);
			my$perPho  = substr((($AA_ind_Pho[$p]/$$AA_Pho_freq[$i])*100),0,4);
			my$logratio= substr(log2($perPho/($perNonp+0.01)),0,4);
			print A "$perPho,";
			print B "$logratio,";
		}
		print A "\n";
		print B "\n";
	}
	close(A);
	close(B);
}

sub UniqueNonPmotif
{
	my($Motif1,$Motif2) = @_;
	my%Filter;
	my@M_New = ();
	foreach my$r(@{$Motif1})
	{
		$Filter{$r} = 1;
	}
	for(my$i=0;$i <scalar@{$Motif2};$i++)
	{
		if(!exists$Filter{$$Motif2[$i]})
		{
			push(@M_New,$$Motif2[$i]);
		}
	}
	return(@M_New);
}
	
sub logfraction
{
	my($A,$B) = @_;
	my@Normalize = ();
	for(my$i=0;$i<scalar@{$A};$i++)
	{
		my$r = ($$A[$i]/$$B[$i]);
		my$log = substr(log2($r),0,4);
		push(@Normalize,$log);	
	}
	return(@Normalize);
}


sub log2 
{ 
    my $n = shift; 
    if ($n ==0)
    {
		return $n;
	}
	else
	{
		return log($n) / log(2); 
	}
} 

sub AA_Ind_count
{
	my($data,$r) = @_;
	my@M = (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
	
	for(my$t =0;$t < scalar@{$data};$t++)
	{
		my@w = split("",$$data[$t]);
		for(my$i = 0;$i<scalar@w;$i++)
		{
			if($w[$i] eq $r)
			{
				$M[$i]  += 1;
			}
		}
	}
	return(@M);
}

sub AA_count
{
	my$data = shift;
	my@AA = ("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","W","Y","V");
	my@M = (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
	my@M_normaltize = ();
	for(my$t = 1;$t<scalar@{$data};$t++)
	{
		my@w = split("",$$data[$t]);
		for(my$i = 0;$i<scalar@w;$i++)
		{
			for(my$j=0;$j <scalar@AA;$j++)
			{
				if($w[$i] eq $AA[$j])
				{
					$M[$j]  += 1;
				}
			}
		}
	}
	return(@M);
}

sub Motif_finder
{
	my@d  = @_;
	my@Motif_All = ();
	for(my$i = 1;$i <scalar@d;$i++)
	{
		my@w = ();
		@w = split("!",$d[$i]);
		my@Spnon = Serinefinder($w[9]);
		my@Motif = NonPhosmotif(\@Spnon,$w[9]);
		push(@Motif_All,@Motif);
	}
	return(@Motif_All);
}

sub NonPhosmotif
{
	my($S,$seq) = @_;
	my@motif = ();
	for(my$i=0;$i<scalar@{$S};$i++)
	{
		if(($$S[$i]  > 7))
		{
			my$start = $$S[$i] - 6;
			my$str = substr($seq,$start,13);
			push(@motif,$str);
		}
	}
	return(@motif);
}
	
sub Serinefinder
{
	my$string = shift;
	my@array = split("",$string);
	my@S = ();
	for(my$i=0;$i<scalar@array;$i++)
	{
		if($array[$i] eq "S")
		{
			push(@S,$i);
		}
	}
	return(@S);
}
