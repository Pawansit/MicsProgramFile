#!/usr/bin/perl -w


use strict; 
use warnings; 

###### Input File ##################

system("rm Protein_nodes_400.txt Match_PP_state.csv Match_Anno_state.csv");			##Remove the Intermediate file.
system("cut -f6 78_posi_dat_uniprot-v2.txt| sort -u > Protein_nodes_400.txt");		## Extract the Plasmodb ID for Interaction search in STRING Database.
open(C,"Protein_nodes_400.txt");													## Open the Plasmodb ID file.



####### OutPut File ################
open(D,">CDPK_interaction_Protein_nodes_400.sif");							# This file store the Pair wise intractio information
open(L,">CDPK_search_400.log");												# This file will store the intermediate information	
open(V,">CDPK_interaction_Protein_nodes_400.csv");							# This file store the Node protein functon information 

############# parameters ###########
my$Interaction_cutoff = 400;											# Cutoff for string pair confidance based selection 




my@entry = <C>;
my%annotation;

for(my$i = 0;$i <scalar@entry;$i++)										## Iterate over all Nodes.
{
	my$j = $i + 1;
	$entry[$i] =~ s/^\s+|\s*$//g;
	my($ID1,$ID2) = Interaction_String($entry[$i],$Interaction_cutoff);	# This will return the list of interacting pairs in used cutoff value
	my$count = scalar@{$ID1};
	print L "$entry[$i]\t$count\n";										# Number of interaction pairs per Node.
	
	
	for(my$e =0;$e <scalar@{$ID1};$e++)
	{
		my$first  = `perl ../Plasmo-conversion/PlasmoOld2New.pl ../Plasmo-conversion/PlasmoConversion-v1.txt $$ID1[$e]`;		## Extract the New Plasmodb ID and corresponding function
		my$second = `perl ../Plasmo-conversion/PlasmoOld2New.pl ../Plasmo-conversion/PlasmoConversion-v1.txt $$ID2[$e]`;
		my@l = split("=",$first);
		my@m = split("=",$second);
		my@t = split(/\|/,$l[1]);
		my@s = split(/\|/,$m[1]);
		$t[0] =~ s/^\s|\s*$//g;
		$s[0] =~ s/^\s|\s*$//g;
		$t[1] =~ s/^\s|\s*$//g;											# protein function Annotation
		$s[1] =~ s/^\s|\s*$//g;
		my$PPmatch = PP_matcher($$ID1[$e],$t[0]);						## match entry with Schizont Proteome data
		Annotaion($$ID1[$e],$t[1]);										## Match the prtein type in the protein function
		if(!exists $annotation{$$ID1[$e]})
		{
			$annotation{$$ID1[$e]} = $t[1];
		}
		print D "$$ID1[$e]\tinteracts with\t$$ID2[$e]\n";
	}
}
	print V " \"Name\",\"Description\"\n";
	foreach my$u (keys %annotation){
		print V "\"$u\",\"$annotation{$u}\"\n";
	}
	
 close(C);
 close(D);
 close(L);
 close(V);

system("echo Name,Fuction,state > tmp1");
system("cat tmp1 Match_Anno_state.csv > Match_Anno_state-v1.csv");
system("mv Match_Anno_state-v1.csv Match_Anno_state.csv");

system("echo Name,NewName,PPstate >tmp2");
system("cat tmp2 Match_PP_state.csv > Match_PP_state-v1.csv");
system("mv Match_PP_state-v1.csv Match_PP_state.csv");

system("rm tmp1 tmp2 PlasmoConversion-v1.txt");
####################### Subroutine #####################################
 
 
sub Update_table
{
	open(B,"Match_Anno_state_400.csv");										
	open(C,"Match_PP_state_400.csv");
	open(F,">Annotation_table_400.csv");										# Annotation table for the network annotation
	my@c = <C>;
	print F "Name,PPState,ProteinType,ProteinTypeV\n";
	while(my$line = <B>)
	{
		my@w = split("\t", $line);
		for(my$i=0;$i<scalar@c;$i++)
		{
			my@qs = split(",",$c[$i]);
			
			if($w[0] eq $qs[0])
			{
				$w[2] =~ s/^\s+|\s*$//g;
				$qs[2] =~ s/^\s+|\s*$//g;
				print F "\"$w[0]\",\"$qs[2]\",\"$w[1]\",\"$w[2]\"\n";
			}
		}
	}	
	close(B);
	close(C);
	close(F);
}


 
sub Annotaion															#Annotate by protein Function type
{
	my($pattern_o,$anno) = @_;
	open(O,">>Match_Anno_state_400.csv");									#Network ID match with Plasmodb annotation 
	
	my@A = qw(kinase transport membrane);
	#my@AA = ($anno =~ /$A[1]/g);
	#print "$pattern_o,@AA\n";
	
	if($anno =~ m/$A[0]/g)
	{
		print O "$pattern_o,Kinase,1\n";
	}
	elsif($anno =~ m/$A[1]/g)
	{
		print O "$pattern_o,Transport,2\n";
	}
	elsif($anno =~ m/$A[2]/g)
	{
		print O "$pattern_o,membrane,3\n";
	}
	close(O);
}
 

sub PP_matcher
{
	my($pattern_o,$pattern_n) = @_;
	open(P,"PP-All-unique.csv");										# Schizont stage Proteome Plasmodb list
	open(O,">>Match_PP_state_400.csv");										# Network Plasmodb ID matching with schizont stage proteome ID.
	#print O "Plasmo_o,Plasmo_N,Match\n";
	while(my$line = <P>)
	{
		$line =~ s/^\s+|\s*$//g;
		if($line eq $pattern_n){ print O "$pattern_o,$line,1\n"; return(1);} 
	}
	close(P);
	close(O);
}

sub Interaction_String
{
	my($pattern,$cutoff) =@_;
	open(A,"../../../String-v11/5833.protein.links.full.v11.0.txt");	# String Plasmodium collection version 11
	my@store = ();
	my%pair;
	my@keys = ();
	my@value = ();
	while(my$line = <A>)
	{
		if($line =~ m/$pattern/)
		{
			my@w1 = split(" ",$line);
			if($w1[15] >= $cutoff)
			{
				my$s = substr($w1[0],5,15);
				my$ss = substr($w1[1],5,15);
				if(!exists $pair{$s})
				{
					$pair{$s} = $pattern;
				}
			}
		}
	}
	
	for my$key (keys %pair)
	{
		push(@keys,$key);
		push(@value,$pair{$key});
		#print "$key = $pair{$key}\t$number\n";
	}
	return(\@keys,\@value);
	close(A);
}
