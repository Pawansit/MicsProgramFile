#!/usr/bin/perl -w


use strict;
use warnings;


########################################################################
#
#
#	Script for finding the euqilavent Uniport and NCBI ID in PlasmodB 
#	database. This script use the blastp for the finding the simialrity 
#	based ID mappig. 
#
#	Useage: 
#	perl Database-ID-convertor.pl <Database name> <Input ID filename>
#
#	here Applicable databasesare NCBI or Uniport.
#
#
#
#
#
#
#
#
#
#
#
#
########################################################################



my$InputDatabaseName = $ARGV[0];																			### Database name to be used for the mapping
my$InputDatabaseID = $ARGV[1];																				### ID file


my%mapperID1 = ID_hashing(`cat NCBI2PlasmoDB-mappig.csv`);													## Create the has of the already store Converted ID
my%mapperID2 = ID_hashing(`cat Uniport2PlasmoDB-mappig.csv`);
my%mapperID3 = ID_hashing(`cat PlasmoOld2PlasmoNew.csv`);


open(G,">>NCBI2PlasmoDB-mappig.csv");
open(P,">>Uniport2PlasmoDB-mappig.csv");
open(H,">>PlasmoOld2PlasmoNew.csv");

my@ID  = InputModify(`cat $InputDatabaseID`);																## Remove the white space from the ID and store in the array

for(my$i=0;$i < scalar@ID;$i++)
{
	my($InputNCBI,$plasmodbIDNCBI,$similarityNCBI,$inputUniport,$PlasmodbIDUniPort,$similarityUniport);
	if($InputDatabaseName =~ /NCBI/)																		## If IDs belong to the NCBI
	{
		if(!exists $mapperID1{$ID[$i]})																		## Search in the alerdy store file
		{
			($InputNCBI,$plasmodbIDNCBI,$similarityNCBI)= NCBI2Plasmodb($ID[$i]);
			print G "$InputNCBI,$plasmodbIDNCBI,$similarityNCBI\n";
			print   "$InputNCBI,$plasmodbIDNCBI,$similarityNCBI\n";
		}
		else
		{
			my@we = split(/\./,$mapperID1{$ID[$i]});
			#$we =~ s/^\s|\s*$//g;
			print "$ID[$i] = $we[0]\n";
		}
	}
	if($InputDatabaseName =~ /Uniport/)
	{
		if(!exists $mapperID2{$ID[$i]})
		{
			my($inputUniport,$PlasmodbIDUniPort,$similarityUniport) = UniPort2Plasmodb($ID[$i]);
			print "$inputUniport\t$PlasmodbIDUniPort\t$similarityUniport\n";
			print P "$inputUniport,$PlasmodbIDUniPort,$similarityUniport\n";
		}
		else
		{
			my@we = split(/\./,$mapperID2{$ID[$i]});
			
			#$we =~ s/^\s|\s*$//g;
			print "$ID[$i] = $we[0]\n";
		}
	}
	
	if($InputDatabaseName =~ /Plasmodb/)
	{
		if(!exists $mapperID3{$ID[$i]})
		{
			my($inputplasmodb,$fd ) = PlmasmoOld2PlasmoNew($ID[$i]);
			my@PlasmodbNew = split(/\|/,$fd);
			$PlasmodbNew[1] =~  s/^\s|\s*$//g;
			print "$inputplasmodb\t$PlasmodbNew[0]\t$PlasmodbNew[1]\n";
			print H "$inputplasmodb,$PlasmodbNew[0],$PlasmodbNew[1]\n";
		}
		else
		{
			my$we = $mapperID3{$ID[$i]};
			$we =~ s/^\s|\s*$//g;
			print "$ID[$i] = $we\n";
		}
	}
}

close(G);
close(P);
close(H);



sub PlmasmoOld2PlasmoNew
{
	my($ID) = shift;
	system("./PlasmoOld2New.sh $ID");
	my$f = `grep 'PF3D7' tmp | cut -d'|' -f1,3| tr ">" " "`;
	print "$ID=$f\n";
	system("rm tmp");
	return($ID,$f);
}



sub NCBI2Plasmodb
{
	my($ID) = shift;
	system("wget -O pawan.fa 'https://www.ncbi.nlm.nih.gov/search/api/sequence/$ID/?report=fasta'");
	system("blastp  -query pawan.fa -db /home/pawan/NII-Work/Projects/CellDesigner/Program/Phospho-protein-analysis/PP-Downloaded-dataset/Plasmodb-protein_ID_mapping/backup/Plasmodb -out AlignOut.txt -outfmt 6 -num_alignments >= 1 -num_threads >= 4");
	my$w1  = `cut -f1 AlignOut.txt`;
	my$w11 = `cut -f2 AlignOut.txt`;
	my$w111 = `cut -f3 AlignOut.txt`;
	$w1 =~ s/^\s+|\s*$//;
	$w11 =~ s/^\s+|\s*$//;
	$w111 =~ s/^\s+|\s*$//;
	system("rm pawan.fa AlignOut.txt =");
	return($ID,$w11,$w111);
}



sub UniPort2Plasmodb
{
	my($UniID) = shift;
	
	if($UniID =~/^UPI*/)
	{
		system("wget https://www.uniprot.org/uniparc/$UniID.fasta");
		system("blastp  -query $UniID.fasta -db /home/pawan/NII-Work/Projects/CellDesigner/Program/Phospho-protein-analysis/PP-Downloaded-dataset/Plasmodb-protein_ID_mapping/backup/Plasmodb -out AlignOut.txt -outfmt 6 -num_alignments >= 1 -num_threads >= 4");
		system("rm $UniID.fasta");
		my@w1 = `cut -f1 AlignOut.txt`;
		my@w11 = `cut -f2 AlignOut.txt`;
		my@w111 = `cut -f3 AlignOut.txt`;
		$w1[0] =~ s/^\s|\s*$//;
		my@w3 = split(/\-/, $w11[0]);
		return($w1[0],$w3[0],$w111[0]);
	}
	else
	{
		system("wget https://www.uniprot.org/uniprot/$UniID.fasta");
		system("blastp  -query $UniID.fasta -db /home/pawan/NII-Work/Projects/CellDesigner/Program/Phospho-protein-analysis/PP-Downloaded-dataset/Plasmodb-protein_ID_mapping/backup/Plasmodb -out AlignOut.txt -outfmt 6 -num_alignments >= 1 -num_threads >= 4");
		system("rm $UniID.fasta");
		my@w1 = `cut -f1 AlignOut.txt`;
		my@w11 = `cut -f2 AlignOut.txt`;
		my@w111 = `cut -f3 AlignOut.txt`;
		my@w2 = split(/\|/, $w1[0]);
		my@w3 = split(/\-/, $w11[0]);
		return($w2[1],$w3[0],$w111[0]);
	}
}

sub InputModify
{
	my@ID = @_;
	my@Out = ();
	foreach my$r (@ID){
		$r =~ s/^\s|\s*$//;
		push(@Out,$r);
	}
	return(@Out);
}	

sub ID_hashing
{
	my(@data) = @_;
	my%match;
	foreach my$r (@data){
		my@q = split(",",$r);
		$match{$q[0]} = $q[1];
	}
	return(%match);
}
