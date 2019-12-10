#!/usr/bin/perl -w 

##############################################################
#    	Author: Pawan Kumar
#		Use: Convert the Json Output file from the Cytoscape to the CellDesigner input file
#
#############################################################




use strict;
use warnings;

################### Read Json Output File from the Cytoscape and write Nodes and Edges Information in file ###########
my$JsonReader = "PKA-new.cyjs";

system("Rscript parseJSON.R $JsonReader");


########################### Read Nodes and Edges information ###########
my@file1 = split(/\./,$JsonReader);

open(A,$file1[0]."_Nodes.csv");
open(B,$file1[0]."_Edges.csv");

my@Nodes = <A>;
my@Edges = <B>;

###### Read Nodes and Edges information and create input for Cell Designer  routine ######

my($Idss,$Xcoor, $Ycoor, $Sourceprotein, $targetprotein) = NodesandEdges(\@Nodes,\@Edges);
	
my@IDs 				= Sequence("sa",$Idss);   						
my@Species 			= Sequence("s",$Idss); 							
my@x 				= @{$Xcoor};
my@y 				= @{$Ycoor};
my@protein 			= Sequence("pr",$Idss);       					
my$Num_Nodes        = scalar@$Idss;
my@annotation 		= repetition("A",$Num_Nodes);
my$Num_reaction 	= scalar@$Sourceprotein;
my@reaction 		= repetition("re",$Num_reaction); 				
my@reactant 		= Sequence("s",$Sourceprotein);					
my@reactant_alias 	= Sequence("sa",$Sourceprotein);
my@reactant_Cpmd 	= repetition("CDMT0000",$Num_reaction); 
my@product  		= Sequence("s",$targetprotein); 
my@product_alias 	= Sequence("sa",$targetprotein);
my@product_Cpmd 	= repetition("CDMT0000",$Num_reaction); 



########## Write the final xml output in a file ###########

CellDesigner(\@IDs,\@Species,\@x,\@y,\@protein,\@annotation,\@reaction,\@reactant,\@reactant_alias,\@reactant_Cpmd,\@product,\@product_alias,\@product_Cpmd,$file1[0]);


############ SubRoutines ######################

sub NodesandEdges
{
	my($Nodes,$Edges) = @_;
	my@ID = ();
	my@X = ();
	my@Y = ();
	my@Source = ();
	my@target = ();
	my$j = 0;
	for(my$i=1;$i<scalar@{$Nodes};$i++)
	{
		my@w = split(",",$$Nodes[$i]);
		$w[1] =~ s/"//g;
		$w[2] =~ s/"//g;
		$w[3] =~ s/"//g;
		$w[3] =~ s/\n//g;
		$ID[$j] = $w[1];
		$X[$j]  = $w[2];
		$Y[$j]  = $w[3];
		$j++;
	}
	my@Xnew = translation(@X);
	my@Ynew = translation(@Y);
	
	my$z = 0;
	for(my$j = 1;$j <scalar@{$Edges};$j++)
	{
		my@M = split(",",$$Edges[$j]);
		$M[1] =~ s/"//g;
		$M[2] =~ s/"//g;
		$M[2] =~ s/\n//g;
		$Source[$z] = $M[1];
		$target[$z] = $M[2];
		$z++;
	}
	return(\@ID,\@Xnew,\@Ynew,\@Source,\@target);
}
	

sub CellDesigner
{
	my($ID,$Species,$x,$y,$protein,$annotation,$reaction,$reactant,$reactant_alias,$reactant_Cpmd,$product,$product_alias,$product_Cpmd,$fileinput) = @_;
	
	my$filestring = $fileinput.".xml";
	open(OUT,">$filestring");

	print OUT '<?xml version="1.0" encoding="UTF-8"?>',"\n";
	print OUT '<sbml xmlns="http://www.sbml.org/sbml/level2/version4" xmlns:celldesigner="http://www.sbml.org/2001/ns/celldesigner" level="2" version="4">',"\n";
	print OUT '<model metaid="Demo1" id="Demo1">',"\n";
	
	
	my$z = 80.0;
	my$h = 40.0;
	print OUT  "\<annotation\>\n";
	print OUT  "\<celldesigner:extension\>\n";
	print OUT  "\<celldesigner:modelVersion\>4.0\</celldesigner:modelVersion\>\n";
	print OUT  "\<celldesigner:modelDisplay sizeX=\"3000\" sizeY=\"2000\"\/\>\n";
	print OUT  "\<celldesigner:listOfCompartmentAliases\/\>\n";
	print OUT  "\<celldesigner:listOfComplexSpeciesAliases\/\>\n";
	print OUT  "\<celldesigner:listOfSpeciesAliases\>\n";
	
	for(my$i=0;$i<scalar@{$ID};$i++)
	{
		print OUT  "\<celldesigner:speciesAlias id=\"$$ID[$i]\" species=\"$$Species[$i]\" \>\n";
		print OUT  "\<celldesigner:activity>inactive</celldesigner:activity\>\n";
		print OUT  "\<celldesigner:bounds x=\"$$x[$i]\" y=\"$$y[$i]\" w=\"$z\" h=\"$h\"\/\>\n";
		print OUT  "\<celldesigner:font size=\"12\"\/\>\n";
		print OUT  "\<celldesigner:view state=\"usual\"\/\>\n";
		print OUT  "\<celldesigner:usualView\>\n";
		print OUT  "\<celldesigner:innerPosition x=\"0.0\" y=\"0.0\"\/\>\n";
		print OUT  "\<celldesigner:boxSize width=\"80.0\" height=\"40.0\"\/\>\n";
		print OUT  "\<celldesigner:singleLine width=\"1.0\"\/\>\n";
		print OUT  "\<celldesigner:paint color=\"ffccffcc\" scheme=\"Color\"\/\>\n";
		print OUT  "\<\/celldesigner:usualView\>\n";
		print OUT  "\<celldesigner:briefView\>\n";
		print OUT  "\<celldesigner:innerPosition x=\"0.0\" y=\"0.0\"\/\>\n";
		print OUT  "\<celldesigner:boxSize width=\"80.0\" height=\"60.0\"\/\>\n";
		print OUT  "\<celldesigner:singleLine width=\"0.0\"\/\>\n";
		print OUT  "\<celldesigner:paint color=\"3fff0000\" scheme=\"Color\"\/\>\n";
		print OUT  "\</celldesigner:briefView\>\n";
		print OUT  "\<celldesigner:info state=\"empty\" angle=\"-1.5707963267948966\"\/\>\n";
		print OUT  "\<\/celldesigner:speciesAlias\>\n";
	}
	print OUT  "\<\/celldesigner:listOfSpeciesAliases\>\n";
	print OUT  "\<celldesigner:listOfGroups\/\>\n";
	print OUT  "\<celldesigner:listOfProteins\>\n";
	for(my$i = 0;$i < scalar@{$ID};$i++)
	{
		print OUT  "\<celldesigner:protein id=\"$$protein[$i]\" name=\"$$Species[$i]\" type=\"GENERIC\"\/\>\n";
	}
	print OUT  "\<\/celldesigner:listOfProteins\>\n";
	print OUT  "\<celldesigner:listOfGenes\/\>\n";
	print OUT  "\<celldesigner:listOfRNAs\/\>\n";
	print OUT  "\<celldesigner:listOfAntisenseRNAs\/\>\n";
	print OUT  "\<celldesigner:listOfLayers\/\>\n";
	print OUT  "\<celldesigner:listOfBlockDiagrams\/\>\n";
	print OUT  "\<\/celldesigner:extension\>\n";
	print OUT  "\<\/annotation\>\n";
	
	units();
	compartment();
	print OUT  "\<listOfSpecies\>\n";
	
	for(my$i = 0;$i < scalar@{$ID};$i++)
	{
		if(length($$annotation[$i]) >= 5 )
		{
		  print OUT  "\<species metaid=\"$$Species[$i]\" id=\"$$Species[$i]\" name=\"$$Species[$i]\" compartment=\"default\" initialAmount=\"0\"\>\n";
		  print OUT  "\<notes\>\n";
		  print OUT   '<html xmlns="http://www.w3.org/1999/xhtml">',"\n";
		  print OUT  "\<head\>\n";
		  print OUT  "\<title\/\>\n";
		  print OUT  "\<\/head\>\n";
		  print OUT  "\<body\>Link to String v10:\n";
		  print OUT  "\<a href=\"http:\/\/version10.string-db.org\/interactions\/P62344\" target=\"#blank\"\>P62344\<\/a\>\n";
		  print OUT  "\<\/body\>\n";
		  print OUT  "\<\/html\>\n";
		  print OUT  "\<\/notes\>\n";
		  print OUT  "\<annotation\>\n";
		  print OUT  "\<celldesigner:extension\>\n";
		  print OUT  "\<celldesigner:positionToCompartment\>inside\<\/celldesigner:positionToCompartment\>\n";
		  print OUT  "\<celldesigner:speciesIdentity\>\n";
		  print OUT  "\<celldesigner:class\>PROTEIN\<\/celldesigner:class\>\n";
		  print OUT  "\<celldesigner:proteinReference\>$$protein[$i]\<\/celldesigner:proteinReference\>\n";
		  print OUT  "\<\/celldesigner:speciesIdentity\>\n";
		  print OUT  "\<\/celldesigner:extension\>\n";
		  print OUT  "\<\/annotation\>\n";
		  print OUT  "\<\/species\>\n";
		}
		else
		{
		  print OUT  "\<species metaid=\"$$Species[$i]\" id=\"$$Species[$i]\" name=\"$$Species[$i]\" compartment=\"default\" initialAmount=\"0\"\>\n";
		  print OUT  "\<annotation\>\n";
		  print OUT  "\<celldesigner:extension\>\n";
		  print OUT  "\<celldesigner:positionToCompartment\>inside\<\/celldesigner:positionToCompartment\>\n";
		  print OUT  "\<celldesigner:speciesIdentity\>\n";
		  print OUT  "\<celldesigner:class\>PROTEIN\<\/celldesigner:class\>\n";
		  print OUT  "\<celldesigner:proteinReference\>$$protein[$i]\<\/celldesigner:proteinReference\>\n";
		  print OUT  "\<\/celldesigner:speciesIdentity\>\n";
		  print OUT  "\<\/celldesigner:extension\>\n";
		  print OUT  "\<\/annotation\>\n";
		  print OUT  "\<\/species\>\n";
		}
	}
	print OUT  "\<\/listOfSpecies\>\n";
	print OUT  "\<listOfReactions\>\n";
	

	
	
	for(my$j=0;$j<scalar@{$reactant};$j++)
	{
		print OUT  "\<reaction metaid=\"$$reaction[$j]\" id=\"$$reaction[$j]\" reversible=\"false\"\>\n";
		print OUT  "\<annotation\>\n";
		print OUT  "\<celldesigner:extension\>\n";
		print OUT  "\<celldesigner:reactionType\>STATE_TRANSITION\<\/celldesigner:reactionType\>\n";
		print OUT  "\<celldesigner:baseReactants\>\n";
		print OUT  "\<celldesigner:baseReactant species=\"$$reactant[$j]\" alias=\"$$reactant_alias[$j]\"\>\n";  
		print OUT  "\<celldesigner:linkAnchor position=\"E\"\/\>\n";
		print OUT  "\<\/celldesigner:baseReactant\>\n";
		print OUT  "\<\/celldesigner:baseReactants\>\n";
		print OUT  "\<celldesigner:baseProducts\>\n";
		print OUT  "\<celldesigner:baseProduct species=\"$$product[$j]\" alias=\"$$product_alias[$j]\"\>\n";
		print OUT  "\<celldesigner:linkAnchor position=\"W\"\/\>\n";
		print OUT  "\<\/celldesigner:baseProduct\>\n";
		print OUT  "\<\/celldesigner:baseProducts\>\n";
		print OUT  "\<celldesigner:connectScheme connectPolicy=\"direct\" rectangleIndex=\"0\"\>\n";
		print OUT  "\<celldesigner:listOfLineDirection\>\n";
		print OUT  "\<celldesigner:lineDirection index=\"0\" value=\"unknown\"\/\>\n";
		print OUT  "\<\/celldesigner:listOfLineDirection\>\n";
		print OUT  "\<\/celldesigner:connectScheme\>\n";
		print OUT  "\<celldesigner:line width=\"1.0\" color=\"ff000000\"\/\>\n";
		print OUT  "\<\/celldesigner:extension\>\n";
		print OUT  "\<\/annotation\>\n";
		print OUT  "\<listOfReactants\>\n";
		print OUT  "\<speciesReference metaid=\"$$reactant_Cpmd[$j]\" species=\"$$reactant[$j]\"\>\n";
		print OUT  "\<annotation\>\n";
		print OUT  "\<celldesigner:extension\>\n";
		print OUT  "\<celldesigner:alias\>$$reactant_alias[$j]\<\/celldesigner:alias\>\n";
		print OUT  "\<\/celldesigner:extension\>\n";
		print OUT  "\<\/annotation\>\n";
		print OUT  "\<\/speciesReference\>\n";
		print OUT  "\<\/listOfReactants\>\n";
		print OUT  "\<listOfProducts\>\n";
		print OUT  "\<speciesReference metaid=\"$$product_Cpmd[$j]\" species=\"$$product[$j]\"\>\n";
		print OUT  "\<annotation\>\n";
		print OUT  "\<celldesigner:extension\>\n";
		print OUT  "\<celldesigner:alias\>$$product_alias[$j]\<\/celldesigner:alias\>\n";
		print OUT  "\<\/celldesigner:extension\>\n";
		print OUT  "\<\/annotation\>\n";
		print OUT  "\<\/speciesReference\>\n";
		print OUT  "\<\/listOfProducts\>\n";
		print OUT  "\<\/reaction\>\n";
	}
	print OUT  "\<\/listOfReactions\>\n";
	print OUT  "\</model\>\n";
	print OUT  "\<\/sbml\>\n";
	
}

sub units
{
	print OUT "\<listOfUnitDefinitions\>\n";
	print OUT '<unitDefinition metaid="substance" id="substance" name="substance">',"\n";
	print OUT "\<listOfUnits\>\n";
	print OUT '<unit metaid="CDMT00001" kind="mole"/>',"\n";
	print OUT "\</listOfUnits\>\n";
	print OUT "\<\/unitDefinition\>\n";
	print OUT '<unitDefinition metaid="volume" id="volume" name="volume">',"\n";
	print OUT "\<listOfUnits\>\n";
	print OUT '<unit metaid="CDMT00002" kind="litre"/>',"\n";
	print OUT '</listOfUnits>',"\n";
	print OUT '</unitDefinition>',"\n";
	print OUT '<unitDefinition metaid="area" id="area" name="area">',"\n";
	print OUT '<listOfUnits>',"\n";
	print OUT '<unit metaid="CDMT00003" kind="metre" exponent="2"/>',"\n";
	print OUT '</listOfUnits>',"\n";
	print OUT '</unitDefinition>',"\n";
	print OUT '<unitDefinition metaid="length" id="length" name="length">',"\n";
	print OUT '<listOfUnits>',"\n";
	print OUT '<unit metaid="CDMT00004" kind="metre"/>',"\n";
	print OUT '</listOfUnits>',"\n";
	print OUT '</unitDefinition>',"\n";
	print OUT '<unitDefinition metaid="time" id="time" name="time">',"\n";
	print OUT '<listOfUnits>',"\n";
	print OUT '<unit metaid="CDMT00005" kind="second"/>',"\n";
	print OUT '</listOfUnits>',"\n";
	print OUT '</unitDefinition>',"\n";
	print OUT '</listOfUnitDefinitions>',"\n";	
}

sub compartment
{
	print OUT  "\<listOfCompartments\>\n";
	print OUT  '<compartment metaid="default" id="default" size="1" units="volume"/>',"\n";
	print OUT  "</listOfCompartments>\n";	
}


	
	
sub repetition
{
	my($pattern,$number) = @_;
	my@array = ();
	for(my$i = 0;$i<$number;$i++)
	{
		my$r = $pattern.$i;
		push(@array,$r);
	}
	return @array;
}
	
sub Sequence
{
	my($pattern,$Number) = @_;
	my@output = ();
	for(my$i=0;$i<scalar@{$Number};$i++)
	{
		my$r = $pattern.$$Number[$i];
		push(@output,$r);
	}
	return(@output);
}
	
	

sub translation
{
	my(@X)  = @_;
	my@W    = ();
	my@NewX = ();
	foreach (@X){
		push(@W,abs($_));
	}
	my$max = MAX(@W);
	my$Transtale = $max +50;
	
	foreach (@X){
		my$x = $_+ $Transtale;
		push(@NewX,$x);
	}
	return(@NewX);
}

sub MAX
{
	my@X = @_;
	my@m = ();
	foreach ( sort { $a <=> $b } @X )
	{
		push(@m,$_);
	}
	return(pop(@m))
}
	
