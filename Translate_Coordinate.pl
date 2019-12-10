#!/usr/bin/perl -w


####################################################
#
#		Author: Pawan Kumar
#		Translate the Network-Nodes coordiate (x,y) to the (++) space
#		which will be useful for CellDesigner input creation.
#
####################################################




use strict;
use warnings;


open(A,"PKA-All_Nodes.csv");

my@ID = ();
my@X = ();
my@Y = ();
my$j = 0;
my@Nodes = <A>;

for(my$i=1;$i<scalar@Nodes;$i++)
{
	my@w = split(",",$Nodes[$i]);
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


sub translation
{
	my(@X)  = @_;
	my@W    = ();
	my@NewX = ();
	foreach (@X){
		push(@W,abs($_));
	}
	my$max = MAX(@W);
	my$Transtale = $max +10;
	
	foreach (@X){
		my$x = $_+ $Transtale;
		print "$_ $x\n";
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
	
