#!/usr/bin/perl -w

use strict;
use warnings;

my$file = $ARGV[0];
my@data = `cat $file`;


my$entry = $ARGV[1];

open(OUT,">>PlasmoConversion-v1.txt");

my%match;
foreach my$r (@data){
	my@q = split("=",$r);
	$match{$q[0]} = $q[1];
}



if(!exists $match{$entry})
{
	system("./PlasmoOld2New.sh $entry");
	my@f = `grep 'PF3D7' tmp | cut -d'|' -f1,3| tr ">" " "`;
	#$f =~ s/^\s|\s*$//g;
	print OUT "$entry=@f\n";
	print "$entry=@f\n";
	
}
else
{
	my$we = $match{$entry};
	$we =~ s/^\s|\s*$//g;
	print "$entry = $we\n";
}





