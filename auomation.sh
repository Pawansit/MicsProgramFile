#!/usr/bin/bash

while read line
do
	perl PlasmoOld2New.pl PlasmoConversion-v1.txt $line
	
done < "col.txt"
