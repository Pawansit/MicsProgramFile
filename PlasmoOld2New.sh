#!/usr/bin/bash

echo $1

curl https://plasmodb.org/cgi-bin/geneSrt?project_id=PlasmoDB\&ids=$1\&type=genomic\&upstreamAnchor=Start\&upstreamSign=minus\&upstreamOffset=10\&downstreamAnchor=End\&downstreamSign=plus\&downstreamOffset=2000 >tmp
#grep 'PF3D7' tmp | cut -d'|' -f1| tr ">" " "


