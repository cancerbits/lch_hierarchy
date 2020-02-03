#!/bin/bash

N="$1"
if [ "$N" == "" ] ; then
        N=`date '+%Y-%m-%d-%H:%M:%S'`
fi

echo "Execute complete analysis:"
echo $N

$TOOLS/Rscript $CODEBASE/lch/src/main.R $N
