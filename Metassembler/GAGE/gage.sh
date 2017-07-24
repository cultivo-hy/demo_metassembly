#!/bin/bash

GAGE_DIR="__GAGE_DIR"

ref=$1
ctg=$2
scf=$3

USAGE='USAGE: '$(basename $0)' <ref> <ctg> <scf>'

if [ ! $ref ]; then echo -e $USAGE'\n'; exit 1; fi
if [ ! $ctg ]; then echo -e $USAGE'\n'; exit 1; fi
if [ ! $scf ]; then echo -e $USAGE'\n'; exit 1; fi

$GAGE_DIR/getCorrectnessStats.sh $ref $ctg $scf 
