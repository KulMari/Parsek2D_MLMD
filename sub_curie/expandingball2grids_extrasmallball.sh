#!/bin/bash
case='expandingball2grids_extrasmallball'
inputfile=$HOME'/Parsek2D-multilevel/inputfiles/inputfile.'$case
mkdir -p $SCRATCHDIR'/'$case'/restart'
mkdir -p $WORKDIR'/'$case'/restart'
#MSUB -r $case # Request name
#MSUB -n 128 # Number of tasks to use
#MSUB -T 86400 # Elapsed time limit in seconds
#MSUB -o $SCRATCHDIR'/'$case'/'%I.out # Standard output. %I is the job id
#MSUB -e $SCRATCHDIR'/'$case'/'%I.err # Error output. %I is the job id
#MSUB -q xlarge # The queue choice


cp $HOME'/Parsek2D-multilevel/ParsekEM' $WORKDIR'/'$case
cp $inputfile $WORKDIR'/'$case

BRIDGE_MSUB_PWD=$WORKDIR'/'$case
export BRIDGE_MSUB_PWD

set -x
cd ${BRIDGE_MSUB_PWD}
ccc_mprun ./ParsekEM inputfile.$case
