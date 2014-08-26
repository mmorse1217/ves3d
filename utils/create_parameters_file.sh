#!/bin/bash
#

executable=ParabolicFlow.exe

initFile=$1
: ${initFile:=$executable.in}

outFile=$2
: ${outFile:=$executable.out}

centFile=$3
: ${centFile:=$executable.cntrs.in}

echo $initFile
echo $outFile
echo $centFile

# saveData
# initFile : precomputed/dumbbell_12_double.txt
# outFile  : ParabolicFlow.out
# centFile : ParabolicFlow.cntrs.in

# bendingModulus  : 1e-2
# bgFlowParam     : 0.01
# errorFactor     : 1
# nSurfs          : 1
# positionIterMax : 24
# positionRestart : 1
# positionTol     : 1e-6
# repMaxIter      : 200
# repTimeStep     : 1e-1
# repTol          : 1e-5
# saveStride      : .5
# scheme          : SemiImplicit
# shOrder         : 12
# singularStokes  : Direct
# tensionIterMax  : 24
# tensionRestart  : 1
# tensionTol      : 1e-6
# timeHorizon     : 100
# timeStep        : 1e-1
