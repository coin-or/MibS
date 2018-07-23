#!/usr/bin/env bash

cd ../..

PWD=`pwd`

MibSPATH=$PWD
RESULTSPATH=$PWD/output
SCRIPTPATH=$PWD/scripts

cd $SCRIPTPATH/analyze

rm -r -f tmp

figNum=5
echo Plotting Figure $figNum

#noHeuristics
for testSet in dataMIBLP-XU dataIBLP-DEN dataIBLP-FIS
do
    if [ $testSet == dataMIBLP-XU ]
    then
	appendix=MiblpXu
    elif
    then
	appendix=IblpDen
    else
	appendix=IblpFis
    fi

    CURRENTOUTPUTPATH=$RESULTSPATH/method4Output/$testSet
    LISTPATH=$SCRIPTPATH/analyze/usefulFiles/r1LeqR2List$appendix

    for file in ${CURRENTINSTPATH}/*.out
    do
	instance_name=`basename ${file%.*}`
	if grep -r "${instance_name}.aux" "$LISTPATH" > /dev/null
	then
	    cp $CURRENTINSTPATH/${instance_name}.out tmp
	fi
    done

    CURRENTOUTPUTPATH=$RESULTSPATH/method7Output/$testSet
    LISTPATH=$SCRIPTPATH/analyze/usefulFiles/r1GeR2List$appendix

    for file in ${CURRENTINSTPATH}/*.out
    do
	instance_name=`basename ${file%.*}`
	if grep -r "${instance_name}.aux" "$LISTPATH" > /dev/null
	then
	    cp $CURRENTINSTPATH/${instance_name}.out tmp
	fi
    done
done

cd tmp

cat *.out > outNoHeuristics

cd ..

cp tmp/outNoHeuristic parse

rm -r -f tmp

cd parse

awk -f Timeparse-mibs.awk outNoHeuristic

mv time.summary noHeuristics

cp noHeuristics ../performance

rm -f noHeuristics

rm -f outNoHeuristic

cd ..

#impObjectiveCut, secondLevelPriority and weightedSums  
for i in {12..14}
do
    mkdir tmp
    for testSet in dataMIBLP-XU dataIBLP-DEN dataIBLP-FIS
    do
	name=$RESULTSPATH/method${i}Output/$testSet
	cp $name/*.out tmp

    done

    cd tmp

    cat *.out > outMethod$i

    cd ..

    cp tmp/outMethod$i parse

    rm -r -f tmp

    cd parse

    awk -f Timeparse-mibs.awk outMethod$i

    if [ $i == 12 ]
    then
	name=impObjectiveCut
    elif [ $i == 13 ]
    then
	name=secondLevelPriority
    else
	name=weightedSums
    fi

    mv time.summary $name

    cp $name ../performance

    rm -r -f $name

    cd ..

done

cd performance

python perf.py -c 1 --legend --x-limit=5 noHeuristics impObjectiveCut secondLevelPriority weightedSums > figure5.eps

rm noHeuristics
rm impObjectiveCut
rm secondLevelPriority
rm weightedSums

	

    
	
