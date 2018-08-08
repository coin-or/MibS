#!/usr/bin/env bash

cd ../..

PWD=`pwd`

MibSPATH=$PWD
RESULTSPATH=$PWD/output
SCRIPTPATH=$PWD/scripts
USEFULPATH=$SCRIPTPATH/analyze/usefulFiles
instanceListDir=$USEFULPATH/finalInstanceList

cd $SCRIPTPATH/analyze

rm -r -f tmp

mkdir tmp

figNum=5
echo Plotting Figure $figNum

#noHeuristics
for testSet in dataMIBLP-XU dataIBLP-DEN dataIBLP-FIS
do
    if [ $testSet == dataMIBLP-XU ]
    then
	appendix=MiblpXu
    elif [ $testSet == dataIBLP-DEN ]
    then
	appendix=IblpDen
    else
	appendix=IblpFis
    fi

    CURRENTOUTPUTPATH=$RESULTSPATH/method4Output/$testSet
    LISTPATH=$SCRIPTPATH/analyze/usefulFiles/r1LeqR2List$appendix
    solvedInstDir=$USEFULPATH/finalInstanceList

    for file in ${CURRENTOUTPUTPATH}/*.out
    do
	instance_name=`basename ${file%.*}`
	if grep -r "${instance_name}.aux" "$LISTPATH" > /dev/null
	then
	    if grep -r "${instance_name}.aux" "$solvedInstDir" > /dev/null
	    then
		cp $CURRENTOUTPUTPATH/${instance_name}.out tmp
	    fi
	fi
    done

    CURRENTOUTPUTPATH=$RESULTSPATH/method7Output/$testSet
    LISTPATH=$SCRIPTPATH/analyze/usefulFiles/r1GeR2List$appendix

    for file in ${CURRENTOUTPUTPATH}/*.out
    do
	instance_name=`basename ${file%.*}`
	if grep -r "${instance_name}.aux" "$LISTPATH" > /dev/null
	then
	    if grep -r "${instance_name}.aux" "$solvedInstDir" > /dev/null
	    then
		cp $CURRENTOUTPUTPATH/${instance_name}.out tmp
	    fi
	fi
    done
done

cd tmp

cat *.out > outNoHeuristics

cd ..

cp tmp/outNoHeuristics parse

rm -r -f tmp

cd parse

awk -f Timeparse-mibs.awk outNoHeuristics

mv time.summary noHeuristics

cp noHeuristics ../performance

rm -f noHeuristics

rm -f outNoHeuristics

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

    rm -r -f outMethod$i

    cd ..

done

cd performance

python perf.py -c 1 --legend --x-limit=5 noHeuristics impObjectiveCut secondLevelPriority weightedSums > figure5.eps

rm noHeuristics
rm impObjectiveCut
rm secondLevelPriority
rm weightedSums

	

    
	
