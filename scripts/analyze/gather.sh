#!/usr/bin/env bash

cd ../..

PWD=`pwd`

MibSPATH=$PWD
RESULTSPATH=$PWD/output
SCRIPTPATH=$PWD/scripts
INSTPATH=$mibs_dir/testSets/BilevelLib/general

cd $SCRIPTPATH/analyze

rm -r -f usefulFiles
mkdir usefulFiles
USEFULPATH=$SCRIPTPATH/analyze/usefulFiles

rm -f $SCRIPTPATH/analyze/performance/figure*
rm -r -f tmp

i=1
for testSet in dataMIBLP-XU dataIBLP-DEN dataIBLP-FIS
do
    mkdir tmp
    name=$RESULTSPATH'/method'$i'Output/'$testSet
    cp $name'/'*'.out' tmp

    cd tmp

    cat *.out > outMethod$i

    cd ..

    cp tmp/outMethod$i parse

    rm -r -f tmp

    cd parse

    awk -f InstanceSeparator-mibs.awk outMethod$i

    rm -f outMethod$i

    if [ $testSet == dataMIBLP-XU ]
    then
	appendix=MiblpXu
	touch r1LeqR2List.summary
    elif [ $testSet == dataIBLP-DEN ]
    then
	appendix=IblpDen
    else
	appendix=IblpFis
    fi
    name=r1LeqR2List$appendix
    mv r1LeqR2List.summary $name
    cp $name $USEFULPATH
    rm -f $name
    name=r1GeR2List$appendix
    mv r1GeR2List.summary $name
    cp $name $USEFULPATH
    rm -f $name
    cd ..
done

for figNum in '2a' '2b' '3a' '3b' '4a' '4b'
do
    echo Plotting Figure $figNum 
    if [ $figNum == '2a' ]
    then
	for i in {1..11}
	do
	    mkdir tmp
	    name=$RESULTSPATH'/method'$i'Output/dataIBLP-DEN'
            cp $name'/'*'.out' tmp

            name=$RESULTSPATH'/method'$i'Output/dataIBLP-FIS'
            cp $name'/'*'.out' tmp

            cd tmp

            cat *.out > outMethod$i

            cd ..

            cp tmp/outMethod$i parse

            rm -r -f tmp

            cd parse

            awk -f Timeparse-mibs.awk outMethod$i

            mv time.summary timeMethod$i.summary

	    cp timeMethod$i.summary ../refinement

	    rm timeMethod$i.summary

            rm outMethod$i

	    cd ..
	done

        cd refinement

        python refine.py 71

        for i in {1..11}
	do
	    rm timeMethod$i.summary
        done

        for i in {1..5}
        do
	    cp finalTimeMethod$i ../performance
        done

        for i in {1..11}
	do
	    rm finalTimeMethod$i
        done

        cd ../performance

        mv finalTimeMethod1 whenLInt-LInt
        mv finalTimeMethod2 whenLFixed-LFixed
        mv finalTimeMethod3 whenXYInt-LFixed
        mv finalTimeMethod4 whenXYIntOrLFixed-LFixed
        mv finalTimeMethod5 whenLInt-LFixed

        python perf.py -c 1 --legend --x-limit=12 whenLInt-LInt whenLInt-LFixed whenLFixed-LFixed whenXYInt-LFixed whenXYIntOrLFixed-LFixed > figure2a.eps

        rm when*

        cd ..
    elif [ $figNum == '2b' ]
    then
	for i in {1..11}
	do
	    mkdir tmp
	    name=$RESULTSPATH'/method'$i'Output/dataMIBLP-XU'

	    cp $name'/'*'.out' tmp

	    cd tmp

	    cat *.out > outMethod$i

	    cd ..

	    cp tmp/outMethod$i parse

	    rm -r -f tmp

	    cd parse

	    awk -f Timeparse-mibs.awk outMethod$i

	    mv time.summary timeMethod$i.summary

 	    cp timeMethod$i.summary ../refinement

	    rm timeMethod$i.summary

	    rm outMethod$i

	    cd ..
	done

        cd refinement

        python refine.py 100

        for i in {1..11}
	do
	    rm timeMethod$i.summary
	done

        for i in {1..5}
	do
	    cp finalTimeMethod$i ../performance
	done

	for i in {1..11}
	do
	    rm finalTimeMethod$i
	done

	cd ../performance

        mv finalTimeMethod1 whenLInt-LInt
        mv finalTimeMethod2 whenLFixed-LFixed
        mv finalTimeMethod3 whenXYInt-LFixed
        mv finalTimeMethod4 whenXYIntOrLFixed-LFixed
        mv finalTimeMethod5 whenLInt-LFixed

        python perf.py -c 1 --legend --x-limit=6.9 whenLInt-LInt whenLInt-LFixed whenLFixed-LFixed whenXYInt-LFixed whenXYIntOrLFixed-LFixed > figure2b.eps

        rm when*

        cd ..
    elif [ $figNum == '3a' ] || [ $figNum == '3b' ]
    then
	for i in {1..11}
	do
	    mkdir tmp
	    name=$RESULTSPATH'/method'$i'Output/dataIBLP-DEN'
	    cp $name'/'*'.out' tmp

	    name=$RESULTSPATH'/method'$i'Output/dataIBLP-FIS'
	    cp $name'/'*'.out' tmp

	    name=$RESULTSPATH'/method'$i'Output/dataMIBLP-XU'
	    cp $name'/'*'.out' tmp

	    cd tmp

	    cat *.out > outMethod$i

	    cd ..

	    cp tmp/outMethod$i parse

	    rm -r -f tmp

	    cd parse

	    if [ $figNum == '3a' ]
	    then
		awk -f Timeparse-mibsFig3a.awk outMethod$i
	    fi

	    if [ $figNum == '3b' ]
	    then
		awk -f Timeparse-mibsFig3b.awk outMethod$i
	    fi

	    mv time.summary timeMethod$i.summary

	    cp timeMethod$i.summary ../refinement

	    rm timeMethod$i.summary

	    rm outMethod$i

	    cd ..
	done

        cd refinement
    
        if [ $figNum == '3a' ]
	then
	    python refine.py 41
	fi

	if [ $figNum == '3b' ]
	then
	    python refine.py 130
	fi

	for i in {1..11}
	do
	    rm timeMethod$i.summary
	done

        cp finalTimeMethod4 ../performance
        cp finalTimeMethod7 ../performance

	for i in {1..11}
	do
	    rm finalTimeMethod$i
	done

        cd ../performance

        mv finalTimeMethod4 linkingBranchingStrategy
        mv finalTimeMethod7 fractionalBranchingStrategy

        if [ $figNum == '3a' ]
	then
	    python perf.py -c 1 --legend --x-limit=25 fractionalBranchingStrategy linkingBranchingStrategy > figure3a.eps
	fi

        if [ $figNum == '3b' ]
	then
	    python perf.py -c 1 --legend --x-limit=20 fractionalBranchingStrategy linkingBranchingStrategy > figure3b.eps
	fi

        rm *Strategy

        cd ..

    elif [ $figNum == '4a' ] || [ $figNum == '4b' ]
    then
	for i in {1..11}
	do
	    mkdir tmp
	    name=$RESULTSPATH'/method'$i'Output/dataIBLP-DEN'
	    cp $name'/'*'.out' tmp

	    name=$RESULTSPATH'/method'$i'Output/dataIBLP-FIS'
	    cp $name'/'*'.out' tmp

	    name=$RESULTSPATH'/method'$i'Output/dataMIBLP-XU'
	    cp $name'/'*'.out' tmp

	    cd tmp

	    cat *.out > outMethod$i

	    cd ..

	    cp tmp/outMethod$i parse

	    rm -r -f tmp

	    cd parse

	    awk -f Timeparse-mibs.awk outMethod$i

	    mv time.summary timeMethod$i.summary

	    cp timeMethod$i.summary ../refinement

	    rm timeMethod$i.summary

	    rm outMethod$i

	    cd ..
	done

        cd refinement

        python refine.py 171

        for i in {1..11}
	do
	    rm timeMethod$i.summary
	done

	if [ $figNum == '4a' ]
	then
	    cp finalTimeMethod6 ../performance
	    cp finalTimeMethod7 ../performance
	    cp finalTimeMethod9 ../performance
	    cp finalTimeMethod11 ../performance
	fi

	if [ $figNum == '4b' ]
	then
	    cp finalTimeMethod3 ../performance
	    cp finalTimeMethod4 ../performance
	    cp finalTimeMethod8 ../performance
	    cp finalTimeMethod10 ../performance
	fi

        for i in {1..11}
	do
	    rm finalTimeMethod$i
	done

        cd ../performance

	if [ $figNum == '4a' ]
	then
	    mv finalTimeMethod11 withoutPoolWhenXYInt-LFixed
	    mv finalTimeMethod6 withPoolWhenXYInt-LFixed
	    mv finalTimeMethod9 withoutPoolWhenXYIntOrLFixed-LFixed
	    mv finalTimeMethod7 withPoolWhenXYIntOrLFixed-LFixed	
	    python perf.py -c 1 --legend --x-limit=11 withoutPoolWhenXYInt-LFixed withPoolWhenXYInt-LFixed withoutPoolWhenXYIntOrLFixed-LFixed withPoolWhenXYIntOrLFixed-LFixed> figure4a.eps
	fi

	if [ $figNum == '4b' ]
	then
	    mv finalTimeMethod10 withoutPoolWhenXYInt-LFixed
	    mv finalTimeMethod3 withPoolWhenXYInt-LFixed
	    mv finalTimeMethod8 withoutPoolWhenXYIntOrLFixed-LFixed
	    mv finalTimeMethod4 withPoolWhenXYIntOrLFixed-LFixed
            python perf.py -c 1 --legend --x-limit=12 withoutPoolWhenXYInt-LFixed withPoolWhenXYInt-LFixed withoutPoolWhenXYIntOrLFixed-LFixed withPoolWhenXYIntOrLFixed-LFixed> figure4b.eps
	fi

        rm with*

        cd ..
    fi
done

#gathering results to see which instances should be solved for
# investigation of the impact of heuristics
for i in {1..11}
do
    mkdir tmp
    name=$RESULTSPATH'/method'$i'Output/dataIBLP-DEN'
    cp $name'/'*'.out' tmp

    name=$RESULTSPATH'/method'$i'Output/dataIBLP-FIS'
    cp $name'/'*'.out' tmp

    name=$RESULTSPATH'/method'$i'Output/dataMIBLP-XU'
    cp $name'/'*'.out' tmp

    cd tmp

    cat *.out > outMethod$i

    cd ..

    cp tmp/outMethod$i parse

    rm -r -f tmp

    cd parse

    awk -f Timeparse-mibs.awk outMethod$i

    if [ $i == 11 ]
    then
	awk -f Instanceparse_mibs.awk outMethod$i
	cp instanceList.summary ../refinement
	rm instanceList.summary
    fi
    
    mv time.summary timeMethod$i.summary

    cp timeMethod$i.summary ../refinement

    rm timeMethod$i.summary

    rm outMethod$i

    cd ..
done

cd refinement

python refine.py 171 instanceList.summary

rm instanceList.summary

for i in {1..11}
do
    rm timeMethod$i.summary
    rm finalTimeMethod$i
done

cp finalInstanceList ../usefulFiles

rm -f finalInstanceList

cd ..

#Table 2
rm -r -f tmp

for i in 4 7 8 9
do
    echo method$i
    mkdir tmp
    name=$RESULTSPATH'/method'$i'Output/dataIBLP-DEN'
    cp $name'/'*'.out' tmp

    name=$RESULTSPATH'/method'$i'Output/dataIBLP-FIS'
    cp $name'/'*'.out' tmp

    name=$RESULTSPATH'/method'$i'Output/dataMIBLP-XU'
    cp $name'/'*'.out' tmp

    cd tmp

    cat *.out > outMethod$i

    cd ..

    cp tmp/outMethod$i parse

    rm -r -f tmp

    cd parse

    awk -f Timeparse-mibs.awk outMethod$i
    awk -f Instanceparse_mibs.awk outMethod$i
    awk -f VFNumparse_mibs.awk outMethod$i
    awk -f UBNumparse_mibs.awk outMethod$i
    awk -f VFTimeparse_mibs.awk outMethod$i
    awk -f UBTimeparse_mibs.awk outMethod$i

    cp time.summary ../table
    cp instanceList.summary ../table
    cp VFNum.summary ../table
    cp VFTime.summary ../table
    cp UBNum.summary ../table
    cp UBTime.summary ../table

    rm -f outMethod$i
    rm -f time.summary
    rm -f instanceList.summary
    rm -f VFNum.summary
    rm -f VFTime.summary
    rm -f UBNum.summary
    rm -f UBTime.summary

    cd ../table

    echo run python

    python mipTime.py 171 $i $USEFULPATH

    cd ..

done

cd table

python makeTable2.py
    
rm -f *.summary

rm -f tmpOutputMethod*




