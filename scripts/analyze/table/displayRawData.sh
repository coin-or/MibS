#!/usr/bin/env bash

cd ../..

PWD=`pwd`

MibSPATH=$PWD
RESULTSPATH=$PWD/output
INSTPATH=$mibs_dir/testSets/BilevelLib/general
SCRIPTPATH=$PWD/scripts
USEFULPATH=$SCRIPTPATH/analyze/usefulFiles
usedInstDir=$USEFULPATH/finalInstanceList

cd $SCRIPTPATH/analyze

rm -r -f tmp

for tableNum in {3..9}
do
    echo Plotting Table ${tableNum}

    if [ $tableNum == 3 ]
    then
	dataSet=(dataIBLP-DEN dataIBLP-FIS)
	methodSet=(1 2 3 4 5)
    elif [ $tableNum == 4 ]
    then
        dataSet=(dataMIBLP-XU)
	methodSet=(1 2 3 4 5)
    elif [ $tableNum == 5 ]
    then
	dataSet=(dataIBLP-DEN dataIBLP-FIS)
	methodSet=(4 7)
    elif [ $tableNum == 6 ]
    then
	dataSet=(dataMIBLP-XU dataIBLP-DEN dataIBLP-FIS)
	methodSet=(4 7)
    elif [ $tableNum == 7 ]
    then
	dataSet=(dataMIBLP-XU dataIBLP-DEN dataIBLP-FIS)
	methodSet=(6 7 9 11)
    elif [ $tableNum == 8 ]
    then
	dataSet=(dataMIBLP-XU dataIBLP-DEN dataIBLP-FIS)
	methodSet=(3 4 8 10)
    elif [ $tableNum == 9 ]
    then
	dataSet=(dataMIBLP-XU dataIBLP-DEN dataIBLP-FIS)
	#15: noHeuristics
	methodSet=(15 12 13 14)
    fi

    
    for testSet in "${dataSet[@]}"
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

	count=0

	for i in "${methodSet[@]}"
	do
	    rm -r -f tmp
	    mkdir tmp
	    
	    if [ $i != 15 ]
	    then
		
		CURRENTOUTPUTPATH=$RESULTSPATH/method${i}Output/$testSet

		if [ $tableNum == 5 ]
		then
		    LISTPATH=$USEFULPATH/r1LeqR2List$appendix
		elif [ $tableNum == 6 ]
		then
		    LISTPATH=$USEFULPATH/r1GeR2List$appendix
		fi

		for file in ${CURRENTOUTPUTPATH}/*.out
		do
		    instance_name=`basename ${file%.*}`
		    if grep -r "${instance_name}.aux" "$usedInstDir" > /dev/null
		    then
			if [ $tableNum == 5 -o  $tableNum == 6 ]
			then
			    if grep -r "${instance_name}.aux" "$LISTPATH" > /dev/null
			    then
				cp $CURRENTOUTPUTPATH/${instance_name}.out tmp
			    fi
			else
			    cp $CURRENTOUTPUTPATH/${instance_name}.out tmp
			fi
		    fi
		done
	    else
		CURRENTOUTPUTPATHLeq=$RESULTSPATH/method4Output/$testSet
		CURRENTOUTPUTPATHGe=$RESULTSPATH/method7Output/$testSet
		LISTPATHLeq=$USEFULPATH/r1LeqR2List$appendix
		LISTPATHGe=$USEFULPATH/r1GeR2List$appendix
		for file in ${CURRENTOUTPUTPATHLeq}/*.out
		do
		    instance_name=`basename ${file%.*}`
		    if grep -r "${instance_name}.aux" "$usedInstDir" > /dev/null
		    then
			if grep -r "${instance_name}.aux" "$LISTPATHLeq" > /dev/null
			then
			    cp $CURRENTOUTPUTPATHLeq/${instance_name}.out tmp
			fi
		    fi
		done

		for file in ${CURRENTOUTPUTPATHGe}/*.out
		do
		    instance_name=`basename ${file%.*}`
		    if grep -r "${instance_name}.aux" "$usedInstDir" > /dev/null
		    then
			if grep -r "${instance_name}.aux" "$LISTPATHGe" > /dev/null
			then
			    cp $CURRENTOUTPUTPATHGe/${instance_name}.out tmp
			fi
		    fi
		done
	    fi

	    cd tmp

	    name=out${appendix}Method${i}
	    cat *.out > $name

	    cd ..

	    cp tmp/$name parse

	    rm -r -f tmp

	    cd parse

	    if [ $count == 0 ]
	    then
		awk -f Instanceparse_mibs.awk $name
		mv instanceList.summary instanceList${appendix}.summary
		cp instanceList${appendix}.summary ../table
		rm -f instanceList${appendix}.summary
	    fi
	    awk -f Timeparse-mibs.awk $name
	    awk -f Nodeparse_mibs.awk $name
	    awk -f Costparse_mibs.awk $name
	    if [ $tableNum == 5 -o  $tableNum == 6 ]
	    then
		if [ $count == 0 ]
		then
		    awk -f IntVarNum-mibs.awk $name
		    mv uLIntVarNum.summary uLIntVarNum${appendix}.summary
		    mv lLIntVarNum.summary lLIntVarNum${appendix}.summary
		    cp uLIntVarNum${appendix}.summary ../table
		    cp lLIntVarNum${appendix}.summary ../table
		    rm -f uLIntVarNum${appendix}.summary
		    rm -f lLIntVarNum${appendix}.summary
		fi
	    fi
	    count=1

	    mv time.summary time${appendix}Method${i}.summary
	    cp time${appendix}Method${i}.summary ../table
	    rm -f time${appendix}Method${i}.summary

	    mv node.summary node${appendix}Method${i}.summary
	    cp node${appendix}Method${i}.summary ../table
	    rm -f node${appendix}Method${i}.summary

	    mv cost.summary cost${appendix}Method${i}.summary
	    cp cost${appendix}Method${i}.summary ../table
	    rm -f cost${appendix}Method${i}.summary

	    rm -f $name

	    cd ..
	done
    done

    cd table

    python makeRawDataTable.py table${tableNum}.csv $INSTPATH

    for testSet in "${dataSet[@]}"
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
	rm -f instanceList${appendix}.summary
	rm -f uLIntVarNum${appendix}.summary
	rm -f lLIntVarNum${appendix}.summary
	for i in "${methodSet[@]}"
	do
	    rm -f time${appendix}Method${i}.summary
	    rm -f node${appendix}Method${i}.summary
	    rm -f cost${appendix}Method${i}.summary
	done
    done
    
		   

    cd ..
    
done
