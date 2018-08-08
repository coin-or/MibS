#!/usr/bin/env bash

# Parse arguments
function parse_args {
    echo "Script run with the following arguments:"
    for arg in "$@"
    do
	echo $arg
	option=
	option_arg=
	case $arg in
	    *=*)
		option=`expr "x$arg" : 'x\(.*\)=[^=]*'`
		option_arg=`expr "x$arg" : 'x[^=]*=\(.*\)'`
		# with bash, one could also do it in the following way:
		# option=${arg%%=*}    # remove longest suffix matching =*
		# option_arg=${arg#*=} # remove shortest prefix matching *=
		case $option in
		    --build-dir)
			if [ "x$option_arg" != x ]; then
			    case $option_arg in
				[\\/$]* | ?:[\\/]* | NONE | '' )
				    build_dir=$option_arg
				    ;;
				*)
				    build_dir=$PWD/$option_arg
				    ;;
			    esac
			else
			    echo "No path provided for --build-dir"
			    exit 3
			fi
			;;
		    --mibs-dir)
			if [ "x$option_arg" != x ]; then
			    case $option_arg in
				[\\/$]* | ?:[\\/]* | NONE | '' )
				    mibs_dir=$option_arg
				    ;;
				*)
				    mibs_dir=$PWD/$option_arg
				    ;;
			    esac
			else
			    echo "No path provided for --mibs-dir"
			    exit 3
			fi
			;;
		    *)
			echo "Unrecognized command...exiting"
			exit 3
			;;
		esac
		;;
	    *)
		echo "Unrecognized command...exiting"
                exit 3
		;;
	esac
    done
}

cd ../..

PWD=`pwd`

build_dir=$PWD/build
mibs_dir=$PWD

parse_args "$@"

MibSPATH=$mibs_dir
EXECUTABLE=$build_dir/MibS/src
RESULTSPATH=$mibs_dir/output
INSTPATH=$mibs_dir/testSets/BilevelLib/general
SCRIPTPATH=$mibs_dir/scripts

cd $MibSPATH

rm -r -f testSets
rm -r -f output

#creating the directories for stroing the outputs
mkdir -p output

cd $RESULTSPATH

for i in {1..14}
do
    cd $RESULTSPATH
    
    name=method${i}Output
    mkdir -p $name
    cd $name

    name=dataIBLP-DEN
    mkdir -p $name

    name=dataIBLP-FIS
    mkdir -p $name

    name=dataMIBLP-XU
    mkdir -p $name
done
    
#storing the test sets
cd $MibSPATH
mkdir -p testSets
cd testSets
git clone https://github.com/tkralphs/BilevelLib.git

cd $EXECUTABLE

for i in {1..11}
do
    rm -f mibs.par

    name=$SCRIPTPATH/parameters/method${i}/mibs.par
    cp $name $EXECUTABLE

    #running the instances in MIBLP-XU
    CURRENTINSTPATH=$INSTPATH/MIBLP-XU

    for file in ${CURRENTINSTPATH}/*.mps
    do
	instance_name=`basename ${file%.*}`
	./mibs -param mibs.par -Alps_instance $CURRENTINSTPATH/${instance_name}.mps -MibS_auxiliaryInfoFile $CURRENTINSTPATH/${instance_name}.aux -MibS_useIntersectionCut 1 -MibS_intersectionCutType 3 | tee $RESULTSPATH/method${i}Output/dataMIBLP-XU/${instance_name}.out
    done

    #running the instances in IBLP-DEN
    CURRENTINSTPATH=$INSTPATH/RANDOM/RAND_BILEVEL

    for file in ${CURRENTINSTPATH}/*.mps
    do
	instance_name=`basename ${file%.*}`
	./mibs -param mibs.par -Alps_instance $CURRENTINSTPATH/${instance_name}.mps -MibS_auxiliaryInfoFile $CURRENTINSTPATH/${instance_name}.aux | tee $RESULTSPATH/method${i}Output/dataIBLP-DEN/${instance_name}.out
    done

    #running the instances in IBLP-FIS
    CURRENTINSTPATH=$INSTPATH/IBLP-FIS

    for instance_name in lseu-0.100000 lseu-0.900000 p0033-0.100000 p0033-0.500000 p0033-0.900000 p0201-0.900000 p0548-0.100000 p0548-0.500000 p0548-0.900000 p2756-0.100000 p2756-0.500000 p2756-0.900000 seymour-0.100000 seymour-0.500000 seymour-0.900000 stein27-0.100000 stein27-0.500000 stein27-0.900000 stein45-0.100000 stein45-0.500000 stein45-0.900000
    do
	./mibs -param mibs.par -Alps_instance $CURRENTINSTPATH/${instance_name}.mps -MibS_auxiliaryInfoFile $CURRENTINSTPATH/${instance_name}.aux | tee $RESULTSPATH/method${i}Output/dataIBLP-FIS/${instance_name}.out
    done
done

cd $SCRIPTPATH/analyze

$SCRIPTPATH/analyze/gather.sh

#experiments for impact of the heuristics
cd $EXECUTABLE

USEFULPATH=$SCRIPTPATH/analyze/usefulFiles
instanceListDir=$USEFULPATH/finalInstanceList

for i in {12..14}
do
    rm -f mibs.par

    name=$SCRIPTPATH/parameters/method${i}/mibs.par
    cp $name $EXECUTABLE

    #running the instances in MIBLP-XU
    CURRENTINSTPATH=$INSTPATH/MIBLP-XU

    for file in ${CURRENTINSTPATH}/*.mps
    do
	instance_name=`basename ${file%.*}`
	if grep -r "${instance_name}.aux" "$instanceListDir" > /dev/null
	then
	    ./mibs -param mibs.par -Alps_instance $CURRENTINSTPATH/${instance_name}.mps -MibS_auxiliaryInfoFile $CURRENTINSTPATH/${instance_name}.aux -MibS_useIntersectionCut 1 -MibS_intersectionCutType 3 | tee $RESULTSPATH/method${i}Output/dataMIBLP-XU/${instance_name}.out
	fi
    done

    #running the instances in IBLP-DEN
    CURRENTINSTPATH=$INSTPATH/RANDOM/RAND_BILEVEL

    for file in ${CURRENTINSTPATH}/*.mps
    do
	instance_name=`basename ${file%.*}`
	if grep -r "${instance_name}.aux" "$instanceListDir" > /dev/null
	then
	    ./mibs -param mibs.par -Alps_instance $CURRENTINSTPATH/${instance_name}.mps -MibS_auxiliaryInfoFile $CURRENTINSTPATH/${instance_name}.aux | tee $RESULTSPATH/method${i}Output/dataIBLP-DEN/${instance_name}.out
	fi
    done

    #running the instances in IBLP-FIS
    CURRENTINSTPATH=$INSTPATH/IBLP-FIS

    for instance_name in lseu-0.100000 lseu-0.900000 p0033-0.100000 p0033-0.500000 p0033-0.900000 p0201-0.900000 p0548-0.100000 p0548-0.500000 p0548-0.900000 p2756-0.100000 p2756-0.500000 p2756-0.900000 seymour-0.100000 seymour-0.500000 seymour-0.900000 stein27-0.100000 stein27-0.500000 stein27-0.900000 stein45-0.100000 stein45-0.500000 stein45-0.900000
    do
	if grep -r "${instance_name}.aux" "$instanceListDir" > /dev/null
	then	    
	    ./mibs -param mibs.par -Alps_instance $CURRENTINSTPATH/${instance_name}.mps -MibS_auxiliaryInfoFile $CURRENTINSTPATH/${instance_name}.aux | tee $RESULTSPATH/method${i}Output/dataIBLP-FIS/${instance_name}.out
	fi
    done
done

cd $SCRIPTPATH/analyze

$SCRIPTPATH/analyze/gatherHeuristic.sh

cd $SCRIPTPATH/analyze

$SCRIPTPATH/analyze/table/displayRawData.sh
