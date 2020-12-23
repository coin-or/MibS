#!/usr/bin/env bash

# modified from run.sh

cd ../..

PWD=`pwd`

build_dir=$PWD/../build-mibs
mibs_dir=$PWD

#parse_args "$@"

MibSPATH=$mibs_dir
EXECUTABLE=$build_dir/bin
RESULTSPATH=$mibs_dir/output
INSTPATH=$mibs_dir/testSets/BilevelLib/general
SCRIPTPATH=$mibs_dir/scripts
PARAMPATH=$mibs_dir/testSets/BilevelLib

cd $MibSPATH

#rm -r -f testSets
#rm -r -f output

#creating the directories for stroing the outputs
#mkdir -p output

cd $RESULTSPATH

#for i in {1..14}
# gaps="10 25 50"
# gaps="20 30 40 60 70 80 90"
gaps="30"
for i in $gaps
do
    cd $RESULTSPATH
    
    # name=BR${i}Output
    name=BR${i}Output
    mkdir -p $name
    cd $name

    name=dataIBLP-DEN
    mkdir -p $name

    # name=dataIBLP-FIS
    # mkdir -p $name

    # name=dataMIBLP-XU
    # mkdir -p $name
done
    
#storing the test sets
cd $MibSPATH
# mkdir -p testSets
# cd testSets
#git clone https://github.com/tkralphs/BilevelLib.git

cd $EXECUTABLE

#for i in {1..11}
for i in $gaps
do
    # rm -f mibs.par

    # name=$SCRIPTPATH/parameters/method${i}/mibs.par
    # cp $name $EXECUTABLE

    #running the instances in MIBLP-XU
    # CURRENTINSTPATH=$INSTPATH/MIBLP-XU

    # for file in ${CURRENTINSTPATH}/*.mps
    # do
	# instance_name=`basename ${file%.*}`
	# ./mibs -param mibs.par -Alps_instance $CURRENTINSTPATH/${instance_name}.mps -MibS_auxiliaryInfoFile $CURRENTINSTPATH/${instance_name}.aux -MibS_useIntersectionCut 1 -MibS_intersectionCutType 3 | tee $RESULTSPATH/method${i}Output/dataMIBLP-XU/${instance_name}.out
    # done
	
    #running the instances in IBLP-DEN
    #all params set to default with time limit 3600
    CURRENTINSTPATH=$INSTPATH/RANDOM/RAND_BILEVEL

    for file in ${CURRENTINSTPATH}/*.mps
    do
	instance_name=`basename ${file%.*}`
    # echo `pwd`
    echo ${i} $CURRENTINSTPATH/${instance_name}.mps
	# ./mibs -param mibs.par -Alps_instance $CURRENTINSTPATH/${instance_name}.mps -MibS_auxiliaryInfoFile $CURRENTINSTPATH/${instance_name}.aux | tee $RESULTSPATH/method${i}Output/dataIBLP-DEN/${instance_name}.out
    ./mibs -gp ${i} -param $PARAMPATH/mibs.par -Alps_instance $CURRENTINSTPATH/${instance_name}.mps -MibS_auxiliaryInfoFile $CURRENTINSTPATH/${instance_name}.aux | tee $RESULTSPATH/BR${i}Output/dataIBLP-DEN/${instance_name}.out
    done
    
    #running the instances in IBLP-FIS
    # CURRENTINSTPATH=$INSTPATH/IBLP-FIS

    # for instance_name in lseu-0.100000 lseu-0.900000 p0033-0.100000 p0033-0.500000 p0033-0.900000 p0201-0.900000 p0548-0.100000 p0548-0.500000 p0548-0.900000 p2756-0.100000 p2756-0.500000 p2756-0.900000 seymour-0.100000 seymour-0.500000 seymour-0.900000 stein27-0.100000 stein27-0.500000 stein27-0.900000 stein45-0.100000 stein45-0.500000 stein45-0.900000
    # do
	# ./mibs -param mibs.par -Alps_instance $CURRENTINSTPATH/${instance_name}.mps -MibS_auxiliaryInfoFile $CURRENTINSTPATH/${instance_name}.aux | tee $RESULTSPATH/method${i}Output/dataIBLP-FIS/${instance_name}.out
    # done
	
done

# use python script to analyze output files
cd $SCRIPTPATH/analyze
#$SCRIPTPATH/analyze/gather.sh
