#!/bin/sh
#PBS -N SquareRootTest
#PBS -q batch

# script to test sqrt.py on quark


echo " "
echo    "//////////////////////////////////////////////////////////////////////"
echo    "//"
echo -n "// Starting nde simulation run on ";hostname
echo    "//"
echo    "////////////////////////////////////////////////////////////////////// "
echo " "
echo It is now `date` at the start.
echo " "

# define subdirectories needed to run code.
WORKAREA=/tmp/adil
GMNAREA=/home/aa9pb/adil

echo "The work area is ${WORKAREA}/${CASE}. The GMNAREA is ${GMNAREA}."
# templates for inputs inserted with sed in submit_eodSim.pl
case=;
number=;
source_dir=;
remote_dir=;

# print out the inputs inserted by the submit script (submit_gmnsim.pl)
echo "Parameters:"
echo "   case="$case
echo "   number="$number
echo "   source_dir="$source_dir
echo "   remote_dir="$remote_dir

# Which node am I using?
echo I was given the following slots to run:
cat ${PBS_NODEFILE}
echo " ";

# delete old directories on remote nodes and create them anew for this run.
rm -r ${WORKAREA}/$case
mkdir -p ${WORKAREA}/$case

# move the appropriate input data to the remote node.
cp $GMNAREA/sqrt.py  ${WORKAREA}/$case
echo "Local files:";
echo " ";
cd ${WORKAREA}/$case
pwd
ls -la

echo " "
echo -n "=============== Start of nde simulation run on ";date
echo " "
startTime=$(date +"%s")

# run gemc. creates e_pi_n.ev for reconstruction input.

echo ;echo "Simulate events with gemc. ******************************"; echo
bash
#scl enable python33 'python sqrt.py $number'
python sqrt.py $number
echo ;echo "Done with gemc!. *********************************"; echo
echo " ";

echo It is now `date`.


# Copy output files back to the master.
if [ "$case" -lt 10 ]
    then
    echo "The work area is ${WORKAREA}/$case. The run area is ${GMNAREA}."
    ls -la
    #cp $WORKAREA/$case/e_pi_n.ev $GMNAREA/e_pi_n000$case.ev
    #cp $WORKAREA/$case/e_pi_n_Rec.0.evio $GMNAREA/e_pi_n_Rec_000$case.ev
    cp $WORKAREA/$case $GMNAREA/e_pi_n000$case.ev
    elif  [ $case -lt 100 ]
    then
    echo "The work area is ${WORKAREA}/$case."
    ls -la
    cp $WORKAREA/$case $GMNAREA/e_pi_n000$case.ev

    elif [ $case -lt 1000 ]
    then
    echo "The work area is ${WORKAREA}/$case."
    ls -la
    cp $WORKAREA/$case $GMNAREA/e_pi_n000$case.ev

    else
    echo "The work area is ${WORKAREA}/$case."
    ls -la
    cp $WORKAREA/$case $GMNAREA/e_pi_n000$case.ev

fi


echo " ";
echo It is now `date` at the end of the run.
echo " "
stopTime=$(date +"%s")
diff=$(($stopTime-$startTime))

echo "Case $case with events took a total of $(($diff / 60)) minutes and $(($diff % 60)) seconds."

