#!/usr/bin/perl
#
# script to simulate the gmn experiment on the Richmond cluster from Advanced
# Clustering. it copies the results back to the master. 
#                                          - gpg (6/10/15)
#
# usage: submit_gmnsim.pl
#
# to generate the cases edit run_gemc_on_node.sh.template

# number of nodes to use
$segments=10;

# where to put stuff
$log_dir="/home/aa9pb/Parallel/logs";

# do housekeeping on the master to get rid of old files from previous runs.
#
print "\n";
use POSIX qw(strftime);
print strftime("Starting nde simulation run on %a, %d %b %Y %H:%M:%S", localtime(time())) , "\n\n";
print "Clean things up from previous runs.\n\n";
#system("rm /home/aa9pb/Parallel/test/*");
system("rm /home/aa9pb/Parallel/logs/*");
system("rm /home/aa9pb/Parallel/inputScripts/*.sh");

# housekeeping on the remote nodes to get rid of old files from previous runs.
print "Clean up old stuff on the cluster.\n";
system("act -g nodes -x physics05 rm -r /tmp/adil");
system("act -g nodes -x physics05 mkdir -p /tmp/adil");
#

# loop over each case and submit.
print "Loop over jobs.\n";
for ($counter=0; $counter <= $segments; $counter++) {

# build the name of log file for this case.
    if ($counter < 10) {
	$logfile ="$log_dir/ndesim_log000$counter";
    } elsif ($counter > 9 && ($counter < 100)) {
	$logfile ="$log_dir/ndesim_log00$counter";
    } elsif (($counter > 99) && ($counter < 1000)) {
	$logfile ="$log_dir/ndesim_log0$counter";
    } elsif (($counter > 999)) {
	$logfile ="$log_dir/nedsim_log$counter";
    }
    # create a file with all the inputs needed to run this job on the remote node.                                                                                                            
    system("sed '24s%case=%case=$counter%' run_Parallel.sh.template > step1.sh");
    system("sed '25s%segments=%segments=$segments%' step1.sh > step2.sh");
    system("sed '26s%source_dir=%source_dir=/home/aa9pb/Parallel%' step2.sh > step3.sh");
    system("sed '27s%remote_dir=%remote_dir=/tmp/adil/$counter%' step3.sh > run_Parallel.sh");
    # label this input file with the case number and store it with all the others. This prevents the input script
    # from getting overwritten by later submissions.
    system("cp run_Parallel.sh ./run_Parallel$counter.sh");
#
# run ndesim script to simulate clas12.  ************************************************************************
#
    #print "qsub -l nodes=1:ppn=1 -l walltime=32:00:00 -l mem=8GB -j oe -o $logfile inputScripts/#run_ndesim_on_node$counter.sh \n";
    system("qsub -l nodes=1:ppn=24 -l walltime=72:00:00 -l mem=18GB -j oe -o $logfile ~/Parallel/run_Parallel$counter.sh \n");
    if ($counter <=  ($njobs-1)) {
	print "Job $counter submitted. Sleep for 1 second.\n\n";
	sleep 0.5;
    }
} # end of for loop.

print "\nEnd of submissions.\n";
