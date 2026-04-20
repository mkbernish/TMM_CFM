#!/bin/csh -f
#  Ninv_stp_Fe.cmd
#
#  SGE job for Ninv_stp_Fe built Sat Jan  7 20:42:14 PST 2012
#
#  The following items pertain to this script
#  Use current working directory
#$ -cwd
#  input           = /dev/null
#  output          = /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe.joblog.$JOB_ID
#$ -o /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe.joblog.$JOB_ID
#  error           = Merged with joblog
#$ -j y
#  The following items pertain to the user program
#  user program    = /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe.m
#  arguments       = 
#  program input   = Specified by user program
#  program output  = Specified by user program
#  Resources requested
# 
#$ -l h_data=1024M,h_rt=24:00:00
#
#  Name of application for log
#$ -v QQAPP=matlab
#  Email address to notify
#$ -M tweber@mail
#  Notify at beginning and end of job
#$ -m bea
#  Job is not rerunable
#$ -r n
#
# Initialization for serial execution
#
  unalias *
  set qqversion = 
  set qqapp     = "matlab serial"
  set qqmtasks  = 1
  set qqidir    = /u/home/campus/tweber/MITgcm/Ninv_Fe
  set qqjob     = Ninv_stp_Fe
  set qqodir    = /u/home/campus/tweber/MITgcm/Ninv_Fe
  cd     /u/home/campus/tweber/MITgcm/Ninv_Fe
  source /u/local/bin/qq.sge/qr.runtime
  if ($status != 0) exit (1)
#
  echo "SGE job for Ninv_stp_Fe built Sat Jan  7 20:42:14 PST 2012"
  echo ""
  echo "  Ninv_stp_Fe directory:"
  echo "    "/u/home/campus/tweber/MITgcm/Ninv_Fe
  echo "  Submitted to SGE:"
  echo "    "$qqsubmit
  echo "  SCRATCH directory:"
  echo "    "$qqscratch
#
  echo ""
  echo "Ninv_stp_Fe started on:   "` hostname -s `
  echo "Ninv_stp_Fe started at:   "` date `
  echo ""
#
# Run the user program
#
  source /u/local/Modules/default/init/modules.csh
  module load matlab
  setenv LM_LICENSE_FILE /u/local/licenses/license.matlab
  set path = ( $path /sbin )
#
  echo "mcc -m -R -nodisplay,-singleCompThread Ninv_stp_Fe.m"
requeue:
  /u/local/apps/matlab/7.11/bin/mcc -m -R -nodisplay,-singleCompThread Ninv_stp_Fe.m
  set rc = $status
  # waiting for bluearc per ppk
  sleep 60
  #
  if( `grep -c 'Maximum number of users'                  /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe.joblog.$JOB_ID` > 0 ||\
      `grep -c 'Licensed number of users already reached' /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe.joblog.$JOB_ID` > 0 ||\
      `grep -c 'License checkout failed'                  /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe.joblog.$JOB_ID` > 0 ||\
      `grep -c 'Could not check out a Compiler license'   /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe.joblog.$JOB_ID` > 0 ) then
    head -n 13 /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe.joblog.$JOB_ID > /u/home/campus/tweber/MITgcm/Ninv_Fe/.$$
    mv /u/home/campus/tweber/MITgcm/Ninv_Fe/.$$ /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe.joblog.$JOB_ID
    echo "------------ waiting for a license. retrying mcc command."
    sleep 90
    goto requeue
  endif # waiting for license
#
  echo ""
  if( $rc != 0 || ! -e /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe ) then
    echo "============================================================"
    echo "matlab.queue mcc step failed with status $rc"
    echo "============================================================"
    #
    echo ============================================================ >> /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe.output.$JOB_ID
    echo matlab.queue mcc step failed with status $rc >> /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe.output.$JOB_ID
    echo ============================================================ >> /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe.output.$JOB_ID
    # stop here if mcc failed
  else
    # execute mcc-compiled executable
    chmod +x /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe
    echo Ninv_stp_Fe  \>\& Ninv_stp_Fe.output.$JOB_ID
    /usr/bin/time /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe  >& /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe.output.$JOB_ID
    set rc = $status
    #
    if( $rc != 0 ) then
      rm -f /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe
      echo "matlab.queue execute mcc-compiled Ninv_stp_Fe step failed with status $rc"
      echo "retrying with matlab executable..."
      #
      echo matlab.queue execute mcc-compiled Ninv_stp_Fe step failed with status $rc >> /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe.output.$JOB_ID
      echo retrying with matlab executable... >> /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe.output.$JOB_ID
      #
      set qqargs = ( -nojvm -nodisplay -nosplash  )
      if( 1 == 1 ) set qqargs = ( -singleCompThread $qqargs )
      #
      # run with matlab
      echo matlab $qqargs -r Ninv_stp_Fe -logfile /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe.output.$JOB_ID
requeue2:
      /usr/bin/time /u/local/apps/matlab/7.11/bin/matlab $qqargs -r Ninv_stp_Fe -logfile /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe.output.$JOB_ID
      #
      # waiting for bluearc per ppk
      sleep 60
      #
      if( `grep -c 'Maximum number of users'                  /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe.joblog.$JOB_ID` > 0 ||\
          `grep -c 'Licensed number of users already reached' /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe.joblog.$JOB_ID` > 0 ||\
          `grep -c 'License checkout failed'                  /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe.joblog.$JOB_ID` > 0 ) then
        head -n 17 /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe.joblog.$JOB_ID > /u/home/campus/tweber/MITgcm/Ninv_Fe/.$$
        mv /u/home/campus/tweber/MITgcm/Ninv_Fe/.$$ /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe.joblog.$JOB_ID
        echo "------------ waiting for a license. retrying matlab executable."
        sleep 90
        goto requeue2
      endif # waiting for license
    endif # execute mcc-executable failed
  endif # mcc step succeeded
#
  echo ""
  echo "Ninv_stp_Fe finished at:  "` date `
#
# Cleanup after serial execution
#
  rm -f /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe.prj
  rm -f /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe_main.c
  rm -f /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe_mcc_component_data.c
  rm -f /u/home/campus/tweber/MITgcm/Ninv_Fe/mccExcludedFiles.log
  rm -f /u/home/campus/tweber/MITgcm/Ninv_Fe/readme.txt
  rm -f /u/home/campus/tweber/MITgcm/Ninv_Fe/run_Ninv_stp_Fe.sh

  source /u/local/bin/qq.sge/qr.runtime
#
  echo "-------- /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe.joblog.$JOB_ID --------" >> /u/local/apps/queue.logs/matlab.log.serial
 if (`wc -l /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe.joblog.$JOB_ID  | awk '{print $1}'` >= 1000) then
        head -50 /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe.joblog.$JOB_ID >> /u/local/apps/queue.logs/matlab.log.serial
        echo " "  >> /u/local/apps/queue.logs/matlab.log.serial
        tail -10 /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe.joblog.$JOB_ID >> /u/local/apps/queue.logs/matlab.log.serial
  else
        cat /u/home/campus/tweber/MITgcm/Ninv_Fe/Ninv_stp_Fe.joblog.$JOB_ID >> /u/local/apps/queue.logs/matlab.log.serial
  endif
  exit (0)
