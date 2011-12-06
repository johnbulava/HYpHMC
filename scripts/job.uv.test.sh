# !/bin/bash
#PBS -N PPHMCAlpha
#PBS -j oe
#PBS -l walltime=00:30:00
#PBS -l nodes=1:ppn=16:uv
#PBS -l feature=uv
#PBS -q testq

JOB_SCRIPT_NAME=job.uv.sh

cd $PBS_O_WORKDIR
EXEC=pHMC
PATH=$PATH:/home/b/bepkalla/scripts

# LOAD MODULES
module load gcc/4.3.4
module load intel.compiler/11.0.083
#module load mvapich2/1.4.0-gcc
#module load mpt

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sw/numerics/fftw-3.2.2/gcc-sse4.1/lib:/sw/local/lib:/sw/numerics/gsl-gcc/1.13/lib:/home/b/bepkalla/libs

OLD_CHECK=$(checkDSimulation.py)
echo "$OLD_CHECK"
PRE_STATUS=$(echo "$OLD_CHECK" | grep -oE "STATUS=[0-9,]+" | sed -e "s/STATUS=//")

echo "START: " `date`
START_TIME_STAMP=$(date "+%s")

./${EXEC}

END_TIME_STAMP=$(date "+%s")
echo "END: " `date`

CHECK=$(checkDSimulation.py)
echo "$CHECK"

DELTA_TS=$(echo "$END_TIME_STAMP - $START_TIME_STAMP" | bc)
DELTA_TH=$(echo "scale=1;${DELTA_TS}/3600" | bc)
echo "TOTAL RUN TIME: $DELTA_TH HOURS"

POST_STATUS=$(echo "$CHECK" | grep -oE "STATUS=[0-9,]+" | sed -e "s/STATUS=//")
STATUS_MAX_INDEX=$(echo "$POST_STATUS" | grep -o "," | wc -l)
STATUS_COUNT=$(echo "$STATUS_MAX_INDEX + 1" | bc)
IDX_ARR=$(seq 1 "$STATUS_COUNT")
PERFORMANCE_STR=$(printf "PERFORMANCE24h=")
for i in ${IDX_ARR}; do
   OLD_STEPS=$(echo "$PRE_STATUS" | cut -d "," -f${i})
   NEW_STEPS=$(echo "$POST_STATUS" | cut -d "," -f${i})

   TMP=$(echo "scale=0;24*60*60*(${NEW_STEPS} - ${OLD_STEPS})/${DELTA_TS}" |bc)
   if [ "$i" -eq 1 ]; then
      PERFORMANCE_STR=$(printf "%s%d:%1.0f" ${PERFORMANCE_STR} $i ${TMP})
   else
      PERFORMANCE_STR=$(printf "%s,%d:%1.0f" ${PERFORMANCE_STR} $i ${TMP})
   fi
done
echo "$PERFORMANCE_STR"

Q_FLAG=1
if [ "$Q_FLAG" -eq 1 ]; then
    PROCEED_FLAG=$(echo "$CHECK" | grep -oE "PROCEED=[0-9]+" | grep -oE "[0-9]")
    if [ "$PROCEED_FLAG" -eq 1 -a "$DELTA_TH" -gt 1 ]; then
	echo "Resubmitting job."
	msub -A bep00024 ${JOB_SCRIPT_NAME}
    else
	echo "AT LEAST ONE SUB RUN CANNOT PROCEED"
    fi
fi
