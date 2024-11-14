#!/bin/bash

# Define basic settings
script="stress01.py"
base=`echo ${script} | rev | cut -f2- -d. | rev`
workdir="/tmp/tlankau"
startdir=`pwd`
homedir=`pwd`
pid=$SLURM_JOB_ID
if [ -z ${pid} ]
  then
    # echo "Create pid for shell"
    pid=`echo $$`
fi
slout="slurm-${pid}.out"
homedir="${homedir}/Test-${pid}"

# Report basic settings
echo -n "Runnig ${script} on "
hostname
echo "  name base: ${base}"
echo "  startdir : ${startdir}"
echo "  workdir  : ${workdir}"
echo "  homedir  : ${homedir}"
echo "  pid      : ${pid}"
echo "  slurm out: ${slout}"
echo

# Test for workdir and create if it doesn't exists
if test -d ${workdir}
then
  echo "${workdir} found"
  echo "Delete old version of ${workdir}"
  rm -rf ${workdir}
  echo "Build a new work directory"
  mkdir -p ${workdir}
else
  echo "${workdir} NOT found -> create"
  mkdir -p ${workdir}
fi

# Copy files and run job
echo "Copy ${script} to ${workdir}"
echo
cp ${script} ${workdir}
cd ${workdir}
fdate=$(TZ="Asia/Taipei" date '+%d-%m-%Y %H:%M:%S')
echo "$fdate Job starts"
# Avoid problems with the CCDC license
. /ceph/sharedfs/work/NTHU_tyang/scripts/csd/activate.sh
python ${script} | tee ${base}.log
fdate=$(TZ="Asia/Taipei" date '+%d-%m-%Y %H:%M:%S')
echo "$fdate Job ends"
echo

# clean up (don't waste space on /tmp)
echo "Move Results from ${workdir} to ${homedir}"
if test -d ${homedir}
then
  echo "${homedir} found -> delete old, create new"
  rm -rf ${homedir}
else
  echo "${homedir} NOT found -> create new (by mv command)"
fi
echo "Move ${workdir} to ${homedir}"
mv ${workdir} ${homedir}
echo "Copy output from the SLURM demon to ${workdir}"
cd ${startdir}
if test -f ${slout}
then
  mv ${slout} ${homedir}
fi
echo "Free space on \tmp and delete ${workdir}"
rm -rf ${workdir}

exit 0