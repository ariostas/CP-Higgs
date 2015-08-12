
 root_script=$1
      sample=$2
   conf_file=$3
   cross_sec=$4
     eosflag=$5
      

cd /afs/cern.ch/user/a/ariostas/CP-Higgs/Delphes/Selection

h=`basename $0`
echo "Script:    $h"
echo "Arguments: $*"

# some basic printing
echo " "; echo "${h}: Show who and where we are";
echo " "
echo " user executing: "`id`;
echo " running on    : "`hostname`;
echo " executing in  : "`pwd`;
echo " submitted from: $HOSTNAME";
echo " ";

# initialize the CMSSW environment
echo " "; echo "${h}: Initialize CMSSW (in $CMSSW_BASE)"; echo " "
workDir=`pwd`
cd   $CMSSW_BASE
eval `scram runtime -sh`
cd -

while read line
do
  array=($line)
    if [ "${array[0]}" != "#" ]; then 
	root -b -l -q "${root_script}+(\"${sample}\", \"${array[0]}\", ${cross_sec}, ${eosflag})"

    fi
done < ${conf_file}

rm ${conf_file}

# get the return code from the root job
status=`echo $?`
echo "${h}: Status - $status"

exit $status
