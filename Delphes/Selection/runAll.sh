rm -r LSFJOB_*

rm -r /afs/cern.ch/work/a/ariostas/public/CP-Higgs_Samples_small/*

while read line #loop over lines in ${conf_file}
do
	array=($line)
	if [ "${array[0]}" != "#" ]; then 

	        ./run.sh select.C ${array[0]} ${array[1]} ${array[2]}

	fi
done < xsec.txt
