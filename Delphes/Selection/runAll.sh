
rm -r /afs/cern.ch/work/a/ariostas/public/CP-Higgs_Samples_small/*

while read line #loop over lines in ${conf_file}
do
	array=($line)
	if [ "${array[0]}" != "#" ]; then 

            root -b -l -q "select.C+(\"${array[0]}\", ${array[1]})"

	fi
done < xsec.txt
