rm -r LSFJOB_*

rm -r /afs/cern.ch/work/a/ariostas/public/CP-Higgs_Samples_small/*

while read line #loop over lines in ${conf_file}
do
	array=($line)
	if [ "${array[0]}" != "#" ]; then 

        cd /afs/cern.ch/work/a/ariostas/public/CP-Higgs_Samples/
        ls ${array[0]}* | grep root > "../CP-Higgs_Samples_small/${array[0]}.txt"
        cd -
        cd /afs/cern.ch/work/a/ariostas/public/CP-Higgs_Samples_small/
        split -l 1 "${array[0]}.txt" ${array[0]}
        rm "${array[0]}.txt"
        cd -
        for file in /afs/cern.ch/work/a/ariostas/public/CP-Higgs_Samples_small/${array[0]}*
        do
	        bsub -q 8nh -W 200 -J ${array[0]} run.sh select.C ${array[0]} $file ${array[1]} ${array[2]}
        done

	fi
done < xsec.txt
