rm -r LSFJOB_*

while read line #loop over lines in ${conf_file}
do
	array=($line)
	if [ "${array[0]}" != "#" ]; then 

	mkdir /afs/cern.ch/work/a/ariostas/private/CP-Higgs_temp/${array[0]}
	rm -r /afs/cern.ch/work/a/ariostas/private/CP-Higgs_temp/${array[0]}/*
        cd /afs/cern.ch/work/a/ariostas/private/CP-Higgs_temp/${array[0]}
        /afs/cern.ch/project/eos/installation/0.2.41/bin/eos.select ls /store/user/arapyan/higgs_cp_samples/${array[0]}/Bacon | grep root > "files.txt"
        split -l 100 "files.txt"
        rm "files.txt"
        cd -
        for file in /afs/cern.ch/work/a/ariostas/private/CP-Higgs_temp/${array[0]}/x*
        do
	        bsub -q 8nh -W 200 -J ${array[0]} run.sh select.C ${array[0]} $file ${array[1]} ${array[2]}
        done

	fi
done < xsec.txt
