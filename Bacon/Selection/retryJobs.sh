while read line
do
	array=($line)
	if [ "${array[0]}" != "#" ]; then 

        for file in /afs/cern.ch/work/a/ariostas/private/CP-Higgs_temp/${array[0]}/x*
        do
	        bsub -q 8nh -W 200 -J ${array[0]} run.sh select.C ${array[0]} $file ${array[1]} ${array[2]}
        done

	fi
done < xsec.txt