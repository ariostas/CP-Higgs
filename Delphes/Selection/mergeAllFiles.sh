while read line
do
    array=($line)
  
    if [ "${array[0]}" != "#" ]; then 

        count=0

        for file in /afs/cern.ch/work/a/ariostas/public/CP-Higgs_Samples_small/out_${array[0]}*
        do
            file_array[$count]=$file
            let "count+=1"
        done

        hadd -f /afs/cern.ch/work/a/ariostas/public/CP-Higgs_Samples_small/${array[0]}_temp.root ${file_array[@]}

        unset file_array

        #rm /afs/cern.ch/work/a/ariostas/public/CP-Higgs_Samples_small/out_${array[0]}*

        root -l -b -q "cleanUpMergedFiles.C+(\"${array[0]}\")"

        rm /afs/cern.ch/work/a/ariostas/public/CP-Higgs_Samples_small/${array[0]}_temp.root

    fi

done < xsec.txt
