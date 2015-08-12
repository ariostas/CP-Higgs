while read line
do
  array=($line)
  
   if [ "${array[0]}" != "#" ]; then 

	./mergeSelFiles.sh ${array[0]} /afs/cern.ch/work/a/ariostas/private/CP-Higgs_temp/${array[0]}/ /afs/cern.ch/work/a/ariostas/public/CP-Higgs/   

    fi

done < xsec.txt
