#!/bin/bash -l
#SBATCH -t 1:00

bamlist=../Raw_data/gbamlist.txt
out=../Raw_data/GUA_bamlist.txt
info=../Raw_data/plate_info.txt

grep 'GUA' $info > temp.txt
#grep 'ENA' $info >> temp.txt

rm $out

wc=$(wc -l temp.txt | awk '{print $1}')

x=1
while [ $x -le $wc ]
do
        string="sed -n ${x}p temp.txt"
        str=$($string)

        var=$(echo $str | awk -F"\t" '{print $1, $2}')
        set -- $var
        c1=$1
        c2=$2
	c2="${c2}_"

	grep $c1 ${bamlist} | grep $c2 >> ${out}

        x=$(( $x + 1 ))

done

rm temp.txt
