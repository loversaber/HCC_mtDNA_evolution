sortbam=/home/bioinfo/liuqi/tumor/wholexome
#sortbam=/home/bioinfo/liuqi/tumor/wholexome/double2/sorted_bam
samtoolspath=/pub/tool/samtools-1.6/bin/samtools
bwapath=/pub/tool/bwa-0.7.10/bwa
mito2=/home/bioinfo/liuqi/reference/mito_shifted.fa
mito2fai=/home/liuqi/work/reference/mito_shifted.fa.fai
samtoolspath2=/pub/tool/samtools-1.6/bin/samtools
script=/home/bioinfo/liuqi/reference
repair=/home/bioinfo/liuqi/tool/bbmap/repair.sh
parafly=/pub/tool/trinityrnaseq-2.0.6/trinity-plugins/parafly/bin/ParaFly 
#$bwapath index -a bwtsw $mito2
#$samtoolspath faidx $mito2

for i in `find $sortbam -maxdepth 2 -name "*.sorted.bam"`
do
a=${i%%.sorted.bam}
b=${a##*/}
#echo $a
echo $b
#:<<zhushi1
echo "MERGE UNMAPPED"
echo "$samtoolspath index $i" >> p1IndexBam.txt
echo "$samtoolspath view -bh -F 2306 $i -o $b".unmap.bam"" >> p1ExtractBam.txt
echo "$samtoolspath view -bh -f 2 -F 2304 $i MT -o $b".MT.bam"" >> p1ExtractBam.txt
echo "$samtoolspath merge $b".merge.bam" $b".unmap.bam" $b".MT.bam"" >> p1MergeBam.txt
echo "$samtoolspath sort -n -@ 12 -T $b"_tmp" -o $b".merge.sortn.bam" $b".merge.bam"" >> p1SortnBam.txt

#samtools sort -n $b".merged.nodup.bam" $b".merged.nodup.sort"

echo "bedtools bamtofastq -i $b".merge.sortn.bam"  -fq $b".1.fq" -fq2 $b".2.fq"" >> p2Bwa_fq.txt
echo "bwa mem -R '@RG\tID:MT\tSM:'"$b"'\tLB:biancheng\tPL:ILLUMINA' -t 12 $mito2 $b".1.fq" $b".2.fq"|$samtoolspath view -@ 12 -Shu -|$samtoolspath sort -n -@ 12 -o $b".mem.sortn.bam" -" >> p2Bwa_mem.txt

echo "$samtoolspath view -b -f 2 -F 2304 $b".mem.sortn.bam" > $b".mem.p.bam"" >> p3EndView.txt
echo "$samtoolspath sort -@ 12 -o $b".mem.p.sort.bam" $b".mem.p.bam"" >> p3Sort.txt

echo "bash -c 'perl $script/pileup2baseindel.pl -i <($samtoolspath mpileup -B -Q 30 -d 1000000 -L 10000 -f $mito2 $b".mem.p.sort.bam") -bq 30 -prefix $b".mem"'" >>p4TxtGet.txt

done

cycleParafly(){
$parafly -c ${1} -CPU 10 ${1}
#$parafly -c ${1} -CPU 10 ${2}
}

paraflytxt=(p1IndexBam.txt p1ExtractBam.txt p1MergeBam.txt p1SortnBam.txt p2Bwa_fq.txt p2Bwa_mem.txt p3EndView.txt p3Sort.txt p4TxtGet.txt)

for parafil in ${paraflytxt[@]}
do 
cycleParafly $parafil
done

ls *.mem1.txt > fils.txt
python3 /home/bioinfo/liuqi/tumor/wholexome/end/plus500.py fils.txt

for i in `ls *.mem1.txt`
do
a=${i%%.mem1.txt}
b=${a##*/}
echo $b
mitofil=$b"_mito.fa"
$bwapath index -a bwtsw $mitofil
$samtoolspath faidx $mitofil

#echo "bedtools bamtofastq -i $b".merge.sortn.bam"  -fq $b".1.fq" -fq2 $b".2.fq"" >> map_p2Bwa_fq.txt
echo "bwa mem -R '@RG\tID:MT\tSM:'"$b"'\tLB:biancheng\tPL:ILLUMINA' -t 12 $mitofil $b".1.fq" $b".2.fq"|$samtoolspath view -@ 12 -Shu -|$samtoolspath sort -n -@ 12 -o $b"_2.mem.sortn.bam" -" >> map2_Bwamem.txt

echo "$samtoolspath view -b -f 2 -F 2304 $b"_2.mem.sortn.bam" > $b"_2.mem.p.bam"" >> map2_EndView.txt
echo "$samtoolspath sort -@ 12 -o $b"_2.mem.p.sort.bam" $b"_2.mem.p.bam"" >> map2_Sort.txt

echo "bash -c 'perl $script/pileup2baseindel.pl -i <($samtoolspath mpileup -B -Q 30 -d 1000000 -L 10000 -f $mitofil $b"_2.mem.p.sort.bam") -bq 30 -prefix $b"_2.mem"'" >>map2_p4TxtGet.txt

done

paraflytxt2=(map2_Bwamem.txt map2_EndView.txt map2_Sort.txt map2_p4TxtGet.txt)

for parafil in ${paraflytxt2[@]}
do
cycleParafly $parafil
done

ls *_2.mem1.txt > fils_2.txt
python3 /home/bioinfo/liuqi/tumor/wholexome/end/plus500.py fils_2.txt

arr=(0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1)
for((i=0;i<=${#arr[@]}-1;i++));do
(cat $script/analysis2.sas; find . -maxdepth 1 -name "*_plus.txt" |sed 's/\.\/\(.*\).txt/%all(\1,'"${arr[$i]}"',3);/') > t.sas
pwd=`pwd`;sed 's|THEPATH|'"$pwd"'|g' t.sas -i
/pub/tool/SASHome/sas t.sas

array=(302 303 304 305 306 307 308 309 310 311 312 313 314 315 316 513 514 515 516 517 518 519 520 521 522 523 524 525 526 566 567 568 569 570 571 572 573 3106 3107 16181 16182 16183 16184 16185 16186 16187 16188 16189 16190 16191 16192 16193 16194 16519) #only human
for a in `ls *.xls`; do file=${a%.xls}; awk 'BEGIN{while(getline<"'$a'"){if($27=="Smarker"||$27==1||($28==1&&$29<0.05||$29=="<.0001")){n++;list[n]=$0}}}{listb[$1]=1}END{for(x=1;x<=n;x++){split(list[x],arr,"\t");if(!(arr[2] in listb)){print list[x]}}}' RS=' ' <<< ${array[@]} > $file\_filtered.xls; done

if [ ${arr[$i]} != ${arr[-1]} ]; then
mkdir result_${arr[$i]}; mv *.xls result_${arr[$i]}; cp *_plus.txt result_${arr[$i]}; rm -f t.sas t.log t.lst
else
mkdir result_${arr[$i]}; mv *.xls result_${arr[$i]}; cp *_plus.txt result_${arr[$i]}; rm -f t.sas t.log t.lst
fi

done
