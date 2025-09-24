########break point, use CSCC3 as a example
cd /data2/home/lijie/CSCC/cell_ranger/CSCC3/outs

samtools view -q 5 possorted_bam.bam|cut -f 3,4,5,7,8,9,17,18,19,20 > possorted_bam_mate_length.txt
awk -F "\t" '$6>500'  possorted_bam_mate_length.txt > possorted_bam_mate_longer.txt
grep -v '=' possorted_bam_mate_length.txt > possorted_bam_mate_differ_chr.txt
rm possorted_bam_mate_length.txt

samtools view -q 5 -f 1 possorted_bam.bam|cut -f 1,9,17,18,19,20 > possorted_bam_mate_strand.txt
sort -k 1 possorted_bam_mate_strand.txt > temp.txt

awk '$2!=0' temp.txt|cut -f 1|uniq -d > temp_paired_end.txt
sort -k 1b,1 temp_paired_end.txt >a.txt
sort -k 1b,1 temp.txt >b.txt
join a.txt b.txt > same_strand_1.txt

awk '{a[$1]+=$2}END{for(i in a){print i,a[i]}}' same_strand_1.txt > same_strand_2.txt
sort -k 1 same_strand_2.txt > same_strand_2_sort.txt
join same_strand_1.txt same_strand_2_sort.txt|awk '$7!=0' > same_strand.txt
rm possorted_bam_mate_strand.txt temp.txt temp_paired_end.txt a.txt b.txt same_strand_1.txt same_strand_2.txt same_strand_2_sort.txt
