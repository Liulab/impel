########################################################
###Basic process of DNA methylation and RNA sequencing #
########################################################

##################
##DNA methylation#
##################
#algnment bash
vim bismark_tmp.sh
#
path=/rawData/chennan/mnt/S1435-L4-L7-56004-kunyuan
bismark --parallel 7 -p 4 --non_directional --genome /data1/ref/hg38/GATK_bundle/ \
-1 $path/HXH-FR180601001-D01-L03-M_S1_L006_R1_001.fastq.gz \
-2 $path/HXH-FR180601001-D01-L03-M_S1_L006_R2_001.fastq.gz 
#
for i in {1..80} ; do sed "s/HXH-FR180601001-D01-L03-M_S1/HXH-FR18060100${i}-D01-L03-M_S${i}/g" bismark_tmp.sh >> bismark_r1.sh; done
nohup sh bismark_r1.sh > bismark_r1.log 2>&1 &

#sort bash
vim sort_tmp1.sh
#
samtools sort -n -@ 8 -o HXH-FR180601001-D01-L03-M_S1.sort.bam ../1.bismark/HXH-FR180601001-D01-L03-M_S1_L006_R1_001_bismark_bt2_pe.bam 
#
for i in {1..80} ; do sed "s/HXH-FR180601001-D01-L03-M_S1/HXH-FR18060100${i}-D01-L03-M_S${i}/g" sort_tmp1.sh >> sort_r1.sh; done
nohup sh sort_r1.sh > sort_r1.log 2>&1 &

#deduplicated bash
vim deduplicated_tmp1.sh
#
deduplicate_bismark -p --bam  HXH-FR180601001-D01-L03-M_S1.sort.bam
#
nohup sh deduplicated_tmp1.sh > deduplicated_tmp1.log 2>&1 &
for i in {1..80} ; do sed "s/HXH-FR180601001-D01-L03-M_S1/HXH-FR18060100${i}-D01-L03-M_S${i}/g" deduplicated_tmp1.sh >> deduplicated_r1.sh; done
nohup sh deduplicated_r1.sh > deduplicated_r1.log 2>&1 &

#bismark_methylation_extractor bash
vim extract_tmp.sh
###########################
bismark_methylation_extractor -p --ignore_r2 2 \
--comprehensive --parallel 10 --scaffolds \
--gzip --bedGraph --zero_based --buffer_size 100G \
--genome_folder /data1/ref/hg38/GATK_bundle/ ../1.bismark/HXH-FR180601001-D01-L03-M_S1_L006_R1_001_bismark_bt2_pe.bam
###########################
for i in {1..80} ; do sed "s/HXH-FR180601001-D01-L03-M_S1/HXH-FR18060100${i}-D01-L03-M_S${i}/g" extract_tmp.sh >> extract.sh; done
nohup sh extract0.sh > extract0.log 2>&1 &

#################
##RNA sequecing #
#################
#fastQC
fastqc -o QC -t 20 ${i}.clean.fq.gz 

#hisat2
path=/rawData/transcript/X101SC19091569-Z01-J008/2.cleandata
hisat2 -p 8 --dta -t --rna-strandness RF -x /data1/ref/hg38/hisatIndex/genome/genome \
-1 ${i}_1.clean.fq.gz \
-2 ${i}_2.clean.fq.gz \
-S ${i}.sam
samtools view -@ 8 -bS ${i}.sam > ${i}.bam
samtools sort -@ 8 -o ${i}.sort.bam ${i}.bam

#stringtie
stringtie ${i}.sort.bam -p 8 -G /data1/ref/hg38/gencode.v28.annotation.gtf -l ${i} > ${i}.gtf
stringtie --merge -p 8 -G /data1/ref/hg38/gencode.v28.annotation.gtf mergelist.txt >stringtie_merged.gtf
stringtie ${i}.sort.bam -e -p 8 -G stringtie_merged.gtf -l ${i} -A ${i}.tab > ${i}.gtf
 