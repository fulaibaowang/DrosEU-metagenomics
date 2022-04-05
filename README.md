# DrosEU-metagenomics
The bioinformatics pipeline for the metagenomics analyses of the DrosEU (https://droseu.net/) data from 2014-2016.
2014 raw sequences data are available at https://www.ncbi.nlm.nih.gov/search/all/?term=PRJNA388788

## A) QC, trim, map 

### 1) run FASTQC to find the overrepresented sequences (adapters)

```bash
fastqc -t 12  *.gz
```


### 2) use R to merge all overrepresented sequences into a single file

```R
#https://www.biostars.org/p/321827/
#https://www.biostars.org/p/366687/
#install.packages('fastqcr')
library(fastqcr)

# Aggregating Multiple FastQC Reports into a Data Frame 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Demo QC directory containing zipped FASTQC reports
# copy folder under fastqcr pacakage folder
qc.dir <- system.file("droseu_output/2014", package = "fastqcr")
list.files(qc.dir)
qc <- qc_aggregate(qc.dir)
qc

false_overseq=qc[qc$module=='Overrepresented sequences'&qc$status!='PASS',]$sample
# "Dme7_DSFP0597-w_HVNKCCCXX_L7_2.fq"
# "Dme7_DSFP0597-w_HVNKCCCXX_L7_2_fastqc.zip"
# "Dme7_DSFP0597-w_HVNKCCCXX_L7_2.fastqc.zip"
# paste(substr(false_overseq,1,nchar(false_overseq)-3),"_fastqc.zip",sep='')

#library("magrittr")
false_overseq2=c()
for (i in false_overseq) {
  #print(i)}
  isubstr=substr(i,nchar(i)-5,nchar(i))
  if (isubstr=='.fq.gz') {
    i=paste(substr(i,1,nchar(i)-6),"_fastqc.zip",sep='')
  }
  else {i=paste(i,"_fastqc.zip",sep='')}
  false_overseq2=append(false_overseq2,i)
}
false_overseq2

Sequence=c()
Source=c()
percentage=c()

for (i in false_overseq2) {
  qc.file <- system.file("droseu_output/2014", i,  package = "fastqcr")
  #print(qc.file)
  a=fastqcr::qc_read(qc.file)$overrepresented_sequences$Sequence
  b=fastqcr::qc_read(qc.file)$overrepresented_sequences$`Possible Source`
  c=fastqcr::qc_read(qc.file)$overrepresented_sequences$Percentage
  Sequence=append(Sequence,a)
  Source=append(Source,b)
  percentage=append(percentage,c)
      }
#fastqcr::qc_read(qc.file)$overrepresented_sequences %>%
#  dplyr::mutate(name=paste(">",1:n(),"-",Count,sep=""),fa=paste(name,Sequence,sep="\n")) %>%
#  dplyr::pull(fa) %>% 
#  readr::write_lines("overrepresented.fa")


df= cbind(Sequence,Source)
df
unique(df)
write.csv(unique(df),"DrosEU_adapters.fa.csv")
```

### 3) trim with BBduk (remove adapters, minimum length > 75bp, BQ > 18

The basic command looks like this:

```bash
bbduk.sh -Xmx1g in1=xx_R1.fastq.gz_q18.fq.gz in2=xx_R2.fastq.gz_q18.fq.gz out1=xx_R1.fastq.gz.out out2=10_R2.fastq.gz.out ref=DrosEU_adapters.fa ktrim=r k=23 mink=11 hdist=1 overwrite=t tbo=t tpe=t minlength=75 qtrim=rl trimq=18

```

I generate a .sh file for each pair of fastq files.
```bash
unset listA
unset listB
listA=(*_R1.*fq.gz)
listB=(*_R2.*fq.gz)
n=${#listA[@]}
for i in $(seq $n);do echo -e '#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --cpus-per-task=4\n#SBATCH --time=00:39:59\n#SBATCH --output=cut_'${listA[$i]}'.txt\ncd /pfs/work7/workspace/scratch/fr_yw50-restore_DrosEU-0/2014/qtrim_reads/\n/pfs/work7/workspace/scratch/fr_yw50-restore_DrosEU-0/bbmap/bbduk.sh -Xmx1g in1='${listA[$i]}' in2='${listB[$i]}' outu1='${listA[$i]}'_cleanr.fq.gz outu2='${listB[$i]}'_cleanr.fq.gz ref=DrosEU_adapters2014.fa ktrim=r k=23 mink=11 hdist=1 overwrite=t tbo=t tpe=t minlength=75\n/pfs/work7/workspace/scratch/fr_yw50-restore_DrosEU-0/bbmap/bbduk.sh -Xmx1g in1='${listA[$i]}'_cleanr.fq.gz in2='${listB[$i]}'_cleanr.fq.gz outu1='${listA[$i]}'_cleanrl.fq.gz outu2='${listB[$i]}'_cleanrl.fq.gz ref=DrosEU_adapters.fa ktrim=r k=23 mink=11 hdist=1 overwrite=t tbo=t tpe=t minlength=75' > ${listA[$i]}_bbduk.sh1 ;done
```

And then I could submit all .sh files in the cluster(Slurm) simultaneously.

```bash
for f in *sh1; do sbatch -p single $f;done 
```

### 4) mapping on fly genome

Mapping on D.melanogaster and D.simulans genome (flies.fna.gz). We only need unmapped reads (microbe reads) as output.
The basic command looks like this:

```bash
cd to/folder/after/trimming
listA=(*_R1.*fq.gz)
listB=(*_R2.*fq.gz)
bbmap.sh -Xmx30g ref=flies.fna.gz in1='${listA[$i]}' in2='${listB[$i]}' outu1='${listA[$i]}'_unmapped.fq.gz outu2='${listB[$i]}'_unmapped.fq.gz
```

In cluster I again generate .sh file for each pair of fastq files.

```bash
unset listA
unset listB
listA=(*_R1.*cl*.fq.gz)
listB=(*_R2.*cl*.fq.gz)
echo ${#listA[*]} 
n=${#listA[@]}
for i in $(seq $n);do echo -e '#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --cpus-per-task=16\n#SBATCH --time=7:59:59\n#SBATCH --output=bbmap_run_'${listA[$i]}'.txt
date\ncd /pfs/work7/workspace/scratch/fr_yw50-restore_DrosEU-0/2014/qtrim_reads
/pfs/work7/workspace/scratch/fr_yw50-restore_DrosEU-0/bbmap/bbmap.sh -Xmx30g ref=flies.fna.gz in1='${listA[$i]}' in2='${listB[$i]}' outu1='${listA[$i]}'_unmapped.fq.gz outu2='${listB[$i]}'_unmapped.fq.gz\ndate' > ${listA[$i]}.sh2 ;done

for f in 1*sh2; do sbatch -p single $f;done 
```
## B) short reads analysis with Diamond and Megan

### 1) reads annotation with Diamond blastx
The input files are reads that are not mapped on fly genome.

The basic command looks like this:
```bash
#build nr.dmnd from ncbi database, this only needs to be done once
diamond makedb --in nr.gz --db nr

diamond blastx --query 'xx.unmapped.fq.gz' --db nr.dmnd --daa xx.unmapped.fq.gz.daa
```

In cluster I again generate .sh file for each fastq file and submit them simultaneously.

```bash
unset listA
listA=(*unmapped.fq.gz)
n=${#listA[@]}
for i in $(seq $n); do echo -e '#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --cpus-per-task=16\n#SBATCH --time=5:19:59\n#SBATCH --output=diamond_'${listA[$i]}'.txt\ncd /pfs/work7/workspace/scratch/fr_yw50-restore_DrosEU-0/2014\ndate
/home/fr/fr_fr/fr_yw50/diamond/diamond blastx --query '${listA[$i]}' --db /home/fr/fr_fr/fr_yw50/diamond/nr.dmnd --daa daa/'${listA[$i]}'.daa\ndate' > ${listA[$i]}_daa.sh4 ;done
for f in sh4; do sbatch -p single $f;done  
```

### 2) convert daa files to rma6 files
The input files are R1-R2 pairs of daa files.
The basic command looks like this:
```bash
daa2rma -p true -ps 1 -i daa/1_R1.unmapped.fq.gz.daa daa/1_R2.unmapped.fq.gz.daa -o rma/ -a2t /diamond/prot_acc2tax-Jul2019X1.abin -a2eggnog /acc2eggnog-Jul2019X.abin
```
In cluster I again generate .sh file for each fastq file and submit them simultaneously.

```bash
unset listA
unset listB
listA=(*_R1.*.daa)
listB=(*_R2.*.daa)
n=${#listA[@]}
for i in $(seq $n); do echo -e '#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --cpus-per-task=16\n#SBATCH --time=10:19:59\n#SBATCH --output=megan_'${listA[$i]}'.txt\ncd /pfs/work7/workspace/scratch/fr_yw50-restore_DrosEU-0/2014/daa\ndate
/home/fr/fr_fr/fr_yw50/diamond/megan/tools/daa2rma -p true -ps 1 -i '${listA[$i]}' '${listB[$i]}' -o /pfs/work7/workspace/scratch/fr_yw50-restore_DrosEU-0/2014/rma/ -a2t /home/fr/fr_fr/fr_yw50/diamond/prot_acc2tax-Jul2019X1.abin -a2eggnog /home/fr/fr_fr/fr_yw50/diamond/acc2eggnog-Jul2019X.abin -a2interpro2go /home/fr/fr_fr/fr_yw50/diamond/acc2interpro-Jul2019X.abin -a2seed /home/fr/fr_fr/fr_yw50/diamond/acc2seed-May2015XX.abin\ndate' > ${listA[$i]}_rma.sh6 ;done
```
The ouput rma6 files can be analyzed in MEGAN6.

## C) Mapping short reads on assembly 

### 1) mapping on assembly with bowtie2
Mapping the unmapped (microbe) reads on the assembly, the assembly file final.contigs.fa was generated from Megahit.

The basic command looks like this:
```bash
#build library for the assembly, this only needs to be done once
bowtie2-build ./final.contigs.fa mapping/contig

bowtie2 --threads 8 -x mapping/contig -1 R1_clean_unmapped.fq -2 R2_clean_unmapped.fq -S Sample.sam
```

In cluster I again generate .sh file for each pair of fastq files.

```bash
unset listA
unset listB
listA=(*R1*.gz_unmapped.fq)
listB=(*R2*.gz_unmapped.fq)
n=${#listA[@]}
for i in $(seq $n); do bowtie2 --threads 16 -x mapping/contig -1 ${listA[$i]} -2 ${listB[$i]} -S Sample_${listA[$i]:0:2}.sam; done  >  ${listA[$i]}.sh3; done

for f in 1*sh2; do sbatch -p single $f;done 
```

### 5)convert sam to bam, sort and index bam

```bash
for x in *.sam;
do samtools view -F 4 -bS $x > $x-RAW.bam;
done

#sort
for x in *RAW.bam;
do samtools sort $x > $x-sort.bam;
done

#index
for x in *RAW.bam;
do samtools sort $x > $x-sort.bam;
done

#calculate coverage
for x in *.sam;
do /pfs/work7/workspace/scratch/fr_yw50-restore_DrosEU-0/bbmap/pileup.sh -Xmx8g in=$x out=$x.cov.txt;
done
```

## D) Binning approach with Maxbin2, Metabat2, Concoct, Das tool 

### 1) Binning with Maxbin2
Maxbin2 needs an abundance file for each sample. I first prepare for the abundance file.

```bash
for x in *.cov.txt;
do awk '{print $1"\t"$2}' $x | grep -v '^#' > $x.abundance.txt;
done
```
And then make a abundancelist.txt containing all sample names like followings:
Sample_201401.sam.cov.txt.abundance.txt
Sample_201402.sam.cov.txt.abundance.txt
Sample_201403.sam.cov.txt.abundance.txt
...
Sample_201670.sam.cov.txt.abundance.txt

Now we can run maxbin2
```bash
run_MaxBin.pl -thread 14 -contig final.contigs.fa -out maxbin1500 -abund_list abundancelist.txt -min_contig_length 1500
```

### 2) Binning with Concoct
run Concoct in Bioconda and prepare the assembly
```bash
conda activate concoct_env
conda install mkl #important
cut_up_fasta.py final.contigs.fa -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
```

Make a coverage_table.tsv
```bash
#this step  requires sorted and indexed bam 
conda activate concoct_env
concoct_coverage_table.py contigs_10K.bed *sort.bam > coverage_table.tsv
```

run Concoct binning
```bash
conda activate concoct_env
concoct --composition_file contigs_10K.fa --coverage_file coverage_table.all.tsv -b concoct_output/ -l 1500 -t 40
merge_cutup_clustering.py concoct_output/clustering_gt1500.csv > concoct_output/clustering_merged.csv
mkdir concoct_output/fasta_bins
extract_fasta_bins.py final.contigs.fa concoct_output/clustering_merged.csv --output_path concoct_output/fasta_bins
```

### 3) Binning with Metabat2
```bash
cd metabat2
./bin/jgi_summarize_bam_contig_depths --outputDepth depth.txt path/to/*.bam
./bin/metabat2 -i final.contigs.fa -a depth.txt -o bins_dir/bin -t 40 -m 1500
```

### 4) Binning with Das tool using the results from Maxbin2, Metabat2, Concoct
```bash
conda activate dastool
#Preparation of input files
#Converting MaxBin fasta output into tab separated scaffolds2bin file:
Fasta_to_Scaffolds2Bin.sh -i /concoct_output/fasta_bins -e fa > /DAStool/concoct.scaffolds2bin.tsv
Fasta_to_Scaffolds2Bin.sh -i /metabat2/bins_dir -e fa > /DAStool/metabat.scaffolds2bin.tsv
Fasta_to_Scaffolds2Bin.sh -i /maxbin_outfasta_1500/ -e fasta > /DAStool/maxbin.scaffolds2bin.tsv

cd DAStool/
DAS_Tool -i concoct.scaffolds2bin.tsv,metabat.scaffolds2bin.tsv,maxbin.scaffolds2bin.tsv \
	 -l concoct,metabat,maxbin \
         -c //megahit/final.contigs.fa \
         -o sample_output/ \
         --threads 24 \
	 --score_threshold 0.6 \
	 --write_bins 1 \
	 --search_engine diamond 
```

## E) Assess the bins (MAGs)
### 1) using CAT/BAT for taxonomy classification for the bins by DAStool
```bash
/CAT-master/CAT_pack/CAT bins -b /path/to/DASTool_bins/ -d /CAT_prepare_20210107/2021-01-07_CAT_database -t /CAT_prepare_20210107/2021-01-07_taxonomy -s .fa
/home/fr/fr_fr/fr_yw50/downloads/CAT-master/CAT_pack/CAT add_names -i out.BAT.bin2classification.txt -o out.BAT.bin2classification.taxname.out -t /CAT_prepare_20210107/2021-01-07_taxonomy
```
### 2) using busco to get a figure of the completeness of MAGs
```bash
conda activate busco
cd /DASTool_bins/
for x in *.fa; do busco -i $x -l bacteria_odb10 -o "${x//./}"out -m genome; done
cp  *out/short_summary.*.txt BUSCO_summaries/.
generate_plot.py -wd ./BUSCO_summaries
```
