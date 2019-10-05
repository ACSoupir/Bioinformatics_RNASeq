---
title: "RNA Seq for STAT736"
author: "Alex Soupir"
date: "9/10/2019"
output: 
  html_document:
    keep_md: true
---



# **Background**

For STAT736-Fall-2019, we are analyzing the RNA-Seq from the publication [Genome-wide analysis of p53 transcriptional programs in B cells upon exposure to genotoxic stress in vivo.](https://www.ncbi.nlm.nih.gov/pubmed/26372730?dopt=Abstract) We are only using the [sequences](https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP061386) *B cells from spleen* and not the *non-B cells from spleen* from the SRA Run Selector on NCBI.

The mice were exposed to whole-body ionizing radiation and sequences were extracted from both Bcells and non-B cells from the spleens of the mice. Two genotypes of mice were used: mice with p53 knocked out and the wild-type C57/Bl6. There were 4 different group combinations including the 2 different genotypes; each genotype was subjected to the ionizing radiation as well as control/mock.

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>Treatment groups of the mice that were either controls or treated with ionizing radiation to determine reaction of p53.</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> Genotype </th>
   <th style="text-align:left;"> Treatment </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Group 1 </td>
   <td style="text-align:left;"> p53 </td>
   <td style="text-align:left;"> Mock </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Group 2 </td>
   <td style="text-align:left;"> C57/Bl6 </td>
   <td style="text-align:left;"> Mock </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Group 3 </td>
   <td style="text-align:left;"> p53 </td>
   <td style="text-align:left;"> IR </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Group 5 </td>
   <td style="text-align:left;"> C57/Bl6 </td>
   <td style="text-align:left;"> IR </td>
  </tr>
</tbody>
</table>

The pipeline used in this analysis used **conda** on South Dakota State University's High Performance Computing cluster to run the programs **FastQC**, **Trimmomatic**, and **Tophat**. 

This is different than previous RNA-Seq analyses where I used my workstation pc with **Ubuntu 18.04** to run **FastQC**, **Trimmomatic**, **HiSat2**, **HTSeq**, and **DESeq2** locally. Also, the previous RNA-Seq alayses were of Soybean with treatment combinations of mycorrhizae and rhizobia inoculation.

# **Cleaning Data**

## Programs used?

+ FastQC
+ Trimmomatic 0.39
+ Bowtie 2.2.5.0
+ Tophat 2.1.1
+ STAR
+ Cufflinks
+ Kraken2
+ MultiQC
+ featureCount

### Picking the right node

To find a node that we can use on our own, we need to see which nodes are already allocated to jobs and which ones are idle. To do this, we can run **sinfo**. We want to pick one of the nodes that are marked 'idle' so we get the whole thing and we aren't interrupting someone elses job. For the sake of this exercise, lets work on **big-mem**.

Once a node that is idle has been found, you can *ssh* into it by typing **ssh -X big-mem00#** where # is the node number.

```bash
ssh -X big-mem005
```

Once on the node, the modules will have to be pulled from the shared folder again, otherwise we will be left with very basic ones. **NOTE:** if running programs that are in a personal folder such as miniconda (these examples), it is not necessary to add the other modules.


```bash
module use /cm/shared/modulefiles_local/
```

After loading the modules you can use it just as you would any other command line.

### Creating slurm scripts

When running on a cluster, it can sometime be difficult to find open nodes with the resources needed to run the jobs that we have. Making a slurm script is really easy. Fist we make a new file with the *touch* command.


```bash
touch commands.slurm
```

Now in our directory we have the file **commands.slurm** which we can edit to hold our code in. We can edit it with the *vi* command.


```bash
vi commands.slurm
```

We have a few things that we need to put in the file header so slurm knows what to do with our commands.


```bash
#!/bin/bash

#SBATCH --job-name=example
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --output=job-%j-%N.log
#SBATCH --partition=bigmem
#SBATCH --time=10:00:00
```

When we break this down, we see *--job-name* which is what we will see when we look at whats running later, *--nodes* is the number of nodes we have, *--ntasks-per-node* is the number of cores that we are requesting to have allocated, *--output* is the output log file of the job (here it names the output file with the job number and the node that we used), *--partition* here is requesting a big-mem node but **compute** can also be used, and finally *--time* is how long we are requesting the allocation for. 

If the time runs out before the job is done I believe that it just kills the job even if not finished so we need to think a little about how much time to set. If the time is set too low, the job is killed and if the time is set too long, we may face issues with getting the node allocated to us.

To submit a job we can use *sbatch commands.slurm* and then we have the job ID. To check the status of our submission we use *sbatch* and then it shows all of the submitted jobs and how long they have been running along with the name that we set in the script.

## Acquiring sequences

To download the sequences from the sequence read archive (SRA), the SRA Toolkit was used. The downloading of the files took a very long time, so this was left to run over night. The **--gzip** was used to keep the files a relatively small, although this can be left out to download uncompressed files, and **--split-files** was used to split the forward read from the reverse read for paired end read trimming through **Trimmomatic**.


```bash
~/tools/sratoolkit.2.9.6-1-centos_linux64/bin/fastq-dump --gzip --split-files SRR2121770
```

This is an example of the single file, but the above code needed to be ran for all of the following SRA numbers:

+ SRR2121770
+ SRR2121771
+ SRR2121774
+ SRR2121775
+ SRR2121778
+ SRR2121779
+ SRR2121780
+ SRR2121781
+ SRR2121786
+ SRR2121787
+ SRR2121788
+ SRR2121789

The results from downloading with **--split-files** gives 2 files per SRR, as mentioned before, one forward and one reverse. The suffix of the split files is one with **\_1.fastq.gz** and another with **\_2.fastq.gz**.

## FastQC

FastQC can be run on all of the read files by using the wild card (\*) as in **\*.fastq.gz**. This prevents the need to hard code each individual read file into a FastQC command, which saves a lot of time since there are 24 read files in total for these 12 samples.


```bash
~/miniconda2/bin/fastqc *.fastq.gz
```

The output from running FastQC is a zipped folder and an HTML file for each of the **\fastq.gz** files in the folder. The HTML document looks something like this:

![FastQC of Raw SRR2121770_1.fastq.gz Read](./TopFastQCRaw.PNG)

This is just the top of the file, and every category under the **Summary** heading has a graph that shows how the read quality looks for that particular metric. These reports can give insight into whether the reads are of decent quality or if the quality is poor. 

The raw reads we have here all passed for adapter content and sequence length distribution and everything failed per base sequence content. SRR2121770, SRR2121771, SRR2121774, SRR2121775, SRR2121788, SRR2121781-2, and SRR2121789-1 were fairly decent quality reads. SRR2121778, SRR2121779, SRR2121780, SRR2121786, SRR2121787, SRR2121781-1, and SRR2121789-2 were of fairly lower quality (failing 3 or more in both reads. All of them failed both per base sequence quality and per tile sequence quality. 

### MultiQC

First lets install **multiqc** with conda. The command for this is ***conda install \-c bioconda multiqc***.

When that is finished, we can run MultiQC in the folder with the QC files (they should be moved into a folder alone so things don't get cluttered later on in the analyses).


```bash
~/miniconda/bin/multiqc .
```

When MultiQC is finished running, there will be a new folder called **multiqc_data** where the summaries are stored. Now lets go back up a level where our raw data folder and fastqc folder is and make a new folder for all of our MultiQC data. We will copy the FastQC output from MultiQC to this new folder.


```bash
mkdir MultiQC_All

cp RawQC/multiqc_data/multiqc_fastqc.txt MultiQC_All/
```

## Trimming with Trimmomatic

Conda was used again to run Trimmomatic. This isn't as easy as using the wildcard like with FastQC because each output has to be personalized for the read files that are input into Trimmomatic. Also, we have to make sure that the adapter sequences are in the same folder that we are running so we can refer to them easily when calling the Trimmomatic program. In this case, we are using the TruSeq3-PE-2.fa adapter sequences For example:


```bash
~/miniconda2/bin/trimmomatic PE SRR2121770_1.fastq.gz SRR2121770_2.fastq.gz 770_fp.fq.gz 770_fu.fq.gz 770_rp.fq.gz 770_ru.fq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2keepBothReads LEADING:3 TRALING:3 MINLEN:36 &
```

This would be repeated for each of the pairs (12 in total). We are trimming paired-end reads with the TruSeq3-PE-2 adapters. We are chopping off the first and last 3 bases and if we end up with a sequence less than 36 bases, we get rid of it. We want to make sure that there are enough bases in a read to work with. These parameters can be tweaked for possibly better end results with less being discarded. When Trimmomatic is finished running, it will out put the total number of reads, the total number from both the forward and reverse reads that are kept, the number of only forward reads kept, the number of only reverse reads kept, and the number of discarded reads.

The highest number of reads dropped was from trimming SRR2121786, where 20.55% dropped. Most reads were between 5% and 10% dropped. SRR2121786, SRR2121787, and SRR2121779 had sequence drops greater than 15%.

The trimmed reads can be analyzed again with FastQC to see how well the trimming worked to make the file better quality. After running FastQC on the trimmed files we see that the quality of those that were really bad quality were improved. There were a few different metrics throughout all of the files that bounced from a warning before the failing, or from passing before to a warning, and so forth, overall creating better quality read files.

## Alignment

### Using Tophat

Tophat can be installed using the same **conda install** (***conda install \-c bioconda tophat***). When this is finished installing, then we will need to get the mouse genome from the [Johns Hopkins Univeristy Center for computational BIology](http://ccb.jhu.edu/software/tophat/igenomes.shtml). The version of the mouse genome that I am using here is the [NCBI build37.2](ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/NCBI/build37.2/Mus_musculus_NCBI_build37.2.tar.gz). Instead of downloading this from the website and having to move it to the cluser, I will just download it using wget into the folder that has the raw reads, trimmed reads, and the FastQC files.


```bash
wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/NCBI/build37.2/Mus_musculus_NCBI_build37.2.tar.gz
```

This will take a long time to download because the file is a little less than 16GB zipped.

We notice here that we have a zipped tar file. To make this file easier to use, lets unzip it.


```bash
tar zxvf Mus_muculus_NCBI_build37.2.tar.gz
```

Since Tophat is requiring **\*.bt21** files (large index) and the files downloaded for the genome above are only small index files, we have to create a large index using **bowtie2-build**. For this, lets navigate to the WholeGenomeFasta folder within the extracted folder and then run **bowtie2-build**.


```bash
~/miniconda2/bin/bowtie2-build --large-index genome.fa genome
```

This process took about 26 minutes to run. Now lets copy the index files to a folder close to our reads so we can access them easier, rather than having to refer to the longer path where we build them. After they are copied to a new folder closer to our working directory, I went ahead and unzipped the trimmed read files to try and make the Tophat faster but it turned out not to work. The multicore call with **\-p** didn't use more cores than 1 until **bowtie2-align-s**, then 20 cores were used.


```bash
~/miniconda2/bin/tophat --no-converage-search -p 20 -G Mus_musculus/NCBI/build37.2/Annotation/Archives/archive-2015-07-17-14-32-40/Genes/genes.gtf -0 770_thout ./Index/genome 770_fp.fq.gz 770_rp.fq.gz 770_fu.fq 770_ru.fq
```

This run took almost 3 hours to complete.. Running with 80 cores rather than 20 cores took just 4 minutes less, so the whole process must be limited by a single core and the core clock speed. The process does use close to 8,000% at its peak so there is a benefit to multicore, just isn't very scalable.

### Using STAR

STAR can be installed the same way as the previous programs with **conda install** (***conda install \-c bioconda star***). In order to run STAR, we need to creaate indices just like with tophat, but STAR has this built in. I'm going to be using the same genome and GTF file as previously downloaded, but Dr. Ge uses a different zipped genome from the *gencode* database.


```bash
~/miniconda2/bin/STAR \
--runThreadN 80 \
--runMode genomeGenerate \
--genomeDir starIndex \
--genomeFastaFiles Index/genome.fa \ #same when we made the bowtie indices 
--sjdbGTFfile Mus_musculus/NCBI/build37.2/Annotation/Archives/archive-2015-07-17-14-32-40/Genes/genes.gtf
```

With the index files made, we can start aligning with STAR. It's important here than we only pick the paired end reads and not use all of the reads. Tophat is able to use all 4 reads but STAR doesn't allow that, so we need to make sure that we feed in the large files from trimming.


```bash
~/miniconda2/bin/STAR --runThreadN 80 --genomeDir starIndex --readFilesIn 770_fp.fq 770_rp.fq --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix 2121770 --outSAMtype BAM SortedByCoordinate
```

## Assembling transcripts with Cufflinks

Once STAR is done running, we can assemble the transcripts with Cufflinks. This can also be installed with **conda install** (***conda install \-c bioconda cufflinks***).


```bash
~/miniconda2/bin/cufflinks -p 20 -o SRR2121771_clout --library-type fr-firststrand 2121770Aligned.sortedByCoord.out.bam
```

## Checking for Contamination

### PhiX contamination

Now we will look at what kind of contamination we are looking at. When samples are sequenced with Illumina, a PhiX control is run along side them. This control is for cluster generation, sequencing, alignment, and calibration for cross-talk matrix generation. We will use Bowtie to create a file to determine the PhiX contamination level.


```bash
~/miniconda2/bin/bowtie2 -p 20 -x PhiX/Illumina/RTA/Sequence/Bowtie2Index/genome \
-1 TrimmedReads/770_fp.fq  -2 TrimmedReads/770_rp.fq -S phix.sam &> PhiXout/SRR2121770_phix.out
```

When the job is done running, the output file will show how much PhiX contamination we have. For example, lookin at the **SRR2121770_phix.out** created above, we see that 0.11% of the reads aligned with PhiX. The lower this value the better.

### rRNA Sequences

To retreive the rRNA sequences for mouse, we need to search the taxonomy database on NCBI for *Mus musculus*. Click on *Mus musculus* on the next page, and then the top *Mus musculus* at the head of the list. Now, select the top subtree link in the **Nucleotide** database. Select rRNA sequences on the left side of the page and download full list just downloading with Send > Complete Record > File > FASTA > Create File. Drag the file using WinSCP to the raw folder on the cluster and rename it to rRNA.fa.

We are going to need to install **bwa** with conda in order to get the alignments to work. This can be done with ***conda install -c bioconda bwa***. Following this, we will need to make indixes for the rRNA that we downloaded. To make this more clean, lets make a directory for the rRNA sequences that we downloaded and the indices that we make.


```bash
mkdir rRNA
```

Then we move the **rRNA.fa** to the new **rRNA** folder with WinSCP and then we can run the bwa.


```bash
time ~/miniconda2/bin/bwa mem -t 20 rRNA/rRNA.fa TrimmedReads/770_fp.fq TrimmedReads/770_rp.fq > rnaAlign/770_rna.sam
```

When we are done creating the new *\*.sam* files for all of the forward/reverse read combinations, we can use samtools to convert the *\*.sam* file to *\*.bam* files which are essentially the same file just that *sam* is easier for us to look at while *bam* is binary. Samtools can be installed with ***conda install \-c bioconda samtools***.


```bash
~/miniconda2/bin/samtools view -@ 10 -bS -o rnaAlign/770_rna.bam rnaAlign/770.sam
```

Now in the rnaAlign folder we have our sam and bam file for each of the libraries. Lets create an output file with *flagstat*.


```bash
~/miniconda2/bin/samtools flagstat -@ 10 rnaAlign/770_rna.out
```

Wihtin this file we will be able to see the summary of our alignments to the rRNA file that we downloaded from NCBI.


```bash
#From 770_rna.out
205559289 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
85 + 0 supplementary
0 + 0 duplicates
4265179 + 0 mapped (2.07% : N/A)
205559204 + 0 paired in sequencing
102779602 + 0 read1
102779602 + 0 read2
4151684 + 0 properly paired (2.02% : N/A)
4187608 + 0 with itself and mate mapped
77486 + 0 singletons (0.04% : N/A)
4222 + 0 with mate mapped to a different chr
1026 + 0 with mate mapped to a different chr (mapQ>=5)
```

### Bacterial contamination

In order to find out the contamination, we need to install Kraken2 with ***conda install \-c bioconda kraken2*** and download a pre-built database containing bacteria, archaea, and viral sequences. The database we are going to download only contains about 5% of k-mers from the original database (but directions are sort of lacking to build an entirely new database). More information can be found at https://ccb.jhu.edu/software/kraken/ for the pre-built databases.

Using the code in the next chunk will download the 8GB database and then extract the files so we can use them with the **Kraken2** program. Lets do this in the main project folder.


```bash
wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken2_v2_8GB_201904_UPDATE.tgz

tar xzf minikraken2_v2_8GB_201904_UPDATE.tgz
```

Now lets make a directory for the output.


```bash
mkdir krakenOut
```

We can call **Kraken** with the extracted database folder and point it to the location of out paired end reads from trimming and to the output folder that we just created for the outputs.


```bash
~/miniconda2/bin/kraken2 --db minikraken2_v2_8GB_201904_UPDATE/ --output krakenOut/770.out --threads 10 --paired TrimmedReads/770_fp.fq TrimmedReads/770_rp.fq
```

When Kraken is done running, it will print out the number (and percentage) of reads that were classified. In this case, we have used 102779602 sequences, of which 19142843 sequences were classified (18.63%) and 83636759 sequences were unclassified (81.37%). My interpretation of this is that 18.63% of the reads are possibly from microbial cell contamination.

## Counting Transcripts

Since we have the **bam** files from the alignments of the different samples, we can count the features for each and get the transcipt counts using featureCounts form ***conda install \-c bioconda/label/cf201901 subread***. The genome and annotations that we previously downloaded were from genome **mm9** so we have to specify to *featureCounts* what we want to actually count. FeatureCounts defaults to using **gene_id** which our output bam files don't have described correctly for *featureCounts* to read them. This is a single line of code because we can use a wildcard to run through all of the **bam** files. 


```bash
#Move to the Star Alignment output folder for a working directory
cd StarOut

~/miniconda2/bin/featureCounts -a /gpfs/scratch/alex.soupir/Mus/raw/Mus_musculus/NCBI/build37.2/Annotation/Archives/archive-2015-07-17-14-32-40/Genes/genes.gtf -g 'transcript_id' -o readCounds.txt *bam
```

With the files that we are working with, this will take between 3.5 minutes to 5 minutes per **bam** file. The output will be a file that can be imported into excel and saved as csv which we then can work with in R.

### Final QC of cleaning the data

Lets look at the data that we have collected from all of the MultiQC runs that we had with initial FastQC, Trimmomatic, STAR alignment, PhiX contamination, rRNA contamination, and the final feature counts.


```r
qc = read.csv('Whole Data QC.csv', header=FALSE, na.strings="")
kable(qc) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"))
```

<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> V1 </th>
   <th style="text-align:left;"> V2 </th>
   <th style="text-align:left;"> V3 </th>
   <th style="text-align:left;"> V4 </th>
   <th style="text-align:left;"> V5 </th>
   <th style="text-align:left;"> V6 </th>
   <th style="text-align:left;"> V7 </th>
   <th style="text-align:left;"> V8 </th>
   <th style="text-align:left;"> V9 </th>
   <th style="text-align:left;"> V10 </th>
   <th style="text-align:left;"> V11 </th>
   <th style="text-align:left;"> V12 </th>
   <th style="text-align:left;"> V13 </th>
   <th style="text-align:left;"> V14 </th>
   <th style="text-align:left;"> V15 </th>
   <th style="text-align:left;"> V16 </th>
   <th style="text-align:left;"> V17 </th>
   <th style="text-align:left;"> V18 </th>
   <th style="text-align:left;"> V19 </th>
   <th style="text-align:left;"> V20 </th>
   <th style="text-align:left;"> V21 </th>
   <th style="text-align:left;"> V22 </th>
   <th style="text-align:left;"> V23 </th>
   <th style="text-align:left;"> V24 </th>
   <th style="text-align:left;"> V25 </th>
   <th style="text-align:left;"> V26 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Sample </td>
   <td style="text-align:left;"> SRR2121770_1 </td>
   <td style="text-align:left;"> SRR2121770_2 </td>
   <td style="text-align:left;"> SRR2121771_1 </td>
   <td style="text-align:left;"> SRR2121771_2 </td>
   <td style="text-align:left;"> SRR2121774_1 </td>
   <td style="text-align:left;"> SRR2121774_2 </td>
   <td style="text-align:left;"> SRR2121775_1 </td>
   <td style="text-align:left;"> SRR2121775_2 </td>
   <td style="text-align:left;"> SRR2121778_1 </td>
   <td style="text-align:left;"> SRR2121778_2 </td>
   <td style="text-align:left;"> SRR2121779_1 </td>
   <td style="text-align:left;"> SRR2121779_2 </td>
   <td style="text-align:left;"> SRR2121780_1 </td>
   <td style="text-align:left;"> SRR2121780_2 </td>
   <td style="text-align:left;"> SRR2121781_1 </td>
   <td style="text-align:left;"> SRR2121781_2 </td>
   <td style="text-align:left;"> SRR2121786_1 </td>
   <td style="text-align:left;"> SRR2121786_2 </td>
   <td style="text-align:left;"> SRR2121787_1 </td>
   <td style="text-align:left;"> SRR2121787_2 </td>
   <td style="text-align:left;"> SRR2121788_1 </td>
   <td style="text-align:left;"> SRR2121788_2 </td>
   <td style="text-align:left;"> SRR2121789_1 </td>
   <td style="text-align:left;"> SRR2121789_2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FastQC </td>
   <td style="text-align:left;"> adapter_content </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Sequences flagged as poor quality </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> sequence_duplication_levels </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> warn </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> avg_sequence_length </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Encoding </td>
   <td style="text-align:left;"> Sanger / Illumina 1.9 </td>
   <td style="text-align:left;"> Sanger / Illumina 1.9 </td>
   <td style="text-align:left;"> Sanger / Illumina 1.9 </td>
   <td style="text-align:left;"> Sanger / Illumina 1.9 </td>
   <td style="text-align:left;"> Sanger / Illumina 1.9 </td>
   <td style="text-align:left;"> Sanger / Illumina 1.9 </td>
   <td style="text-align:left;"> Sanger / Illumina 1.9 </td>
   <td style="text-align:left;"> Sanger / Illumina 1.9 </td>
   <td style="text-align:left;"> Sanger / Illumina 1.9 </td>
   <td style="text-align:left;"> Sanger / Illumina 1.9 </td>
   <td style="text-align:left;"> Sanger / Illumina 1.9 </td>
   <td style="text-align:left;"> Sanger / Illumina 1.9 </td>
   <td style="text-align:left;"> Sanger / Illumina 1.9 </td>
   <td style="text-align:left;"> Sanger / Illumina 1.9 </td>
   <td style="text-align:left;"> Sanger / Illumina 1.9 </td>
   <td style="text-align:left;"> Sanger / Illumina 1.9 </td>
   <td style="text-align:left;"> Sanger / Illumina 1.9 </td>
   <td style="text-align:left;"> Sanger / Illumina 1.9 </td>
   <td style="text-align:left;"> Sanger / Illumina 1.9 </td>
   <td style="text-align:left;"> Sanger / Illumina 1.9 </td>
   <td style="text-align:left;"> Sanger / Illumina 1.9 </td>
   <td style="text-align:left;"> Sanger / Illumina 1.9 </td>
   <td style="text-align:left;"> Sanger / Illumina 1.9 </td>
   <td style="text-align:left;"> Sanger / Illumina 1.9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> per_base_sequence_quality </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> fail </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> sequence_length_distribution </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Sequence length </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:left;"> 51 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> File type </td>
   <td style="text-align:left;"> Conventional base calls </td>
   <td style="text-align:left;"> Conventional base calls </td>
   <td style="text-align:left;"> Conventional base calls </td>
   <td style="text-align:left;"> Conventional base calls </td>
   <td style="text-align:left;"> Conventional base calls </td>
   <td style="text-align:left;"> Conventional base calls </td>
   <td style="text-align:left;"> Conventional base calls </td>
   <td style="text-align:left;"> Conventional base calls </td>
   <td style="text-align:left;"> Conventional base calls </td>
   <td style="text-align:left;"> Conventional base calls </td>
   <td style="text-align:left;"> Conventional base calls </td>
   <td style="text-align:left;"> Conventional base calls </td>
   <td style="text-align:left;"> Conventional base calls </td>
   <td style="text-align:left;"> Conventional base calls </td>
   <td style="text-align:left;"> Conventional base calls </td>
   <td style="text-align:left;"> Conventional base calls </td>
   <td style="text-align:left;"> Conventional base calls </td>
   <td style="text-align:left;"> Conventional base calls </td>
   <td style="text-align:left;"> Conventional base calls </td>
   <td style="text-align:left;"> Conventional base calls </td>
   <td style="text-align:left;"> Conventional base calls </td>
   <td style="text-align:left;"> Conventional base calls </td>
   <td style="text-align:left;"> Conventional base calls </td>
   <td style="text-align:left;"> Conventional base calls </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> basic_statistics </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> per_sequence_gc_content </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Total Sequences </td>
   <td style="text-align:left;"> 118323219 </td>
   <td style="text-align:left;"> 118323219 </td>
   <td style="text-align:left;"> 103127231 </td>
   <td style="text-align:left;"> 103127231 </td>
   <td style="text-align:left;"> 91225885 </td>
   <td style="text-align:left;"> 91225885 </td>
   <td style="text-align:left;"> 108661623 </td>
   <td style="text-align:left;"> 108661623 </td>
   <td style="text-align:left;"> 136766553 </td>
   <td style="text-align:left;"> 136766553 </td>
   <td style="text-align:left;"> 124265478 </td>
   <td style="text-align:left;"> 124265478 </td>
   <td style="text-align:left;"> 99685598 </td>
   <td style="text-align:left;"> 99685598 </td>
   <td style="text-align:left;"> 106360532 </td>
   <td style="text-align:left;"> 106360532 </td>
   <td style="text-align:left;"> 99988292 </td>
   <td style="text-align:left;"> 99988292 </td>
   <td style="text-align:left;"> 75410310 </td>
   <td style="text-align:left;"> 75410310 </td>
   <td style="text-align:left;"> 82750449 </td>
   <td style="text-align:left;"> 82750449 </td>
   <td style="text-align:left;"> 92454136 </td>
   <td style="text-align:left;"> 92454136 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> per_base_n_content </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> fail </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> per_base_sequence_content </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> overrepresented_sequences </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> warn </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> %GC </td>
   <td style="text-align:left;"> 47 </td>
   <td style="text-align:left;"> 47 </td>
   <td style="text-align:left;"> 46 </td>
   <td style="text-align:left;"> 46 </td>
   <td style="text-align:left;"> 47 </td>
   <td style="text-align:left;"> 47 </td>
   <td style="text-align:left;"> 47 </td>
   <td style="text-align:left;"> 47 </td>
   <td style="text-align:left;"> 46 </td>
   <td style="text-align:left;"> 47 </td>
   <td style="text-align:left;"> 46 </td>
   <td style="text-align:left;"> 46 </td>
   <td style="text-align:left;"> 46 </td>
   <td style="text-align:left;"> 46 </td>
   <td style="text-align:left;"> 46 </td>
   <td style="text-align:left;"> 46 </td>
   <td style="text-align:left;"> 47 </td>
   <td style="text-align:left;"> 47 </td>
   <td style="text-align:left;"> 46 </td>
   <td style="text-align:left;"> 46 </td>
   <td style="text-align:left;"> 47 </td>
   <td style="text-align:left;"> 47 </td>
   <td style="text-align:left;"> 47 </td>
   <td style="text-align:left;"> 48 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> total_deduplicated_percentage </td>
   <td style="text-align:left;"> 52.20958673 </td>
   <td style="text-align:left;"> 55.28743506 </td>
   <td style="text-align:left;"> 90.47297502 </td>
   <td style="text-align:left;"> 97.22714143 </td>
   <td style="text-align:left;"> 91.11411337 </td>
   <td style="text-align:left;"> 97.47431986 </td>
   <td style="text-align:left;"> 72.31423101 </td>
   <td style="text-align:left;"> 95.57923295 </td>
   <td style="text-align:left;"> 65.27895305 </td>
   <td style="text-align:left;"> 83.85617925 </td>
   <td style="text-align:left;"> 73.08810256 </td>
   <td style="text-align:left;"> 79.50169478 </td>
   <td style="text-align:left;"> 92.46127029 </td>
   <td style="text-align:left;"> 94.09779679 </td>
   <td style="text-align:left;"> 54.09833163 </td>
   <td style="text-align:left;"> 52.4832411 </td>
   <td style="text-align:left;"> 68.88309902 </td>
   <td style="text-align:left;"> 74.63550918 </td>
   <td style="text-align:left;"> 93.00930552 </td>
   <td style="text-align:left;"> 94.9666922 </td>
   <td style="text-align:left;"> 43.56442006 </td>
   <td style="text-align:left;"> 44.95245995 </td>
   <td style="text-align:left;"> 42.6643141 </td>
   <td style="text-align:left;"> 57.36052687 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Filename </td>
   <td style="text-align:left;"> SRR2121770_1.fastq.gz </td>
   <td style="text-align:left;"> SRR2121770_2.fastq.gz </td>
   <td style="text-align:left;"> SRR2121771_1.fastq.gz </td>
   <td style="text-align:left;"> SRR2121771_2.fastq.gz </td>
   <td style="text-align:left;"> SRR2121774_1.fastq.gz </td>
   <td style="text-align:left;"> SRR2121774_2.fastq.gz </td>
   <td style="text-align:left;"> SRR2121775_1.fastq.gz </td>
   <td style="text-align:left;"> SRR2121775_2.fastq.gz </td>
   <td style="text-align:left;"> SRR2121778_1.fastq.gz </td>
   <td style="text-align:left;"> SRR2121778_2.fastq.gz </td>
   <td style="text-align:left;"> SRR2121779_1.fastq.gz </td>
   <td style="text-align:left;"> SRR2121779_2.fastq.gz </td>
   <td style="text-align:left;"> SRR2121780_1.fastq.gz </td>
   <td style="text-align:left;"> SRR2121780_2.fastq.gz </td>
   <td style="text-align:left;"> SRR2121781_1.fastq.gz </td>
   <td style="text-align:left;"> SRR2121781_2.fastq.gz </td>
   <td style="text-align:left;"> SRR2121786_1.fastq.gz </td>
   <td style="text-align:left;"> SRR2121786_2.fastq.gz </td>
   <td style="text-align:left;"> SRR2121787_1.fastq.gz </td>
   <td style="text-align:left;"> SRR2121787_2.fastq.gz </td>
   <td style="text-align:left;"> SRR2121788_1.fastq.gz </td>
   <td style="text-align:left;"> SRR2121788_2.fastq.gz </td>
   <td style="text-align:left;"> SRR2121789_1.fastq.gz </td>
   <td style="text-align:left;"> SRR2121789_2.fastq.gz </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> per_tile_sequence_quality </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> warn </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> fail </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> per_sequence_quality_scores </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> fail </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
   <td style="text-align:left;"> pass </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Trimming </td>
   <td style="text-align:left;"> surviving </td>
   <td style="text-align:left;"> 102779602 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 88509413 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 77240029 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 94179097 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 108337271 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 92924182 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 76172416 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 95488961 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 69906417 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 54427039 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 73771492 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 82764590 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> surviving_pct </td>
   <td style="text-align:left;"> 86.86 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 85.83 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 84.67 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 86.67 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 79.21 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 74.78 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 76.41 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 89.78 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 69.91 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 72.17 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 89.15 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 89.52 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> reverse_only_surviving_pct </td>
   <td style="text-align:left;"> 4.48 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2.52 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2.53 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2.29 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2.11 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2.44 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2.52 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.8 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2.59 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2.65 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.85 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.77 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> forward_only_surviving </td>
   <td style="text-align:left;"> 3202267 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 5041955 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4667583 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4832339 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 9375417 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 7916665 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 6350522 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3894412 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 6944977 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 5349296 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3489058 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3797776 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> dropped </td>
   <td style="text-align:left;"> 7036792 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 6973727 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 7013028 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 7160952 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 16166923 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 20393774 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 14647690 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 5067412 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 20551506 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 13638339 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3962243 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4252397 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> dropped_pct </td>
   <td style="text-align:left;"> 5.95 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 6.76 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 7.69 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 6.59 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 11.82 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 16.41 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 14.69 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4.76 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 20.55 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 18.09 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4.79 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4.6 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> forward_only_surviving_pct </td>
   <td style="text-align:left;"> 2.71 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4.89 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 5.12 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4.45 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 6.86 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 6.37 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 6.37 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3.66 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 6.95 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 7.09 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4.22 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4.11 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> input_read_pairs </td>
   <td style="text-align:left;"> 118323219 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 103127231 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 91225885 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 108661623 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 136766553 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 124265478 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 99685598 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 106360532 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 99988292 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 75410310 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 82750449 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 92454136 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> reverse_only_surviving </td>
   <td style="text-align:left;"> 5304558 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2602136 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2305245 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2489235 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2886942 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3030857 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2514970 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1909747 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2585392 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1995636 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1527656 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1639373 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> STAR Alignment </td>
   <td style="text-align:left;"> uniquely_mapped_percent </td>
   <td style="text-align:left;"> 81.73 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 82.35 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 81.67 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 83.1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 83.43 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 84.55 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 82.48 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 86.82 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 80.75 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 80.27 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 80.81 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 81.62 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> num_splices </td>
   <td style="text-align:left;"> 8669080 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 8000968 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 10613774 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 11755566 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 10864591 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 8361644 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 6866296 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 7991948 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 6183074 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4335214 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 7285454 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 8581902 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> num_GCAG_splices </td>
   <td style="text-align:left;"> 90012 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 83679 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 116149 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 124779 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 116447 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 86763 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 74268 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 80836 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 65893 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 44765 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 80392 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 97422 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> insertion_length </td>
   <td style="text-align:left;"> 1.11 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.12 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.17 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.13 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.11 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.11 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.13 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.16 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.15 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.13 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.15 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> deletion_length </td>
   <td style="text-align:left;"> 1.35 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.33 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.37 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.26 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.36 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.35 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.35 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.22 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.25 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> unmapped_tooshort_percent </td>
   <td style="text-align:left;"> 4.08 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4.21 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3.42 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4.09 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3.51 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3.87 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3.66 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3.45 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 5.83 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 6.42 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 5.61 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4.53 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> avg_mapped_read_length </td>
   <td style="text-align:left;"> 100.73 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 101.11 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 101.1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 101.14 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 100.6 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 100.51 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 100.51 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 100.81 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 100.41 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 100.33 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 101.18 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 101.22 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> deletion_rate </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.01 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.01 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.01 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.01 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.01 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.01 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.01 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.01 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> mismatch_rate </td>
   <td style="text-align:left;"> 0.23 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.29 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.29 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.17 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.17 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.22 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.22 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.17 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.26 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.28 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.16 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.15 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> avg_input_read_length </td>
   <td style="text-align:left;"> 101 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 101 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 101 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 101 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 100 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 100 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 100 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 101 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 100 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 100 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 101 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 101 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> num_ATAC_splices </td>
   <td style="text-align:left;"> 8039 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 7497 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 10782 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 11677 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 9722 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 7455 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 6186 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 7544 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 5492 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3832 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 6602 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 7593 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> num_annotated_splices </td>
   <td style="text-align:left;"> 8581504 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 7919864 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 10513157 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 11643271 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 10755797 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 8273449 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 6799294 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 7913548 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 6103902 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4260450 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 7202217 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 8487777 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> num_GTAG_splices </td>
   <td style="text-align:left;"> 8571029 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 7909792 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 10486843 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 11619110 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 10738422 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 8267426 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 6785842 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 7903568 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 6111689 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4286617 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 7198460 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 8476887 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> uniquely_mapped </td>
   <td style="text-align:left;"> 83997598 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 72886249 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 63080382 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 78265974 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 90382084 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 78565331 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 62830380 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 82902461 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 56450144 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 43685904 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 59617353 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 67556146 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> multimapped_toomany </td>
   <td style="text-align:left;"> 1852495 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1405409 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1242903 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1061446 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1736807 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1362362 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1374461 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 846193 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1009216 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 955877 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 703236 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 890704 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> unmapped_mismatches </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> unmapped_mismatches_percent </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> total_reads </td>
   <td style="text-align:left;"> 102779602 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 88509413 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 77240029 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 94179097 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 108337271 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 92924182 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 76172416 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 95488961 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 69906417 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 54427039 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 73771492 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 82764590 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> unmapped_other </td>
   <td style="text-align:left;"> 432340 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 370957 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 254880 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 377196 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 476580 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 454868 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 335571 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 468726 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 370629 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 359051 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 317451 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 405340 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> insertion_rate </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> unmapped_other_percent </td>
   <td style="text-align:left;"> 0.42 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.42 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.33 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.4 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.44 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.49 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.44 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.49 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.53 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.66 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.43 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.49 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> multimapped_percent </td>
   <td style="text-align:left;"> 11.96 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 11.44 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 12.97 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 11.27 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 11.02 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 9.63 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 11.61 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 8.35 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 11.44 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 10.9 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 12.19 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 12.28 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> multimapped </td>
   <td style="text-align:left;"> 12297290 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 10128397 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 10020381 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 10617654 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 11939994 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 8949092 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 8840667 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 7971370 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 7999509 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 5933625 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 8991826 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 10165077 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> num_noncanonical_splices </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> unmapped_tooshort </td>
   <td style="text-align:left;"> 4199879 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3718401 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2641483 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3856827 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3801806 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3592529 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2791337 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3300211 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4076919 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3492582 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4141626 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3747323 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> multimapped_toomany_percent </td>
   <td style="text-align:left;"> 1.8 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.59 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.61 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.13 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.6 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.47 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.8 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.89 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.44 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.76 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.95 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.08 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PhiX Contamination </td>
   <td style="text-align:left;"> overall_alignment_rate </td>
   <td style="text-align:left;"> 0.11 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.13 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.18 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.12 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.12 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.15 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.18 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.16 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.2 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> paired_aligned_mate_multi_halved </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> paired_aligned_none </td>
   <td style="text-align:left;"> 102674701 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 88400163 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 77108986 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 93998943 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 108212089 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 92820803 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 76103405 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 95350795 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 69784944 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 54344793 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 73629819 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 82605166 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> paired_aligned_mate_multi </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> paired_aligned_discord_one </td>
   <td style="text-align:left;"> 5081 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 5514 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 6603 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 9196 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 6325 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 5222 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3475 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 7190 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 6424 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4218 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 7405 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 8139 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> paired_aligned_mate_one </td>
   <td style="text-align:left;"> 2761 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2466 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2825 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3974 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2693 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2277 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1515 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2999 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2696 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1961 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2913 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3585 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> paired_total </td>
   <td style="text-align:left;"> 102779602 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 88509413 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 77240029 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 94179097 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 108337271 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 92924182 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 76172416 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 95488961 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 69906417 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 54427039 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 73771492 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 82764590 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> paired_aligned_mate_none_halved </td>
   <td style="text-align:left;"> 102668239.5 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 88393416 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 77100970.5 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 93987760 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 108204417.5 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 92814442.5 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 76099172.5 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 95342105.5 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 69777172 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 54339594.5 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 73620957.5 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 82595234.5 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> paired_aligned_multi </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> paired_aligned_one </td>
   <td style="text-align:left;"> 104901 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 109250 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 131043 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 180154 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 125182 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 103379 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 69011 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 138166 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 121473 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 82246 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 141673 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 159424 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> total_reads </td>
   <td style="text-align:left;"> 102779602 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 88509413 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 77240029 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 94179097 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 108337271 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 92924182 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 76172416 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 95488961 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 69906417 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 54427039 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 73771492 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 82764590 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> paired_aligned_mate_one_halved </td>
   <td style="text-align:left;"> 1380.5 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1233 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1412.5 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1987 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1346.5 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1138.5 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 757.5 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1499.5 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1348 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 980.5 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1456.5 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1792.5 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> paired_aligned_mate_none </td>
   <td style="text-align:left;"> 205336479 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 176786832 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 154201941 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 187975520 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 216408835 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 185628885 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 152198345 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 190684211 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 139554344 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 108679189 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 147241915 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 165190469 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rRNA Contamination </td>
   <td style="text-align:left;"> mapped_passed </td>
   <td style="text-align:left;"> 4265172 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3371478 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1366577 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4031209 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 5947728 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 5093417 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4109793 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 5042021 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 7103552 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4879843 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 5080844 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 7013490 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> duplicates_failed </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> secondary_passed </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> paired in sequencing_failed </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> duplicates_passed </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> read2_passed </td>
   <td style="text-align:left;"> 102779602 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 88509413 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 77240029 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 94179097 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 108337271 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 92924182 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 76172416 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 95488961 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 69906417 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 54427039 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 73771492 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 82764590 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> read1_failed </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> read1_passed </td>
   <td style="text-align:left;"> 102779602 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 88509413 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 77240029 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 94179097 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 108337271 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 92924182 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 76172416 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 95488961 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 69906417 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 54427039 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 73771492 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 82764590 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> with mate mapped to a different chr_passed </td>
   <td style="text-align:left;"> 4240 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3672 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 692 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2742 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 7850 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 5990 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4782 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 5840 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 5348 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2434 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2250 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3646 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> total_failed </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> properly paired_passed_pct </td>
   <td style="text-align:left;"> 2.02 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.84 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.86 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2.07 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2.7 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2.69 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2.65 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2.59 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4.98 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4.39 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3.36 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4.15 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> singletons_passed </td>
   <td style="text-align:left;"> 77493 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 85229 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 21005 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 88202 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 48177 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 47851 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 41299 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 53429 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 60929 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 55645 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 58038 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 57600 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> supplementary_passed </td>
   <td style="text-align:left;"> 85 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 65 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 16 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 23 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 43 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 44 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 44 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 72 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 19 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 18 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 58 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 50 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> singletons_passed_pct </td>
   <td style="text-align:left;"> 0.04 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.05 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.01 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.05 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.02 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.03 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.03 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.03 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.04 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.05 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.04 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.03 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> mapped_failed_pct </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> mapped_passed_pct </td>
   <td style="text-align:left;"> 2.07 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1.9 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.88 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2.14 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2.75 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2.74 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2.7 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2.64 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 5.08 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4.48 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3.44 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4.24 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> supplementary_failed </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> with itself and mate mapped_failed </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> mapped_failed </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> total_passed </td>
   <td style="text-align:left;"> 205559289 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 177018891 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 154480074 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 188358217 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 216674585 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 185848408 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 152344876 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 190977994 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 139812853 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 108854096 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 147543042 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 165529230 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> properly paired_failed </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> flagstat_total </td>
   <td style="text-align:left;"> 205559289 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 177018891 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 154480074 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 188358217 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 216674585 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 185848408 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 152344876 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 190977994 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 139812853 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 108854096 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 147543042 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 165529230 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> with mate mapped to a different chr (mapQ &gt;= 5)_passed </td>
   <td style="text-align:left;"> 1026 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 854 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 118 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 442 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1093 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 940 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 801 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1206 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 964 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 487 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 823 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 951 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> properly paired_failed_pct </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> with mate mapped to a different chr (mapQ &gt;= 5)_failed </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> with itself and mate mapped_passed </td>
   <td style="text-align:left;"> 4187594 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3286184 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1345556 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3942984 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 5899508 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 5045522 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4068450 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4988520 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 7042604 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4824180 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 5022748 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 6955840 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> read2_failed </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> with mate mapped to a different chr_failed </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> properly paired_passed </td>
   <td style="text-align:left;"> 4151444 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3256926 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1334312 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3903816 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 5846848 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 5001956 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4035832 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4945218 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 6955726 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4775228 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4958480 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 6869584 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> paired in sequencing_passed </td>
   <td style="text-align:left;"> 205559204 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 177018826 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 154480058 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 188358194 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 216674542 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 185848364 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 152344832 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 190977922 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 139812834 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 108854078 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 147542984 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 165529180 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> singletons_failed </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> secondary_failed </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> singletons_failed_pct </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> nan </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Feature Counts </td>
   <td style="text-align:left;"> Unassigned_Ambiguity </td>
   <td style="text-align:left;"> 5553184 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4245795 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4120211 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4873854 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 5225069 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4143370 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3730746 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3744259 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3208052 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 2154031 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3756137 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 4158307 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Unassigned_MappingQuality </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> percent_assigned </td>
   <td style="text-align:left;"> 38.92490502 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 41.43212911 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 45.0933741 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 46.82554793 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 45.79057473 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 42.0896936 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 43.45341423 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 41.01237187 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 46.15918368 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 44.71355742 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 45.65729376 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 46.21716065 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Unassigned_Nonjunction </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Unassigned_Duplicate </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Unassigned_Chimera </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Unassigned_Unmapped </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Assigned </td>
   <td style="text-align:left;"> 92630635 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 83632332 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 81770639 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 99643999 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 112432428 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 87231124 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 76094464 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 86121025 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 71761506 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 53877618 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 75554738 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 86661345 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Unassigned_MultiMapping </td>
   <td style="text-align:left;"> 69977472 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 56081314 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 55175498 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 56266414 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 64772012 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 50119894 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 49456582 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 44182994 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 42564996 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 33123230 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 46247590 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 52396728 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Unassigned_Overlapping_Length </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Unassigned_Secondary </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Unassigned_FragmentLength </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Unassigned_NoFeatures </td>
   <td style="text-align:left;"> 69811377 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 57894371 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 40269914 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 52014095 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 63106671 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 65756168 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 45835550 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 75939638 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 37930730 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 31340159 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 39923831 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 44292640 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Total </td>
   <td style="text-align:left;"> 237972668 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 201853812 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 181336262 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 212798362 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 245536180 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 207250556 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 175117342 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 209987916 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 155465284 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 120495038 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 165482296 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 187509020 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
</tbody>
</table>

# **Differentially Expressed Sequence Identification**

Programs Used

+ R
+ RStudio

Packages

+ readr (***install.packages('readr')***)
+ limma (***BiocManager::install('limma')***)
+ DESeq2 (***BiocManager::install('DESeq2')***)
+ dplyr (***install.packages("dplyr")***)
+ ggplot2 (***install.packages("ggplot2")***)
+ gplots (***install.packages("gplots")***)
+ Annotations (***BiocManager::install('AnnotationDbi')***)
+ org.Hs.eg.db (***BiocManager::install('org.Hs.eg.db')***)
  + This is for Human
+ org.Mm.eg.db (***BiocManager::install('org.Mm.eg.db')***)
  + This is for Mouse
+ ggrepel (***install.packages("ggrepel")***)
+ ReportingTools (***BiocManager::install('ReportingTools')***)
+ GO.db (***BiocManager::install('GO.db')***)
+ GOstats (***BiocManager::install('GOstats')***)
+ pathview (***BiocManager::install('pathview')***)
+ gage (***BiocManager::install('gage')***)
+ gageData (***BiocManager::install('gageData')***)
+ select (***BiocManager::install('Select')***)

With these, you most certainly will have to step through each and install extra things when you start calling the packages. Take it step by step to ensure that each dependency is installed. 

## Analyzing Reads Counts


```r
library(limma)
library(DESeq2)
library(dplyr)
library(readr)

countData = read_csv("readCounts.csv", skip = 1)
```

This gives us our dataframe from out featureCounts program, but if we look at the data we see that featureCounts added some extra information that characterizes each gene_id.


```r
head(countData)
```

```
## # A tibble: 6 x 18
##   Geneid Chr   Start End   Strand Length SRR2121770Align~ SRR2121771Align~
##   <chr>  <chr> <chr> <chr> <chr>   <dbl>            <dbl>            <dbl>
## 1 ENSMU~ chr1  3073~ 3074~ +        1070                0                0
## 2 ENSMU~ chr1  3102~ 3102~ +         110                0                0
## 3 ENSMU~ chr1~ 3205~ 3207~ -;-;-~   6094                0                0
## 4 ENSMU~ chr1  3252~ 3253~ +         480                0                0
## 5 ENSMU~ chr1  3365~ 3368~ -        2819                0                0
## 6 ENSMU~ chr1  3375~ 3377~ -        2233                0                0
## # ... with 10 more variables:
## #   SRR2121774Aligned.sortedByCoord.out.bam <dbl>,
## #   SRR2121775Aligned.sortedByCoord.out.bam <dbl>,
## #   SRR2121778Aligned.sortedByCoord.out.bam <dbl>,
## #   SRR2121779Aligned.sortedByCoord.out.bam <dbl>,
## #   SRR2121780Aligned.sortedByCoord.out.bam <dbl>,
## #   SRR2121781Aligned.sortedByCoord.out.bam <dbl>,
## #   SRR2121786Aligned.sortedByCoord.out.bam <dbl>,
## #   SRR2121787Aligned.sortedByCoord.out.bam <dbl>,
## #   SRR2121788Aligned.sortedByCoord.out.bam <dbl>,
## #   SRR2121789Aligned.sortedByCoord.out.bam <dbl>
```

We also need to set out row names to the gene_id. We will do some data frame manipulation and then look at the data again.


```r
countData = as.data.frame(countData)
rownames(countData) = countData$Geneid
countData = countData[,-c(1:6)]
head(countData)
```

```
##                      SRR2121770Aligned.sortedByCoord.out.bam
## ENSMUSG00000102693.1                                       0
## ENSMUSG00000064842.1                                       0
## ENSMUSG00000051951.5                                       0
## ENSMUSG00000102851.1                                       0
## ENSMUSG00000103377.1                                       0
## ENSMUSG00000104017.1                                       0
##                      SRR2121771Aligned.sortedByCoord.out.bam
## ENSMUSG00000102693.1                                       0
## ENSMUSG00000064842.1                                       0
## ENSMUSG00000051951.5                                       0
## ENSMUSG00000102851.1                                       0
## ENSMUSG00000103377.1                                       0
## ENSMUSG00000104017.1                                       0
##                      SRR2121774Aligned.sortedByCoord.out.bam
## ENSMUSG00000102693.1                                       0
## ENSMUSG00000064842.1                                       0
## ENSMUSG00000051951.5                                       0
## ENSMUSG00000102851.1                                       0
## ENSMUSG00000103377.1                                       0
## ENSMUSG00000104017.1                                       0
##                      SRR2121775Aligned.sortedByCoord.out.bam
## ENSMUSG00000102693.1                                       0
## ENSMUSG00000064842.1                                       0
## ENSMUSG00000051951.5                                       4
## ENSMUSG00000102851.1                                       0
## ENSMUSG00000103377.1                                       2
## ENSMUSG00000104017.1                                       0
##                      SRR2121778Aligned.sortedByCoord.out.bam
## ENSMUSG00000102693.1                                       0
## ENSMUSG00000064842.1                                       0
## ENSMUSG00000051951.5                                       0
## ENSMUSG00000102851.1                                       0
## ENSMUSG00000103377.1                                       0
## ENSMUSG00000104017.1                                       0
##                      SRR2121779Aligned.sortedByCoord.out.bam
## ENSMUSG00000102693.1                                       0
## ENSMUSG00000064842.1                                       0
## ENSMUSG00000051951.5                                       0
## ENSMUSG00000102851.1                                       0
## ENSMUSG00000103377.1                                       0
## ENSMUSG00000104017.1                                       0
##                      SRR2121780Aligned.sortedByCoord.out.bam
## ENSMUSG00000102693.1                                       0
## ENSMUSG00000064842.1                                       0
## ENSMUSG00000051951.5                                       0
## ENSMUSG00000102851.1                                       0
## ENSMUSG00000103377.1                                       0
## ENSMUSG00000104017.1                                       0
##                      SRR2121781Aligned.sortedByCoord.out.bam
## ENSMUSG00000102693.1                                       0
## ENSMUSG00000064842.1                                       0
## ENSMUSG00000051951.5                                       0
## ENSMUSG00000102851.1                                       0
## ENSMUSG00000103377.1                                       0
## ENSMUSG00000104017.1                                       0
##                      SRR2121786Aligned.sortedByCoord.out.bam
## ENSMUSG00000102693.1                                       4
## ENSMUSG00000064842.1                                       0
## ENSMUSG00000051951.5                                       4
## ENSMUSG00000102851.1                                       2
## ENSMUSG00000103377.1                                       6
## ENSMUSG00000104017.1                                       0
##                      SRR2121787Aligned.sortedByCoord.out.bam
## ENSMUSG00000102693.1                                       2
## ENSMUSG00000064842.1                                       0
## ENSMUSG00000051951.5                                       2
## ENSMUSG00000102851.1                                       0
## ENSMUSG00000103377.1                                      14
## ENSMUSG00000104017.1                                       0
##                      SRR2121788Aligned.sortedByCoord.out.bam
## ENSMUSG00000102693.1                                       0
## ENSMUSG00000064842.1                                       0
## ENSMUSG00000051951.5                                       0
## ENSMUSG00000102851.1                                       0
## ENSMUSG00000103377.1                                       0
## ENSMUSG00000104017.1                                       0
##                      SRR2121789Aligned.sortedByCoord.out.bam
## ENSMUSG00000102693.1                                       0
## ENSMUSG00000064842.1                                       0
## ENSMUSG00000051951.5                                       0
## ENSMUSG00000102851.1                                       0
## ENSMUSG00000103377.1                                       0
## ENSMUSG00000104017.1                                       0
```

### Quick Data Exploration


```r
dim(countData)
```

```
## [1] 55487    12
```


```r
summary(countData)
```

```
##  SRR2121770Aligned.sortedByCoord.out.bam
##  Min.   :      0.0                      
##  1st Qu.:      0.0                      
##  Median :      0.0                      
##  Mean   :    273.8                      
##  3rd Qu.:     25.0                      
##  Max.   :1085454.0                      
##  SRR2121771Aligned.sortedByCoord.out.bam
##  Min.   :     0.0                       
##  1st Qu.:     0.0                       
##  Median :     0.0                       
##  Mean   :   235.7                       
##  3rd Qu.:    20.0                       
##  Max.   :834016.0                       
##  SRR2121774Aligned.sortedByCoord.out.bam
##  Min.   :     0.0                       
##  1st Qu.:     0.0                       
##  Median :     0.0                       
##  Mean   :   229.6                       
##  3rd Qu.:    15.0                       
##  Max.   :783768.0                       
##  SRR2121775Aligned.sortedByCoord.out.bam
##  Min.   :      0.0                      
##  1st Qu.:      0.0                      
##  Median :      0.0                      
##  Mean   :    269.6                      
##  3rd Qu.:     24.0                      
##  Max.   :1032589.0                      
##  SRR2121778Aligned.sortedByCoord.out.bam
##  Min.   :      0.0                      
##  1st Qu.:      0.0                      
##  Median :      0.0                      
##  Mean   :    286.6                      
##  3rd Qu.:     25.0                      
##  Max.   :1233905.0                      
##  SRR2121779Aligned.sortedByCoord.out.bam
##  Min.   :      0.0                      
##  1st Qu.:      0.0                      
##  Median :      0.0                      
##  Mean   :    257.7                      
##  3rd Qu.:     24.0                      
##  Max.   :1232421.0                      
##  SRR2121780Aligned.sortedByCoord.out.bam
##  Min.   :     0.0                       
##  1st Qu.:     0.0                       
##  Median :     0.0                       
##  Mean   :   193.7                       
##  3rd Qu.:    17.0                       
##  Max.   :805463.0                       
##  SRR2121781Aligned.sortedByCoord.out.bam
##  Min.   :      0.0                      
##  1st Qu.:      0.0                      
##  Median :      0.0                      
##  Mean   :    268.2                      
##  3rd Qu.:     26.0                      
##  Max.   :1210716.0                      
##  SRR2121786Aligned.sortedByCoord.out.bam
##  Min.   :     0.0                       
##  1st Qu.:     0.0                       
##  Median :     3.0                       
##  Mean   :   160.2                       
##  3rd Qu.:    21.0                       
##  Max.   :641911.0                       
##  SRR2121787Aligned.sortedByCoord.out.bam
##  Min.   :     0.0                       
##  1st Qu.:     0.0                       
##  Median :     4.0                       
##  Mean   :   119.4                       
##  3rd Qu.:    21.0                       
##  Max.   :516752.0                       
##  SRR2121788Aligned.sortedByCoord.out.bam
##  Min.   :     0.0                       
##  1st Qu.:     0.0                       
##  Median :     0.0                       
##  Mean   :   184.3                       
##  3rd Qu.:    16.0                       
##  Max.   :706363.0                       
##  SRR2121789Aligned.sortedByCoord.out.bam
##  Min.   :     0.0                       
##  1st Qu.:     0.0                       
##  Median :     1.0                       
##  Mean   :   211.1                       
##  3rd Qu.:    18.0                       
##  Max.   :756086.0
```

Let's also go ahead and change the names to describe out data a little better.


```r
columns = c('Trp53m_mock_1', 'Trp53m_mock_2', 'Trp53m_4h7Gy_1', 'Trp53m_4h7Gy_2', 'Trp53p_mock_1', 'Trp53p_mock_2', 'Trp53p_mock_3', 'Trp53p_mock_4', 'Trp53p_4h7Gy_1', 'Trp53p_4h7Gy_2', 'Trp53p_4h7Gy_3', 'Trp53p_4h7Gy_4')
colnames(countData) = columns
```

Here we have to make sure that we convert the *+/+* and *-/-* to characters. These characters


```r
barplot( colSums(countData)/1e6, col="green",las=3,main="Total read counts (millions)", ylab="Total reads in millions")
```

![](Mouse_RNA_Seq_p53_genotoxic_files/figure-html/unnamed-chunk-35-1.png)<!-- -->


```r
hist(countData[,1], br=200)
```

![](Mouse_RNA_Seq_p53_genotoxic_files/figure-html/unnamed-chunk-36-1.png)<!-- -->

We can see that our count data is highly skewed to the right. This is a great case for using **log** transformation!


```r
logCountData = log2(1+countData)
par(mfrow = c(1, 2))  # two columns
hist(logCountData[,1])
boxplot(logCountData,las=3)
```

![](Mouse_RNA_Seq_p53_genotoxic_files/figure-html/unnamed-chunk-37-1.png)<!-- -->


```r
x <- logCountData
myColors = rainbow(dim(x)[2])
plot(density(x[,1]),col = myColors[1], lwd=2,
     xlab="Expresson values", ylab="Density", main= "Distribution of transformed data",
     ylim=c(0, max(density(x[,1])$y)+.02 ) )
  
for( i in 2:dim(x)[2] )
lines(density(x[,i]),col=myColors[i], lwd=2)
legend("topright", cex=1.1,colnames(x), lty=rep(1,dim(x)[2]), col=myColors )	
```

![](Mouse_RNA_Seq_p53_genotoxic_files/figure-html/unnamed-chunk-38-1.png)<!-- -->


```r
plot(logCountData[,1],logCountData[,2])
```

![](Mouse_RNA_Seq_p53_genotoxic_files/figure-html/unnamed-chunk-39-1.png)<!-- -->

### Filtering, Normalization, and Trasformation using DESeq2


```r
detectGroups <- function (x){  # x are col names
  tem <- gsub("[0-9]*$","",x) # Remove all numbers from end
  #tem = gsub("_Rep|_rep|_REP","",tem)
  tem <- gsub("_$","",tem); # remove "_" from end
  tem <- gsub("_Rep$","",tem); # remove "_Rep" from end
  tem <- gsub("_rep$","",tem); # remove "_rep" from end
  tem <- gsub("_REP$","",tem)  # remove "_REP" from end
  return( tem )
}

groups = as.character ( detectGroups( colnames( countData ) ) )
groups
```

```
##  [1] "Trp53m_mock"  "Trp53m_mock"  "Trp53m_4h7Gy" "Trp53m_4h7Gy"
##  [5] "Trp53p_mock"  "Trp53p_mock"  "Trp53p_mock"  "Trp53p_mock" 
##  [9] "Trp53p_4h7Gy" "Trp53p_4h7Gy" "Trp53p_4h7Gy" "Trp53p_4h7Gy"
```


```r
colData = cbind(colnames(countData), groups )
colData
```

```
##                        groups        
##  [1,] "Trp53m_mock_1"  "Trp53m_mock" 
##  [2,] "Trp53m_mock_2"  "Trp53m_mock" 
##  [3,] "Trp53m_4h7Gy_1" "Trp53m_4h7Gy"
##  [4,] "Trp53m_4h7Gy_2" "Trp53m_4h7Gy"
##  [5,] "Trp53p_mock_1"  "Trp53p_mock" 
##  [6,] "Trp53p_mock_2"  "Trp53p_mock" 
##  [7,] "Trp53p_mock_3"  "Trp53p_mock" 
##  [8,] "Trp53p_mock_4"  "Trp53p_mock" 
##  [9,] "Trp53p_4h7Gy_1" "Trp53p_4h7Gy"
## [10,] "Trp53p_4h7Gy_2" "Trp53p_4h7Gy"
## [11,] "Trp53p_4h7Gy_3" "Trp53p_4h7Gy"
## [12,] "Trp53p_4h7Gy_4" "Trp53p_4h7Gy"
```


```r
colData = as.data.frame(cbind(colnames(countData), groups ))
sampleID = as.factor( c(1:2,1:2,1:4,1:4) )
colData = cbind(data.frame(colData),sampleID)

colData
```

```
##                V1       groups sampleID
## 1   Trp53m_mock_1  Trp53m_mock        1
## 2   Trp53m_mock_2  Trp53m_mock        2
## 3  Trp53m_4h7Gy_1 Trp53m_4h7Gy        1
## 4  Trp53m_4h7Gy_2 Trp53m_4h7Gy        2
## 5   Trp53p_mock_1  Trp53p_mock        1
## 6   Trp53p_mock_2  Trp53p_mock        2
## 7   Trp53p_mock_3  Trp53p_mock        3
## 8   Trp53p_mock_4  Trp53p_mock        4
## 9  Trp53p_4h7Gy_1 Trp53p_4h7Gy        1
## 10 Trp53p_4h7Gy_2 Trp53p_4h7Gy        2
## 11 Trp53p_4h7Gy_3 Trp53p_4h7Gy        3
## 12 Trp53p_4h7Gy_4 Trp53p_4h7Gy        4
```

```r
str(colData)
```

```
## 'data.frame':	12 obs. of  3 variables:
##  $ V1      : Factor w/ 12 levels "Trp53m_4h7Gy_1",..: 3 4 1 2 9 10 11 12 5 6 ...
##  $ groups  : Factor w/ 4 levels "Trp53m_4h7Gy",..: 2 2 1 1 4 4 4 4 3 3 ...
##  $ sampleID: Factor w/ 4 levels "1","2","3","4": 1 2 1 2 1 2 3 4 1 2 ...
```

Creating a DESeq Dataset


```r
dds = DESeqDataSetFromMatrix(countData=countData,
                             colData=colData,
                             design= ~ sampleID+groups)   # note that the study design is changed.
dds = DESeq(dds)  # main function
nrow(dds)
```

```
## [1] 55487
```

Filtering


```r
dds <- dds[ rowSums(counts(dds)) > 5, ]
nrow(dds)
```

```
## [1] 34139
```

Regularized log transformation


```r
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)
```

```
##                      Trp53m_mock_1 Trp53m_mock_2 Trp53m_4h7Gy_1
## ENSMUSG00000102693.1    -0.7468115    -0.7465790     -0.7462666
## ENSMUSG00000051951.5    -0.1371545    -0.1368755     -0.1364994
## ENSMUSG00000103377.1     1.0412975     1.0422670      1.0435633
##                      Trp53m_4h7Gy_2 Trp53p_mock_1 Trp53p_mock_2
## ENSMUSG00000102693.1     -0.7467248    -0.7468175    -0.7467610
## ENSMUSG00000051951.5     -0.0655419    -0.1371616    -0.1370939
## ENSMUSG00000103377.1      1.0866746     1.0412726     1.0415087
##                      Trp53p_mock_3 Trp53p_mock_4 Trp53p_4h7Gy_1
## ENSMUSG00000102693.1    -0.7462584    -0.7468557    -0.62955781
## ENSMUSG00000051951.5    -0.1364894    -0.1372074    -0.03185913
## ENSMUSG00000103377.1     1.0435975     1.0411128     1.22758742
##                      Trp53p_4h7Gy_2 Trp53p_4h7Gy_3 Trp53p_4h7Gy_4
## ENSMUSG00000102693.1    -0.68345459     -0.7458901     -0.7462950
## ENSMUSG00000051951.5    -0.06916313     -0.1363120     -0.1365337
## ENSMUSG00000103377.1     1.49681594      1.0442045      1.0434458
```

Variance Stabilizing Transformation


```r
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
```

```
##                      Trp53m_mock_1 Trp53m_mock_2 Trp53m_4h7Gy_1
## ENSMUSG00000102693.1      6.914084      6.914084       6.914084
## ENSMUSG00000051951.5      6.914084      6.914084       6.914084
## ENSMUSG00000103377.1      6.914084      6.914084       6.914084
##                      Trp53m_4h7Gy_2 Trp53p_mock_1 Trp53p_mock_2
## ENSMUSG00000102693.1       6.914084      6.914084      6.914084
## ENSMUSG00000051951.5       7.152884      6.914084      6.914084
## ENSMUSG00000103377.1       7.083037      6.914084      6.914084
##                      Trp53p_mock_3 Trp53p_mock_4 Trp53p_4h7Gy_1
## ENSMUSG00000102693.1      6.914084      6.914084       7.211668
## ENSMUSG00000051951.5      6.914084      6.914084       7.211668
## ENSMUSG00000103377.1      6.914084      6.914084       7.278228
##                      Trp53p_4h7Gy_2 Trp53p_4h7Gy_3 Trp53p_4h7Gy_4
## ENSMUSG00000102693.1       7.152550       6.914084       6.914084
## ENSMUSG00000051951.5       7.152550       6.914084       6.914084
## ENSMUSG00000103377.1       7.540786       6.914084       6.914084
```

For the log2 approach, we need to first estimate size factors to account for sequencing depth, and then specify normalized=TRUE. Sequencing depth correction is done automatically for the rlog and the vst.


```r
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
```

```
##  Trp53m_mock_1  Trp53m_mock_2 Trp53m_4h7Gy_1 Trp53m_4h7Gy_2  Trp53p_mock_1 
##      1.2892566      1.0903014      0.8973654      1.2078255      1.2952327 
##  Trp53p_mock_2  Trp53p_mock_3  Trp53p_mock_4 Trp53p_4h7Gy_1 Trp53p_4h7Gy_2 
##      1.2406213      0.8931001      1.3347016      0.7767888      0.6056051 
## Trp53p_4h7Gy_3 Trp53p_4h7Gy_4 
##      0.8227361      0.9123249
```

We will first look at the log transformed data


```r
slog <- log2(counts(dds, normalized=TRUE)+1)
head(slog)
```

```
##                       Trp53m_mock_1 Trp53m_mock_2 Trp53m_4h7Gy_1
## ENSMUSG00000102693.1       0.000000      0.000000              0
## ENSMUSG00000051951.5       0.000000      0.000000              0
## ENSMUSG00000103377.1       0.000000      0.000000              0
## ENSMUSG00000103201.1       0.000000      0.000000              0
## ENSMUSG00000102592.1       0.000000      0.000000              0
## ENSMUSG00000025900.12      1.734188      1.503021              0
##                       Trp53m_4h7Gy_2 Trp53p_mock_1 Trp53p_mock_2
## ENSMUSG00000102693.1        0.000000     0.0000000             0
## ENSMUSG00000051951.5        2.108269     0.0000000             0
## ENSMUSG00000103377.1        1.409184     0.0000000             0
## ENSMUSG00000103201.1        0.000000     0.0000000             0
## ENSMUSG00000102592.1        0.000000     0.0000000             0
## ENSMUSG00000025900.12       0.000000     0.8254291             0
##                       Trp53p_mock_3 Trp53p_mock_4 Trp53p_4h7Gy_1
## ENSMUSG00000102693.1              0      0.000000       2.620447
## ENSMUSG00000051951.5              0      0.000000       2.620447
## ENSMUSG00000103377.1              0      0.000000       3.125007
## ENSMUSG00000103201.1              0      0.000000       2.281566
## ENSMUSG00000102592.1              0      0.000000       3.922280
## ENSMUSG00000025900.12             0      2.246759       3.922280
##                       Trp53p_4h7Gy_2 Trp53p_4h7Gy_3 Trp53p_4h7Gy_4
## ENSMUSG00000102693.1        2.105169              0       0.000000
## ENSMUSG00000051951.5        2.105169              0       0.000000
## ENSMUSG00000103377.1        4.592001              0       0.000000
## ENSMUSG00000103201.1        3.447241              0       0.000000
## ENSMUSG00000102592.1        2.926941              0       1.674552
## ENSMUSG00000025900.12       4.777149              0       1.674552
```


```r
par(mfrow = c(1, 3))  # 3 columns
plot(slog[,1],slog[,2])
plot(assay(rld)[,1],assay(rld)[,2])
plot(assay(vsd)[,1],assay(vsd)[,2])
```

![](Mouse_RNA_Seq_p53_genotoxic_files/figure-html/unnamed-chunk-49-1.png)<!-- -->


```r
par(mfrow = c(1, 3))  # 3 columns
slog <- log2(counts(dds, normalized=TRUE)+1)
plot(slog[,1],slog[,2])
slog <- log2(counts(dds, normalized=TRUE)+4)
plot(slog[,1],slog[,2], xlim=c(0,20))
slog <- log2(counts(dds, normalized=TRUE)+20)
plot(slog[,1],slog[,2],xlim=c(0,20))
```

![](Mouse_RNA_Seq_p53_genotoxic_files/figure-html/unnamed-chunk-50-1.png)<!-- -->


```r
library("dplyr")
library("ggplot2")

df <- bind_rows(
  as_data_frame(slog[,1:2]) %>%
         mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))
```

```
## Warning: `as_data_frame()` is deprecated, use `as_tibble()` (but mind the new semantics).
## This warning is displayed once per session.
```

```r
colnames(df)[1:2] <- c("x", "y")  

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)
```

![](Mouse_RNA_Seq_p53_genotoxic_files/figure-html/unnamed-chunk-51-1.png)<!-- -->

### Exploratory Data Analysis

PCA plot


```r
plotPCA(rld, intgroup = c("groups")) + theme(aspect.ratio=1)
```

![](Mouse_RNA_Seq_p53_genotoxic_files/figure-html/unnamed-chunk-52-1.png)<!-- -->

A prettier PCA plot created with GGPlot2


```r
pca.object <- prcomp(t(assay(rld))) # PCA 
pcaData = as.data.frame(pca.object$x[,1:2]); 
pcaData = cbind(pcaData,detectGroups(colnames(assay(rld)) ))
colnames(pcaData) = c("PC1", "PC2", "Type")
percentVar=round(100*summary(pca.object)$importance[2,1:2],0)
#plot
p=ggplot(pcaData, aes(PC1, PC2, color=Type, shape = Type)) + geom_point(size=5) 
p=p+xlab(paste0("PC1: ",percentVar[1],"% variance")) 
p=p+ylab(paste0("PC2: ",percentVar[2],"% variance")) 
p=p+ggtitle("Principal component analysis (PCA)")+coord_fixed(ratio=1.0)+ 
    theme(plot.title = element_text(size = 16,hjust = 0.5)) + theme(aspect.ratio=1) +
    theme(axis.text.x = element_text( size = 16),
    axis.text.y = element_text( size = 16),
    axis.title.x = element_text( size = 16),
    axis.title.y = element_text( size = 16) ) +
  theme(legend.text=element_text(size=16))
print(p)
```

![](Mouse_RNA_Seq_p53_genotoxic_files/figure-html/unnamed-chunk-53-1.png)<!-- -->

Multidimensional Scaling Plot


```r
dist2 <- function(x, ...)   # distance function = 1-PCC (Pearson's correlation coefficient)
  as.dist(1-cor(t(x), method="pearson"))

fit = cmdscale( dist2(t(assay(rld))) , eig=T, k=2)
mdsData <- as.data.frame(fit$points[,1:2]); 
mdsData <- cbind(mdsData,detectGroups(colnames(assay(rld))) )
colnames(mdsData) = c("x1", "x2", "Type")
	
p<-ggplot(mdsData, aes(x1, x2, color=Type, shape = Type)) + geom_point(size=5) 
p=p+xlab("Dimension 1") 
p=p+ylab("Dimension 2") 
p=p+ggtitle("Multidimensional scaling (MDS)")+ coord_fixed(ratio=1.)+ 
     theme(plot.title = element_text(hjust = 0.5)) + theme(aspect.ratio=1) +
	 	 theme(axis.text.x = element_text( size = 16),
        axis.text.y = element_text( size = 16),
        axis.title.x = element_text( size = 16),
        axis.title.y = element_text( size = 16) ) +
	   theme(legend.text=element_text(size=16))
print(p)
```

![](Mouse_RNA_Seq_p53_genotoxic_files/figure-html/unnamed-chunk-54-1.png)<!-- -->

Creating a heatmap


```r
library(gplots)

hclust2 <- function(x, method="average", ...)  # average linkage in hierarchical clustering
  hclust(x, method=method, ...)

n=50 # number of top genes by standard deviation

x = assay(rld)
if(n>dim(x)[1]) n = dim(x)[1] # max	as data

x = x[order(apply(x,1,sd),decreasing=TRUE),]  # sort genes by standard deviation

x = x[1:n,]   # only keep the n genes

# this will cutoff very large values, which could skew the color 
x=as.matrix(x[1:n,])-apply(x[1:n,],1,mean)
cutoff = median(unlist(x)) + 4*sd (unlist(x)) 
x[x>cutoff] <- cutoff
cutoff = median(unlist(x)) - 4*sd (unlist(x)) 
x[x< cutoff] <- cutoff
	
groups = detectGroups(colnames(x) )
groups.colors = rainbow(length(unique(groups) ) )


	lmat = rbind(c(5,4),c(0,1),c(3,2))
	lwid = c(1.5,4)
	lhei = c(1,.2,4)

heatmap.2(x, distfun = dist2,hclustfun=hclust2,
	 col=greenred(75), density.info="none", trace="none", scale="none", keysize=.5
	,key=T, symkey=F
	,ColSideColors=groups.colors[ as.factor(groups)]
	,margins=c(8,8)
	,cexRow=1
	,srtCol=45
	,cexCol=1.  # size of font for sample names
	,lmat = lmat, lwid = lwid, lhei = lhei
	)
```

![](Mouse_RNA_Seq_p53_genotoxic_files/figure-html/heatmap-1.png)<!-- -->

K-means clustering of genes



### Identification of Differentially Expressed Genes


```r
dds <- DESeq(dds)
res <- results(dds)
res
```

```
## log2 fold change (MLE): groups Trp53p mock vs Trp53m 4h7Gy 
## Wald test p-value: groups Trp53p mock vs Trp53m 4h7Gy 
## DataFrame with 34139 rows and 6 columns
##                               baseMean     log2FoldChange            lfcSE
##                              <numeric>          <numeric>        <numeric>
## ENSMUSG00000102693.1 0.704323893476708  0.944373860101389 5.04011012102995
## ENSMUSG00000051951.5 0.980301949412825  -1.93606050718611 4.90304157802561
## ENSMUSG00000103377.1  2.70811247007369 0.0164654841466683 4.92883452955993
## ENSMUSG00000103201.1  1.14745830269913   1.21394868463973 5.02701308026115
## ENSMUSG00000102592.1   1.9131690413215  0.642969613528605 5.02743169014771
## ...                                ...                ...              ...
## ENSMUSG00000102053.1   5.4606917322249  -5.79608463544145 3.08302312189078
## ENSMUSG00000095852.7  5.67546279438922  -5.99404661143399 4.83967404452342
## ENSMUSG00000100145.1  3.36574617296549  -5.16394977182676  4.3023162124544
## ENSMUSG00000099405.1  0.94707084987755  -3.02778143921305 4.95303595408583
## ENSMUSG00000091987.8  1.02158581031203  -2.48887433546013 4.96776834749799
##                                     stat             pvalue
##                                <numeric>          <numeric>
## ENSMUSG00000102693.1   0.187371671932518  0.851369231767202
## ENSMUSG00000051951.5  -0.394869281929633   0.69293936113219
## ENSMUSG00000103377.1 0.00334064453734835  0.997334556258204
## ENSMUSG00000103201.1   0.241485085727421  0.809179174368776
## ENSMUSG00000102592.1   0.127892262522162  0.898234234374561
## ...                                  ...                ...
## ENSMUSG00000102053.1   -1.88000037829323 0.0601080263659187
## ENSMUSG00000095852.7   -1.23852279229773  0.215522276221424
## ENSMUSG00000100145.1   -1.20027202019184  0.230033712629537
## ENSMUSG00000099405.1  -0.611298094195216  0.541002252458883
## ENSMUSG00000091987.8  -0.501004507731051  0.616367950479261
##                                   padj
##                              <numeric>
## ENSMUSG00000102693.1                NA
## ENSMUSG00000051951.5                NA
## ENSMUSG00000103377.1                NA
## ENSMUSG00000103201.1                NA
## ENSMUSG00000102592.1                NA
## ...                                ...
## ENSMUSG00000102053.1 0.143688035730787
## ENSMUSG00000095852.7 0.376681084450195
## ENSMUSG00000100145.1                NA
## ENSMUSG00000099405.1                NA
## ENSMUSG00000091987.8                NA
```

DESeq2 uses the Benjamini-Hochberg (BH) adjustment (Benjamini and Hochberg 1995) as implemented in the base R p.adjust function


```r
res <- results(dds, alpha = 0.5, lfcThreshold=.01)
summary(res)
```

```
## 
## out of 34139 with nonzero total read count
## adjusted p-value < 0.5
## LFC > 0.01 (up)    : 6203, 18%
## LFC < -0.01 (down) : 5915, 17%
## outliers [1]       : 0, 0%
## low counts [2]     : 11914, 35%
## (mean count < 3)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

Now lets sort genes by fold change


```r
res <- res[order(abs( res$log2FoldChange), decreasing=TRUE),]
head(res)
```

```
## log2 fold change (MLE): groups Trp53p mock vs Trp53m 4h7Gy 
## Wald test p-value: groups Trp53p mock vs Trp53m 4h7Gy 
## DataFrame with 6 rows and 6 columns
##                               baseMean    log2FoldChange            lfcSE
##                              <numeric>         <numeric>        <numeric>
## ENSMUSG00000044022.5  1.04775840423632 -22.8083207137305 5.01253450505715
## ENSMUSG00000045871.5  1.61728198380297 -22.7627770976699 4.12039782777052
## ENSMUSG00000110839.1 0.865416000692208 -22.6593426505353 4.97255213974549
## ENSMUSG00000025527.9  1.69240483905238 -22.3509710830054 4.24986810302321
## ENSMUSG00000104025.1 0.728198191952891 -22.1641354425874 4.88426986117945
## ENSMUSG00000104283.1 0.728198191952891 -22.1641354425874 4.88426986117945
##                                   stat               pvalue      padj
##                              <numeric>            <numeric> <numeric>
## ENSMUSG00000044022.5 -4.54826209988765 5.40907435855941e-06        NA
## ENSMUSG00000045871.5 -5.52198550934124 3.35190139219697e-08        NA
## ENSMUSG00000110839.1 -4.55487283270488 5.24172787416592e-06        NA
## ENSMUSG00000025527.9 -5.25686222287999 1.46534006583066e-07        NA
## ENSMUSG00000104025.1 -4.53581314551642 5.73819675582155e-06        NA
## ENSMUSG00000104283.1 -4.53581314551642 5.73819675582155e-06        NA
```

MA Plot


```r
DESeq2::plotMA(res,  ylim = c(-5, 5))
```

![](Mouse_RNA_Seq_p53_genotoxic_files/figure-html/unnamed-chunk-58-1.png)<!-- -->

Volcano plot


```r
library(dplyr)
res1 = as.data.frame(res)
# add a new column using the mutate function in dplyr
res1 = mutate(res1, sig=ifelse(res1$padj<0.05, "FDR<0.05", "Not Sig"))
res1[which(abs(res1$log2FoldChange)<0.5),'sig'] <- "Not Sig"


p = ggplot(res1, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(col=sig)) +
  scale_color_manual(values=c("red", "black"))
p
```

```
## Warning: Removed 11914 rows containing missing values (geom_point).
```

![](Mouse_RNA_Seq_p53_genotoxic_files/figure-html/unnamed-chunk-59-1.png)<!-- -->

### Gene Annotations

Plot counts of top gene


```r
topGene <- rownames(res)[1]
plotCounts(dds, gene = topGene, intgroup=c("groups"))
```

![](Mouse_RNA_Seq_p53_genotoxic_files/figure-html/unnamed-chunk-60-1.png)<!-- -->

Here we see an interesting point or our normalized counts under the Trp53p_4h7Gy group that seems to be extremely high, while the other 3 replicates are around 0.5. I cannot get this portion to work, however. I get an error stating that None of the keys entered are valid keys for 'SYMBOL'. 

Let's look at the keys that we have to work with for our ***res*** file from DESeq2


```r
head(row.names(res))
```

```
## [1] "ENSMUSG00000044022.5" "ENSMUSG00000045871.5" "ENSMUSG00000110839.1"
## [4] "ENSMUSG00000025527.9" "ENSMUSG00000104025.1" "ENSMUSG00000104283.1"
```

Now we need to find the same key in the Mm database.


```r
library(AnnotationDbi)
library(org.Mm.eg.db)

columns(org.Mm.eg.db)
```

```
##  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT" 
##  [5] "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"    
##  [9] "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"       
## [13] "IPI"          "MGI"          "ONTOLOGY"     "ONTOLOGYALL" 
## [17] "PATH"         "PFAM"         "PMID"         "PROSITE"     
## [21] "REFSEQ"       "SYMBOL"       "UNIGENE"      "UNIPROT"
```

```r
#key = gsub("\\..*","", row.names(res))
res$symbol <- gsub("\\..*","", row.names(res))
```


```r
key.en = keys(org.Mm.eg.db, keytype="ENSEMBL")
head(key.en)
```

```
## [1] "ENSMUSG00000030359" "ENSMUSG00000020804" "ENSMUSG00000025375"
## [4] "ENSMUSG00000015243" "ENSMUSG00000028125" "ENSMUSG00000026944"
```

These are ENSEMBL symbols, so we need to designate that when looking for the genes that we have. 


```r
library(biomaRt)
mart = useMart("ENSEMBL_MART_ENSEMBL")

mart = useDataset('mmusculus_gene_ensembl', mart)

ens = row.names(res)
ensLookup <- gsub("\\.[0-9]*$", "", ens)
#ensLookup

annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "gene_biotype", "external_gene_name","entrezgene_id"),
  filter="ensembl_gene_id",
  values=ensLookup,
  uniqueRows=TRUE)

annotLookup <- data.frame(
  ens[match(annotLookup$ensembl_gene_id, ensLookup)],
  annotLookup)

colnames(annotLookup) <- c(
  "original_id",
  c("ensembl_gene_id", "gene_biotype", "external_gene_name", "entrezgene_id"))

head(annotLookup)
```

```
##             original_id    ensembl_gene_id   gene_biotype
## 1 ENSMUSG00000001366.14 ENSMUSG00000001366 protein_coding
## 2  ENSMUSG00000002007.5 ENSMUSG00000002007 protein_coding
## 3  ENSMUSG00000005864.7 ENSMUSG00000005864 protein_coding
## 4  ENSMUSG00000006651.8 ENSMUSG00000006651 protein_coding
## 5 ENSMUSG00000010376.15 ENSMUSG00000010376 protein_coding
## 6  ENSMUSG00000011463.5 ENSMUSG00000011463 protein_coding
##   external_gene_name entrezgene_id
## 1              Fbxo9         71538
## 2              Srpk3         56504
## 3              Cnga2         12789
## 4              Aplp1         11803
## 5              Nedd8         18002
## 6               Cpb1         76703
```

```r
head(res)
```

```
## log2 fold change (MLE): groups Trp53p mock vs Trp53m 4h7Gy 
## Wald test p-value: groups Trp53p mock vs Trp53m 4h7Gy 
## DataFrame with 6 rows and 7 columns
##                               baseMean    log2FoldChange            lfcSE
##                              <numeric>         <numeric>        <numeric>
## ENSMUSG00000044022.5  1.04775840423632 -22.8083207137305 5.01253450505715
## ENSMUSG00000045871.5  1.61728198380297 -22.7627770976699 4.12039782777052
## ENSMUSG00000110839.1 0.865416000692208 -22.6593426505353 4.97255213974549
## ENSMUSG00000025527.9  1.69240483905238 -22.3509710830054 4.24986810302321
## ENSMUSG00000104025.1 0.728198191952891 -22.1641354425874 4.88426986117945
## ENSMUSG00000104283.1 0.728198191952891 -22.1641354425874 4.88426986117945
##                                   stat               pvalue      padj
##                              <numeric>            <numeric> <numeric>
## ENSMUSG00000044022.5 -4.54826209988765 5.40907435855941e-06        NA
## ENSMUSG00000045871.5 -5.52198550934124 3.35190139219697e-08        NA
## ENSMUSG00000110839.1 -4.55487283270488 5.24172787416592e-06        NA
## ENSMUSG00000025527.9 -5.25686222287999 1.46534006583066e-07        NA
## ENSMUSG00000104025.1 -4.53581314551642 5.73819675582155e-06        NA
## ENSMUSG00000104283.1 -4.53581314551642 5.73819675582155e-06        NA
##                                  symbol
##                             <character>
## ENSMUSG00000044022.5 ENSMUSG00000044022
## ENSMUSG00000045871.5 ENSMUSG00000045871
## ENSMUSG00000110839.1 ENSMUSG00000110839
## ENSMUSG00000025527.9 ENSMUSG00000025527
## ENSMUSG00000104025.1 ENSMUSG00000104025
## ENSMUSG00000104283.1 ENSMUSG00000104283
```


```r
resTemp = as.data.frame(res)
res.2 = merge(resTemp, annotLookup, by.x = "symbol", by.y= "ensembl_gene_id", all.x=TRUE)

write.csv(res.2, file = "results.csv")
head(res.2)
```

```
##               symbol   baseMean log2FoldChange     lfcSE       stat
## 1 ENSMUSG00000000001  22.128060      0.3044189 0.6588361  0.4468773
## 2 ENSMUSG00000000003   1.518339      2.0359587 3.4948357  0.5797007
## 3 ENSMUSG00000000028 263.555301     -1.0296405 0.3126001 -3.2618043
## 4 ENSMUSG00000000031   2.905246      0.7980923 3.3530677  0.2350362
## 5 ENSMUSG00000000037   2.983225      3.7704963 3.5584777  1.0567711
## 6 ENSMUSG00000000049   2.482652      0.5360606 2.2102653  0.2380079
##        pvalue       padj           original_id   gene_biotype
## 1 0.654963647 0.83874584  ENSMUSG00000000001.4 protein_coding
## 2 0.562116480         NA ENSMUSG00000000003.15 protein_coding
## 3 0.001107055 0.00613115 ENSMUSG00000000028.15 protein_coding
## 4 0.814180627         NA ENSMUSG00000000031.16         lncRNA
## 5 0.290616087         NA ENSMUSG00000000037.17 protein_coding
## 6 0.811874960         NA ENSMUSG00000000049.11 protein_coding
##   external_gene_name entrezgene_id
## 1              Gnai3         14679
## 2               Pbsn         54192
## 3              Cdc45         12544
## 4                H19            NA
## 5              Scml2        107815
## 6               Apoh         11818
```

Let's make a file with just the genes with an adjusted p-value < 0.5


```r
resSig = as.data.frame(subset(res.2,padj<0.5) )
resSig = resSig[order(resSig$log2FoldChange,decreasing=TRUE),]
head(resSig)
```

```
##                   symbol baseMean log2FoldChange    lfcSE     stat
## 34210 ENSMUSG00000118465 3.344910       17.47401 2.981323 5.857805
## 12264 ENSMUSG00000045053 3.614708       13.91114 3.876613 3.585898
## 3419  ENSMUSG00000022454 3.543100       13.42247 4.011089 3.343847
## 5308  ENSMUSG00000026705 5.590915       11.85766 3.354682 3.531681
## 14096 ENSMUSG00000053391 3.306057       10.84705 5.026153 2.156133
## 9290  ENSMUSG00000035401 5.954097       10.56069 3.251278 3.245091
##             pvalue         padj           original_id   gene_biotype
## 34210 4.690250e-09 8.378504e-08  ENSMUSG00000118465.1 protein_coding
## 12264 3.359204e-04 2.186188e-03  ENSMUSG00000045053.2 protein_coding
## 3419  8.262539e-04 4.741413e-03 ENSMUSG00000022454.17 protein_coding
## 5308  4.129276e-04 2.601280e-03 ENSMUSG00000026705.16 protein_coding
## 14096 3.107329e-02 1.001892e-01  ENSMUSG00000053391.5 protein_coding
## 9290  1.174131e-03 6.441636e-03  ENSMUSG00000035401.9 protein_coding
##       external_gene_name entrezgene_id
## 34210         CT867961.1     100040448
## 12264              Kcng3        225030
## 3419               Nell2         54003
## 5308              Klhl20        226541
## 14096            Olfr211        258914
## 9290                Emsy        233545
```

```r
print(paste("Number of significant genes:", nrow(resSig), sep=" "))
```

```
## [1] "Number of significant genes: 12145"
```

```r
write.csv(resSig,"SigGenes.csv")
```

Here is a volcano plot that shows the symbol that we created at each point.


```r
library(dplyr)
res1 = as.data.frame(subset(res.2, !is.na(res.2$padj)))
# add a new column using the mutate function in dplyr
res1 = mutate(res1, sig=ifelse(res1$padj<0.5, "FDR<0.05", "Not Sig"))
res1[which(abs(res1$log2FoldChange)<1),'sig'] <- "Not Sig"
p = ggplot(res1, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(col=sig)) +
  scale_color_manual(values=c("red", "black"))
p+geom_text(data=filter(res1, padj<0.05), aes(label=symbol))
```

![](Mouse_RNA_Seq_p53_genotoxic_files/figure-html/unnamed-chunk-67-1.png)<!-- -->


```r
library(dplyr)
# Install ggrepel package if needed
# install.packages("devtools")
# devtools::install_github("slowkow/ggrepel")
library(ggrepel) # "repels" overlapping text

p+geom_text_repel(data=filter(res1,  abs(log2FoldChange)>5 | padj < 1e-100 ), aes(label=symbol))
```

![](Mouse_RNA_Seq_p53_genotoxic_files/figure-html/unnamed-chunk-68-1.png)<!-- -->


```r
library("ReportingTools")
```

```
## 
```

```
## Registered S3 method overwritten by 'GGally':
##   method from   
##   +.gg   ggplot2
```

```
## 
```

```
## Registered S3 method overwritten by 'R.oo':
##   method        from       
##   throw.default R.methodsS3
```

```r
htmlRep <- HTMLReport(shortName="report", title="My report",
                      reportDirectory="./report")
publish(resSig, htmlRep)
url <- finish(htmlRep)
#browseURL(url)
```

















































































