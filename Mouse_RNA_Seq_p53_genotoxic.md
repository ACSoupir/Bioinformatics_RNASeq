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
qc = read.csv('Whole Data QC.csv', header=TRUE, na.strings="")
kable(qc) %>%
  kable_styling() %>%
  scroll_box(width = "1000px", height = "800px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:800px; overflow-x: scroll; width:1000px; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> X </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Sample </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121770_1 </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121770_2 </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121771_1 </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121771_2 </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121774_1 </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121774_2 </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121775_1 </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121775_2 </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121778_1 </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121778_2 </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121779_1 </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121779_2 </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121780_1 </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121780_2 </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121781_1 </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121781_2 </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121786_1 </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121786_2 </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121787_1 </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121787_2 </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121788_1 </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121788_2 </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121789_1 </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121789_2 </th>
  </tr>
 </thead>
<tbody>
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
</table></div>

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

When the count file is completed, we can import it into R and start working with it to determine differentially expressed genes. First we will import it into R


```r
library(limma)
library(DESeq2)
library(dplyr)
library(readr)

countData = read_csv("readCounts.csv", skip = 1)
```

This gives us our dataframe from out featureCounts program, but if we look at the data we see that featureCounts added some extra information that characterizes each gene_id.


```r
kable(head(countData)) %>%
  kable_styling() %>%
  scroll_box(width = "1000px", height = "320px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:320px; overflow-x: scroll; width:1000px; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Geneid </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Chr </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Start </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> End </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Strand </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Length </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121770Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121771Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121774Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121775Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121778Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121779Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121780Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121781Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121786Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121787Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121788Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121789Aligned.sortedByCoord.out.bam </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000102693.1 </td>
   <td style="text-align:left;"> chr1 </td>
   <td style="text-align:left;"> 3073253 </td>
   <td style="text-align:left;"> 3074322 </td>
   <td style="text-align:left;"> + </td>
   <td style="text-align:right;"> 1070 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000064842.1 </td>
   <td style="text-align:left;"> chr1 </td>
   <td style="text-align:left;"> 3102016 </td>
   <td style="text-align:left;"> 3102125 </td>
   <td style="text-align:left;"> + </td>
   <td style="text-align:right;"> 110 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000051951.5 </td>
   <td style="text-align:left;"> chr1;chr1;chr1;chr1;chr1;chr1;chr1 </td>
   <td style="text-align:left;"> 3205901;3206523;3213439;3213609;3214482;3421702;3670552 </td>
   <td style="text-align:left;"> 3207317;3207317;3215632;3216344;3216968;3421901;3671498 </td>
   <td style="text-align:left;"> -;-;-;-;-;-;- </td>
   <td style="text-align:right;"> 6094 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000102851.1 </td>
   <td style="text-align:left;"> chr1 </td>
   <td style="text-align:left;"> 3252757 </td>
   <td style="text-align:left;"> 3253236 </td>
   <td style="text-align:left;"> + </td>
   <td style="text-align:right;"> 480 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000103377.1 </td>
   <td style="text-align:left;"> chr1 </td>
   <td style="text-align:left;"> 3365731 </td>
   <td style="text-align:left;"> 3368549 </td>
   <td style="text-align:left;"> - </td>
   <td style="text-align:right;"> 2819 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000104017.1 </td>
   <td style="text-align:left;"> chr1 </td>
   <td style="text-align:left;"> 3375556 </td>
   <td style="text-align:left;"> 3377788 </td>
   <td style="text-align:left;"> - </td>
   <td style="text-align:right;"> 2233 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
</tbody>
</table></div>

We also need to set out row names to the gene_id. We will do some data frame manipulation and then look at the data again.


```r
countData = as.data.frame(countData)
rownames(countData) = countData$Geneid
countData = countData[,-c(1:6)]
kable(head(countData)) %>%
  kable_styling() %>%
  scroll_box(width = "1000px", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:1000px; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121770Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121771Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121774Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121775Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121778Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121779Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121780Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121781Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121786Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121787Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121788Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121789Aligned.sortedByCoord.out.bam </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000102693.1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000064842.1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000051951.5 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000102851.1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000103377.1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000104017.1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
</tbody>
</table></div>

### Quick Data Exploration


```r
dim(countData)
```

```
## [1] 55487    12
```


```r
kable(summary(countData)) %>%
  kable_styling() %>%
  scroll_box(width = "1000px", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:1000px; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121770Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121771Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121774Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121775Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121778Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121779Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121780Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121781Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121786Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121787Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121788Aligned.sortedByCoord.out.bam </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SRR2121789Aligned.sortedByCoord.out.bam </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> Min.   :      0.0 </td>
   <td style="text-align:left;"> Min.   :     0.0 </td>
   <td style="text-align:left;"> Min.   :     0.0 </td>
   <td style="text-align:left;"> Min.   :      0.0 </td>
   <td style="text-align:left;"> Min.   :      0.0 </td>
   <td style="text-align:left;"> Min.   :      0.0 </td>
   <td style="text-align:left;"> Min.   :     0.0 </td>
   <td style="text-align:left;"> Min.   :      0.0 </td>
   <td style="text-align:left;"> Min.   :     0.0 </td>
   <td style="text-align:left;"> Min.   :     0.0 </td>
   <td style="text-align:left;"> Min.   :     0.0 </td>
   <td style="text-align:left;"> Min.   :     0.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> 1st Qu.:      0.0 </td>
   <td style="text-align:left;"> 1st Qu.:     0.0 </td>
   <td style="text-align:left;"> 1st Qu.:     0.0 </td>
   <td style="text-align:left;"> 1st Qu.:      0.0 </td>
   <td style="text-align:left;"> 1st Qu.:      0.0 </td>
   <td style="text-align:left;"> 1st Qu.:      0.0 </td>
   <td style="text-align:left;"> 1st Qu.:     0.0 </td>
   <td style="text-align:left;"> 1st Qu.:      0.0 </td>
   <td style="text-align:left;"> 1st Qu.:     0.0 </td>
   <td style="text-align:left;"> 1st Qu.:     0.0 </td>
   <td style="text-align:left;"> 1st Qu.:     0.0 </td>
   <td style="text-align:left;"> 1st Qu.:     0.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> Median :      0.0 </td>
   <td style="text-align:left;"> Median :     0.0 </td>
   <td style="text-align:left;"> Median :     0.0 </td>
   <td style="text-align:left;"> Median :      0.0 </td>
   <td style="text-align:left;"> Median :      0.0 </td>
   <td style="text-align:left;"> Median :      0.0 </td>
   <td style="text-align:left;"> Median :     0.0 </td>
   <td style="text-align:left;"> Median :      0.0 </td>
   <td style="text-align:left;"> Median :     3.0 </td>
   <td style="text-align:left;"> Median :     4.0 </td>
   <td style="text-align:left;"> Median :     0.0 </td>
   <td style="text-align:left;"> Median :     1.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> Mean   :    273.8 </td>
   <td style="text-align:left;"> Mean   :   235.7 </td>
   <td style="text-align:left;"> Mean   :   229.6 </td>
   <td style="text-align:left;"> Mean   :    269.6 </td>
   <td style="text-align:left;"> Mean   :    286.6 </td>
   <td style="text-align:left;"> Mean   :    257.7 </td>
   <td style="text-align:left;"> Mean   :   193.7 </td>
   <td style="text-align:left;"> Mean   :    268.2 </td>
   <td style="text-align:left;"> Mean   :   160.2 </td>
   <td style="text-align:left;"> Mean   :   119.4 </td>
   <td style="text-align:left;"> Mean   :   184.3 </td>
   <td style="text-align:left;"> Mean   :   211.1 </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> 3rd Qu.:     25.0 </td>
   <td style="text-align:left;"> 3rd Qu.:    20.0 </td>
   <td style="text-align:left;"> 3rd Qu.:    15.0 </td>
   <td style="text-align:left;"> 3rd Qu.:     24.0 </td>
   <td style="text-align:left;"> 3rd Qu.:     25.0 </td>
   <td style="text-align:left;"> 3rd Qu.:     24.0 </td>
   <td style="text-align:left;"> 3rd Qu.:    17.0 </td>
   <td style="text-align:left;"> 3rd Qu.:     26.0 </td>
   <td style="text-align:left;"> 3rd Qu.:    21.0 </td>
   <td style="text-align:left;"> 3rd Qu.:    21.0 </td>
   <td style="text-align:left;"> 3rd Qu.:    16.0 </td>
   <td style="text-align:left;"> 3rd Qu.:    18.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> Max.   :1085454.0 </td>
   <td style="text-align:left;"> Max.   :834016.0 </td>
   <td style="text-align:left;"> Max.   :783768.0 </td>
   <td style="text-align:left;"> Max.   :1032589.0 </td>
   <td style="text-align:left;"> Max.   :1233905.0 </td>
   <td style="text-align:left;"> Max.   :1232421.0 </td>
   <td style="text-align:left;"> Max.   :805463.0 </td>
   <td style="text-align:left;"> Max.   :1210716.0 </td>
   <td style="text-align:left;"> Max.   :641911.0 </td>
   <td style="text-align:left;"> Max.   :516752.0 </td>
   <td style="text-align:left;"> Max.   :706363.0 </td>
   <td style="text-align:left;"> Max.   :756086.0 </td>
  </tr>
</tbody>
</table></div>

Let's also go ahead and change the names to describe out data a little better.


```r
columns = c('Trp53m_mock_1', 'Trp53m_mock_2', 'Trp53m_4h7Gy_1', 'Trp53m_4h7Gy_2', 'Trp53p_mock_1', 'Trp53p_mock_2', 'Trp53p_mock_3', 'Trp53p_mock_4', 'Trp53p_4h7Gy_1', 'Trp53p_4h7Gy_2', 'Trp53p_4h7Gy_3', 'Trp53p_4h7Gy_4')
colnames(countData) = columns
```

Here we have to make sure that we convert the *+/+* and *-/-* to characters. These characters


```r
par(mar=c(8,4,4,1)+0.1)
barplot( colSums(countData)/1e6, col="green",las=3,main="Total read counts (millions)", ylab="Total read counts in millions")
```

![](Mouse_RNA_Seq_p53_genotoxic_files/figure-html/unnamed-chunk-35-1.png)<!-- -->


```r
hist(countData[,1], br=200, xlab="Number of Reads Counts per Feature", main="Histogram of Read Counts for Trp53-/- Mock")
```

![](Mouse_RNA_Seq_p53_genotoxic_files/figure-html/unnamed-chunk-36-1.png)<!-- -->

We can see that our count data is highly skewed to the right. This is a great case for using **log** transformation!


```r
logCountData = log2(1+countData)
par(mfrow = c(1, 2), mar=c(8,4,4,1))  # two columns
hist(logCountData[,1], main="Histogram of Log Read Counts", xlab="Log transformed counts")
boxplot(logCountData,las=3, main="Boxplot of Log Read Counts")
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
plot(logCountData[,1],logCountData[,2], xlab="Trp53-/- mock replication 1", ylab="Trp53-/- mock replication 2")
```

![](Mouse_RNA_Seq_p53_genotoxic_files/figure-html/unnamed-chunk-39-1.png)<!-- -->

### Filtering, Normalization, and Trasformation using DESeq2

We have to make the experiment design into a small dataframe so we can tell DESeq how we want to analyze the data. Here will will make a small table that has the rep names that we changed the column names to previously, and then a column for which columns are Trp53+/+ or Trp53-/-, and which columns were control mice and which columns were treated with ionizing radiation.


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
p53 = c("m", "m", "m", "m",
        "p", "p", "p", "p", "p", "p", "p", "p")
treatment = c("control", "control", "IR", "IR", "control", "control", "control", "control",
              "IR", "IR", "IR", "IR")
```


```r
colData = cbind(colnames(countData), p53 )
colData
```

```
##                        p53
##  [1,] "Trp53m_mock_1"  "m"
##  [2,] "Trp53m_mock_2"  "m"
##  [3,] "Trp53m_4h7Gy_1" "m"
##  [4,] "Trp53m_4h7Gy_2" "m"
##  [5,] "Trp53p_mock_1"  "p"
##  [6,] "Trp53p_mock_2"  "p"
##  [7,] "Trp53p_mock_3"  "p"
##  [8,] "Trp53p_mock_4"  "p"
##  [9,] "Trp53p_4h7Gy_1" "p"
## [10,] "Trp53p_4h7Gy_2" "p"
## [11,] "Trp53p_4h7Gy_3" "p"
## [12,] "Trp53p_4h7Gy_4" "p"
```


```r
colData = as.data.frame(cbind(colnames(countData), p53, treatment))
colData
```

```
##                V1 p53 treatment
## 1   Trp53m_mock_1   m   control
## 2   Trp53m_mock_2   m   control
## 3  Trp53m_4h7Gy_1   m        IR
## 4  Trp53m_4h7Gy_2   m        IR
## 5   Trp53p_mock_1   p   control
## 6   Trp53p_mock_2   p   control
## 7   Trp53p_mock_3   p   control
## 8   Trp53p_mock_4   p   control
## 9  Trp53p_4h7Gy_1   p        IR
## 10 Trp53p_4h7Gy_2   p        IR
## 11 Trp53p_4h7Gy_3   p        IR
## 12 Trp53p_4h7Gy_4   p        IR
```

```r
str(colData)
```

```
## 'data.frame':	12 obs. of  3 variables:
##  $ V1       : Factor w/ 12 levels "Trp53m_4h7Gy_1",..: 3 4 1 2 9 10 11 12 5 6 ...
##  $ p53      : Factor w/ 2 levels "m","p": 1 1 1 1 2 2 2 2 2 2 ...
##  $ treatment: Factor w/ 2 levels "control","IR": 1 1 2 2 1 1 1 1 2 2 ...
```

Creating a DESeq Dataset


```r
dds = DESeqDataSetFromMatrix(countData=countData,
                             colData=colData,
                             design= ~ p53+treatment+p53*treatment)   # note that the study design is changed.
dds = DESeq(dds)  # main function
nrow(dds)
```

```
## [1] 55487
```

Filtering: we will only keep rows that have a sum count between all samples greater than 5. This will remove most of the genes that mostly have "0" counts.


```r
dds <- dds[ rowSums(counts(dds)) > 5, ]
nrow(dds)
```

```
## [1] 34139
```

Regularized log transformation - used for clustering


```r
rld <- rlog(dds, blind = FALSE)
kable(head(assay(rld), 6)) %>%
  kable_styling() %>%
  scroll_box(width = "1000px", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:1000px; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53m_mock_1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53m_mock_2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53m_4h7Gy_1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53m_4h7Gy_2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53p_mock_1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53p_mock_2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53p_mock_3 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53p_mock_4 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53p_4h7Gy_1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53p_4h7Gy_2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53p_4h7Gy_3 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53p_4h7Gy_4 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000102693.1 </td>
   <td style="text-align:right;"> -0.7873810 </td>
   <td style="text-align:right;"> -0.7869647 </td>
   <td style="text-align:right;"> -0.7864123 </td>
   <td style="text-align:right;"> -0.7872251 </td>
   <td style="text-align:right;"> -0.7873917 </td>
   <td style="text-align:right;"> -0.7872901 </td>
   <td style="text-align:right;"> -0.7863979 </td>
   <td style="text-align:right;"> -0.7874606 </td>
   <td style="text-align:right;"> -0.6313032 </td>
   <td style="text-align:right;"> -0.7027829 </td>
   <td style="text-align:right;"> -0.7854203 </td>
   <td style="text-align:right;"> -0.7864622 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000051951.5 </td>
   <td style="text-align:right;"> -0.1689760 </td>
   <td style="text-align:right;"> -0.1684855 </td>
   <td style="text-align:right;"> -0.1678308 </td>
   <td style="text-align:right;"> -0.0729011 </td>
   <td style="text-align:right;"> -0.1689887 </td>
   <td style="text-align:right;"> -0.1688691 </td>
   <td style="text-align:right;"> -0.1678136 </td>
   <td style="text-align:right;"> -0.1690696 </td>
   <td style="text-align:right;"> -0.0298983 </td>
   <td style="text-align:right;"> -0.0789877 </td>
   <td style="text-align:right;"> -0.1675075 </td>
   <td style="text-align:right;"> -0.1678901 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000103377.1 </td>
   <td style="text-align:right;"> 0.9367497 </td>
   <td style="text-align:right;"> 0.9385253 </td>
   <td style="text-align:right;"> 0.9408618 </td>
   <td style="text-align:right;"> 1.0007595 </td>
   <td style="text-align:right;"> 0.9367037 </td>
   <td style="text-align:right;"> 0.9371385 </td>
   <td style="text-align:right;"> 0.9409228 </td>
   <td style="text-align:right;"> 0.9364087 </td>
   <td style="text-align:right;"> 1.1895736 </td>
   <td style="text-align:right;"> 1.5284052 </td>
   <td style="text-align:right;"> 0.9420021 </td>
   <td style="text-align:right;"> 0.9406517 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000103201.1 </td>
   <td style="text-align:right;"> -0.0694969 </td>
   <td style="text-align:right;"> -0.0689413 </td>
   <td style="text-align:right;"> -0.0682009 </td>
   <td style="text-align:right;"> -0.0692892 </td>
   <td style="text-align:right;"> -0.0695112 </td>
   <td style="text-align:right;"> -0.0693758 </td>
   <td style="text-align:right;"> -0.0681815 </td>
   <td style="text-align:right;"> -0.0696029 </td>
   <td style="text-align:right;"> 0.0406074 </td>
   <td style="text-align:right;"> 0.1775249 </td>
   <td style="text-align:right;"> -0.0678357 </td>
   <td style="text-align:right;"> -0.0682680 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000102592.1 </td>
   <td style="text-align:right;"> 0.6079270 </td>
   <td style="text-align:right;"> 0.6090832 </td>
   <td style="text-align:right;"> 0.6106138 </td>
   <td style="text-align:right;"> 0.6083601 </td>
   <td style="text-align:right;"> 0.6078971 </td>
   <td style="text-align:right;"> 0.6081797 </td>
   <td style="text-align:right;"> 0.6106538 </td>
   <td style="text-align:right;"> 0.6077056 </td>
   <td style="text-align:right;"> 0.9928379 </td>
   <td style="text-align:right;"> 0.8039315 </td>
   <td style="text-align:right;"> 0.6113644 </td>
   <td style="text-align:right;"> 0.6829781 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000025900.12 </td>
   <td style="text-align:right;"> 1.7329730 </td>
   <td style="text-align:right;"> 1.7166717 </td>
   <td style="text-align:right;"> 1.6536060 </td>
   <td style="text-align:right;"> 1.6483113 </td>
   <td style="text-align:right;"> 1.6766519 </td>
   <td style="text-align:right;"> 1.6478854 </td>
   <td style="text-align:right;"> 1.6536997 </td>
   <td style="text-align:right;"> 1.7812828 </td>
   <td style="text-align:right;"> 2.0563917 </td>
   <td style="text-align:right;"> 2.2880207 </td>
   <td style="text-align:right;"> 1.6553593 </td>
   <td style="text-align:right;"> 1.7301108 </td>
  </tr>
</tbody>
</table></div>

Variance Stabilizing Transformation

# Interactions cause a difference between the lfc betwen pooled data, e.g. p53+/+ (control and IR) and p53-/- (control and IR)


```r
vsd <- vst(dds, blind = FALSE)
kable(head(assay(vsd), 6)) %>%
  kable_styling() %>%
  scroll_box(width = "1000px", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:1000px; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53m_mock_1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53m_mock_2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53m_4h7Gy_1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53m_4h7Gy_2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53p_mock_1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53p_mock_2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53p_mock_3 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53p_mock_4 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53p_4h7Gy_1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53p_4h7Gy_2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53p_4h7Gy_3 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53p_4h7Gy_4 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000102693.1 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.784352 </td>
   <td style="text-align:right;"> 6.714627 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000051951.5 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.715020 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.784352 </td>
   <td style="text-align:right;"> 6.714627 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000103377.1 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.632583 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.862784 </td>
   <td style="text-align:right;"> 7.171182 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000103201.1 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.737470 </td>
   <td style="text-align:right;"> 6.919246 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000102592.1 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 7.013160 </td>
   <td style="text-align:right;"> 6.830644 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.662579 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000025900.12 </td>
   <td style="text-align:right;"> 6.669513 </td>
   <td style="text-align:right;"> 6.643045 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.569347 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.732883 </td>
   <td style="text-align:right;"> 7.013160 </td>
   <td style="text-align:right;"> 7.220955 </td>
   <td style="text-align:right;"> 6.433022 </td>
   <td style="text-align:right;"> 6.662579 </td>
  </tr>
</tbody>
</table></div>

For the log2 approach, we need to first estimate size factors to account for sequencing depth, and then specify normalized=TRUE. Sequencing depth correction is done automatically for the rlog and the vst.

Size Factor


```r
dds <- estimateSizeFactors(dds)
kable(sizeFactors(dds)) %>%
  kable_styling() %>%
  scroll_box(width = "300px", height = "520px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:520px; overflow-x: scroll; width:300px; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Trp53m_mock_1 </td>
   <td style="text-align:right;"> 1.2892566 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Trp53m_mock_2 </td>
   <td style="text-align:right;"> 1.0903014 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Trp53m_4h7Gy_1 </td>
   <td style="text-align:right;"> 0.8973654 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Trp53m_4h7Gy_2 </td>
   <td style="text-align:right;"> 1.2078255 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Trp53p_mock_1 </td>
   <td style="text-align:right;"> 1.2952327 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Trp53p_mock_2 </td>
   <td style="text-align:right;"> 1.2406213 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Trp53p_mock_3 </td>
   <td style="text-align:right;"> 0.8931001 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Trp53p_mock_4 </td>
   <td style="text-align:right;"> 1.3347016 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Trp53p_4h7Gy_1 </td>
   <td style="text-align:right;"> 0.7767888 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Trp53p_4h7Gy_2 </td>
   <td style="text-align:right;"> 0.6056051 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Trp53p_4h7Gy_3 </td>
   <td style="text-align:right;"> 0.8227361 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Trp53p_4h7Gy_4 </td>
   <td style="text-align:right;"> 0.9123249 </td>
  </tr>
</tbody>
</table></div>

We will first look at the log transformed data


```r
slog <- log2(counts(dds, normalized=TRUE)+1)
kable(head(slog)) %>%
  kable_styling() %>%
  scroll_box(width = "1000px", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:1000px; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53m_mock_1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53m_mock_2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53m_4h7Gy_1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53m_4h7Gy_2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53p_mock_1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53p_mock_2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53p_mock_3 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53p_mock_4 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53p_4h7Gy_1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53p_4h7Gy_2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53p_4h7Gy_3 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Trp53p_4h7Gy_4 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000102693.1 </td>
   <td style="text-align:right;"> 0.000000 </td>
   <td style="text-align:right;"> 0.000000 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.000000 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.000000 </td>
   <td style="text-align:right;"> 2.620447 </td>
   <td style="text-align:right;"> 2.105169 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000051951.5 </td>
   <td style="text-align:right;"> 0.000000 </td>
   <td style="text-align:right;"> 0.000000 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 2.108269 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.000000 </td>
   <td style="text-align:right;"> 2.620447 </td>
   <td style="text-align:right;"> 2.105169 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000103377.1 </td>
   <td style="text-align:right;"> 0.000000 </td>
   <td style="text-align:right;"> 0.000000 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1.409184 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.000000 </td>
   <td style="text-align:right;"> 3.125008 </td>
   <td style="text-align:right;"> 4.592001 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000103201.1 </td>
   <td style="text-align:right;"> 0.000000 </td>
   <td style="text-align:right;"> 0.000000 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.000000 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.000000 </td>
   <td style="text-align:right;"> 2.281566 </td>
   <td style="text-align:right;"> 3.447241 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000102592.1 </td>
   <td style="text-align:right;"> 0.000000 </td>
   <td style="text-align:right;"> 0.000000 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.000000 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.000000 </td>
   <td style="text-align:right;"> 3.922280 </td>
   <td style="text-align:right;"> 2.926941 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1.674552 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000025900.12 </td>
   <td style="text-align:right;"> 1.734188 </td>
   <td style="text-align:right;"> 1.503021 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.000000 </td>
   <td style="text-align:right;"> 0.8254291 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 2.246759 </td>
   <td style="text-align:right;"> 3.922280 </td>
   <td style="text-align:right;"> 4.777149 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1.674552 </td>
  </tr>
</tbody>
</table></div>


```r
par(mfrow = c(1, 3))  # 3 columns
plot(slog[,1],slog[,2])
plot(assay(rld)[,1],assay(rld)[,2])
plot(assay(vsd)[,1],assay(vsd)[,2])
```

![](Mouse_RNA_Seq_p53_genotoxic_files/figure-html/unnamed-chunk-49-1.png)<!-- -->

As the log transformation constant increases, the information of the data is lost.


```r
par(mfrow = c(1, 3))  # 3 columns
slog <- log2(counts(dds, normalized=TRUE)+1)
plot(slog[,1],slog[,2])
slog <- log2(counts(dds, normalized=TRUE)+4)
plot(slog[,1],slog[,2], xlim=c(0,20))
slog <- log2(counts(dds, normalized=TRUE)+20)
plot(slog[,1],slog[,2], xlim=c(0,20))
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
plotPCA(rld, intgroup = c("p53", "treatment")) + theme(aspect.ratio=1)
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

n=100 # number of top genes by standard deviation

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
	,margins=c(8,12)
	,cexRow=1
	,srtCol=45
	,cexCol=1.  # size of font for sample names
	,lmat = lmat, lwid = lwid, lhei = lhei
	)
```

![](Mouse_RNA_Seq_p53_genotoxic_files/figure-html/heatmap-1.png)<!-- -->

### Identification of Differentially Expressed Genes


```r
dds <- DESeq(dds)
res <- results(dds)

kable(head(res)) %>%
  kable_styling() %>%
  scroll_box(width = "1000px", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:1000px; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> baseMean </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> log2FoldChange </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> lfcSE </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> stat </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> pvalue </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> padj </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000102693.1 </td>
   <td style="text-align:right;"> 0.7043239 </td>
   <td style="text-align:right;"> 3.5623773 </td>
   <td style="text-align:right;"> 6.556337 </td>
   <td style="text-align:right;"> 0.5433487 </td>
   <td style="text-align:right;"> 0.5868897 </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000051951.5 </td>
   <td style="text-align:right;"> 0.9803019 </td>
   <td style="text-align:right;"> 0.3244198 </td>
   <td style="text-align:right;"> 6.442842 </td>
   <td style="text-align:right;"> 0.0503535 </td>
   <td style="text-align:right;"> 0.9598407 </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000103377.1 </td>
   <td style="text-align:right;"> 2.7081125 </td>
   <td style="text-align:right;"> 3.1757639 </td>
   <td style="text-align:right;"> 5.419836 </td>
   <td style="text-align:right;"> 0.5859520 </td>
   <td style="text-align:right;"> 0.5579078 </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000103201.1 </td>
   <td style="text-align:right;"> 1.1474583 </td>
   <td style="text-align:right;"> 4.2647949 </td>
   <td style="text-align:right;"> 6.546806 </td>
   <td style="text-align:right;"> 0.6514314 </td>
   <td style="text-align:right;"> 0.5147681 </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000102592.1 </td>
   <td style="text-align:right;"> 1.9131690 </td>
   <td style="text-align:right;"> 5.0084508 </td>
   <td style="text-align:right;"> 5.437877 </td>
   <td style="text-align:right;"> 0.9210305 </td>
   <td style="text-align:right;"> 0.3570345 </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000025900.12 </td>
   <td style="text-align:right;"> 4.2877013 </td>
   <td style="text-align:right;"> 6.7470978 </td>
   <td style="text-align:right;"> 3.112343 </td>
   <td style="text-align:right;"> 2.1678514 </td>
   <td style="text-align:right;"> 0.0301700 </td>
   <td style="text-align:right;"> 0.1442115 </td>
  </tr>
</tbody>
</table></div>

DESeq2 uses the Benjamini-Hochberg (BH) adjustment (Benjamini and Hochberg 1995) as implemented in the base R p.adjust function


```r
res <- results(dds, alpha = 0.5, lfcThreshold=0.01)
summary(res)
```

```
## 
## out of 34139 with nonzero total read count
## adjusted p-value < 0.5
## LFC > 0.01 (up)    : 5794, 17%
## LFC < -0.01 (down) : 4094, 12%
## outliers [1]       : 3, 0.0088%
## low counts [2]     : 8605, 25%
## (mean count < 2)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

Now lets sort genes by fold change


```r
res <- res[order(abs( res$log2FoldChange), decreasing=TRUE),]
kable(head(res)) %>%
  kable_styling() %>%
  scroll_box(width = "1000px", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:1000px; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> baseMean </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> log2FoldChange </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> lfcSE </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> stat </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> pvalue </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> padj </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000107277.1 </td>
   <td style="text-align:right;"> 178.134530 </td>
   <td style="text-align:right;"> 12.82857 </td>
   <td style="text-align:right;"> 1.906777 </td>
   <td style="text-align:right;"> 6.722639 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000036281.13 </td>
   <td style="text-align:right;"> 17.898231 </td>
   <td style="text-align:right;"> 11.90953 </td>
   <td style="text-align:right;"> 1.997989 </td>
   <td style="text-align:right;"> 5.955753 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.0000002 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000034189.5 </td>
   <td style="text-align:right;"> 21.363195 </td>
   <td style="text-align:right;"> 11.19200 </td>
   <td style="text-align:right;"> 2.064077 </td>
   <td style="text-align:right;"> 5.417433 </td>
   <td style="text-align:right;"> 0.0000001 </td>
   <td style="text-align:right;"> 0.0000035 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000111241.1 </td>
   <td style="text-align:right;"> 26.215073 </td>
   <td style="text-align:right;"> 10.91199 </td>
   <td style="text-align:right;"> 1.885793 </td>
   <td style="text-align:right;"> 5.781118 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.0000005 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000027479.14 </td>
   <td style="text-align:right;"> 4.844616 </td>
   <td style="text-align:right;"> 10.82957 </td>
   <td style="text-align:right;"> 3.264449 </td>
   <td style="text-align:right;"> 3.314364 </td>
   <td style="text-align:right;"> 0.0009185 </td>
   <td style="text-align:right;"> 0.0130936 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000014932.15 </td>
   <td style="text-align:right;"> 77.883535 </td>
   <td style="text-align:right;"> 10.65416 </td>
   <td style="text-align:right;"> 1.686891 </td>
   <td style="text-align:right;"> 6.309929 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.0000000 </td>
  </tr>
</tbody>
</table></div>

MA Plot

Show the significant genes. The lower the average read counts for all samples and the higher the variation between the samples, the less significant those genes are. 


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
## Warning: Removed 8608 rows containing missing values (geom_point).
```

![](Mouse_RNA_Seq_p53_genotoxic_files/figure-html/unnamed-chunk-59-1.png)<!-- -->

### Gene Annotations

Plot counts of top gene


```r
topGene <- rownames(res)[1]
plotCounts(dds, gene = topGene, intgroup=c("p53", "treatment"))
```

![](Mouse_RNA_Seq_p53_genotoxic_files/figure-html/unnamed-chunk-60-1.png)<!-- -->

Here we see an interesting point or our normalized counts under the Trp53p_4h7Gy group that seems to be extremely high, while the other 3 replicates are around 0.5. I cannot get this portion to work, however. I get an error stating that None of the keys entered are valid keys for 'SYMBOL'. 

Let's look at the keys that we have to work with for our ***res*** file from DESeq2


```r
head(row.names(res))
```

```
## [1] "ENSMUSG00000107277.1"  "ENSMUSG00000036281.13" "ENSMUSG00000034189.5" 
## [4] "ENSMUSG00000111241.1"  "ENSMUSG00000027479.14" "ENSMUSG00000014932.15"
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
#res$symbol <- gsub(" ","",row.names(res)) 
```


```r
message("Ensembl IDs")
```

```
## Ensembl IDs
```

```r
key.en = keys(org.Mm.eg.db, keytype="ENSEMBL")
head(key.en)
```

```
## [1] "ENSMUSG00000030359" "ENSMUSG00000020804" "ENSMUSG00000025375"
## [4] "ENSMUSG00000015243" "ENSMUSG00000028125" "ENSMUSG00000026944"
```

```r
cat("\n\n")
```

```r
message("SYMBOL names")
```

```
## SYMBOL names
```

```r
key.sy = keys(org.Mm.eg.db, keytype="SYMBOL")
head(key.sy)
```

```
## [1] "Pzp"   "Aanat" "Aatk"  "Abca1" "Abca4" "Abca2"
```

These are ENSEMBL symbols, so we need to designate that when looking for the genes that we have. 


```r
res$ensembl <- gsub("\\..*","", row.names(res))
res$entrez <- mapIds(org.Mm.eg.db,
                     keys= res$ensembl,
                     column="ENTREZID",
                     keytype="ENSEMBL", #Out ID is ENSMBL
                     multiVals="first")
```

```
## 'select()' returned 1:many mapping between keys and columns
```

```r
res$symbol <- mapIds(org.Mm.eg.db,
                     keys= res$ensembl,
                     column="SYMBOL",
                     keytype="ENSEMBL", #Out ID is ENSMBL
                     multiVals="first")
```

```
## 'select()' returned 1:many mapping between keys and columns
```

```r
write.csv(res, file = "results.csv")
```


```r
kable(head(res)) %>%
  kable_styling() %>%
  scroll_box(width = "1000px", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:1000px; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> baseMean </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> log2FoldChange </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> lfcSE </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> stat </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> pvalue </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> padj </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> symbol </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> ensembl </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> entrez </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000107277.1 </td>
   <td style="text-align:right;"> 178.134530 </td>
   <td style="text-align:right;"> 12.82857 </td>
   <td style="text-align:right;"> 1.906777 </td>
   <td style="text-align:right;"> 6.722639 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ENSMUSG00000107277 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000036281.13 </td>
   <td style="text-align:right;"> 17.898231 </td>
   <td style="text-align:right;"> 11.90953 </td>
   <td style="text-align:right;"> 1.997989 </td>
   <td style="text-align:right;"> 5.955753 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.0000002 </td>
   <td style="text-align:left;"> Snapc4 </td>
   <td style="text-align:left;"> ENSMUSG00000036281 </td>
   <td style="text-align:left;"> 227644 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000034189.5 </td>
   <td style="text-align:right;"> 21.363195 </td>
   <td style="text-align:right;"> 11.19200 </td>
   <td style="text-align:right;"> 2.064077 </td>
   <td style="text-align:right;"> 5.417433 </td>
   <td style="text-align:right;"> 0.0000001 </td>
   <td style="text-align:right;"> 0.0000035 </td>
   <td style="text-align:left;"> Hsdl1 </td>
   <td style="text-align:left;"> ENSMUSG00000034189 </td>
   <td style="text-align:left;"> 72552 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000111241.1 </td>
   <td style="text-align:right;"> 26.215073 </td>
   <td style="text-align:right;"> 10.91199 </td>
   <td style="text-align:right;"> 1.885793 </td>
   <td style="text-align:right;"> 5.781118 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.0000005 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> ENSMUSG00000111241 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000027479.14 </td>
   <td style="text-align:right;"> 4.844616 </td>
   <td style="text-align:right;"> 10.82957 </td>
   <td style="text-align:right;"> 3.264449 </td>
   <td style="text-align:right;"> 3.314364 </td>
   <td style="text-align:right;"> 0.0009185 </td>
   <td style="text-align:right;"> 0.0130936 </td>
   <td style="text-align:left;"> Mapre1 </td>
   <td style="text-align:left;"> ENSMUSG00000027479 </td>
   <td style="text-align:left;"> 13589 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENSMUSG00000014932.15 </td>
   <td style="text-align:right;"> 77.883535 </td>
   <td style="text-align:right;"> 10.65416 </td>
   <td style="text-align:right;"> 1.686891 </td>
   <td style="text-align:right;"> 6.309929 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:left;"> Yes1 </td>
   <td style="text-align:left;"> ENSMUSG00000014932 </td>
   <td style="text-align:left;"> 22612 </td>
  </tr>
</tbody>
</table></div>

Let's make a file with just the genes with an adjusted p-value < 0.5


```r
resSig = as.data.frame(subset(res,padj<0.5) )
resSig = resSig[order(resSig$log2FoldChange,decreasing=TRUE),]
head(resSig)
```

```
##                         baseMean log2FoldChange    lfcSE     stat
## ENSMUSG00000107277.1  178.134530       12.82857 1.906777 6.722639
## ENSMUSG00000036281.13  17.898231       11.90953 1.997989 5.955753
## ENSMUSG00000034189.5   21.363195       11.19200 2.064077 5.417433
## ENSMUSG00000111241.1   26.215073       10.91199 1.885793 5.781118
## ENSMUSG00000027479.14   4.844616       10.82957 3.264449 3.314364
## ENSMUSG00000014932.15  77.883535       10.65416 1.686891 6.309929
##                             pvalue         padj symbol            ensembl
## ENSMUSG00000107277.1  1.784620e-11 2.190535e-09   <NA> ENSMUSG00000107277
## ENSMUSG00000036281.13 2.588774e-09 2.052608e-07 Snapc4 ENSMUSG00000036281
## ENSMUSG00000034189.5  6.046073e-08 3.540420e-06  Hsdl1 ENSMUSG00000034189
## ENSMUSG00000111241.1  7.420598e-09 5.248069e-07   <NA> ENSMUSG00000111241
## ENSMUSG00000027479.14 9.185180e-04 1.309363e-02 Mapre1 ENSMUSG00000027479
## ENSMUSG00000014932.15 2.791627e-10 2.795021e-08   Yes1 ENSMUSG00000014932
##                       entrez
## ENSMUSG00000107277.1    <NA>
## ENSMUSG00000036281.13 227644
## ENSMUSG00000034189.5   72552
## ENSMUSG00000111241.1    <NA>
## ENSMUSG00000027479.14  13589
## ENSMUSG00000014932.15  22612
```

```r
write.csv(resSig,"SigGenes.csv")
```

Here is a volcano plot that shows the symbol that we created at each point.


```r
library(dplyr)
res1 = as.data.frame(res)
# add a new column using the mutate function in dplyr
res1 = mutate(res1, sig=ifelse(res1$padj<0.5, "FDR<0.05", "Not Sig"))
res1[which(abs(res1$log2FoldChange)<1),'sig'] <- "Not Sig"
p = ggplot(res1, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(col=sig)) +
  scale_color_manual(values=c("red", "black"))
p+geom_text(data=filter(res1, padj<1e-50), aes(label=symbol))
```

```
## Warning: Removed 7357 rows containing missing values (geom_point).
```

![](Mouse_RNA_Seq_p53_genotoxic_files/figure-html/unnamed-chunk-67-1.png)<!-- -->





# 7. GO Enrichment analysis using GOstats

Here we will do a GO Enrichment analysis for genes that have a decreased fold-change of 5 or more


```r
library(GO.db)
library(GOstats)
selectedGenes = unique(resSig[resSig$log2FoldChange>5,'entrez'])  # upregulated genes
universeGenes =  unique( mapIds(org.Mm.eg.db,
                     keys= res$ensembl,
                     column="ENTREZID",
                     keytype="ENSEMBL", #Out ID is ENSMBL
                     multiVals="first")
                    )


 hgCutoff <- 0.001
 params <- new("GOHyperGParams",
     geneIds=selectedGenes,
     universeGeneIds=universeGenes,
     annotation="org.Mm.eg.db",
     ontology="BP",
     pvalueCutoff=hgCutoff,
     conditional=FALSE,
     testDirection="over")

hgOver <- hyperGTest(params)
summary(hgOver)[1:10,]
```

```
##        GOBPID       Pvalue OddsRatio   ExpCount Count Size
## 1  GO:0019236 2.643566e-07  5.062236   4.423515    18   78
## 2  GO:0016049 1.730359e-05  2.028026  26.144107    49  461
## 3  GO:0048588 2.580339e-05  2.442482  13.951085    31  246
## 4  GO:1990138 3.035698e-05  2.725414  10.208111    25  180
## 5  GO:0050789 3.546769e-05  1.302786 562.863914   624 9925
## 6  GO:0001558 5.334949e-05  2.045100  22.174286    42  391
## 7  GO:0032879 5.962358e-05  1.395721 147.223650   191 2596
## 8  GO:0060560 1.073049e-04  2.275476  14.348068    30  253
## 9  GO:0061564 1.204488e-04  1.884894  26.711225    47  471
## 10 GO:0050794 1.250445e-04  1.272536 527.135525   584 9295
##                                              Term
## 1                           response to pheromone
## 2                                     cell growth
## 3                       developmental cell growth
## 4                     neuron projection extension
## 5                regulation of biological process
## 6                       regulation of cell growth
## 7                      regulation of localization
## 8  developmental growth involved in morphogenesis
## 9                                axon development
## 10                 regulation of cellular process
```


```r
summary(hgOver)[1:10,c("GOBPID","Pvalue","Term")]
```

```
##        GOBPID       Pvalue                                           Term
## 1  GO:0019236 2.643566e-07                          response to pheromone
## 2  GO:0016049 1.730359e-05                                    cell growth
## 3  GO:0048588 2.580339e-05                      developmental cell growth
## 4  GO:1990138 3.035698e-05                    neuron projection extension
## 5  GO:0050789 3.546769e-05               regulation of biological process
## 6  GO:0001558 5.334949e-05                      regulation of cell growth
## 7  GO:0032879 5.962358e-05                     regulation of localization
## 8  GO:0060560 1.073049e-04 developmental growth involved in morphogenesis
## 9  GO:0061564 1.204488e-04                               axon development
## 10 GO:0050794 1.250445e-04                 regulation of cellular process
```



```r
params1 <- params
ontology(params1) <- "CC"
hgOver <- hyperGTest(params1)
summary(hgOver)[1:10,c("GOCCID","Pvalue","Term")]
```

```
##        GOCCID       Pvalue                                        Term
## 1  GO:0031226 2.194852e-10      intrinsic component of plasma membrane
## 2  GO:0005887 2.497976e-10       integral component of plasma membrane
## 3  GO:0030424 3.300807e-07                                        axon
## 4  GO:0098978 7.089583e-07                       glutamatergic synapse
## 5  GO:0033267 1.725126e-06                                   axon part
## 6  GO:0044459 1.872288e-06                        plasma membrane part
## 7  GO:0099055 3.042517e-06 integral component of postsynaptic membrane
## 8  GO:0099699 4.152767e-06     integral component of synaptic membrane
## 9  GO:0045211 6.111879e-06                       postsynaptic membrane
## 10 GO:0016021 7.169574e-06              integral component of membrane
```


```r
params1 <- params
ontology(params1) <- "MF"
hgOver <- hyperGTest(params1)
summary(hgOver)[1:10,c("GOMFID","Pvalue","Term")]
```

```
##        GOMFID       Pvalue                                         Term
## 1  GO:0016503 4.685420e-07                  pheromone receptor activity
## 2  GO:0005550 1.234794e-05                            pheromone binding
## 3  GO:0005549 2.293309e-05                              odorant binding
## 4  GO:0004888 3.280256e-05    transmembrane signaling receptor activity
## 5  GO:0038023 5.222020e-05                  signaling receptor activity
## 6  GO:0060089 8.100909e-05                molecular transducer activity
## 7  GO:0005488 6.931421e-04                                      binding
## 8  GO:0015267 8.237763e-04                             channel activity
## 9  GO:0022803 8.237763e-04   passive transmembrane transporter activity
## 10 GO:0046873 8.631819e-04 metal ion transmembrane transporter activity
```

## GO Enrichment analysis  of downregulated genes

Next we will have a look at the genes that are upregulated by a fold-change of 5 or greater. 


```r
selectedGenes = unique(resSig[resSig$log2FoldChange<5,'entrez'])  # upregulated genes

 params <- new("GOHyperGParams",
     geneIds=selectedGenes,
     universeGeneIds=universeGenes,
     annotation="org.Mm.eg.db",
     ontology="BP",
     pvalueCutoff=hgCutoff,
     conditional=FALSE,
     testDirection="over")

hgOver <- hyperGTest(params)
summary(hgOver)[1:10,c("GOBPID","Pvalue","Term")]
```

```
##        GOBPID       Pvalue                                          Term
## 1  GO:0016043 9.317561e-12               cellular component organization
## 2  GO:0008152 1.757957e-11                             metabolic process
## 3  GO:0071840 4.299006e-11 cellular component organization or biogenesis
## 4  GO:0044237 7.869800e-11                    cellular metabolic process
## 5  GO:0044238 1.013489e-09                     primary metabolic process
## 6  GO:0006807 1.956109e-09           nitrogen compound metabolic process
## 7  GO:0071704 2.363963e-09           organic substance metabolic process
## 8  GO:0006996 2.779804e-09                        organelle organization
## 9  GO:0033036 1.835538e-08                    macromolecule localization
## 10 GO:0051179 2.179269e-08                                  localization
```

# 8. Pathway analysis using expression data


```r
# bioconductor packages
# source("https://bioconductor.org/biocLite.R");
# biocLite(c("pathview","gage","gageData"))
library(pathview) 
library(gage) 
```

## Prepare data


```r
foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
head(foldchanges)
```

```
##     <NA>   227644    72552     <NA>    13589    22612 
## 12.82857 11.90953 11.19200 10.91199 10.82957 10.65416
```


```r
library(gageData)
data(go.sets.mm)
data(go.subs.mm)
gobpsets = go.sets.mm[go.subs.mm$BP]
gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)
#lapply(gobpres, head)
message("Greater")
```

```
## Greater
```

```r
kable(head(gobpres$greater)) %>%
  kable_styling() %>%
  scroll_box(width = "1000px", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:1000px; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.geomean </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> stat.mean </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.val </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> q.val </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> set.size </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> exp1 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> GO:0019236 response to pheromone </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 7.637046 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 76 </td>
   <td style="text-align:right;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0007420 brain development </td>
   <td style="text-align:right;"> 0.0012778 </td>
   <td style="text-align:right;"> 3.026272 </td>
   <td style="text-align:right;"> 0.0012778 </td>
   <td style="text-align:right;"> 0.7925518 </td>
   <td style="text-align:right;"> 400 </td>
   <td style="text-align:right;"> 0.0012778 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0045664 regulation of neuron differentiation </td>
   <td style="text-align:right;"> 0.0021154 </td>
   <td style="text-align:right;"> 2.869812 </td>
   <td style="text-align:right;"> 0.0021154 </td>
   <td style="text-align:right;"> 0.7925518 </td>
   <td style="text-align:right;"> 355 </td>
   <td style="text-align:right;"> 0.0021154 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0007156 homophilic cell adhesion </td>
   <td style="text-align:right;"> 0.0023079 </td>
   <td style="text-align:right;"> 2.879945 </td>
   <td style="text-align:right;"> 0.0023079 </td>
   <td style="text-align:right;"> 0.7925518 </td>
   <td style="text-align:right;"> 70 </td>
   <td style="text-align:right;"> 0.0023079 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0001944 vasculature development </td>
   <td style="text-align:right;"> 0.0030643 </td>
   <td style="text-align:right;"> 2.746850 </td>
   <td style="text-align:right;"> 0.0030643 </td>
   <td style="text-align:right;"> 0.7925518 </td>
   <td style="text-align:right;"> 487 </td>
   <td style="text-align:right;"> 0.0030643 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0033555 multicellular organismal response to stress </td>
   <td style="text-align:right;"> 0.0038427 </td>
   <td style="text-align:right;"> 2.719970 </td>
   <td style="text-align:right;"> 0.0038427 </td>
   <td style="text-align:right;"> 0.7925518 </td>
   <td style="text-align:right;"> 53 </td>
   <td style="text-align:right;"> 0.0038427 </td>
  </tr>
</tbody>
</table></div>

```r
message("Less")
```

```
## Less
```

```r
kable(head(gobpres$less)) %>%
  kable_styling() %>%
  scroll_box(width = "1000px", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:1000px; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.geomean </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> stat.mean </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.val </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> q.val </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> set.size </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> exp1 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> GO:0007283 spermatogenesis </td>
   <td style="text-align:right;"> 0.0094832 </td>
   <td style="text-align:right;"> -2.352919 </td>
   <td style="text-align:right;"> 0.0094832 </td>
   <td style="text-align:right;"> 0.992104 </td>
   <td style="text-align:right;"> 292 </td>
   <td style="text-align:right;"> 0.0094832 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0048232 male gamete generation </td>
   <td style="text-align:right;"> 0.0100757 </td>
   <td style="text-align:right;"> -2.330049 </td>
   <td style="text-align:right;"> 0.0100757 </td>
   <td style="text-align:right;"> 0.992104 </td>
   <td style="text-align:right;"> 293 </td>
   <td style="text-align:right;"> 0.0100757 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0007266 Rho protein signal transduction </td>
   <td style="text-align:right;"> 0.0180131 </td>
   <td style="text-align:right;"> -2.109476 </td>
   <td style="text-align:right;"> 0.0180131 </td>
   <td style="text-align:right;"> 0.992104 </td>
   <td style="text-align:right;"> 117 </td>
   <td style="text-align:right;"> 0.0180131 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0007094 mitotic cell cycle spindle assembly checkpoint </td>
   <td style="text-align:right;"> 0.0263331 </td>
   <td style="text-align:right;"> -2.065183 </td>
   <td style="text-align:right;"> 0.0263331 </td>
   <td style="text-align:right;"> 0.992104 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.0263331 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0016441 posttranscriptional gene silencing </td>
   <td style="text-align:right;"> 0.0274536 </td>
   <td style="text-align:right;"> -1.957280 </td>
   <td style="text-align:right;"> 0.0274536 </td>
   <td style="text-align:right;"> 0.992104 </td>
   <td style="text-align:right;"> 35 </td>
   <td style="text-align:right;"> 0.0274536 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0007276 gamete generation </td>
   <td style="text-align:right;"> 0.0283921 </td>
   <td style="text-align:right;"> -1.907876 </td>
   <td style="text-align:right;"> 0.0283921 </td>
   <td style="text-align:right;"> 0.992104 </td>
   <td style="text-align:right;"> 388 </td>
   <td style="text-align:right;"> 0.0283921 </td>
  </tr>
</tbody>
</table></div>

```r
message("Stats")
```

```
## Stats
```

```r
kable(head(gobpres$stats)) %>%
  kable_styling() %>%
  scroll_box(width = "1000px", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:1000px; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> stat.mean </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> exp1 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> GO:0019236 response to pheromone </td>
   <td style="text-align:right;"> 7.637046 </td>
   <td style="text-align:right;"> 7.637046 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0007420 brain development </td>
   <td style="text-align:right;"> 3.026272 </td>
   <td style="text-align:right;"> 3.026272 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0045664 regulation of neuron differentiation </td>
   <td style="text-align:right;"> 2.869812 </td>
   <td style="text-align:right;"> 2.869812 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0007156 homophilic cell adhesion </td>
   <td style="text-align:right;"> 2.879945 </td>
   <td style="text-align:right;"> 2.879945 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0001944 vasculature development </td>
   <td style="text-align:right;"> 2.746850 </td>
   <td style="text-align:right;"> 2.746850 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0033555 multicellular organismal response to stress </td>
   <td style="text-align:right;"> 2.719970 </td>
   <td style="text-align:right;"> 2.719970 </td>
  </tr>
</tbody>
</table></div>

## KEGG pathways


```r
library(gageData)
data(kegg.sets.mm)
data(sigmet.idx.mm)
kegg.sets.mm = kegg.sets.mm[sigmet.idx.mm]
#head(kegg.sets.mm, 3)

message("Greater")
```

```
## Greater
```

```r
kable(head(kegg.sets.mm$greater)) %>%
  kable_styling() %>%
  scroll_box(width = "1000px", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:1000px; "><table class="table" style="margin-left: auto; margin-right: auto;">
<tbody>
  <tr>

  </tr>
</tbody>
</table></div>

```r
message("Less")
```

```
## Less
```

```r
kable(head(kegg.sets.mm$less)) %>%
  kable_styling() %>%
  scroll_box(width = "1000px", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:1000px; "><table class="table" style="margin-left: auto; margin-right: auto;">
<tbody>
  <tr>

  </tr>
</tbody>
</table></div>

```r
message("Stats")
```

```
## Stats
```

```r
kable(head(kegg.sets.mm$stats)) %>%
  kable_styling() %>%
  scroll_box(width = "1000px", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:1000px; "><table class="table" style="margin-left: auto; margin-right: auto;">
<tbody>
  <tr>

  </tr>
</tbody>
</table></div>

```r
# Get the results
keggres = gage(foldchanges, gsets=kegg.sets.mm, same.dir=TRUE)

# Look at both up (greater), down (less), and statatistics.
#lapply(keggres, head, n=10)
message("Greater")
```

```
## Greater
```

```r
kable(head(keggres$greater)) %>%
  kable_styling() %>%
  scroll_box(width = "1000px", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:1000px; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.geomean </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> stat.mean </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.val </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> q.val </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> set.size </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> exp1 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> mmu04360 Axon guidance </td>
   <td style="text-align:right;"> 0.0285970 </td>
   <td style="text-align:right;"> 1.910570 </td>
   <td style="text-align:right;"> 0.0285970 </td>
   <td style="text-align:right;"> 0.8131848 </td>
   <td style="text-align:right;"> 128 </td>
   <td style="text-align:right;"> 0.0285970 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mmu03410 Base excision repair </td>
   <td style="text-align:right;"> 0.0355168 </td>
   <td style="text-align:right;"> 1.842982 </td>
   <td style="text-align:right;"> 0.0355168 </td>
   <td style="text-align:right;"> 0.8131848 </td>
   <td style="text-align:right;"> 28 </td>
   <td style="text-align:right;"> 0.0355168 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mmu00591 Linoleic acid metabolism </td>
   <td style="text-align:right;"> 0.0537749 </td>
   <td style="text-align:right;"> 1.628039 </td>
   <td style="text-align:right;"> 0.0537749 </td>
   <td style="text-align:right;"> 0.8131848 </td>
   <td style="text-align:right;"> 40 </td>
   <td style="text-align:right;"> 0.0537749 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mmu04974 Protein digestion and absorption </td>
   <td style="text-align:right;"> 0.0596166 </td>
   <td style="text-align:right;"> 1.567001 </td>
   <td style="text-align:right;"> 0.0596166 </td>
   <td style="text-align:right;"> 0.8131848 </td>
   <td style="text-align:right;"> 76 </td>
   <td style="text-align:right;"> 0.0596166 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mmu01040 Biosynthesis of unsaturated fatty acids </td>
   <td style="text-align:right;"> 0.0662476 </td>
   <td style="text-align:right;"> 1.531834 </td>
   <td style="text-align:right;"> 0.0662476 </td>
   <td style="text-align:right;"> 0.8131848 </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 0.0662476 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mmu03022 Basal transcription factors </td>
   <td style="text-align:right;"> 0.0674085 </td>
   <td style="text-align:right;"> 1.514546 </td>
   <td style="text-align:right;"> 0.0674085 </td>
   <td style="text-align:right;"> 0.8131848 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 0.0674085 </td>
  </tr>
</tbody>
</table></div>

```r
message("Less")
```

```
## Less
```

```r
kable(head(keggres$less)) %>%
  kable_styling() %>%
  scroll_box(width = "1000px", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:1000px; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.geomean </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> stat.mean </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.val </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> q.val </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> set.size </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> exp1 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> mmu03020 RNA polymerase </td>
   <td style="text-align:right;"> 0.0095316 </td>
   <td style="text-align:right;"> -2.445444 </td>
   <td style="text-align:right;"> 0.0095316 </td>
   <td style="text-align:right;"> 0.9111323 </td>
   <td style="text-align:right;"> 23 </td>
   <td style="text-align:right;"> 0.0095316 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mmu03450 Non-homologous end-joining </td>
   <td style="text-align:right;"> 0.0341877 </td>
   <td style="text-align:right;"> -1.999438 </td>
   <td style="text-align:right;"> 0.0341877 </td>
   <td style="text-align:right;"> 0.9111323 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 0.0341877 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mmu00190 Oxidative phosphorylation </td>
   <td style="text-align:right;"> 0.0669600 </td>
   <td style="text-align:right;"> -1.505167 </td>
   <td style="text-align:right;"> 0.0669600 </td>
   <td style="text-align:right;"> 0.9111323 </td>
   <td style="text-align:right;"> 99 </td>
   <td style="text-align:right;"> 0.0669600 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mmu03420 Nucleotide excision repair </td>
   <td style="text-align:right;"> 0.0747360 </td>
   <td style="text-align:right;"> -1.456591 </td>
   <td style="text-align:right;"> 0.0747360 </td>
   <td style="text-align:right;"> 0.9111323 </td>
   <td style="text-align:right;"> 39 </td>
   <td style="text-align:right;"> 0.0747360 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mmu00563 Glycosylphosphatidylinositol(GPI)-anchor biosynthesis </td>
   <td style="text-align:right;"> 0.0757982 </td>
   <td style="text-align:right;"> -1.464592 </td>
   <td style="text-align:right;"> 0.0757982 </td>
   <td style="text-align:right;"> 0.9111323 </td>
   <td style="text-align:right;"> 25 </td>
   <td style="text-align:right;"> 0.0757982 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mmu00240 Pyrimidine metabolism </td>
   <td style="text-align:right;"> 0.1116636 </td>
   <td style="text-align:right;"> -1.222140 </td>
   <td style="text-align:right;"> 0.1116636 </td>
   <td style="text-align:right;"> 0.9111323 </td>
   <td style="text-align:right;"> 87 </td>
   <td style="text-align:right;"> 0.1116636 </td>
  </tr>
</tbody>
</table></div>

```r
message("Stats")
```

```
## Stats
```

```r
kable(head(keggres$stats)) %>%
  kable_styling() %>%
  scroll_box(width = "1000px", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:1000px; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> stat.mean </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> exp1 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> mmu04360 Axon guidance </td>
   <td style="text-align:right;"> 1.910570 </td>
   <td style="text-align:right;"> 1.910570 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mmu03410 Base excision repair </td>
   <td style="text-align:right;"> 1.842982 </td>
   <td style="text-align:right;"> 1.842982 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mmu00591 Linoleic acid metabolism </td>
   <td style="text-align:right;"> 1.628039 </td>
   <td style="text-align:right;"> 1.628039 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mmu04974 Protein digestion and absorption </td>
   <td style="text-align:right;"> 1.567001 </td>
   <td style="text-align:right;"> 1.567001 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mmu01040 Biosynthesis of unsaturated fatty acids </td>
   <td style="text-align:right;"> 1.531834 </td>
   <td style="text-align:right;"> 1.531834 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mmu03022 Basal transcription factors </td>
   <td style="text-align:right;"> 1.514546 </td>
   <td style="text-align:right;"> 1.514546 </td>
  </tr>
</tbody>
</table></div>


```r
# Get the pathways
keggrespathways = data.frame(id=rownames(keggres$less), keggres$less) %>% 
  tbl_df() %>% 
  filter(row_number()<=5) %>% 
  .$id %>% 
  as.character()
keggrespathways
```

```
## [1] "mmu03020 RNA polymerase"                                       
## [2] "mmu03450 Non-homologous end-joining"                           
## [3] "mmu00190 Oxidative phosphorylation"                            
## [4] "mmu03420 Nucleotide excision repair"                           
## [5] "mmu00563 Glycosylphosphatidylinositol(GPI)-anchor biosynthesis"
```


```r
# Get the IDs.
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids
```

```
## [1] "mmu03020" "mmu03450" "mmu00190" "mmu03420" "mmu00563"
```


```r
# Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="mmu", new.signature=FALSE)

# plot multiple pathways (plots saved to disk and returns a throwaway list object)
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="mmu"))
```

## Pathway and regulation of genes for Oxidative phosphorylation.

![Oxidative Phosphorylation](./Mouse_RNA_Seq_p53_genotoxic_files/figure-html/mmu00190.pathview.png)

## Pathway and regulation of genes for Glycosylphosphatidylinositol(GPI)-anchor biosynthesis.

![Glycosylphosphatidylinositol(GPI)-anchor biosynthesis](./Mouse_RNA_Seq_p53_genotoxic_files/figure-html/mmu00563.pathview.png)

## Pathway and regulation of genes for RNA polymerase.

![RNA polymerase](./Mouse_RNA_Seq_p53_genotoxic_files/figure-html/mmu03020.pathview.png)

## Pathway and regulation of genes for Nucleotide excision repair.

![Nucleotide excision repair](./Mouse_RNA_Seq_p53_genotoxic_files/figure-html/mmu03420.pathview.png)

## Pathway and regulation of genes for Non-homologous end-joining.

![Non-homologous end-joining](./Mouse_RNA_Seq_p53_genotoxic_files/figure-html/mmu03450.pathview.png)








































































































