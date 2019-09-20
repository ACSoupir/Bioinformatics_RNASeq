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

# **Analysis**

## Programs used?

+ FastQC
+ Trimmomatic 0.39
+ Bowtie 2.2.5.0
+ Tophat 2.1.1
+ STAR
+ Cufflinks

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









