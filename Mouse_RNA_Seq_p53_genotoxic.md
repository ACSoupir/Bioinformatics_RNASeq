---
title: "RNA Seq for STAT736"
author: "Alex Soupir"
date: "9/10/2019"
output: 
  html_document:
    keep_md: true
---



## Background

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

## Analysis

### Acquiring sequences

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

### FastQC

<<<<<<< HEAD
FastQC can be run over all of the **\.fastq.gz** files by using the * wildcard. This prevents the need to hard code each individual read file, which is nice because there are 24 (12 forward and 12 reverse).
=======
<<<<<<< HEAD
FastQC can be run over all of the **\.fastq.gz** files by using the **\* ** wildcard. This prevents the need to hard code each individual read file, which is nice because there are 24 (12 forward and 12 reverse).
=======
FastQC can be run over all of the **\.fastq.gz** files by using the **\* ** wildcard.
>>>>>>> 71cc2bd27dd47a95f86144834e29602051004210
>>>>>>> 63272c13fe3eb1b3232a480bcc8e886d43ba2acd


```bash
~/miniconda2/bin/fastqc *.fastq.gz
```

<<<<<<< HEAD
=======
<<<<<<< HEAD
The output of running FastQC is a zip folder and an html that can be viewed in a browser. [This](./SRR2121770_1_fastqc.html) is an example 

=======
>>>>>>> 71cc2bd27dd47a95f86144834e29602051004210
>>>>>>> 63272c13fe3eb1b3232a480bcc8e886d43ba2acd



xczvz












































