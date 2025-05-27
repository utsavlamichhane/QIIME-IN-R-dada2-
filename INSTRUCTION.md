# DADA2 Pipeline Tutorial in R

Here I will guide you step by step to run dada2 command in your system using R (preferably R studio). 

You will learn how to:

1. Inspect and filter raw reads  
2. Learn error rates  
3. Dereplicate and infer sequence variants  
4. Merge paired reads  
5. Build a sequence table and remove chimeras  
6. Track reads through the pipeline  
7. (Optional) Assign taxonomy

All example FASTQ files live in a folder called `fastq_files/`. We’ll use three samples (6 FASTQ files total):

fastq_files/
├── Sample_134_R1_001.fastq.gz
├── Sample_134_R2_001.fastq.gz
├── Sample_135_R1_001.fastq.gz
├── Sample_135_R2_001.fastq.gz
├── Sample_136_R1_001.fastq.gz
└── Sample_136_R2_001.fastq.gz



Step 0: Download the fastq files

Download the entire directory from the following URL:

```
https://drive.google.com/drive/folders/18kSe9UBJgshaxdN355GLAPBPI8I5IKYn?usp=sharing
```



Step 1: Library installation and loading


```r
BiocManager::install("dada2", force = TRUE)
library(dada2)

```
ok lets confirm that the package version

```r
packageVersion("dada2")
```

Step2: 2. Define File Paths and Sample Names

Go to session at the top of your R-studio.

Click on the "Set Working Directory"

Then click on "Choose Directory" 

Then you will see something like the following in your console:

```
> setwd("~/Desktop/DESKTOP/CLASSES/6_classes_summer_2025/2_TAing/4_my_files/fastq_files")
```
That is your path.

Lets define the path with the following code:

```
path <- "~/Desktop/DESKTOP/CLASSES/6_classes_summer_2025/2_TAing/4_my_files/fastq_files"
```
