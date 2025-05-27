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



# Step 0: Download the fastq files

Download the entire directory from the following URL:

From one drive:
```
https://outlookuga-my.sharepoint.com/:f:/g/personal/ul54354_uga_edu/EjbLTEqvkRVOj4qYO5WFJWwBbCXy_aJDvEtLUZz3s3zUeQ?e=FcuogT
```
Or from google drive: 

```
https://drive.google.com/drive/folders/18kSe9UBJgshaxdN355GLAPBPI8I5IKYn?usp=sharing
```



# Step 1: Library installation and loading


```r
BiocManager::install("dada2", force = TRUE)
library(dada2)

```
ok lets confirm that the package version

```r
packageVersion("dada2")
```

# Step 2: Define File Paths and Sample Names

Go to session at the top of your R-studio.

Click on the "Set Working Directory"

Then click on "Choose Directory" 

Then you will see something like the following in your console:

```
> setwd("~/Desktop/DESKTOP/CLASSES/6_classes_summer_2025/2_TAing/4_my_files/fastq_files")
```
That is your path.

Ok lets see your file path with the following commad:

```
list.files(path)
```

Lets define the path with the following code:

```
path <- "~/Desktop/DESKTOP/CLASSES/6_classes_summer_2025/2_TAing/4_my_files/fastq_files"
```

Ok lets see your file path with the following commad:

```
list.files(path)
```
You should see the fastq fiels with the full name and extension. 

# Step 4: List Forwards and Reverse Reads,  and get the extracted sample name.

We will assign F and R reads with the following code:

```r
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names=TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names=TRUE))
```
And, lets make the sample names as an R object using the following code: 

```r
sample.names <- sapply(strsplit(basename(fnFs), "_"), function(x) paste(x[1], x[2], sep="_"))
```

Now, the real deal to look if your code is actually doing what you intend to do.

# Step 5: Vizualization of the quality plots
```r
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
```

You should see Q-plots in your "Plots" pannel.

or the plot directly pops if you are a Quartz person. 

# Step 6: Filter and Trim

Firstly, lets assign a directory for the filtered outs, using the following command:

```r
filt_path <- file.path(path, "filtered")
if(!dir.exists(filt_path)) dir.create(filt_path)
```

In the above code i have dir.create just in case if the directory does not exist.

Now, we define the filtered outs (file path this time).

Run the following code:

```r
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
```

Now lets run the actual filtering command:

```r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=c(248,239),   
                     maxN=1,                
                     maxEE=c(2,2),          
                     truncQ=20,             
                     compress=TRUE,         
                     multithread=TRUE)      
```
Few things about the parametes used here;

- truncLen=c(248,239): Here, forward reads are cut at 248 bp, reverse reads at 239 bp.
  Any base beyond that position is discarded.
- maxN=1: Whenever a sequencer can’t confidently call a base (A, C, G or T) at a given position, it writes an “N” there. That “N” is literally an ambiguous base—meaning “I don’t know what this is.”
- truncQ=20, thats the Q score we all know.
- compress=TRUE, means it will zip the file.
- multithread=TRUE, parallel processing we all know.


  
