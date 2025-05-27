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

Lets see how the filter worked using the following command:

```r
head(out)
```

There are so many things we can do but lets keep things simple.

# Step 7: Learn Error Rates

Whats error rate ?

In simple mathematical term:

Error rate = (number of mistakes) ÷ (number of opportunities)

If at quality score 20 you saw 10,000 bases, and 100 of them were actually wrong, your observed error rate is 100 / 10 000 = 0.01 (1%).

lets get the error rate and plot it with the following code:

```r
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

###plotting
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
```

# Step 7: Dereplication

What is derep then ?

whenever we us dada2, it does two very important things:

- Keeps one copy of each unique sequence
- Counts how many times that exact sequence occurred

So instead of storing 10,000 identical “ACGT….” reads, you store one “ACGT….” entry plus a count of 10,000.

ok, lets keep moving and do the actual dereplication with the following command: 

```r
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
```
# Step 8: infering from the samples 

In simple words, it is the step where DADA2 examines your filtered reads and “figures out” which unique sequences are genuine biological variants and which are just noise.

What I mean is the result is, for each sample, a clean list of exact sequence variants (ASVs) with their abundances, ready for merging forward + reverse reads and downstream analysis.

Nothing fancy just checking unique things in our object (dada2 object) from dada2 just like in qiime2. 

Run the following code for the sample inference:

```r
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
```
# Step 8: Merge Paired Reads

We do this to get the high-confidence, full-length amplicon sequence variants (ASVs) for each sample, ready to be turned into a sequence table for downstream analysis.

Run the following command to merge the sequences:

```r
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
```

# Step 9: Make Sequence Table & Remove Chimeras

To create the sequence table as an R object run the following code:

```r
seqtab <- makeSequenceTable(mergers)
dim(seqtab)  

seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus",
                                    multithread=TRUE, verbose=TRUE)
```


# Step 10: Track our pipeline

ok, we have been playing a lot with the pipelines lets keep track for our record or for "reviewer 2".

For the tack of the reads in our pipeline run the following code:

```r
getN <- function(x) sum(getUniques(x))
track <- cbind(out,
               sapply(dadaFs, getN),
               sapply(mergers, getN),
               rowSums(seqtab_nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "nonchim")
rownames(track) <- sample.names
print(track)
```


