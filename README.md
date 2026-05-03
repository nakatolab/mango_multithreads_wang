# Modified mango for multi-threads and replicates

## 1.What is the bugs

One biggest problem of using Mango pipeline for analyzing ChIA-PET is that they can not support multi threads (cpu):
> Threads option is currently disabled to due to errors. We are working on a solution !!

Their bugs comes from the bowtie mapping step. Two important point for this bugs:
- Mango requires that Fastq files are perfectly matched
- The sam files should still be in the same order as the fastq files

If you use multi-threads for bowtie alignments, the sorts of reads in sam file would be 'out-of-order'

Another point is that:
- SAM files being produced by Mango are not canonical SAM files (i.e., they have no headers)

So we cannot use samtools to sort.

## 2.How to solve

Instead, I modified the raw script to deal with reads (fastq) from SRA or Encode:
- If your reads name of fastq is like:
  ``` shell
  $ head ENCSR981FNA_1.fastq
  @HWI-ST689:335:C7ETCACXX:8:1101:1270:1964 1:N:0:
  GTNGACGCCCACAGGGGGCAGGGTCTCGCCCTTGTAGCGTAAGGTTCCCCCAACATTGGCCACAGAGCCGTTGATGACGACAGCAGTTGGATAAGATATCG
  +
  ```
  Please use the `mango_encode.R` and the same parameter of the original Mango
  
- If your reads name of fastq is like:
  ``` shell
  $ head SRR372741_1.fastq
  @SRR372741.1 B2GA003:1:2:1142:995 length=36
  NCCTTTTCAAATTTACCTAC
  +SRR372741.1 B2GA003:1:2:1142:995 length=36
  #111253566@@CC@C@C@@
  ```
   Please use the `mango_SRA.R` and the same parameter of the original Mango
   
- If others:
  Please modify the `alignBowtie` function (line 17) in `mango_encode.R` or `mango_SRA.R` to make sure you can obtain the right sorted sam file ~

## 3.Install (for my personal use, same with the original Mango, )

### step1 shell:
``` shell
git clone https://github.com/YuNagaoka/mango_multithreads_wang.git
mv mango_multithreads_wang mango
R CMD INSTALL --no-multiarch --with-keep.source mango
```

### step2 (required package) R:
``` R
install.packages('hash')
install.packages('Rcpp')
install.packages('optparse')
install.packages('readr')
```

### step3 run the Mango

**Pattern A — single replicate (or explicit file paths)**
``` shell
sample=CTCF
FASTQ=$(ls fastq/$sample | cut -d '_' -f 1 | sort -u)

index=bowtie-indexes/genome
gt=genometable.txt
blacklist=ENCODE-Blacklist/hg38-blacklist.v2.bed

mango="Rscript /opt/mango/mango_SRA.R"
mkdir -p mango/$sample

for i in $FASTQ
do
    setsid $mango --stages 1:5 \
           --prefix ${sample}_${i} \
           --outdir $outdir \
           --chromexclude chrM,chrY \
           --bowtieref $index \
           --bedtoolsgenome $gt \
           --fastq1 fastq/$sample/${i}_1.fastq.gz \
           --fastq2 fastq/$sample/${i}_2.fastq.gz \
           --linkerA CGCGATATCTTATCTGACT \
           --singlelinker TRUE \
           --minlength 15 \
           --maxlength 1000 \
           --keepempty TRUE \
           --threads 10 \
           --shortreads FALSE \
           --MACS_qvalue 0.05 \
           --blacklist $blacklist \
           --reportallpairs TRUE
done
```

**Pattern B — multiple replicates (auto-merge with `--fastqdir`)**

・As an additional change, I have added the `--fastqdir` option, which automatically merges multiple replicates found in the fastq directory.

`--fastqdir` behavior:
- **2 or more** `*_1.fastq.gz` / `*_R1.fastq.gz` , `*_2.fastq.gz` / `*_R2.fastq.gz` pairs found in fastq directory
→ files are merged using `cat` command (without decompression) into temporary files `{prefix}_merged_1.fastq.gz` / `{prefix}_merged_2.fastq.gz`, which are automatically deleted after Stage 1 completes.

- **Exactly 1** pair found → that file is used directly (no merge, no temporary copy).
- Output files are always named using `--prefix`, regardless of whether a merge was performed.

Place all replicate FASTQ files in same directory:

```shell
fastq/CTCF/
  rep1_1.fastq.gz  rep1_2.fastq.gz
  rep2_1.fastq.gz  rep2_2.fastq.gz
  ...
```
Then pass the directory with `--fastqdir`.  
`--fastq1` / `--fastq2` are **not** required when `--fastqdir` is used.

``` shell
sample=CTCF

index=bowtie-indexes/genome
gt=genometable.txt
blacklist=ENCODE-Blacklist/hg38-blacklist.v2.bed

mango="Rscript /opt/mango/mango_SRA.R"
mkdir -p mango/$sample

setsid $mango --stages 1:5 \
       --prefix $sample \
       --outdir mango/$sample \
       --chromexclude chrM,chrY \
       --bowtieref $index \
       --bedtoolsgenome $gt \
       --fastqdir fastq/$sample \
       --linkerA CGCGATATCTTATCTGACT \
       --singlelinker TRUE \
       --minlength 15 \
       --maxlength 1000 \
       --keepempty TRUE \
       --threads 10 \
       --shortreads FALSE \
       --MACS_qvalue 0.05 \
       --blacklist $blacklist \
       --reportallpairs TRUE
```

The additional changes made to the code from the original mango are listed in `fixed_error.log`.
