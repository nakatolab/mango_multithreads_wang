# runs mango chia pet analysis pipeline

# Version info
Mangoversion = "1.2.0"

# Load Packages
suppressPackageStartupMessages(library("Rcpp"))
suppressPackageStartupMessages(library("hash"))
suppressPackageStartupMessages(library("mango"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))

print ("Starting mango ChIA PET analysis tool")
Sys.time()
set.seed(1)

alignBowtie <- function(fastq,output,bowtiepath,bowtieref,
                        shortreads,threads,verbose=TRUE,nlines=10000)
{
  
  # choose alignment parameters
  # note- "-m 1" ensures that only uniquely mapped reads are reported.
  bowtievar=paste("-S -v 0 -k 1 --chunkmbs 500 --sam-nohead --mapq 40 -m 1","--threads",threads)
  if (shortreads == FALSE)
  {
    #bowtievar="-S -v 2 -k 1 --sam-nohead --mapq 40 -m 1 --best"
    bowtievar=paste("-S -n 2 -l 50 -k 1 --chunkmbs 500 --sam-nohead --mapq 40 -m 1 --best","--threads",threads)
  }
  
  # determine the illumina score encoding
  illuminascore = findScore(fastq,nlines =nlines) 
  
  # form full command
  bowtiecommand = paste (bowtiepath, bowtieref, fastq, output, illuminascore, bowtievar)
  
  # print command if desirec
  if (verbose ==TRUE)
  {
    print (bowtiecommand)
  }
  
  # execute command
  # Wrap in bash with SIGHUP ignored so that terminal disconnect does not kill
  # bowtie mid-alignment (R resets SIGHUP handling on startup, overriding nohup).
  exitcode = system(paste0("bash -c 'trap \"\" HUP; ", bowtiecommand, "'"))
  if (exitcode != 0) {
    stop(paste("bowtie alignment failed with exit code", exitcode,
               "for input:", fastq))
  }
  # sort SAM by read name (second dot-delimited field) so paired reads align
  # across both SAM files when buildBedpe reads them sequentially
  sortcpu = paste0('--parallel=', threads)
  sort_exitcode = system(paste0("bash -c 'trap \"\" HUP; sort -t . -k 2 -n ",
                                sortcpu, " ", output, " -o ", output, "'"))
  if (sort_exitcode != 0) {
    stop(paste("sort of SAM file failed with exit code", sort_exitcode,
               "for file:", output))
  }

}



##################################### read commandline paramters #####################################

# read in parameters
option_list <- list(
  
  #---------- GENERAL PARAMETERS ----------#
  
  make_option(c("--stages"),  default="1:5",help="stages of the pipeline to execute"),
  make_option(c("--prefix"), default="mango",help="prefix for all output files"),
  make_option(c("--outdir"),  default="NULL",help="out put directory"),
  make_option(c("--argsfile"),  default="NULL",help="optional argument file used to pass in command line paramters"),
  make_option(c("--bowtieref"),   default="NULL",help="genome reference file for bowtie"),
  make_option(c("--bedtoolsgenome"),  default="NULL",help="bedtools genome file"),
  make_option(c("--chrominclude"),  default="NULL",help="comma separated list of chromosomes to use (e.g. chr1,chr2,chr3,...).  Only these chromosomes will be processed"),
  make_option(c("--chromexclude"),  default="NULL",help="comma separated list of chromosomes to exclude (e.g. chrM,chrY). !!chrM should always be excluded due to its extremely short length!!"),
  make_option(c("--bedtoolspath"),  default="NULL",help="full path to bedtools"),
  make_option(c("--macs2path"),  default="NULL",help="full path to macs2path"),
  make_option(c("--bowtiepath"),  default="NULL",help="full path to bowtiepath"),
  make_option(c("--verboseoutput"),  default="FALSE",help="if true output file will have more columns as well as a header row describing the columns"),
  
  #---------- STAGE 1 PARAMETERS ----------#
  
  make_option(c("--fastq1"),  default="NULL",help="fastq read 1 file"),
  make_option(c("--fastq2"),  default="NULL",help="fastq read 2 file"),
  make_option(c("--linkerA"),  default="GTTGGATAAG",help="linker sequence A to look for"),
  make_option(c("--linkerB"),  default="GTTGGAATGT",help="linker sequence B to look for"),
  make_option(c("--singlelinker"),  default="FALSE",help="If this is true Mango will only look for linkerA"),
  make_option(c("--minlength"),  default="15",help="min length of reads after linker trimming"),
  make_option(c("--maxlength"),  default="25",help="max length of reads after linker trimming"),
  make_option(c("--keepempty"),  default="FALSE",help="Should reads with no linker be kept"),
  
  #---------- STAGE 2 PARAMETERS ----------#
  
  make_option(c("--shortreads"),  default="TRUE",help="should bowtie alignments be done using paramter for very short reads (~20 bp)"),
  make_option(c("--threads"),  default="1",help="(!! This option is currently disabled to due to errors.  We are working on a solution !!) number of threads to be used for bowtie alignment. default = 1"),
  
  #---------- STAGE 4 PARAMETERS ----------#
  
  make_option(c("--npets4dist"),  default="1000000",help="the number of PETS to use to plot PET distance distribution (default = 1000000, use -1 for all PETS)"),
  
  #---------- STAGE 4 PARAMETERS ----------#
  
  make_option(c("--MACS_qvalue"),  default="0.05",help="MACS values"),
  make_option(c("--MACS_shiftsize"),  default="NULL",help="MACS shiftize.  NULL allows MACS to determine it"),
  make_option(c("--peakslop"),  default="500",help="Number of basespairs to extend peaks on both sides"),
  make_option(c("--peakinput"),  default="NULL",help="user supplied peaks file"),
  make_option(c("--blacklist"),  default="NULL",help="BED file of regions to remove from MACS peaks"),
  make_option(c("--gsize"),  default="hs",help="mappable genome size or effective genome size for MACS2"),
  
  #---------- STAGE 5 PARAMETERS ----------#
  
  make_option(c("--distcutrangemin"),  default="1000",help="range in which to look for the self-ligation distance"),
  make_option(c("--distcutrangemax"),  default="100000",help="range in which to look for the self-ligation distance"),
  make_option(c("--biascut"),  default="0.05",help="Self ligation bias cutoff"),
  make_option(c("--numofbins"),  default="50",help="number of bins for probability calculations"),
  make_option(c("--corrMethod"),  default="BH",help="multiple hypothesis tersting correction method"),
  make_option(c("--maxinteractingdist"),  default="1000000",help="maximum disance allowed for an interaction"),
  make_option(c("--FDR"),  default="0.05",help="FDR cutoff for interactions"),
  make_option(c("--extendreads"),  default="120",help="how many bp to extend reads towards peak"),
  make_option(c("--minPETS"),  default="2",help="minimum number of PETs required for an interaction (applied after FDR filtering)"),
  make_option(c("--reportallpairs"),  default="FALSE",help="Should all pairs be reported or just significant pairs"),
  make_option(c("--MHT"),  default="all",help="How should mutliple hypothsesis testing be done?  Correct for 'all' possible pairs of loci or only those 'found' with at least 1 PET"),

  #---------- FASTQ MERGE PARAMETERS ----------#

  make_option(c("--fastqdir"),  default="NULL",help="directory containing replicate *_1.fastq.gz / *_2.fastq.gz files.  When 2 or more pairs are found they are concatenated (without decompression) before stage 1; the merged files are removed after stage 1 completes."))

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))

# check dependencies

# first look in PATH
progs = c("bedtools","macs2","bowtie")
Paths = DefinePaths(progs = progs)
bedtoolspath  = Paths[1]
macs2path     = Paths[2]
bowtiepath    = Paths[3]


# next, look in arguments
if (opt["bedtoolspath"] != "NULL")
{
  bedtoolspath = opt["bedtoolspath"]
}
if (opt["macs2path"] != "NULL")
{
  macs2path = opt["macs2path"]
}
if (opt["bowtiepath"] != "NULL")
{
  bowtiepath = opt["bowtiepath"]
}

if (identical(bedtoolspath, "notfound"))
{
  stop("bedtools not found in path.")
}  

if (identical(macs2path, "notfound"))
{
  stop("macs2 not found in path.")
} 

if (identical(bowtiepath, "notfound"))
{
  stop("bowtie not found in path.")
} 

# get software versions
bedtoolsversion = system(paste(bedtoolspath,"--version"),intern=TRUE)[1]
macs2version    = system(paste(macs2path,"--version 2>&1"),intern=TRUE)[1]
bowtieversion   = system(paste(bowtiepath,"--version"),intern=TRUE)[1]

# print out software versions
print (paste("bedtools version:",bedtoolsversion))
print (paste("macs2 version:",macs2version))
print (paste("bowtie version:",bowtieversion))

# break if dependencies not found
Paths = c(bedtoolspath,macs2path,bowtiepath)
i = 0
pathsOK = T
for (p in Paths)
{
  i = i+ 1
  if (p == "notfound")
  {
    pathsOK = F
    print ("! Configuration Error !")
    print (paste("     Path to ",progs[i]," not in PATH or arguments",sep=""))
    print ("     Please add to PATH or arguments and try again")
    print ("")
  }
}
if (pathsOK == F)
{
  stop ("Exiting Mango.R pipeline.  Check system PATH")
}

# correct stages
if (grepl( ":",opt$stages))
{
  stages = strsplit(opt$stages,split=":")[[1]]
  if (length(stages) == 2)
  {
    opt$stages = seq(as.numeric(stages[1]),as.numeric(stages[2]))
  }
  if (length(stages) == 1)
  {
    opt$stages = as.numeric(stages[1])
  }
}

# set basename for output
opt["outname"] = ""
if (opt["outdir"] == "NULL")
{
  opt["outname"] = file.path(getwd(),opt["prefix"])
}
if (opt["outdir"] != "NULL")
{
  opt["outname"] = file.path( opt["outdir"] ,opt["prefix"]) 
}

logfile = paste(as.character(opt["outname"]),".mango.log",sep="")
if (file.exists(logfile) ==TRUE){file.remove(logfile)}
starttime = paste("Analysis start time:" , as.character(Sys.time()))
write(starttime,file=logfile,append=TRUE)

errorlog = paste(as.character(opt["outname"]), ".error.log", sep="")
if (file.exists(errorlog) == TRUE){file.remove(errorlog)}

##################################### read in arguments #####################################

lines = c()
if (as.character(opt$argsfile) != "NULL")
{
  if (file.exists(as.character(opt$argsfile)) == FALSE)
  {
    stop (paste("Exiting Mango.R pipeline.  No file named ",opt$argsfile,sep=""))
  }
  lines = readLines(as.character(opt$argsfile))
}

for (line in lines)
{
  # remove spaces
  line = gsub(pattern=" ",x=line,replace="")
  
  if (line == "")
  {
    next
  }
  else if (strsplit(line,split="")[[1]][1] == "#" )
  {
    next
  }
  else
  {
    lineinfo = strsplit(line,split="#")[[1]][1]
    arginfo  = strsplit(lineinfo,split="=")[[1]]
    opt[arginfo[1]] = arginfo[2]
  }
}
       
resultshash = hash()
stage5_warnings = c()
errors_occurred = c()

##################################### merge replicate fastq.gz files if fastqdir is given #####################################

merged_fastq1 <- NULL
merged_fastq2 <- NULL

if (as.character(opt["fastqdir"]) != "NULL") {
  fastqdir <- as.character(opt["fastqdir"])
  if (!dir.exists(fastqdir)) {
    stop(paste("fastqdir does not exist:", fastqdir))
  }

  files1_num <- sort(list.files(fastqdir, pattern = "_1\\.fastq\\.gz$", full.names = TRUE))
  files1_R   <- sort(list.files(fastqdir, pattern = "_R1\\.fastq\\.gz$", full.names = TRUE))

  if (length(files1_num) > 0 && length(files1_R) > 0) {
    stop(paste("Mixed naming: both *_1.fastq.gz and *_R1.fastq.gz files found in fastqdir:", fastqdir))
  }

  if (length(files1_num) > 0) {
    files1 <- files1_num
    files2_num <- sort(list.files(fastqdir, pattern = "_2\\.fastq\\.gz$", full.names = TRUE))
    files2_R   <- sort(list.files(fastqdir, pattern = "_R2\\.fastq\\.gz$", full.names = TRUE))
    if (length(files2_num) > 0 && length(files2_R) > 0) {
      stop(paste("Mixed naming: both *_2.fastq.gz and *_R2.fastq.gz files found in fastqdir:", fastqdir))
    }
    files2 <- files2_num
    suffix1 <- "_1"; suffix2 <- "_2"
  } else if (length(files1_R) > 0) {
    files1 <- files1_R
    files2_R   <- sort(list.files(fastqdir, pattern = "_R2\\.fastq\\.gz$", full.names = TRUE))
    files2_num <- sort(list.files(fastqdir, pattern = "_2\\.fastq\\.gz$", full.names = TRUE))
    if (length(files2_R) > 0 && length(files2_num) > 0) {
      stop(paste("Mixed naming: both *_R2.fastq.gz and *_2.fastq.gz files found in fastqdir:", fastqdir))
    }
    files2 <- files2_R
    suffix1 <- "_R1"; suffix2 <- "_R2"
  } else {
    stop(paste("No *_1.fastq.gz or *_R1.fastq.gz files found in fastqdir:", fastqdir))
  }

  if (length(files1) != length(files2)) {
    stop(paste("Mismatch: found", length(files1), paste0("*", suffix1, ".fastq.gz"), "files but",
               length(files2), paste0("*", suffix2, ".fastq.gz"), "files in fastqdir:", fastqdir))
  }

  if (length(files1) == 1) {
    print(paste("fastqdir: single replicate found, using directly:", files1[1]))
    opt["fastq1"] <- files1[1]
    opt["fastq2"] <- files2[1]
  } else {
    print(paste("fastqdir:", length(files1), "replicates found, merging with cat (no decompression)"))
    outname_for_merge <- as.character(opt["outname"])
    merged_fastq1 <- paste0(outname_for_merge, "_merged_1.fastq.gz")
    merged_fastq2 <- paste0(outname_for_merge, "_merged_2.fastq.gz")

    cmd1 <- paste("cat", paste(shQuote(files1), collapse = " "), ">", shQuote(merged_fastq1))
    cmd2 <- paste("cat", paste(shQuote(files2), collapse = " "), ">", shQuote(merged_fastq2))

    print(paste("Merging _1 files:", cmd1))
    ret1 <- system(cmd1)
    if (ret1 != 0) stop(paste("Failed to merge *_1.fastq.gz files. Command:", cmd1))

    print(paste("Merging _2 files:", cmd2))
    ret2 <- system(cmd2)
    if (ret2 != 0) {
      if (file.exists(merged_fastq1)) file.remove(merged_fastq1)
      stop(paste("Failed to merge *_2.fastq.gz files. Command:", cmd2))
    }

    print(paste("Merged fastq1:", merged_fastq1))
    print(paste("Merged fastq2:", merged_fastq2))
    opt["fastq1"] <- merged_fastq1
    opt["fastq2"] <- merged_fastq2
  }
}

##################################### parse fastqs #####################################

if (1 %in% opt$stages)
{
  tryCatch({
  checkRequired(opt,c("fastq1","fastq2"))

  # gather arguments
  outname         = as.character(opt["outname"])
  fastq1=as.character(opt["fastq1"])
  fastq2=as.character(opt["fastq2"])
  basename = as.character(opt["outname"])
  minlength = as.numeric(as.character(opt["minlength"]))
  maxlength = as.numeric(as.character(opt["maxlength"]) )
  keepempty=eval(parse(text=as.character(opt["keepempty"])))
  linker1=as.character(opt["linkerA"])
  linker2=as.character(opt["linkerB"])
  singlelinker = as.character(opt["singlelinker"])
  numberlinkers = 2
  if (singlelinker == "TRUE")
  {
    numberlinkers = 1
    print ("singlelinker set to TRUE.  Only looking for one linker sequence")
  }
  
  print ("finding linkers")
  
  parsingresults = parseFastq( fastq1=fastq1,
              fastq2=fastq2,
              basename = basename,
              minlength = minlength,
              maxlength = maxlength, 
              keepempty = keepempty,
              linker1=linker1,
              linker2=linker2,
              numberlinkers=numberlinkers)
  
  # parseFastq returns [total, same, chim, ambi] (R 1-based indexing)
  resultshash[["total PETs"]] = as.numeric(parsingresults[1])
  resultshash[["same PETs"]] = as.numeric(parsingresults[2])
  resultshash[["chimeric PETs"]] = as.numeric(parsingresults[3])
  resultshash[["ambigious PETs"]] = as.numeric(parsingresults[4])
  }, error = function(e) {
    errmsg <- conditionMessage(e)
    cat(paste0("[ERROR] Stage 1 failed: ", errmsg, "\n"), file=stderr())
    errors_occurred <<- c(errors_occurred, paste0("Stage 1: ", errmsg))
    write(paste0("[ERROR] Stage 1: ", errmsg), file=errorlog, append=TRUE)
  })

  # Remove merged fastq.gz files created by --fastqdir to free disk space
  for (mf in c(merged_fastq1, merged_fastq2)) {
    if (!is.null(mf) && file.exists(mf)) {
      file.remove(mf)
      print(paste("Removed merged file:", mf))
    }
  }
}
  
###################################### align reads #####################################

if (2 %in% opt$stages)
{
  tryCatch({
  checkRequired(opt,c("bowtieref"))
  
  # gather arguments
  outname         = as.character(opt["outname"])
  bowtieref       = as.character(opt["bowtieref"])
  shortreads      = as.character(opt["shortreads"])
  threads         = as.character(opt["threads"])
  
  print ("aligning reads")
  # filenames
  fastq1 = paste(outname ,"_1.same.fastq",sep="")
  fastq2 = paste(outname ,"_2.same.fastq",sep="")
  sam1   = paste(outname ,"_1.same.sam",sep="")
  sam2   = paste(outname ,"_2.same.sam",sep="")
  
  # Validate stage-1 FASTQ outputs before attempting alignment
  for (fq in c(fastq1, fastq2)) {
    if (!file.exists(fq)) {
      stop(paste("Stage-1 FASTQ not found:", fq,
                 "- check that stage 1 completed successfully"))
    }
    if (file.info(fq)$size == 0) {
      stop(paste("Stage-1 FASTQ is empty (0 bytes):", fq,
                 "- linker detection may have produced no valid reads"))
    }
    con_check <- file(fq, open = "r")
    first_line <- tryCatch(
      readLines(con_check, n = 1, warn = FALSE),
      error = function(e) character(0),
      finally = close(con_check)
    )
    if (length(first_line) == 0 || !grepl("^@", first_line[1])) {
      stop(paste("Stage-1 FASTQ does not start with '@' (corrupt or not a valid FASTQ):", fq))
    }
  }
  
  # align both ends of each PET
  alignBowtie(fastq=fastq1,output=sam1,bowtiepath=bowtiepath,bowtieref=bowtieref,shortreads,threads)
  alignBowtie(fastq=fastq2,output=sam2,bowtiepath=bowtiepath,bowtieref=bowtieref,shortreads,threads)
  }, error = function(e) {
    errmsg <- conditionMessage(e)
    cat(paste0("[ERROR] Stage 2 failed: ", errmsg, "\n"), file=stderr())
    errors_occurred <<- c(errors_occurred, paste0("Stage 2: ", errmsg))
    write(paste0("[ERROR] Stage 2: ", errmsg), file=errorlog, append=TRUE)
  })
}

##################################### filter reads #####################################

if (3 %in% opt$stages)
{
  tryCatch({
  checkRequired(opt,c("outname"))
  
  # gather arguments
  outname         =  as.character(opt["outname"])
  
  # filenames
  bedpefile          = paste(outname ,".bedpe",sep="")
  bedpefilesort      = paste(outname ,".sort.bedpe",sep="")
  bedpefilesortrmdup = paste(outname ,".sort.rmdup.bedpe",sep="")
  sam1               = paste(outname ,"_1.same.sam",sep="")
  sam2               = paste(outname ,"_2.same.sam",sep="")
  bedpermdupfile     = paste(outname ,".rmdup.bedpe",sep="")
  pdffile            = paste(outname ,".PET.distances.pdf",sep="")
  npets4dist         = as.numeric(as.character(opt["npets4dist"]))
  
  # build bedpe
  print ("building bedpe")
  if (file.exists(bedpefile)){file.remove(bedpefile)}
  buildBedpe(sam1 =sam1, sam2 = sam2, bedpefile = bedpefile);
  
  # split by chromosome and position
  print ("removing duplicate PETs")
  distancesplit = 10000000
  rmdupresults = removeDups(bedpefile,outname,distancesplit)
  
  # add info to results hash for log file
  resultshash[["duplicate PETs"]] = rmdupresults[1]
  resultshash[["nonduplicate PETs"]] = rmdupresults[2]
  resultshash[["interchromosomal PETs"]] = rmdupresults[3]
  resultshash[["intrachromosomal PETs"]] = rmdupresults[4]
  resultshash[["aligned PETs"]] = rmdupresults[5]
  
  # remove temporary files
  for (f in (5:length(rmdupresults)))
  {
    if (file.exists(rmdupresults[f])==TRUE){file.remove(rmdupresults[f])} 
  }
  
  # plot distance distribution
  plotdistancedistribution(bedpefile=bedpermdupfile,pdffile=pdffile,npets=npets4dist)
  
  
#   # split by chrom and sort bedpe
#   print ("sorting bedpe")
#   if (file.exists(bedpefilesort)){file.remove(bedpefilesort)}

#   
#   # this as the way to sort the whol file at once
#   #external_sort(bedpefile, bedpefilesort)
#   
#   # filter duplicates
#   print ("removing PCR duplicates")
#   if (file.exists(bedpefilesortrmdup)){file.remove(bedpefilesortrmdup)}
#   rmdupresults = removeDupBedpe(bedpefilesort,bedpefilesortrmdup,renamePets=TRUE);
#   resultshash[["duplicate PETs"]] = rmdupresults[1]
#   resultshash[["nonduplicate PETs"]] = rmdupresults[2]
#   resultshash[["interchromosomal PETs"]] = rmdupresults[3]
#   resultshash[["intrachromosomal PETs"]] = rmdupresults[4]
  }, error = function(e) {
    errmsg <- conditionMessage(e)
    cat(paste0("[ERROR] Stage 3 failed: ", errmsg, "\n"), file=stderr())
    errors_occurred <<- c(errors_occurred, paste0("Stage 3: ", errmsg))
    write(paste0("[ERROR] Stage 3: ", errmsg), file=errorlog, append=TRUE)
  })
}

##################################### call peaks #####################################

if (4 %in% opt$stages)
{
  tryCatch({
  checkRequired(opt,c("bedtoolsgenome"))
  
  # gather arguments
  outname         = as.character(opt["outname"])
  MACS_qvalue     = as.character(opt["MACS_qvalue"])
  bedtoolsgenome  = as.character(opt["bedtoolsgenome"])
  peakslop        = as.character(opt["peakslop"])
  peakinput       = as.character(opt["peakinput"])
  MACS_shiftsize  = as.character(opt["MACS_shiftsize"])
  blacklist       = as.character(opt["blacklist"])
  gsize           = as.character(opt["gsize"])
  
  # filenames
  bedpefilesortrmdup = paste(outname ,".rmdup.bedpe",sep="")
  tagAlignfile       = paste(outname,".tagAlign",sep="")
  peaksfile          = paste(outname,"_peaks.narrowPeak",sep="")
  peaksfileslop      = paste(outname,"_peaks.slopPeak",sep="")
  
  if (peakinput != "NULL")
  {
    peaksfile = peakinput
  }
  
  if (peakinput == "NULL")
  {
    print ("building tagAlign file")
    # reverse strands for peak calling
    if (file.exists(tagAlignfile)){file.remove(tagAlignfile)}
    buildTagAlign(bedpefilesortrmdup ,tagAlignfile )
    
    # call peaks 
    print ("calling peaks")
    callpeaks(macs2path=macs2path,tagAlignfile,outname,qvalue=MACS_qvalue,
             bedtoolspath=bedtoolspath,bedtoolsgenome=bedtoolsgenome,
             peakslop=peakslop,MACS_shiftsize,gsize=gsize)
  }
  
  # extend and merge peaks according to peakslop
  print ("extending peaks")
  peakcounts = extendpeaks(peaksfile,peaksfileslop,bedtoolspath=bedtoolspath,
             bedtoolsgenome=bedtoolsgenome,peakslop=peakslop,blacklist=blacklist)
  resultshash[["peaks"]] = peakcounts[1]
  resultshash[["mergedpeaks"]] = peakcounts[2]
  }, error = function(e) {
    errmsg <- conditionMessage(e)
    cat(paste0("[ERROR] Stage 4 failed: ", errmsg, "\n"), file=stderr())
    errors_occurred <<- c(errors_occurred, paste0("Stage 4: ", errmsg))
    write(paste0("[ERROR] Stage 4: ", errmsg), file=errorlog, append=TRUE)
  })
}


##################################### new group/score/filter pairs #####################################

if (5 %in% opt$stages)
{
  tryCatch({
  checkRequired(opt,c("outname","bedtoolsgenome"))
  
  # gather arguments
  outname = as.character(opt["outname"])
  distcutrangemin = as.numeric(as.character(opt["distcutrangemin"]))
  distcutrangemax = as.numeric(as.character(opt["distcutrangemax"]))
  bedtoolsgenome    = as.character(opt["bedtoolsgenome"])
  biascut = as.numeric(as.character(opt["biascut"]))
  maxinteractingdist = as.numeric(as.character(opt["maxinteractingdist"]))
  numofbins = as.numeric(as.character(opt["numofbins"]))
  FDR = as.numeric(as.character(opt["FDR"]))
  minPETS = as.numeric(as.character(opt["minPETS"]))
  chrominclude      = as.character(opt["chrominclude"])
  chromexclude      = as.character(opt["chromexclude"])
  reportallpairs    = as.character(opt["reportallpairs"])
  corrMethod = as.character(opt["corrMethod"])
  MHT    = as.character(opt["MHT"])
  extendreads = as.numeric(opt["extendreads"])
  verboseoutput = as.character(opt["verboseoutput"])

  # filenames
  tagAlignfile       = paste(outname,".tagAlign",sep="")
  tagAlignfileExt    = paste(outname ,".tagAlign.extended.bed",sep="")
  temppeakoverlap    = paste(outname ,".temppeakoverlap.bed",sep="")
  peaksfile          = paste(outname ,"_peaks.narrowPeak",sep="")
  peaksfileslop      = paste(outname ,"_peaks.slopPeak",sep="")
  peaksfileslopdepth = paste(outname ,"_peaks.slopPeak.depth",sep="")
  bedpefilesortrmdup = paste(outname ,".rmdup.bedpe",sep="")
  distancefile       = paste(outname ,".distance",sep="")
  distancecutpdf     = paste(outname ,".distance.pdf",sep="")
  modelspdf          = paste(outname ,".models.pdf",sep="")
  allpairsfile       = paste(outname ,".interactions.all.mango",sep="")
  fdrpairsfile       = paste(outname ,".interactions.fdr.mango",sep="")
  
  # counting reads per peak
  print ("counting reads per peak")
  if (file.exists(tagAlignfileExt) ==TRUE){file.remove(tagAlignfileExt)}
  if (file.exists(temppeakoverlap) ==TRUE){file.remove(temppeakoverlap)}
  DeterminePeakDepths(bedtools=bedtoolspath,bedtoolsgenome=bedtoolsgenome,extendreads=extendreads,tagAlignfile=tagAlignfile,
                  tagAlignfileExt=tagAlignfileExt,peaksfileslop=peaksfileslop,temppeakoverlap=temppeakoverlap)
  if (file.exists(tagAlignfileExt) ==TRUE){file.remove(tagAlignfileExt)}
  if (file.exists(temppeakoverlap) ==TRUE){file.remove(temppeakoverlap)}
  
  # build a file of just distances and same / dif
  print ("determining self-ligation distance")
  makeDistanceFile(bedpefilesortrmdup,distancefile,
                   distcutrangemin,
                   distcutrangemax)
  
  # calculate bias and cutoff
  distancecutoff = calcDistBias(distancefile,distancecutpdf=distancecutpdf,
                                range=c(distcutrangemin,distcutrangemax),
                                biascut= biascut)
  # print distancecutoff
  print (paste("self-ligation cutoff =",distancecutoff))
    
  # group PETs into interactions
  print ("grouping PETs into interactions")
  chromosomes = groupPairs(bedpefilesortrmdup=bedpefilesortrmdup,
                           outname=outname,
                           peaksfile=peaksfileslop,
                           bedtoolspath = bedtoolspath,
                           bedtoolsgenome = bedtoolsgenome,
                           extendreads=extendreads,peaksfileslopdepth=peaksfileslopdepth,
                           verbose=FALSE)
  
  # filter out unwanted chromosomes
  originalchroms = chromosomes
  
  # get chromosomes from bedtools
  bedtoolsgenomeinfo = read.table(bedtoolsgenome,header=FALSE,sep="\t")
  chromosomes = bedtoolsgenomeinfo[,1]
  chromosomes = chromosomes[grep("_",chromosomes,invert=TRUE)]
  if(chrominclude[1] != "NULL")
  {
    chromosomes = unlist(strsplit(chrominclude,split=","))
  }
  
  if (chromexclude[1] !=  "NULL")
  {
    chromosomestpremove = unlist(strsplit(chromexclude,split=","))
    chromosomes = chromosomes[-which(chromosomes %in% chromosomestpremove)] 
  }    
  
  print ("modeling PETs based on peak depth and distance")
  #--------------- Gather IAB data ---------------#

  # Helper to write empty interaction output files (used when data are insufficient for modeling)
  write_empty_interactions <- function() {
    putpairs_empty = data.frame(chrom1=character(),start1=integer(),end1=integer(),
                                chrom2=character(),start2=integer(),end2=integer(),
                                name=character(),peak1=character(),peak2=character(),
                                PETs=integer(),distance=numeric(),
                                P_IAB_distance=numeric(),P_combos_distance=numeric(),
                                P_IAB_depth=numeric(),P_combos_depth=numeric(),
                                p_binom=numeric(),P=numeric(),Q=numeric())
    sig_empty = putpairs_empty
    if (verboseoutput == TRUE)
    {
      if (reportallpairs == TRUE) write.table(x=putpairs_empty,file=allpairsfile,quote=FALSE,sep="\t",row.names=FALSE)
      write.table(x=sig_empty,file=fdrpairsfile,quote=FALSE,sep="\t",row.names=FALSE)
    }
    if (verboseoutput == FALSE)
    {
      if (reportallpairs == TRUE) write.table(x=putpairs_empty[,c("chrom1","start1","end1","chrom2","start2","end2","PETs","Q")],file=allpairsfile,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
      write.table(x=sig_empty[,c("chrom1","start1","end1","chrom2","start2","end2","PETs","Q")],file=fdrpairsfile,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
    }
  }

  # gather all putative interactions
  putpairs = combineputativepairs(chromosomes,outname)

  # guard against NULL return when no chromosome pair files were produced
  if (is.null(putpairs)) {
    msg = "No putative interactions found (combineputativepairs returned no data). Writing empty output files."
    warning(msg)
    stage5_warnings = c(stage5_warnings, msg)
    write_empty_interactions()
    resultshash[["putative interactions"]]    = 0
    resultshash[["significant interactions"]] = 0
  } else {

  # calculate interaction distances
  putpairs$distances = abs( (putpairs[,2] + putpairs[,3] ) / 2 - (putpairs[,5] + putpairs[,6] ) / 2  )
  
  # filter out putative interactions that don't fall into distance range
  putpairs = putpairs[which(putpairs$distances < maxinteractingdist & putpairs$distances > distancecutoff),]
  
  # calculate depths
  putpairs$depths = calcDepths(putpairs[,10:11],type="product")
  
  # handle case where there are no putative pairs after filtering
  if (nrow(putpairs) == 0)
  {
    msg = "No putative interactions remain after distance filtering. Writing empty output files."
    warning(msg)
    stage5_warnings = c(stage5_warnings, msg)
    write_empty_interactions()
    resultshash[["putative interactions"]]    = 0
    resultshash[["significant interactions"]] = 0
  } else {

  # reduce numofbins if there are fewer putative pairs than bins
  if (nrow(putpairs) < numofbins)
  {
    msg = paste0("Number of putative pairs (", nrow(putpairs), ") is less than numofbins (", numofbins,
                 "). Reducing numofbins to ", nrow(putpairs), " to avoid errors.")
    warning(msg)
    stage5_warnings = c(stage5_warnings, msg)
    numofbins = nrow(putpairs)
  }

  # Check whether there are enough unique distance values for smooth.spline (requires >= 4)
  n_unique_dist = length(unique(putpairs$distances[!is.na(putpairs$distances)]))
  if (numofbins < 4 || n_unique_dist < 4)
  {
    msg = paste0("Number of putative pairs (", nrow(putpairs), ") yields too few unique distance values (",
                 n_unique_dist, ") for statistical modeling (minimum 4 required). Writing empty output files.")
    warning(msg)
    stage5_warnings = c(stage5_warnings, msg)
    write_empty_interactions()
    resultshash[["putative interactions"]]    = nrow(putpairs)
    resultshash[["significant interactions"]] = 0
  } else {

  totalcombos = 0
  for (reps in (1:2))
  {
    #--------------- Distance Normalization ---------------#
    
    # determine borders to distance bins
    distanceborders = binmaker(putpairs$distances,binmethod="equalocc",numberbins=numofbins)
    
    # model IAB vs distance
    distance_IAB_model = model_chia(x=putpairs$distances,y=putpairs[,12],borders=distanceborders,yvals=TRUE)
    distance_IAB_model_file   = paste(outname ,".distance_IAB_model.",reps, ".text",sep="")
    write.table(distance_IAB_model,file=distance_IAB_model_file,quote = FALSE, sep = "\t",row.names = FALSE,col.names = TRUE)
    dist_IAB_valid = !is.nan(distance_IAB_model[,1]) & !is.na(distance_IAB_model[,1]) & distance_IAB_model[,1] > 0 & is.finite(distance_IAB_model[,3])
    distance_IAB_spline =   smooth.spline(log10(distance_IAB_model[dist_IAB_valid,1]),distance_IAB_model[dist_IAB_valid,3],spar=.75)
    
    #--------------- Depth Normalization ---------------#
    
    # determine borders to depth bins
    depthborders = binmaker(putpairs$depths,binmethod="equalocc",numberbins=numofbins)
    
    # model IAB vs depth
    depth_IAB_model = model_chia(x=putpairs$depths,y=putpairs[,12],borders=depthborders,yvals=TRUE)
    depth_IAB_model_file   = paste(outname ,".depth_IAB_model.",reps, ".text",sep="")
    write.table(depth_IAB_model,file=depth_IAB_model_file,quote = FALSE, sep = "\t",row.names = FALSE,col.names = TRUE)
    depth_IAB_valid = !is.nan(depth_IAB_model[,1]) & !is.na(depth_IAB_model[,1]) & depth_IAB_model[,1] > 0 & is.finite(depth_IAB_model[,3])
    depth_IAB_spline =   smooth.spline(log10(depth_IAB_model[depth_IAB_valid,1]),depth_IAB_model[depth_IAB_valid,3],spar=.75)
    
    # model Combos vs distance
    nbins_dist  = nrow(distance_IAB_model)
    nbins_depth = nrow(depth_IAB_model)
    meanofx_dist  = rep(0,nbins_dist)
    sumofy_dist   = rep(0,nbins_dist)
    pvals_dist    = rep(0,nbins_dist)
    sumofx_dist   = rep(0,nbins_dist)
    countofx_dist = rep(0,nbins_dist)
    meanofx_depth  = rep(0,nbins_depth)
    sumofy_depth   = rep(0,nbins_depth)
    pvals_depth    = rep(0,nbins_depth)
    sumofx_depth   = rep(0,nbins_depth)
    countofx_depth = rep(0,nbins_depth)

    for (chrom in chromosomes)
    {
      # make combos
      combos = makecombos(chrom,outname,mindist=distancecutoff,maxdist=maxinteractingdist)
      
      # calculate distances
      combos$distance = abs( (combos[,2] + combos[,3] ) / 2 - (combos[,5] + combos[,6] ) / 2  )
      
      # calculate depths
      combos$depths = calcDepths(combos[,7:8],type="product")
      
      # model Combo vs distance
      distance_combo_model_chrom = model_chia(x=combos$distance,y=NA,borders=distanceborders,yvals=FALSE)
      
      sumofy_dist   = sumofy_dist   + distance_combo_model_chrom[,2]
      sumofx_dist   = sumofx_dist   + distance_combo_model_chrom[,4]
      countofx_dist = countofx_dist + distance_combo_model_chrom[,5]
       
      # model Combo vs depth
      depth_combo_model_chrom  = model_chia(x=combos$depths,y=NA,borders=depthborders,yvals=FALSE)
      
      sumofy_depth   = sumofy_depth   + depth_combo_model_chrom[,2]
      sumofx_depth   = sumofx_depth   + depth_combo_model_chrom[,4]
      countofx_depth = countofx_depth + depth_combo_model_chrom[,5]
    }
    
    # combine data from all chromosomes
    depth_combo_model    = cbind(sumofx_depth/countofx_depth, sumofy_depth, sumofy_depth/sum(sumofy_depth))
    distance_combo_model = cbind(sumofx_dist /countofx_dist,  sumofy_dist, sumofy_dist/sum(sumofy_dist))
    
    depth_combo_model_file   = paste(outname ,".depth_combo_model.",reps, ".text",sep="")
    write.table(depth_combo_model,file=depth_combo_model_file,quote = FALSE, sep = "\t",row.names = FALSE,col.names = TRUE)
    
    distance_combo_model_file   = paste(outname ,".distance_combo_model.",reps, ".text",sep="")
    write.table(distance_combo_model,file=distance_combo_model_file,quote = FALSE, sep = "\t",row.names = FALSE,col.names = TRUE)
    
    depth_combo_valid    = !is.nan(depth_combo_model[,1]) & !is.na(depth_combo_model[,1]) & depth_combo_model[,1] > 0 & is.finite(depth_combo_model[,3])
    distance_combo_valid = !is.nan(distance_combo_model[,1]) & !is.na(distance_combo_model[,1]) & distance_combo_model[,1] > 0 & is.finite(distance_combo_model[,3])
    depth_combo_spline    =   smooth.spline(log10(depth_combo_model[depth_combo_valid,1]),depth_combo_model[depth_combo_valid,3],spar=.75)
    distance_combo_spline =   smooth.spline(log10(distance_combo_model[distance_combo_valid,1]),distance_combo_model[distance_combo_valid,3],spar=.75)

    if (reps == 2)
    {
      #--------------- Gather IAB data again ---------------#
      
      # gather all putative interactions
      putpairs = combineputativepairs(chromosomes,outname)
      
      # calculate interaction distances
      putpairs$distances = abs( (putpairs[,2] + putpairs[,3] ) / 2 - (putpairs[,5] + putpairs[,6] ) / 2  )
      
      # filter out putative interactions that don't fall into distance range
      putpairs = putpairs[which(putpairs$distances < maxinteractingdist & putpairs$distances > distancecutoff),]
      
      # calculate depths
      putpairs$depths = calcDepths(putpairs[,10:11],type="product")
    }

    #--------------- Score putative interactions ---------------#
    
#    putpairstemfile       = paste(outname ,".putpairs",sep="")
#    write.table(putpairs,file=putpairstemfile,quote = FALSE, sep = "\t",row.names = FALSE)
    
    # Assing the four probabilities
    putpairs$P_IAB_distance    = predict(distance_IAB_spline, log10(putpairs$distances))$y
    putpairs$P_combos_distance = predict(distance_combo_spline,log10(putpairs$distances))$y
    putpairs$P_IAB_depth       = predict(depth_IAB_spline,log10(putpairs$depths))$y
    putpairs$P_combos_depth    = predict(depth_combo_spline,log10(putpairs$depths))$y
    
#     # fix negative values
#     putpairs$P_IAB_distance[which(putpairs$P_IAB_distance <= 0)] = 
#       min(putpairs$P_IAB_distance[which(putpairs$P_IAB_distance > 0)])
#     putpairs$P_combos_distance[which(putpairs$P_combos_distance <= 0)] = 
#       min(putpairs$P_combos_distance[which(putpairs$P_combos_distance > 0)])
#     putpairs$P_IAB_depth[which(putpairs$P_IAB_depth <= 0)] = 
#       min(putpairs$P_IAB_depth[which(putpairs$P_IAB_depth > 0)])
#     putpairs$P_combos_depth[which(putpairs$P_combos_depth <= 0)] = 
#       min(putpairs$P_combos_depth[which(putpairs$P_combos_depth > 0)])

    # cap values to min and max
    putpairs$P_IAB_distance[which(putpairs$P_IAB_distance <= min(distance_IAB_model[,3]))] = min(distance_IAB_model[,3])
    putpairs$P_IAB_distance[which(putpairs$P_IAB_distance >= max(distance_IAB_model[,3]))] = max(distance_IAB_model[,3]) 

    putpairs$P_combos_distance[which(putpairs$P_combos_distance <= min(distance_combo_model[,3]))] = min(distance_combo_model[,3])
    putpairs$P_combos_distance[which(putpairs$P_combos_distance >= max(distance_combo_model[,3]))] = max(distance_combo_model[,3])

    putpairs$P_IAB_depth[which(putpairs$P_IAB_depth <=  min(depth_IAB_model[,3]))] =  min(depth_IAB_model[,3])  
    putpairs$P_IAB_depth[which(putpairs$P_IAB_depth >=  max(depth_IAB_model[,3]))] =  max(depth_IAB_model[,3])
                                                                                          
    putpairs$P_combos_depth[which(putpairs$P_combos_depth <= min(depth_combo_model[,3]))] =  min(depth_combo_model[,3])   
    putpairs$P_combos_depth[which(putpairs$P_combos_depth >= max(depth_combo_model[,3]))] =  max(depth_combo_model[,3])   

    # calculate the binomial probability
    totalcombos =  sum(sumofy_dist)
    putpairs$p_binom           = (putpairs$P_IAB_distance * putpairs$P_IAB_depth) / 
      (putpairs$P_combos_distance * putpairs$P_combos_depth * totalcombos)
   
    # calculate the total IABs
    totalIAB = sum(distance_IAB_model[,2])

    # calculate the final interaction P values
    putpairs$P = apply(cbind(putpairs$V12,rep(totalIAB,nrow(putpairs)),putpairs$p_binom),1,calcP)

    if (reps == 1)
    {
      putpairs[which(putpairs$P < 1/totalcombos ),]$V12 = 0
    }
  }
  
  #--------------- Correct for multiple hypothesis testing ---------------#

  print ("correcting for multiple hypothesis testing")
  n=nrow(putpairs)
  if (MHT == "all")
  {
    n=totalcombos
  }

  putpairs$Q = p.adjust(putpairs$P,method=corrMethod,n=n)

  #--------------- Organize  data ---------------#
  
  pairnames = paste("pair_",(1:nrow(putpairs)),sep="")
  putpairs = cbind(putpairs[,c(1,2,3,4,5,6)],pairnames,putpairs[,c(10,11,12,14,16,17,18,19,20,21,22)])
  names(putpairs) = c("chrom1","start1","end1","chrom2","start2","end2","name",
                      "peak1","peak2","PETs","distance",
                      "P_IAB_distance","P_combos_distance","P_IAB_depth","P_combos_depth",
                      "p_binom","P","Q")
  
  #--------------- Filter interactions ---------------#

  sig = putpairs[which(putpairs$Q < FDR & putpairs$PETs >= minPETS),]

  #--------------- Write outputs ---------------#
  
  print ("writing output files")
  # write results to output
  if (verboseoutput == TRUE)
  {
    if (reportallpairs == TRUE)
    {
      write.table(x=putpairs,file=allpairsfile,quote = FALSE, sep = "\t",row.names = FALSE)
    }
    write.table(x=sig,file=fdrpairsfile,quote = FALSE, sep = "\t",row.names = FALSE)
  }
  if (verboseoutput == FALSE)
  {
    if (reportallpairs == TRUE)
    {
      write.table(x=putpairs[,c("chrom1","start1","end1","chrom2","start2","end2","PETs","Q")],file=allpairsfile,quote = FALSE, sep = "\t",row.names = FALSE,col.names = FALSE)
    }
    write.table(x=sig[,c("chrom1","start1","end1","chrom2","start2","end2","PETs","Q")],file=fdrpairsfile,quote = FALSE, sep = "\t",row.names = FALSE,col.names = FALSE)
  } 
  resultshash[["putative interactions"]]    = nrow(putpairs)
  resultshash[["significant interactions"]] = nrow(sig)

  #--------------- Make plots ---------------#
  
  print ("plotting results")
  
  # plot models  
  pdf(modelspdf)
  par(mfrow=c(2,2))
  par(mgp=c(3,.3,0))
  plot(log10(distance_IAB_model[,1]),   distance_IAB_model[,3],pch=19,col="dodgerblue2",xlab="distance (bp)",ylab="IAB")
  lines(x=log10(distance_IAB_model[,1]),   predict(distance_IAB_spline,log10(distance_IAB_model[,1]))$y)   
  
  plot(log10(distance_combo_model[,1]), distance_combo_model[,3],  pch=19,col="firebrick2",xlab="distance (bp)",ylab="# combos")  
  lines(x=log10(distance_combo_model[,1]),   predict(distance_combo_spline,log10(distance_combo_model[,1]))$y)  
  
  plot(log10(depth_IAB_model[,1])   ,  depth_IAB_model[,3],    pch=19,col="dodgerblue2",xlab="depth (p1 * p2)",ylab="IAB") 
  lines(x=log10(depth_IAB_model[,1]),   predict(depth_IAB_spline,log10(depth_IAB_model[,1]))$y)
  
  plot(log10(depth_combo_model[,1]),   depth_combo_model[,3],  pch=19,col="firebrick2",xlab="depth (p1 * p2)",ylab="# combos") 
  lines(x=log10(depth_combo_model[,1]),   predict(depth_combo_spline,log10(depth_combo_model[,1]))$y)
  dev.off()
  } # end of sufficient pairs block

  } # end of non-empty putpairs block

  } # end of non-NULL putpairs block

  #--------------- Delete temporary files ---------------#
  
  # clean up extra files
  print ("deleting temporary files")
  if (file.exists(distancefile)) file.remove(distancefile)
  for (chrom in originalchroms)
  {
    peaksizecount  = paste(outname,"." ,chrom, "_peaks.count.slopPeak",sep="")
    pairsbedpe     = paste(outname,"." ,chrom, ".pairs.bedpe",sep="")
    bedpefile      = paste(outname,"." ,chrom, ".bedpe",sep="")
    bedfile        = paste(outname,"." ,chrom, ".bed",sep="")
    overlapfile    = paste(outname,"." ,chrom, ".bedNpeak",sep="")
    if (file.exists(peaksizecount)) file.remove(peaksizecount)
    if (file.exists(pairsbedpe)) file.remove(pairsbedpe)
    if (file.exists(bedpefile)) file.remove(bedpefile)
    if (file.exists(bedfile)) file.remove(bedfile)
    if (file.exists(overlapfile)) file.remove(overlapfile)
  }
  }, error = function(e) {
    errmsg <- conditionMessage(e)
    cat(paste0("[ERROR] Stage 5 failed: ", errmsg, "\n"), file=stderr())
    errors_occurred <<- c(errors_occurred, paste0("Stage 5: ", errmsg))
    write(paste0("[ERROR] Stage 5: ", errmsg), file=errorlog, append=TRUE)
  })
} 

##################################### Make Log file #####################################

print ("writing to log file")

stoptime = paste("Analysis end time:" , as.character(Sys.time()))
write(stoptime,file=logfile,append=TRUE)
write("",file=logfile,append=TRUE)

write("Software Versions:",file=logfile,append=TRUE)
write(paste("mango version",":",Mangoversion),file=logfile,append=TRUE)
write(paste("bedtools version",":",bedtoolsversion),file=logfile,append=TRUE)
write(paste("macs2 version",":",macs2version),file=logfile,append=TRUE)
write(paste("bowtie version",":",bowtieversion),file=logfile,append=TRUE)

write("",file=logfile,append=TRUE)
write("Analyzed by Mango using the following parameters:",file=logfile,append=TRUE)

for (key in names(opt))
{
  write(paste( key, ":",opt[key]),file=logfile,append=TRUE)
}
write("",file=logfile,append=TRUE)
write("With the following results:",file=logfile,append=TRUE)
print (resultshash)
for (key in keys(resultshash))
{
  write(paste( key, ":",resultshash[[key]]),file=logfile,append=TRUE)
}
if (length(stage5_warnings) > 0)
{
  write("",file=logfile,append=TRUE)
  write("Warnings:",file=logfile,append=TRUE)
  for (w in stage5_warnings)
  {
    write(paste("[WARNING]", w),file=logfile,append=TRUE)
  }
}
if (length(errors_occurred) > 0)
{
  write("",file=logfile,append=TRUE)
  write("Errors:",file=logfile,append=TRUE)
  for (err in errors_occurred)
  {
    write(paste("[ERROR]", err),file=logfile,append=TRUE)
  }
}

# Write stats.txt summarizing results and any warnings
statsfile = paste(as.character(opt["outname"]), ".stats.txt", sep="")
statslines = c()
statslines = c(statslines, paste("Mango Analysis Summary -", as.character(Sys.time())))
statslines = c(statslines, "")
statslines = c(statslines, "Results:")
for (key in keys(resultshash))
{
  statslines = c(statslines, paste(key, ":", resultshash[[key]]))
}
if (length(stage5_warnings) > 0)
{
  statslines = c(statslines, "")
  statslines = c(statslines, "Warnings:")
  for (w in stage5_warnings)
  {
    statslines = c(statslines, paste("[WARNING]", w))
  }
}
if (length(errors_occurred) > 0)
{
  statslines = c(statslines, "")
  statslines = c(statslines, "Errors:")
  for (err in errors_occurred)
  {
    statslines = c(statslines, paste("[ERROR]", err))
  }
}
writeLines(statslines, con=statsfile)

print("done")

