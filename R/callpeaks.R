
# Define a function that calls peaks using macs2
callpeaks <- function(macs2path,tagAlignfile,peaksfile,qvalue=0.05,
                      bedtoolspath,bedtoolsgenome,peakslop=0,MACS_shiftsize="NULL",gsize="hs")
{
  # call peaks
  shiftsize = ""
  if (MACS_shiftsize != "NULL")
  {
    shiftsize =  paste("--nomodel --extsize",MACS_shiftsize,sep=" ")
  }
  if (MACS_shiftsize == "NULL")
  {
    shiftsize = ""
  }

  command = paste(macs2path," callpeak -t ",tagAlignfile, shiftsize ," -g ",gsize,  " -f BED -n ",peaksfile," -q ",qvalue,sep=" ")
  print (command)
  exitcode = system(command)
  if (exitcode != 0)
  {
    stop(paste("macs2 callpeak failed with exit code", exitcode,
               "- check that the tagAlign file is valid and macs2 is installed correctly"))
  }
  
  # now shorten peak names
  narrowpeakfile <- paste(peaksfile,"_peaks.narrowPeak",sep="")
  peaks <- tryCatch(
    read.table(narrowpeakfile, header=FALSE, sep="\t"),
    error = function(e) {
      if (grepl("no lines available", conditionMessage(e), fixed=TRUE)) {
        warning("MACS2 produced no peaks (empty narrowPeak file). Skipping peak name update.")
        return(NULL)
      }
      stop(e)
    }
  )
  if (!is.null(peaks) && nrow(peaks) > 0) {
    peaks[,4] = paste("peak_",1:nrow(peaks),sep="")
    write.table(peaks,file=narrowpeakfile,quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE,append=FALSE)
  }
  
  # remove unneccesary files 
  listofsuffixes = c("peaks.xls","summits.bed","model.r","","")
  for (suf in listofsuffixes)
  {
    fname = paste(peaksfile,suf,sep="_")
    if (file.exists(fname) == TRUE)
    { 
        file.remove(fname)
    }
  }
}