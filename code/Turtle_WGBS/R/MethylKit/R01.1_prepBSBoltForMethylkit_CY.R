## Data preparation for methylKit
## A. Balard
## 3rd of August 2021

# Adapted James Gilbert 21st July 2022
# Adapted Charley Yen 25th July 2022

# To be run on RStudio on Apocrita OnDemand
# NB: need enough memory to load all files, otherwise it crashes
# Go for around 24 Gb


###### Prep environment ######

# Set directory paths
DIR <- file.path("/data/SBCS-EizaguirreLab/Turtle_WGBS/06_MethylKit/Hatchery_CG_Meth_Calls/BSBolt_mincov1")
DATADIR <- file.path("/data/SBCS-EizaguirreLab/Turtle_WGBS/05_Methylation_Calling/meth_calls/mincov1")
OUTDIR <- file.path(DIR, "formatCG4methylKit")

# Create OUTDIR if it doesn't exist already
# if (!dir.exists(OUTDIR)) { dir.create(OUTDIR)}


###### Import files ######

# Make a list of files to analyse
temp = list.files(path=DATADIR,
                  pattern = ".CGmap.gz",
                  full.names = T)

# Check correct files are included in list
length(temp) # check number of files: should be 58 (hatchlings from no repeated nest -> removed 119-5 and 120-3 for now)
head(temp) # check 1st few files on the list

# Import files into a list
# NB. Takes a while
myfiles = lapply(temp, read.csv, sep="\t", header=F)

# Change the name of the files in myfiles -> just include sample ID
names(myfiles) = lapply(temp, function(x) gsub(pattern = paste0("(", DATADIR, "/)(.*)(.CG.map.*)"), replacement = "\\2", x))

# Check
names(myfiles)
class(myfiles)

# Add column names
new.names=c("chrom", "nucleotide", "position", "context", "sub-context", "methylation_value", "methylated_bases", "all_bases")
myfiles=lapply(myfiles, setNames, new.names)

# Check
head(myfiles)


###### Transform BSBolt output format into MethylKit input format ######

# Create function myrenameFUN() that performs formatting changes
myrenameFUN <- function(BSBDF){
  MKDF=data.frame(chrBase=paste(BSBDF$chrom,BSBDF$position, sep = "."),
                  chr=BSBDF$chrom,
                  base=BSBDF$position,
                  strand=ifelse(BSBDF$nucleotide=="C", yes = "F", no = "R"),
                  coverage=BSBDF$all_bases,
                  freqC=round(BSBDF$methylation_value*100, 2),
                  freqT=round((1-BSBDF$methylation_value)*100,2))
  return(MKDF)
}

# Apply function to every methylation file in list
# NB. Takes a while
myfilesMK=lapply(myfiles, myrenameFUN)


###### Save output files ######

length(myfilesMK)

# Output the transformed files as separate files into OUTDIR
# NB. Takes a while

for(i in 1:length(myfilesMK)){
  write.table(myfilesMK[[i]],
              file=file.path(OUTDIR, paste0(names(myfilesMK)[[i]],".CG4methylkit.txt")),
              quote=FALSE, row.names=FALSE, col.names=TRUE, sep= "\t")
}


# Print minimum all_bases for one file
# print(min(myfilesMK$all_bases))


# Minimum all_bases value for each (loop)
for ( i in names(myfilesMK) ) {
  print(min(myfilesMK[[i]]$all_bases))
}
