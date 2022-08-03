## Data preparation for methylKit
## A. Balard
## 3rd of August 2021

# Adapted James Gilbert 21 July 2022

# Set directory paths
DIR <- file.path("/Users/Chimz/Google Drive/Projects/inProgress/1Turtles/01_Methy/Hatchling_Expt/Data/CG_reduced")
DATADIR <- file.path(DIR, "Raw")
OUTDIR <- file.path(DIR, "formatCG4methylKit")

# create outdir
if (!dir.exists(OUTDIR)) { dir.create(OUTDIR)}

# make a list of files to analyse
temp = list.files(path=DATADIR,
                  pattern = "CG.map.gz",
                  full.names = T)

length(temp) # check number of files

head(temp)

## import files into a list (10 minutes)
myfiles = lapply(temp, read.csv, sep="\t", header=F)

## name the list myfiles
names(myfiles) = lapply(temp, function(x) gsub(pattern = paste0("(", DATADIR, "/)(.*)(.1000n.*)"), replacement = "\\2", x))

names(myfiles)

class(myfiles)

## add column names
new.names=c("chrom", "nucleotide", "position", "context", "sub-context", "methylation_value", "methylated_bases", "all_bases")
myfiles=lapply(myfiles, setNames, new.names)

## Transform BSBolt output format into MethylKit input format:
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

## LONG: make compatible DF
myfilesMK=lapply(myfiles, myrenameFUN)       

## Output the transformed files
length(myfilesMK)

for(i in 1:length(myfilesMK)){
    write.table(myfilesMK[[i]],
                file=file.path(OUTDIR, paste0(names(myfilesMK)[[i]],".CG4methylkit.txt")),
                quote=FALSE, row.names=FALSE, col.names=TRUE, sep= "\t")
}

### Files moved after to:
# OUTDIR
