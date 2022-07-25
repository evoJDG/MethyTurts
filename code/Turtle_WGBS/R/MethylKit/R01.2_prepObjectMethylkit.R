## MethylKit object preparation
## A. Balard
## 25th of August 2021

## Adapted James Gilbert
## 21 July 2022

## NB: change R01.4 loadMethyldata each time this script is amended

## R/4.0.2 to run
library(methylKit)
library(readxl)
library(plyr) # for join (keep row order)
library(dplyr)
library(dendextend) # The package dendextend contains many functions for changing the appearance of a dendrogram and for comparing dendrograms.
library(ggplot2)

## Sources:
## https://www.bioconductor.org/packages/devel/bioc/vignettes/methylKit/inst/doc/methylKit.html
## https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/methylationSeq/Seq_Tutorial.html#load-datasets
## /data/archive/archive-SBCS-EizaguirreLab/RRBS/StickPara_Broject_archive/08CompGenomes_mCextr/04RAnalyses_Methylome/methylbam_peichel...

## load custom functions
# this needs to be in Project dir or full file path will be needed
source("customRfunctions.R")

##### Load prepared dataset (in APOCRITA) #####

DIR <- file.path("/Users/Chimz/Google Drive/Projects/inProgress/1Turtles/01_Methy/Hatchling_Expt/Data/CG_reduced/")
dataPath=file.path(DIR, "formatCG4methylKit")
OUTDIR <- file.path(DIR, "MethylKit_Objects")

temp = list.files(path=dataPath,
                  pattern = ".CG4methylkit.txt",
                  full.names = T)

# exclude the mothers, which have the prefix "SLL"
# grep to print the list elements that contain SLL
# then use the [] subsetting and - to remove those rows
temp <- temp[ -(grep("SLL", temp))]
    

## Add metadata on treatments
metadata <- readxl::read_xlsx(file.path(DIR, "metadata.xlsx"))

# MethylKit only allows numeric treatments. Make a new var - trt_NUM - where shallow = 1 and deep = 2.
metadata <- mutate(metadata,
                   trt_NUM = case_when(
                      Depth == "shallow" ~ 1, 
                      Depth == "deep" ~ 2 
                  )) %>%
  mutate(metadata,
         Mat_ID = as.numeric(substr(metadata$Mama, 4, 6)))


### Make methylkit object
myobj=methylKit::methRead(as.list(temp),
                          mincov=10,
                          sample.id=as.list(metadata$Sample),
                          assembly="Charley's_Loggerhead_Assembly",
                          treatment=metadata$trt_NUM,
                          context="CpG")

###############################
## Filtering and normalising ##
###############################

## Filtering based on coverage:
# It might be useful to filter samples based on coverage. Particularly, if our samples are suffering from PCR bias it would be useful to discard bases with very high read coverage. Furthermore, we would also like to discard bases that have low read coverage, a high enough read coverage will increase the power of the statistical tests. The code below filters a methylRawList and discards bases that have coverage below 10X and also discards the bases that have more than 99.9th percentile of coverage in each sample.
print("Filter")
filtered.myobj=filterByCoverage(myobj, lo.count=10, lo.perc=NULL,
                                hi.count=NULL, hi.perc=99.9)

## normalise the coverage
print("Normalise")
normFil.myobj=normalizeCoverage(filtered.myobj)

#####################
## MERGING SAMPLES ##
#####################
##In order to do further analysis, we will need to get the bases covered in all samples.
# The following function will merge all samples to one object for base-pair locations that are covered in all samples.
# The unite() function will return a methylBase object which will be our main object for all comparative analysis.
# The methylBase object contains methylation information for regions/bases that are covered in all samples.

print("Add CpG present in ALL individuals")
uniteCovALL= methylKit::unite(normFil.myobj, mc.cores=1) # note the number of cores here
uniteCovALL=as(uniteCovALL,"methylBase")

# plot methylation histogram using all normalised, filtered data
getMethylationStats(normFil.myobj[[2]],plot=TRUE,both.strands=FALSE)

## Keep methylated CpG sites observed in at least 50% individuals after filtering and normalising: 
uniteCov_50pc = methylKit::unite(normFil.myobj,
                     min.per.group=as.integer(length(metadata$Sample)/2), mc.cores=1)# try with 8 cores
uniteCov_50pc = as(uniteCov_50pc,"methylBase")


## And save in RDS file for later analyses
save(uniteCovALL, file = file.path(OUTDIR, "uniteCovALL.RDS"))
save(uniteCov_50pc, file = file.path(OUTDIR, "uniteCov_50pc.RDS"))


## Identify differentially methylated sites by maternal_ID and depth
# starting with depth
myDiff_Depth <- calculateDiffMeth(uniteCovALL, metadata) ###### mc.cores can be used for multiple cores

# get all differentially methylated bases (at least 5% difference between depths) ######## difference could be increased - 15 seems more common?
myDiff_Depth_15=getMethylDiff(myDiff_Depth, difference=5, qvalue=0.01)

# define covariates (maternal ID and relocation)
covariates = data.frame(Mat_ID = metadata$Mama, Relocation = metadata$Relocation)

# recalculate differential methylation including covariates
myDiff_Depth_Mat_Reloc <-calculateDiffMeth(uniteCovALL,
                                covariates=covariates) ###### mc.cores can be used for multiple cores

# get all differentially methylated bases (at least 5% difference between depths) ######## difference could be increased - 15 seems more common?
myDiff_Depth_Mat_Reloc_15=getMethylDiff(myDiff_Depth_Mat_Reloc, difference=5, qvalue=0.01)

# save differentially methylated sites comparing 1) depth, then 2) another with maternal ID and relocation as covariates
save(myDiff_Depth_15, file = file.path(OUTDIR, "DMS_Depth.RDS"))
save(myDiff_Depth_Mat_Reloc_15, file = file.path(OUTDIR, "DMS_Depth_Mat_Reloc.RDS"))


##################### Previous tests with ALL numbers of fish 1 to 12:
## we kept for downstream analyses all CpG sites present in at least 1 to 12 individuals per group, or in all individuals:
# print("Unite and store in a list")
# mylist_uniteCov=list()
# for (i in 6:12L){ # done for 1 to 5, then 6 to 12
#     uniteCov=unite(normFil.myobj, min.per.group=i, mc.cores=8)# try with 8 cores
#     uniteCov=as(uniteCov,"methylBase")
#     name=paste0("uniteCov_", as.character(i))
#     mylist_uniteCov[[name]]=uniteCov
# }

## Add CpG present in ALL individuals
# uniteCov=unite(normFil.myobj, mc.cores=8)
# uniteCov=as(uniteCov,"methylBase")
# mylist_uniteCov[["uniteCov_ALL"]]=uniteCov

# CpGALL=length(mylist_uniteCov$uniteCov_ALL$coverage1) # 47238

# Idea: plot number of retained CpG site by nbr of individuals sharing these CpG sites. Preparing file for that:
# print("Make DF")
# CpGDF1_5=data.frame(lapply(mylist_uniteCov, function(x) length(x$coverage1)))
# CpGDF1_5=t(CpGDF1_5)
# CpGDF1_5=data.frame(NbrIndMin=as.numeric(gsub("uniteCov_", "", row.names(CpGDF1_5))),
#                     NbrCpG=CpGDF1_5[,1])
# 
# CpGDF6_12=data.frame(lapply(mylist_uniteCov, function(x) length(x$coverage1)))
# CpGDF6_12=t(CpGDF6_12)
# CpGDF6_12=data.frame(NbrIndMin=as.numeric(gsub("uniteCov_", "", row.names(CpGDF6_12))),
#                      NbrCpG=CpGDF6_12[,1])
# 
# CpGDF=rbind(CpGDF1_5, CpGDF6_12)

# print("Save object for plotting")
# save(CpGDF, file="/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/plots/CpGDF.RDS")

# print("Save plot")
# plotCpGshared <- ggplot(CpGDF, aes(x=NbrIndMin, y=NbrCpG))+
#     geom_smooth(se = F, col = "red")+
#     geom_smooth(method = "lm", se = F, col = "black") +
#     geom_point() +
#     scale_x_continuous("Number of individual fish per treatment group sharing the same methylated CpG sites",
#                        labels = as.character(CpGDF$NbrIndMin), breaks = CpGDF$NbrIndMin)+
#     scale_y_continuous("Number of shared methylated CpG sites") +
#     theme_bw() +
#     geom_hline(yintercept=CpGALL)
# 
# plotCpGshared
# 
# pdf(file="/data/SBCS-EizaguirreLab/Alice/StickParaBroOff/Data/05MethylKit/plots/plotCpGshared.pdf")
# plotCpGshared
# dev.off()
