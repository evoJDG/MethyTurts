##### Calculate global methylation values for hatchling samples
##### and plot by nest depth
## note that we use the residuals of Number of covered sites ~ coverage to account for correlation between them

# JD Gilbert
# 22/07/22
# adapted from Alice's script R_Complete_Analysis.Rmd

devtools::install_github("cardiomoon/ggiraphExtra")
library(predict3d)

# set paths to relevant directories
dataPath=file.path("/Users/Chimz/Google Drive/Projects/inProgress/1Turtles/01_Methy/Hatchling_Expt/Data/CG_reduced/MethylKit_Objects")

# load united (merged) datasets (object called uniteCovALL.RDS etc.)
load(file = file.path(dataPath, "uniteCovALL.RDS")) # sites with required coverage in all samples
load(file = file.path(dataPath, "uniteCov_50pc.RDS")) # sites with required coverage in 50 percent of samples
load(file = file.path(dataPath, "DMS_Depth.RDS")) # sites with differential methylation between depths
load(file = file.path(dataPath, "DMS_Depth_Mat_Reloc.RDS")) # sites with DM by depths when maternal ID and relocation as covariates

######  NOTE the names of the objects in R do not always match the names of the files (doh!)

## Add metadata on treatments
metadata <- readxl::read_xlsx(file.path(DIR, "metadata.xlsx"))
metadata$Trtmt <- as.factor(metadata$Depth)
metadata <- metadata %>%
  mutate(Trtmt = recode(Trtmt, "shallow" = 0, "deep" = 1 )) %>%
  mutate(Trtmt = as.numeric(Trtmt)) %>%
  rename(SampleID = Sample)

# Calculate number of methylated sites, mean coverage, and residuals of methylated sites by covered sites (to account for coverage bias)

mycalcRMS <- function(myUniteCov, myMetaData){
  percMethMat = methylKit::percMethylation(myUniteCov)
  # create a dataframe with all info
  percMethDF = data.frame(SampleID = colnames(percMethMat),
                          Nbr_methCpG = colSums(percMethMat>=70 & !is.na(percMethMat)), ## number of methylated sites
                          Nbr_coveredCpG = colSums(!is.na(percMethMat)), ## number of sites covered in this sample
                          Nbr_NOTcoveredCpG = colSums(is.na(percMethMat)),## number of sites NOT covered in this sample
                          MeanCoverage = colMeans(methylKit::getData(myUniteCov)[,myUniteCov@coverage.index], na.rm = T), ## coverage.index: vector denoting which columns in the data correspond to coverage values
                          OverallPercentageMethylation = colMeans(methylKit::percMethylation(myUniteCov), na.rm = T))
  
  ## RMS in this sample based on covered sites
  percMethDF$RMS_coveredCpG = percMethDF$Nbr_methCpG / percMethDF$Nbr_coveredCpG
  ## merge with original metadata:
  myMetaData = merge(myMetaData, percMethDF)
  # calculate also RMS global, considering CpG covered or not (to compare)
  
  ## this step doesn't work because we don't have M.Seqs_rawReads - I presume read count?
  # myMetaData$RMS_allCpG_coveredOrNot = myMetaData$Nbr_methCpG / (myMetaData$M.Seqs_rawReads*10e6)
  
  
  # calculate residuals of nbr of methCpG by nbr of covered CpG
  myMetaData$res_Nbr_methCpG_Nbr_coveredCpG = residuals(
    lm(myMetaData$Nbr_methCpG ~ myMetaData$Nbr_coveredCpG))
  ## REORDER myMetaData by sample ID
  myMetaData = myMetaData[order(as.numeric(gsub("S", "", myMetaData$SampleID))),]
  return(myMetaData)
}

# calculate the residual methylation ~ coverage
fullMetadata <- mycalcRMS(myUniteCov = uniteCovALL, 
                          myMetaData = metadata)

# does covered correlate with number methylated?
cor.test(fullMetadata$Nbr_coveredCpG,
         fullMetadata$Nbr_methCpG, method = "spearman")

# plot number covered vs number methylated
ggplot(fullMetadata, aes(x=Nbr_coveredCpG, y=Nbr_methCpG))+
  geom_smooth(method = "lm", col="black")+
  geom_point(aes(col=Depth), size = 3)+ scale_color_manual(values = c("grey", "red")) +
  theme_bw() + ggtitle(label = "All CpG sites")

## Check after RMS correction for coverage bias: CORRECTED (p-value = ?? )
cor.test(fullMetadata$Nbr_coveredCpG,
         fullMetadata$RMS_coveredCpG, method = "spearman")
ggplot(fullMetadata, aes(x=Nbr_coveredCpG, y=RMS_coveredCpG))+
  geom_smooth(method = "lm", col="black")+
  geom_point(aes(col=Depth), size = 3)+ scale_color_manual(values = c("grey", "red")) +
  theme_bw() + ggtitle(label = "All CpG sites")

## and with residuals: COMPLETELY CORRECTED p-value = ????
cor.test(fullMetadata$Nbr_coveredCpG,
         fullMetadata$res_Nbr_methCpG_Nbr_coveredCpG, method = "spearman")

ggplot(fullMetadata, aes(x=Nbr_coveredCpG, y=res_Nbr_methCpG_Nbr_coveredCpG))+
  geom_smooth(method = "lm", col="black")+
  geom_point(aes(col=Depth), size = 3)+ scale_color_manual(values = c("grey", "red")) +
  theme_bw() + ggtitle(label = "All CpG sites")

# finally check whether global methylation differs between depths
mod = lm(res_Nbr_methCpG_Nbr_coveredCpG ~ Depth, data = fullMetadata)
summary((mod)) # sex is significant p = 0.000157 ***
anova(mod)

mycols = c("blue","red"); myshape = c(21,22)

ggplot(data = fullMetadata, aes(y = res_Nbr_methCpG_Nbr_coveredCpG, x = Depth, col = Depth)) +
  geom_boxplot() +
  geom_point( alpha = .6) +
  scale_color_manual(name = "Depth",
                     values = mycols,
                     labels = c("Deep", "Shallow"))+
  scale_fill_manual(name = "Depth",
                    values = mycols,
                    labels = c("Deep", "Shallow")) +
  labs( title = "Global Methylation by Nest Depth", y = "Global Methylation\n(Residuals From Model\nNumber Methylated ~ Coverage)") +
  theme_bw() + 
  theme( panel.grid = element_blank(),
         plot.title = element_text(hjust = 0.5))
