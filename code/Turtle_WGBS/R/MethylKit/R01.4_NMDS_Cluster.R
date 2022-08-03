library(vegan)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(MASS)
library(plyr)
library(methylKit)
library(dplyr)
library(dendextend)
library(qualpalr)
# Prepare dendrogram and manhattan plot for initial subsample of methylation data

DIR <- file.path("/Users/Chimz/Google Drive/Projects/inProgress/1Turtles/01_Methy/Hatchling_Expt/Data/CG_reduced/")
dataPath=file.path(DIR, "MethylKit_Objects")
figurePath=file.path(DIR, "Figures")

# load united (merged) datasets (object called uniteCovALL.RDS etc.)
load(file = file.path(dataPath, "uniteCovALL.RDS")) # sites with required coverage in all samples
load(file = file.path(dataPath, "uniteCov_50pc.RDS")) # sites with required coverage in 50 percent of samples
load(file = file.path(dataPath, "DMS_Depth.RDS")) # sites with differential methylation between depths
load(file = file.path(dataPath, "DMS_Depth_Mat_Reloc.RDS")) # sites with DM by depths when maternal ID and relocation as covariates

# print number of rows per object
######  NOTE the names of the objects in R do not always match the names of the files (doh!)
nrow(uniteCovALL)
nrow(uniteCov_50pc)
nrow(myDiff_Depth_15)
nrow(myDiff_Depth_Mat_Reloc_15)

# select differentially methylated rows using rownames from df with DM sites.
# when filtering for sites with DM, the rownames are preserved so they can be used to filter the original dataset
uniteCovDM <- methylKit::select(uniteCovALL, rownames(myDiff_Depth_15))

## Add metadata on treatments
metadata <- readxl::read_xlsx(file.path(DIR, "metadata.xlsx"))
metadata$Trtmt <- as.factor(metadata$Depth)
metadata <- metadata %>%
  mutate(Trtmt = recode(Trtmt, "shallow" = 0, "deep" = 1 )) %>%
  mutate(Trtmt = as.numeric(Trtmt)) %>%
  rename(ID = Sample)

########### NMDS
makePercentMetMat <- function(dataset){
  # creates a matrix containing percent methylation values
  perc.meth=percMethylation(dataset)
  # KOSTAS MBE: "Methylated sites and regions with low variation
  # and a standard deviation below 0.3, that is, noninformative
  # sites across individuals, were excluded from the cluster analyses"
  SD=apply(perc.meth,1, sd, na.rm = TRUE)
  if (length(which(SD<0.3)) >0 ) {
    perc.meth <- perc.meth[-which(SD<0.3),]  
  }
  x=t(perc.meth)
  return(x)
}

## make percent methylation matrix
perc_methy <- makePercentMetMat(uniteCovDM)

#Create NMDS based on bray-curtis distances - metaMDS finds the
NMDS <- metaMDS(comm = perc_methy, distance = "bray", maxit=1000, k = 6)

#extract plotting coordinates
MDS1 = NMDS$points[,1] ; MDS2 = NMDS$points[,2] ; MDS3 = NMDS$points[,3]

NMDS_dt = data.table::data.table(MDS1 = MDS1, MDS2 = MDS2, MDS3 = MDS3,
                                 ID = metadata$ID,
                                 MAT=as.factor(metadata$Mama), 
                                 DEPTH=as.factor(metadata$Trtmt))

dima=1; dimb=2; dimc=3
mycols = c("blue","red"); myshape = c(21,22)

myvar <- "DEPTH"
hulls <- NMDS_dt[, .SD[chull(get(paste0("MDS", dima)), get(paste0("MDS", dimb)))], by = get(myvar) ]

myNMDSplot <- ggplot(NMDS_dt, 
                     aes_string(x=paste0("MDS",dima), y=paste0("MDS",dimb))) +
  geom_polygon(data = hulls, aes_string(fill=myvar), alpha=0.3) +
  scale_color_manual(name = "Depth",
                     values = mycols,
                     labels = c("Deep", "Shallow"))+
  scale_fill_manual(name = "Depth",
                    values = mycols,
                    labels = c("Deep", "Shallow"))+
  scale_shape_manual(name = "Depth",
                     values = myshape,
                     labels = c("Deep", "Shallow")) +
  geom_point(aes_string(fill=myvar, shape=myvar), size = 3, alpha = .6) +
  geom_text_repel(aes(label=gsub("SLL", "", NMDS_dt$MAT)),
                   size = 3,
                   max.overlaps = Inf,
                   box.padding = 0.4)+
  theme_bw() +
  theme(legend.title=element_blank(), legend.position = "top",
        panel.grid = element_blank(),
        text = element_text(size = 18))

myNMDSplot

########### Clustering
## To add color bars:
# Generate color palette
pal = qualpalr::qualpal(55, colorspace=list(h=c(0,360), s=c(0.3,1), l=c(0.2,0.8)))
depth <- factor(metadata$Depth)
col_depth <- c("royalblue1", "indianred1")[1:length(unique(metadata$Depth))][depth]

mat <- factor(metadata$Mama)
col_mat <- pal$hex[ sample(55, 10, replace = F) ][mat]

mydendro <- clusterSamples(uniteCovDM, dist = "correlation", method = "ward", plot = T)
dend = as.dendrogram(mydendro)

## and plot
dend %>% plot(main=paste("Nest Depth and Maternal ID\nMethylation Clustering\n", 
                         "Distance method: correlation; Clustering method: ward.D"),
              ylab = "Height")
colored_bars(cbind(col_depth, col_mat), dend, y_shift = -0.6,
             rowLabels = c("Depth", "Maternal ID"))
# dend %>% rect.dendrogram(k=10, border = 8, lty = 5, lwd = 2)

########### Statistical Analysis
makeDatadistFUN <- function(dataset){
  x=makePercentMetMat(dataset)
  # creates a distance matrix. Method: Bray-Curtis, package vegan
  data.dist = as.matrix((vegdist(x, "bray", upper = FALSE))) 
}

data.dist = makeDatadistFUN(uniteCovDM)

# you can run the permanova to see whether depth is associated with the distance between samples
adonis2(data.dist ~ Depth,
        data = metadata)

# but you may also want to include Maternal ID and relocation in the model
# which you set in the permutations argument ( I haven't read up on this step.... )
perm <- how(nperm = 1000) # 1000 permutations
# define the permutation structure considering maternal ID and relocation
setBlocks(perm) <- with(metadata, Mat_ID, Relocation)
adonis2(data.dist ~ Depth,
        data = metadata, permutations = perm)

# or you could just include them as covariates
# need to look up the difference between setting as permutations and covariates
adonis2(data.dist ~ Depth + Mat_ID + Relocation,
        data = metadata)


