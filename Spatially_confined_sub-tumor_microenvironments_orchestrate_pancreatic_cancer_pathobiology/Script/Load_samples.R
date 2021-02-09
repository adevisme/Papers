library(openxlsx)
library(limma)
library(org.Hs.eg.db)

scaledAnnoHeatmap <- function(mat, ann.row, ann.col, outFile, myColor.ann, clu.col=FALSE, show.col.names=FALSE, show.row.names=FALSE)
{
  require(pheatmap)
  paletteLength <- 25
  myMax <- ceiling(max(abs(mat)))
  
  myBreaks <- seq(-myMax , myMax, length.out=paletteLength)
  myBreaks <- myBreaks[myBreaks != 0]
  myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength-2)
 
  if(nrow(mat)<2) doCluster <- FALSE
  else(doCluster <- TRUE)
  
  pheatmap(mat, color = myColor, breaks = myBreaks, filename = outFile,
           annotation_row = ann.row, annotation_col = ann.col,
           annotation_colors = myColor.ann,
           cluster_cols = clu.col, cluster_rows = doCluster, show_rownames = show.row.names,
           col.names = show.col.names,
           #cellwidth = 12, cellheight = 12
           #height = 9, width = 10,
           border_color = NA
  )
  
}		

ensembl2entrez <- function(ensembl)
{
  entrez <- mget(as.character(ensembl), org.Hs.egENSEMBL2EG, ifnotfound=NA)
  entrez <- lapply(entrez, function(i) return(i[1]))
  return(unlist(entrez))
}

entrez2symbol <- function(entrez)
{
  symbol <- mget(as.character(entrez), org.Hs.egSYMBOL, ifnotfound=NA)
  symbol <- unlist(lapply(symbol, function(i) return(i[1])))
  return(symbol)
}

entrez2name <- function(entrez)
{
  gname <- mget(as.character(entrez), org.Hs.egGENENAME, ifnotfound=NA)
  gname <- unlist(lapply(gname, function(i) return(i[1])))
  return(gname)
}


PathFile = "/home/adevisme/Toronto/Stroma"
setwd(PathFile)
load("Matrix_Info_All_Samples.RData") #MatrixInfo
MatrixInfo <- MatrixInfo[-grep("PCSI_0145", MatrixInfo[,"Sample_ID"]), ]

Paired <- rep(1, dim(MatrixInfo)[1])
MatrixInfo <- cbind(MatrixInfo, Paired)

j=1
for(i in 2:dim(MatrixInfo)[1]){
	if(MatrixInfo[i,1] == MatrixInfo[i-1,1]){
		MatrixInfo[i,6] <- j
		MatrixInfo[i-1,6] <- j
		j = j+1
	}
}

Pathannotation = "/home/adevisme/Toronto/Stroma/doc"
Stroma_anno <- read.xlsx(paste(Pathannotation, "Barbara_PDAC-LCM-proteomes_annotations_2-way_3-way.xlsx", sep="/"), rowNames=TRUE)

Stroma_pos <- match(MatrixInfo[,1], rownames(Stroma_anno))
Stroma_pos <- Stroma_pos[!is.na(Stroma_pos)]

Stroma <- Stroma_anno[Stroma_pos,"LCM.stroma-type_confirmed.3-way"]
Stroma <- c(Stroma, rep("Unknown", (dim(MatrixInfo)[1]-length(Stroma))))

MatrixInfo <- cbind(MatrixInfo, Stroma)
### Removing Intermediate Samples from Matrix Info
#MatrixInfo <- MatrixInfo[-grep("intermediate", MatrixInfo[,"Stroma"]),]



setwd(file.path("/home/adevisme/Toronto/Stroma/Proteome_samples/doc"))
Proteome_Samples <- read.xlsx("RNA_protein_annotation.xlsx", sheet = 1, rowNames = TRUE)
MatrixInfo_Proteome <- MatrixInfo[match(rownames(Proteome_Samples), MatrixInfo[,"Name"]),]
MatrixInfo_Proteome <- MatrixInfo_Proteome[!is.na(MatrixInfo_Proteome[,1]),]
MatrixInfo_Proteome[,"Paired"] <- paste0("pair_", MatrixInfo_Proteome[,"Paired"])

setwd(file.path("/home/adevisme/Toronto/Stroma/Proteome_samples/count/"))
expr.Proteome <- read.table("count.All_Proteome_Samples_WO_PCSI_0145.txt")
expr.Proteome.Tumor <- read.table("count.Proteome_Tumor_Samples_WO_PCSI_0145.txt")
expr.Proteome.Stroma <- read.table("count.Proteome_Stroma_Samples_WO_PCSI_0145.txt")
expr.Proteome.CPM <- read.table("logCPM.all_Proteome_Samples_cpm+1_WO_PCSI_0145.txt")
expr.Proteome.Tumor.CPM <- read.table("logCPM.Proteome_Tumor_Samples_cpm+1_WO_PCSI_0145.txt")
expr.Proteome.Stroma.CPM <- read.table("logCPM.Proteome_Stroma_Samples_cpm+1_WO_PCSI_0145.txt")

