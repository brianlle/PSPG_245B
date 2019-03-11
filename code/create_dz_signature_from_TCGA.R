#this code is to create a disease signature from TCGA data
#samr and rankprod are used (or maybe bayesian); SAMR runs faster than rankprod, but usually gives a less number of genes.

set.seed(97)

#packages
library(GEOquery)
library(Biobase)
library(preprocessCore)

library(gplots)
library(pheatmap)

library("org.Hs.eg.db")

######################
#functions
#if a list of GeneIDs generated, take the first one
take_first_match <- function(geneID){
  if(length(geneID) > 1)
    geneID[1]
  else
    geneID
}

######################

expr_matrix <- read.csv(disease_data_filepath, sep = "\t")
samples <- colnames(expr_matrix)
sample_types <- gsub(".*\\.", "", samples) #extract sample types encoded in last 2 digits

#define cases (tumors) and controls (normal samples)
case <- samples[sample_types == '01'] #primary tumors are '01'
control <- samples[sample_types == '11'] #solid tissue normal are '11'

#narrow down cases by phenotype
#phenotypes <- read.csv("data/BRCA_clinicalMatrix", sep = "\t")
#phenotypes$sampleID <- gsub("-", ".", phenotypes$sampleID) #keep sampleID names consistent; R doesn't like - in column names
#lumA_phenotypes <- phenotypes[phenotypes$PAM50Call_RNAseq == 'LumA', c('sampleID')]
#case <- case[case %in% lumA_phenotypes]

genes <- as.character(expr_matrix$sample)

egm <- data.frame(Symbol = genes, stringsAsFactors = FALSE)
egm$GeneID <- with(egm, mget(x = Symbol, envir = org.Hs.egALIAS2EG, ifnotfound = NA)) #get Entrez gene IDs

#take first Entrez gene ID
egm$GeneID <- sapply(egm$GeneID, take_first_match)

##################################

#create comparison
comparison_frame = subset(expr_matrix, select = c(control,case))
rownames(comparison_frame) = expr_matrix$sample
sample_class=c(rep(0, length(control)), rep(1,length(case)))

#remove probes if half of values are missed
complete_probes <- apply(comparison_frame,1,function(row) {
  ctl_vals <- row[sample_class==0]
  dz_vals <- row[sample_class==1]
  missing_ctl <- length(ctl_vals[is.na(ctl_vals)]) / length(ctl_vals)
  missing_dz <- length(dz_vals[is.na(dz_vals)]) / length(dz_vals)
  ifelse(((missing_ctl > 0.5) || (missing_dz > 0.5)),FALSE,TRUE)
})
comparison_frame <- comparison_frame[complete_probes,]

# deal negative values
#colSummarizeMedian(as.matrix(comparison_frame))
for (i in 1:ncol(comparison_frame)){
  comparison_frame[,i]=comparison_frame[,i]+abs(min(comparison_frame[,i],na.rm=TRUE))+1 #transform negatives to positives, as quantile normailzation will be applied, it would affect 
}

# deal log
#dataset is already logged, as log2(x+1)
values_logged <- TRUE


#normalization
#boxplot(comparison_frame)
comparison_frame1 <- normalize.quantiles(as.matrix(comparison_frame))
row.names(comparison_frame1) <- row.names(comparison_frame)
colnames(comparison_frame1) <- colnames(comparison_frame)
comparison_frame <- comparison_frame1

#collapse probe id using entrez id
comparison_frame <- merge(egm[, c("Symbol", "GeneID")], data.frame(Symbol = rownames(comparison_frame), comparison_frame), by = "Symbol")
#remove probe column
comparison_frame <- comparison_frame[, -1]
#keep only non-NA geneIDs
comparison_frame <- comparison_frame[!is.na(comparison_frame$GeneID),]
comparison_frame <- aggregate(. ~ GeneID, comparison_frame, mean)
rownames(comparison_frame) <- comparison_frame$GeneID
comparison_frame <- comparison_frame[, -1]

if (method_id==2) {
  method <- 'RANKPROD_SINGLE'
} else if (method_id==3){
  method <- "SIGGENES_SAMR"
} else if (method_id==4){
  method <- "SAMR"
}

#probe names
genenames<-rownames(comparison_frame)

if (method == 'RANKPROD_SINGLE') {
  library(RankProd)
  
  # Evaluate the impact of na.rm = TRUE. Alex M seems to think it's OK
  RP_result <- RP(comparison_frame, sample_class,gene.names=genenames, num.perm = 100, logged = T, na.rm = TRUE, plot = FALSE, rand = 123)
  # Leave logged=FALSE because topGene() converts the fold-change incorrectly!
  siggenes <- topGene(RP_result,cutoff=sig_cutoff,method="pfp",logged=FALSE,gene.names=genenames)
  # Normalize the results across methods
  siggenes.result <- list()
  siggenes.result$UP <- siggenes$Table1[,3:5]
  colnames(siggenes.result$UP) <- c("fold.change","q.value","p.value")
  siggenes.result$DOWN <- siggenes$Table2[,3:5]
  colnames(siggenes.result$DOWN) <- c("fold.change","q.value","p.value")
  # Since RANKPROD does the goofy condition 1 / condition 2, inverse the fold-change and convert to log if not logged. 
  if (values_logged) {
    siggenes.result$UP[,"fold.change"] <- -siggenes.result$UP[,"fold.change"]
    siggenes.result$DOWN[,"fold.change"] <- -siggenes.result$DOWN[,"fold.change"]
  } else {
    siggenes.result$UP[,"fold.change"] <- log2(1/siggenes.result$UP[,"fold.change"])
    siggenes.result$DOWN[,"fold.change"] <- log2(1/siggenes.result$DOWN[,"fold.change"])
  }
  siggenes.result$UP=data.frame(siggenes.result$UP,probe=rownames(siggenes.result$UP))
  siggenes.result$DOWN=data.frame(siggenes.result$DOWN,probe=rownames(siggenes.result$DOWN))
  
  
} else if (method == 'SIGGENES_SAMR') {
  # Using only the SAM implementation in SIGGENES, not EBAM
  library(siggenes)
  
  SAM_result <- sam(comparison_frame,sample_class,rand=123,R.unlog=T,gene.names=genenames) #q.version=1
  delta.table <- rbind(c(0,0,0),findDelta(SAM_result,fdr=0.9)) #fdr=0.9, too loose
  siggenes.table <- summary(SAM_result, delta.table[  dim(delta.table)[1],  1] );
  siggenes <- siggenes.table@mat.sig
  
  if( nrow(siggenes) > 0 ) {
    siggenes.result <- list()
    siggenes.result$UP <- siggenes[siggenes$d.value > 0,c("R.fold","q.value","rawp")]
    colnames(siggenes.result$UP) <- c("fold.change","q.value","p.value")
    siggenes.result$DOWN <- siggenes[siggenes$d.value < 0,c("R.fold","q.value","rawp")]
    colnames(siggenes.result$DOWN) <- c("fold.change","q.value","p.value")
    
    
    siggenes.result$UP=data.frame(siggenes.result$UP,probe=rownames(siggenes.result$UP))
    siggenes.result$DOWN=data.frame(siggenes.result$DOWN,probe=rownames(siggenes.result$DOWN))
    
  }
  
} else if (method == 'SAMR') {
  library(samr)
  # class labels for SAMR are 1 or 2
  input_data <- list(
    x=data.matrix(comparison_frame),
    y=sample_class+1,
    geneid=rownames(comparison_frame),
    genenames=rownames(comparison_frame),
    logged2=T
  )
  
  samr.obj <- samr(input_data,resp.type="Two class unpaired",testStatistic="standard",nperms=100)
  delta.table <- samr.compute.delta.table(samr.obj)
  delta.index <- which.max( delta.table[,5] < 0.4 )
  delta=delta.table[delta.index,1]
  
  #replace data with input_data
  siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, input_data, delta.table)
  
  siggenes.result <- list()
  
  if (!is.null(siggenes.table$genes.up)>0){
    UP <- as.data.frame(subset(siggenes.table$genes.up,select=c("Gene ID","Fold Change","q-value(%)")))
    colnames(UP) <- c("probe","fold.change","q.value")
    UP$p.value  <- NA    
    UP$probe <- as.character(UP$probe)
    UP$fold.change <- as.double(as.character(UP$fold.change))
    UP$q.value <- as.double(as.character(UP$q.value))
    siggenes.result$UP <- UP
  }
  if (!is.null(siggenes.table$genes.lo)){
    DOWN <- as.data.frame(subset(siggenes.table$genes.lo,select=c("Gene ID","Fold Change","q-value(%)")))
    colnames(DOWN) <- c("probe","fold.change","q.value")
    DOWN$p.value  <- NA    
    DOWN$probe <- as.character(DOWN$probe)
    DOWN$fold.change <- as.double(as.character(DOWN$fold.change))
    DOWN$q.value <- as.double(as.character(DOWN$q.value))
    siggenes.result$DOWN <- DOWN
  }
  
  
} 

#annotate probe
genes.up <- merge(unique(egm[, c("GeneID", "Symbol")]), siggenes.result$UP, by.x ="GeneID", by.y = "probe", sort=F, all.x=T)
genes.down <- merge( unique(egm[, c("GeneID", "Symbol")]), siggenes.result$DOWN, by.x="GeneID",by.y = "probe", sort=F, all.x=T)

#removed mis-match probes
genes.up <- genes.up[!is.na(genes.up$GeneID) & genes.up$GeneID != "", ]
genes.down <- genes.down[!is.na(genes.down$GeneID) & genes.down$GeneID != "", ]

genes.up[,"up_down"] <- "up"
genes.down[,"up_down"] <- "down"

genes.down$fold.change <- -1 * 1/genes.down$fold.change
genes.sig.up <- subset(genes.up, q.value < q_thresh & abs(fold.change) > fold_change , select=c("GeneID","Symbol","fold.change","q.value","p.value","up_down"))
genes.sig.down <- subset(genes.down, q.value < q_thresh & abs(fold.change) > fold_change, select=c("GeneID","Symbol","fold.change","q.value","p.value","up_down"))

#write disease signature
genes.sig.up <- genes.sig.up[order(-genes.sig.up$fold.change),]
genes.sig.down <- genes.sig.down[order(genes.sig.down$fold.change),]

genes.sig.up <- genes.sig.up[1:min(nrow(genes.sig.up),150),]
genes.sig.down <- genes.sig.down[1:min(nrow(genes.sig.down),150),]
dz_sig <- rbind(genes.sig.up,genes.sig.down)

colnames(dz_sig) <- c("GeneID","Symbol","value","q.value","p.value","up_down")
write.table(dz_sig,dz_sig.file,sep="\t",quote=F,row.names=F,col.names=T)

#################################
#visualize disease signatures
annotation <- data.frame(type = c(rep("case", length(case)), rep("control", length(control))))
rownames(annotation) <- c(case, control)
annotation$type <- as.factor(annotation$type)
annotation <- subset(annotation, select=c("type"))

Var1        <- c("lightblue", "green")
names(Var1) <- c("case", "control")
anno_colors <- list(type = Var1)

my.cols <- bluered(100) 
comparison_frame_subset <- comparison_frame[rownames(comparison_frame) %in% dz_sig$GeneID, ]
pheatmap(scale(comparison_frame_subset), col = my.cols, annotation = annotation,  annotation_colors = anno_colors,
         show_colnames=F, legend=T, show_rownames=F, filename=paste(disease, "/dz_sig_validation.pdf", sep="")
)

