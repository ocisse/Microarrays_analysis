library(AnnotationForge)
library(mouse.db0)
library(org.Mm.eg.db)

# need to do this otherwise won't work because affy text is ill-formated
data <-read.csv("MoGene-2_0-st-v1.na33.mm10.transcript.noheader.csv",header=T,stringsAsFactors=F,sep=",")
sdata <- data[,c(1,9)]
returnRef=function(x){
        refst <- strsplit(x,split="///")[[1]][grep("RefSeq",strsplit(x,split="///")[[1]])[1]]
        refid <- gsub(" ","",strsplit(refst,split="//")[[1]][1])
        return(refid)
}

sdata$refseqids <- sapply(sdata[,2],returnRef)
fdata <- sdata[,-2]
write.table(fdata,"AnnotBuild.txt",
sep="\t",quote=F,row.names=F,col.names=F)

makeDBPackage("MOUSECHIP_DB",
        affy=F,
        prefix="mogene20sttranscriptcluster",
        fileName="AnnotBuild.txt",
        outputDir = ".",
        version="2.11.1",
        baseMapType="refseq",
        manufacturer = "Affymetrix",
        chipName = "Mouse Gene 2.0 ST Array",
        manufacturerUrl = "http://www.affymetrix.com",
        author = "Ousmane Cisse",
        maintainer = "Ousmane Cisse <ousmane.cisse@nih.gov>")
install.packages("mogene20sttranscriptcluster.db", repos=NULL, type="source")

targets_ccmd <- readTargets("affy_targets.txt")
dat <- read.celfiles(filenames=targets_ccmd$FileName)
est <-rma(dat)

annotation(est) <- "mogene20sttranscriptcluster.db"

phenoData(est)$varLabels <- targets_ccmd$Name  
wt <- est_ccmd[, est_ccmd$varLabels %in% c("WT_d35","WT_d120")]
f <- factor(as.character(as.factor(wt$varLabels)))
design <- model.matrix(~f)
fit <- eBayes(lmFit(wt,design))
tbl <- topTable(fit, coef=2)
gene.symbols <- getSYMBOL(rownames(tbl), "mogene20sttranscriptcluster.db")
