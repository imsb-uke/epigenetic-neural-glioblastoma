## To run:
## Rscript process_array.R {idat directory} {output csv path} {path to ref_sample.RData}
## Note: Requires minfi package.
## To install minfi:
## try http:// if https:// URLs are not supported
## source("https://bioconductor.org/biocLite.R")
## biocLite("minfi")

args<-commandArgs(trailingOnly=TRUE)
if (!length(args)==3){
	stop('Usage: Rscript process_array.R {idat directory} {output csv path} {path to ref_sample.RData}')
}


library(minfi)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(IlluminaHumanMethylationEPICv2manifest)

.isRGOrStop <- function(object) {
    if (!is(object, "RGChannelSet")) {
        stop("object is of class '", class(object), "', but needs to be of ",
             "class 'RGChannelSet' or 'RGChannelSetExtended'")
    }
}

.is27k <- function(object) {
    annotation(object)["array"] == "IlluminaHumanMethylation27k"
}
.is450k <- function(object) {
    annotation(object)["array"] == "IlluminaHumanMethylation450k"
}
.isEPIC <- function(object) {
    annotation(object)["array"] == "IlluminaHumanMethylationEPIC"
}
.isAllergy <- function(object) {
    annotation(object)["array"] == "IlluminaHumanMethylationAllergy"
}

.isEPICv2 <- function(object) {
    annotation(object)["array"] == "IlluminaHumanMethylationEPICv2"
}

pmax2 <- function(x, y) {
    pmax(x, y)
}
pmin2 <- function(x, y) {
    pmin(x, y)
}

.getManifestString <- function(annotation) {
    if (length(annotation) == 1) {
        if (annotation == "Unknown") {
            stop("Cannot get Manifest object for an 'Unknown' array")
        }
        return(paste0(annotation, "manifest"))
    }
    if ("array" %in% names(annotation)) {
        if (annotation["array"] == "Unknown") {
            stop("Cannot get Manifest object for an 'Unknown' array")
        }
        return(paste0(annotation["array"], "manifest"))
    }
    stop("unable to get the manifest string for this object")
}

normalize.illumina.control <- function(rgSet, reference = 1) {
    Green <- getGreen(rgSet)
    Red <- getRed(rgSet)

    if (.is450k(rgSet) || .isEPIC(rgSet) || .isAllergy(rgSet) || .isEPICv2(rgSet)) {
        AT.controls <- getControlAddress(
            object = rgSet,
            controlType = c("NORM_A", "NORM_T"))
        CG.controls <- getControlAddress(
            object = rgSet,
            controlType = c("NORM_C", "NORM_G"))
    }
    if (.is27k(rgSet)) {
        AT.controls <- getControlAddress(
            object = rgSet,
            controlType = "Normalization-Red")
        CG.controls <- getControlAddress(
            object = rgSet,
            controlType = "Normalization-Green")
    }

    Green.avg <- colMeans2(Green, rows = match(CG.controls, rownames(Green)))
    Red.avg <- colMeans2(Red, rows = match(AT.controls, rownames(Red)))
    ref <- (Green.avg + Red.avg)[reference] / 2
    if (is.na(ref)) {
        stop("perhaps 'reference' refer to an array that is not present.")
    }
    Green.factor <- ref / Green.avg
    Red.factor <- ref / Red.avg
    Green <- sweep(Green, 2, FUN = "*", Green.factor)
    Red <- sweep(Red, 2, FUN = "*", Red.factor)
    assay(rgSet, "Green") <- Green
    assay(rgSet, "Red") <- Red
    rgSet
}

bgcorrect.illumina <- function(rgSet) {
    .isRGOrStop(rgSet)
    Green <- getGreen(rgSet)
    Red <- getRed(rgSet)
    if (.is450k(rgSet) || .isEPIC(rgSet) || .isAllergy(rgSet) || .isEPICv2(rgSet)) {
        NegControls <- getControlAddress(rgSet, controlType = "NEGATIVE")
    }
    if (.is27k(rgSet)) {
        NegControls <- getControlAddress(rgSet, controlType = "Negative")
    }
    Green.bg <- apply(Green[NegControls, , drop = FALSE], 2, function(xx) {
        sort(as.vector(xx))[31]
    })
    Red.bg <- apply(Red[NegControls, , drop = FALSE], 2, function(xx) {
        sort(as.vector(xx))[31]
    })
    Green <- pmax2(sweep(Green, 2, Green.bg), 0)
    Red <- pmax2(sweep(Red, 2, Red.bg), 0)
    assay(rgSet, "Green") <- Green
    assay(rgSet, "Red") <- Red
    rgSet
}

preprocessIllumina <- function(rgSet, bg.correct = TRUE,
                               normalize = c("controls", "no"), reference = 1) {
    .isRGOrStop(rgSet)
    normalize <- match.arg(normalize)

    if (normalize == "controls") {
        rgSet <- normalize.illumina.control(rgSet, reference = reference)
    }
    if (bg.correct) {
        rgSet <- bgcorrect.illumina(rgSet)
    }
    out <- preprocessRaw(rgSet)
    preprocess <- sprintf(
        "Illumina, bg.correct = %s, normalize = %s, reference = %d",
        bg.correct, normalize, reference)

    out@preprocessMethod <- c(
        rg.norm = preprocess,
        minfi = as.character(packageVersion("minfi")),
        manifest = as.character(
            packageVersion(.getManifestString(rgSet@annotation))))
    out
}

convert_to_epicv1 <- function(betas) {
  base_names <- gsub("_[BT][CO][0-9]+$", "", rownames(betas))
  unique_names <- unique(base_names)
  
  epicv1_betas <- matrix(0, nrow = length(unique_names), ncol = ncol(betas))
  rownames(epicv1_betas) <- unique_names
  colnames(epicv1_betas) <- colnames(betas)
  
  for(i in seq_along(unique_names)) {
    idx <- which(base_names == unique_names[i])
    if(length(idx) == 1) {
      epicv1_betas[i, ] <- betas[idx, ]
    } else {
      epicv1_betas[i, ] <- colMeans(betas[idx, , drop = FALSE], na.rm = TRUE)
    }
  }
  
  return(epicv1_betas)
}

get_betas <- function(grSet,RGChannelSet,p_threshold=0.01,min_beads=3,sex_remove=T){
  # filter by P-value
  betas <- getBeta(grSet)
  pval <- detectionP(RGChannelSet)[rownames(betas),]
  betas[pval>p_threshold] <- NA
  # filter sex chromosomes
  if (sex_remove){
    chrom <- as.character(seqnames(granges(grSet)))
    chrom.sex <- chrom=='chrX' | chrom=='chrY'
    betas <- betas[!chrom.sex,]
  }
  #filter by bead number
  # nbeads <- getNBeads(RGChannelSet)
  # t1 <- getProbeInfo(RGChannelSet,'I')
  # t2 <- getProbeInfo(RGChannelSet,'II')
  # t1.nbeads.u <- nbeads[rownames(nbeads) %in% t1[,2],]; rownames(t1.nbeads.u) <- t1[match(rownames(t1.nbeads.u),t1[,2]),1]
  # t1.nbeads.m <- nbeads[rownames(nbeads) %in% t1[,3],]; rownames(t1.nbeads.m) <- t1[match(rownames(t1.nbeads.m),t1[,3]),1]
  # t1.nbeads.m <- t1.nbeads.m[rownames(t1.nbeads.u),]
  # t1.nbeads.min <- pmin(t1.nbeads.u,t1.nbeads.m)
  # t2.nbeads <- nbeads[rownames(nbeads) %in% t2[,2],]; rownames(t2.nbeads) <- t2[match(rownames(t2.nbeads),t2[,2]),1]
  # nbeads.array <- rbind(t1.nbeads.min,t2.nbeads)
  # cgs <- intersect(rownames(nbeads.array),rownames(betas))
  # nbeads.array <- nbeads.array[cgs,]
  # betas <- betas[cgs,]
  # betas[nbeads.array<min_beads] <- NA
  return(betas)
}

process_array <- function(array_dir, out_path, ref_sample_path=args[3])
{
	array_path <- array_dir
	
	print(paste0("1. process_array script - READING ARRAY - ",Sys.time()))
	rgSet <- read.metharray.exp(array_path,extended=T,force=T)
    rgSet@annotation[["array"]] <- "IlluminaHumanMethylationEPICv2"
    rgSet@annotation[["annotation"]] <- "20a1.hg38"
	
	# Get p-values	
	print(paste0("2. process_array script - P_VALUES CALCULATION - ",Sys.time()))
	# pVal <- detectionP(rgSet, type = "m+u")

	# illumina normalization
	## join with reference sample, normalize, then remove reference sample
	print(paste0("3. process_array script - ILLUMINA - ",Sys.time()))
	# load(ref_sample_path)
	# ref.sample@colData <- rgSet@colData[1,]; rownames(ref.sample@colData)[1] <- 'REF'
	# rgSet <- combineArrays(ref.sample,rgSet)
	MethSet <- preprocessIllumina(rgSet)
	
	# Context removal
	print(paste0("4. process_array script - CONTEXT_REMOVAL - ",Sys.time()))
	MethDrop <- dropMethylationLoci(MethSet, dropRS = TRUE, dropCH = TRUE)

	# Get RatioSet
	print(paste0("5. process_array script - CONVERT TO RATIOSET - ",Sys.time()))
	RSet <- ratioConvert(MethDrop, what = "beta", keepCN = FALSE,offset=100)

	# Get GenomicRatioSet
	print(paste0("6. process_array script - CONVERT TO GENOMICRATIOSET - ",Sys.time()))
    grSet <- mapToGenome(RSet)
	
	# Remove SNPs
	print(paste0("7. process_array script - REMOVE SNPs - ",Sys.time()))
    grSet <- dropLociWithSnps(grSet, snps=c("SBE","CpG"), maf=0.01)
	
	# Get filtered beta values
	print(paste0("8. process_array script - BETA VALUES CALCULATION - ",Sys.time()))
	betas<-getBeta(grSet) #getBeta(grSet,rgSet)
    # betas<-convert_to_epicv1(betas)
	
	print(paste0("13. process_array script - WRITING ARRAY TO CSV - ",Sys.time()))
	# if (length(sampleNames(rgSet)) == 2){
	# 	betas <- as.data.frame(betas)
	# 	colnames(betas)<-sampleNames(rgSet)[2]
	# }
	write.csv(round(betas,5),out_path,quote=FALSE)
	print(paste0("14. process_array script - FINISHED - ",Sys.time()))	
}

process_array(args[1], args[2])
