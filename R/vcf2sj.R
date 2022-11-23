#' annotated VCF to SJ format
#' @param avcf A vcf from VEP.
#' @param tidx The tumor sample index.
#' @param nidx The normal sample index.
#' @export
vcf2sj <- function(avcf, tidx = 1, nidx = 2){
    vcf1 <- expand(readVcf(avcf))
    tid <- colnames(vcf1)[1]
    gann <- lapply(info(vcf1)$CSQ, function(x){
        idx <- grep("RefSeq", x)[1]
        idx <- ifelse(is.na(idx), 1, idx)
        unlist(strsplit(x[idx], split = "\\|"))[1:27]})
    gann <- do.call(rbind, gann)
    cn <- sub(".*: ", "", info(header(vcf1))["CSQ","Description"])
    colnames(gann) <- unlist(strsplit(cn, split = "\\|"))[1:27]

    NM <- gann[, "Feature"]
    AA <- paste0(gann[,"Protein_position"], ":", gann[,"Amino_acids"])

    AD <- geno(vcf1)$AD
    AD_T <- rbind(AD[,tidx,])
    AD_N <- rbind(AD[,nidx,])
    Ref <- as.character(ref(vcf1))
    Mut <- unlist(lapply(alt(vcf1), function(x)as.character(x)[1]))

    Qual <- fixed(vcf1)$FILTER

    data.frame(GeneName=gann[,"SYMBOL"], SJQuality=Qual, Sample = tid, Chr=as.character(seqnames(vcf1)),
               WU_HG19_Pos=start(vcf1), Class=gann[,"Consequence"], AAChange=AA, ProteinGI=NA, mRNA_acc=NM,
               "#Mutant_In_Tumor"=AD_T[,2], "#Total_In_Tumor"=rowSums(AD_T, na.rm=TRUE),
               "#Mutant_In_Normal"=AD_N[,2], "#Total_In_Normal"=rowSums(AD_N, na.rm=TRUE),
               ReferenceAllele=Ref, MutantAllele=Mut, check.names = F)
}

fixIndel <- function(vars, REF){
    Ref <- vars[, "ReferenceAllele"]
    Alt <- vars[, "MutantAllele"]
    Pos <- vars[, "posi"]

    ## in
    idx1 <- vars[,"ReferenceAllele"] == "-"
    rr1 <- GRanges(vars[idx1, "chr"], IRanges(vars[idx1,"posi"], width=1))
    ref1 <- scanFa(REF, param=rr1)
    ref1 <- as.character(ref1)
    Ref[idx1] <- ref1
    Alt[idx1] <- paste0(ref1, Alt[idx1])
    
    ## del
    idx2 <- vars[, "MutantAllele"] == "-"
    pos2 <- vars[idx2, "posi"] - 1
    ref2 <- scanFa(REF, param=GRanges(vars[idx2, "chr"], IRanges(pos2, width=1)))
    ref2 <- as.character(ref2)
    Ref[idx2] <- paste0(ref2, vars[idx2, "ReferenceAllele"])
    Alt[idx2] <- ref2
    Pos[idx2] <- pos2

    vars[, "posi"] <- Pos
    vars[, "ReferenceAllele"] <- Ref
    vars[, "MutantAllele"] <- Alt
    return(vars)
}

#' @export
sj2vcf <- function(v1, nid, ref = "b37", rfile = NULL){
    if(any(grepl("-", v1[,"ReferenceAllele"]) | grepl("-", v1[,"MutantAllele"]))){
        v1 <- fixIndel(v1, rfile)
    }
    r1 <- GRanges(v1$chr, IRanges(v1$posi, v1$posi))
    cdat <- DataFrame(Samples = 1:2, row.names = c(v1$sample[1], nid))
    f1 <- DataFrame(REF = DNAStringSet(v1$ReferenceAllele), ALT = v1$MutantAllele, QUAL = NA_integer_, FILTER = "PASS")
    i1 <- DataFrame(v1[,c("SJQuality", "Gene", "class", "aachange", "mRNA_acc")])
    g1 <- SimpleList(GT = cbind(Tumor = rep("0/1", nrow(v1)), Normal = rep("0/0", nrow(v1))),
                     DP = cbind(Tumor = v1$Total_In_Tumor, Normal = v1$Total_In_Normal),
                     AD = cbind(Tumor = v1$Mutant_In_Tumor, Normal = v1$Mutant_In_Normal))
    g1 <- lapply(g1, function(x){
        colnames(x) <- c(v1$sample[1], nid)
        x})

    h1 <- VCFHeader(reference = ref, samples = c(v1$sample[1], nid),
                    header = DataFrameList(
                        fileformat = DataFrame(Value = "VCFv4.1", row.names = "fileformat"),
                        FORMAT = DataFrame(Number = c(1,1,1),
                                           Type = c("String", "Integer", "Integer"),
                                           Description = c("Genotype", "Depths", "Alt counts"),
                                           row.names = c("GT", "DP", "AD"))))

    vcf <- VCF(rowRanges = r1, colData = cdat, fixed = f1, exptData = list(header = h1),
        info = i1, geno = SimpleList(g1), collapsed = FALSE)
    sort(sortSeqlevels(vcf))
}

#' @export
#' @importFrom Rsamtools scanFa
sj2maf <- function(sj, ref = NULL){
    ## varclass <- rep(NA, nrow(sj))
    varclass <- sj$class
    varclass[sj[,"class"] == "frameshift" & sj[,"ReferenceAllele"]=="-"] <- "Frame_Shift_Ins"
    varclass[sj[,"class"] == "frameshift" & sj[,"MutantAllele"]=="-"] <- "Frame_Shift_Del"
    varclass[sj[,"class"] == "splice"] <- "Splice_Site"
    varclass[sj[,"class"] == "nonsense"] <- "Nonsense_Mutation"
    varclass[sj[,"class"] == "stoploss"] <- "Nonstop_Mutation"
    varclass[sj[,"class"] == "proteinDel"] <- "In_Frame_Del"
    varclass[sj[,"class"] == "proteinIns"] <- "In_Frame_Ins"
    varclass[sj[,"class"] == "missense"] <- "Missense_Mutation"
    varclass[sj[,"class"] == "MNV"] <- "MNV"
    varclass[sj[,"class"] == "silent"] <- "Silent"
    
    vartype <- rep("SNP", nrow(sj))
    idx <- sj[,"ReferenceAllele"]=="-" | (length(sj[,"ReferenceAllele"])==1 & sj[,"MutantAllele"]>1)
    vartype[idx] <- "INS"
    idx <- sj[,"MutantAllele"]=="-" | (length(sj[,"ReferenceAllele"])>1 & sj[,"MutantAllele"]==1)
    vartype[idx] <- "DEL"

    if(any(grepl("-", sj[,"ReferenceAllele"]) | grepl("-", sj[,"MutantAllele"]))){
        sj <- fixIndel(sj, ref)
    }
    ## if(any(grepl("-", sj[, "ReferenceAllele"]))){
    ##     if(is.null(ref))stop("ref required!")
    ##     idx <- grep("-", sj[, "ReferenceAllele"])
    ##     sj1 <- sj[idx,]
    ##     gr <- GRanges(sj1[,"chr"], IRanges(sj1[,"posi"], sj1[, "posi"])) 
    ##     fa1 <- scanFa(ref, gr)
    ##     sj[idx, "ReferenceAllele"] <- as.character(fa1)
    ##     sj[idx, "MutantAllele"] <- paste0(as.character(fa1), sj[idx, "MutantAllele"])
    ## }
    ## if(any(grepl("-", sj[, "MutantAllele"]))){
    ##     if(is.null(ref))stop("ref required!")
    ##     idx <- grep("-", sj[, "MutantAllele"])
    ##     sj1 <- sj[idx,]
    ##     gr <- GRanges(sj1[,"chr"], IRanges(sj1[,"posi"]-1, sj1[, "posi"]-1)) 
    ##     fa1 <- scanFa(ref, gr)
    ##     sj[idx, "ReferenceAllele"] <- paste0(as.character(fa1), sj[idx, "ReferenceAllele"])
    ##     sj[idx, "MutantAllele"] <- as.character(fa1)
    ##     sj[idx, "posi"] <- sj[idx, "posi"]-1
    ## }
    
    maf <- data.frame(Hugo_Symbol = sj[, "Gene"],
                      Entrez_Gene_Id = NA,
                      Center = "RPCCC",
                      NCBI_Build = "37",
                      Chromosome = sj[, "chr"],
                      Start_Position = sj[,"posi"],
                      End_Position = sj[,"posi"],
                      Strand = "+",
                      Variant_Classification = varclass,
                      Variant_Type = vartype,
                      Reference_Allele = sj[,"ReferenceAllele"],
                      Tumor_Seq_Allele1 = sj[,"ReferenceAllele"],
                      Tumor_Seq_Allele2 = sj[,"MutantAllele"],
                      dbSNP_RS = "novel",
                      dbSNP_Val_Status = "",
                      Tumor_Sample_Barcode = sj[,"sample"],
                      Matched_Norm_Sample_Barcode = sub("-.*", "-N", sj[, "sample"]),
                      Protein_position = sj[,"aachange"],
                      Transcript_ID = sj[,"mRNA_acc"])
    maf <- maf[!is.na(maf$Variant_Classification),]
    return(maf)
}

