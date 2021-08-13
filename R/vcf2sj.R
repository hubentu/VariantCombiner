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
