#' Add GT,AD,DP to seurat VCF
#' @param vcf Somatic VCF from seurat
#' @export
seurat_geno <- function(vcf){
    stopifnot(is(vcf, "VCF"))
    if(is(vcf, "CollapsedVCF")){
        vcf <- expand(vcf)
    }

    AR_n <- info(vcf)$AR1
    AR_t <- info(vcf)$AR2
    DP_n <- info(vcf)$DP1
    DP_t <- info(vcf)$DP2
    
    AD_ref <- cbind(round(DP_n * (1-AR_n)), round(DP_t * (1-AR_t)))
    AD_alt <- cbind(round(DP_n * AR_n), round(DP_t * AR_t))
    colnames(AD_ref) <- colnames(AD_alt) <- c("NORMAL", "TUMOR")

    AD <- abind(data.frame(AD_ref[, c("NORMAL", "TUMOR")]),
                data.frame(AD_alt[, c("NORMAL", "TUMOR")]), along = 3)

    f <- DataFrame(Number = c(1, "R", "R"),
                   Type = c("String", "Integer", "Integer"),
                   Description = c("Genotype",
                                   "Read depth",
                                   "Allelic depths for the ref and alt alleles in the order listed"))
    rownames(f) <- c("GT", "DP", "AD")
    vcf@metadata$header@header$FORMAT <- f

    vcf@colData <- DataFrame(row.names=c("NORMAL", "TUMOR"), Samples=c(1,2))
    geno(vcf)$DP <- cbind(NORMAL = DP_n, TUMOR = DP_t)
    geno(vcf)$GT <- cbind(NORMAL = rep("0/0", nrow(vcf)),
                          TUMOR = rep("0/1", nrow(vcf)))
    geno(vcf)$AD <- AD
    return(vcf)
}
