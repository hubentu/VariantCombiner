utils::globalVariables("name")

#' Add GT, AD to strelka2 snv calls
#'
#' @param vcf The SNV VCF object from VariantAnnotation by strelka2.
#' @importFrom tidyr pivot_wider pivot_longer separate
#' @importFrom dplyr left_join
#' @importFrom S4Vectors DataFrame
#' @importFrom methods is
#' @return VCF object.
#' @export
strelka_snv <- function(vcf){
    stopifnot(is(vcf, "VCF"))
    if(is(vcf, "CollapsedVCF")){
        vcf <- expand(vcf)
    }

    if(nrow(vcf) == 0){
        AD_ref <- AD_alt <- data.frame(NORMAL=integer(), TUMOR=integer())
    } else {
        Ref <- as.character(ref(vcf))
        Mut <- unlist(lapply(alt(vcf), function(x)as.character(x)[1]))
        
        AU <- geno(vcf)$AU[,,1, drop=F]
        CU <- geno(vcf)$CU[,,1, drop=F]
        GU <- geno(vcf)$GU[,,1, drop=F]
        TU <- geno(vcf)$TU[,,1, drop=F]
        mcount0 <- data.frame(rn = rownames(vcf), Ref, Mut, A=AU, C=CU, G=GU, T=TU)
        mcount <- pivot_longer(mcount0, cols = -c(1:3))
        mcount <- separate(mcount, name, c("Allele", "TN"), sep = "\\.", extra = "drop")
        AD_ref <- left_join(mcount0, mcount[,-c(2:3)], by = c("rn" = "rn", "Ref" = "Allele"))
        AD_ref <- pivot_wider(AD_ref[, c("rn", "TN", "value")], id_cols = "rn", names_from = "TN")
        AD_alt <- left_join(mcount0, mcount[,-c(2:3)], by = c("rn" = "rn", "Mut" = "Allele"))
        AD_alt <- pivot_wider(AD_alt[, c("rn", "TN", "value")], id_cols = "rn", names_from = "TN")
    }
    
    AD <- abind(data.frame(AD_ref[, c("NORMAL", "TUMOR")]),
                data.frame(AD_alt[, c("NORMAL", "TUMOR")]), along = 3)

    hdf <- header(header(vcf))$FORMAT
    f <- DataFrame(Number = c(1, "R"),
                   Type = c("String", "Integer"),
                   Description = c("Genotype",
                                    "Allelic depths for the ref and alt alleles in the order listed"))
    rownames(f) <- c("GT", "AD")
    hdf <- rbind(f, hdf)
    vcf@metadata$header@header$FORMAT <- hdf
    
    geno(vcf)$GT <- cbind(NORMAL = rep("0/0", nrow(vcf)),
                          TUMOR = rep("0/1", nrow(vcf)))
    geno(vcf)$AD <- AD
    return(vcf)
}


#' Add GT, AD to strelka2 indel calls
#' 
#' @param vcf The Indel VCF object from VariantAnnotation by strelka2.
#' @return VCF object.
#' @export
strelka_indel <- function(vcf){
    stopifnot(is(vcf, "VCF"))
    if(is(vcf, "CollapsedVCF")){
        vcf <- expand(vcf)
    }
    AD_t1 <- geno(vcf)$TIR[,"TUMOR",1]
    AD_t0 <- geno(vcf)$DP[,"TUMOR"] - AD_t1
    AD_n1 <- geno(vcf)$TIR[,"NORMAL",1]
    AD_n0 <- geno(vcf)$DP[,"NORMAL"] - AD_n1

    AD <- abind(cbind(AD_n0, AD_t0),
                cbind(AD_n1, AD_t1), along = 3)
    colnames(AD) <- c("NORMAL", "TUMOR")
    
    hdf <- header(header(vcf))$FORMAT
    f <- DataFrame(Number = c(1, "R"),
                   Type = c("String", "Integer"),
                   Description = c("Genotype",
                                    "Allelic depths for the ref and alt alleles in the order listed"))
    rownames(f) <- c("GT", "AD")
    hdf <- rbind(f, hdf)
    vcf@metadata$header@header$FORMAT <- hdf
    
    geno(vcf)$GT <- cbind(NORMAL = rep("0/0", nrow(vcf)),
                          TUMOR = rep("0/1", nrow(vcf)))
    geno(vcf)$AD <- AD
    return(vcf)
}
