#' SomaticCombiner
#'
#' To merge two VCFs from different somatic callers.
#' @import VariantAnnotation
#' @importFrom abind abind
#' @importFrom GenomicRanges seqnames ranges
#' @importFrom MatrixGenerics rowRanges
#' @importFrom GenomeInfoDb sortSeqlevels
#' @importFrom IRanges DataFrameList
#' @param vcf1 VCF object or vcf file
#' @param vcf2 VCF object or vcf file
#' @param sources The caller of the VCF files.
#' @param GENO The GENO item sources of the mereged VCF. For example
#'     c(GT = 1, AD = 2). The vector names are the GENO names and the
#'     values are the corresponding index.
#' @param id_t The tumor sample ID.
#' @param id_n The normal/control sample ID.
#' @param pass_only If TRUE, only PASS variants will be merged.
#' @return A merged VCF object
#' @importClassesFrom IRanges CompressedList
#' @export

SomaticCombiner <- function(vcf1, vcf2, sources, GENO = c(GT = 1, DP = 1, AD = 1),
                         id_t = "TUMOR", id_n = "NORMAL", pass_only = FALSE){
    ## if(is.character(vcf1)){
    ##     v1 <- expand(readVcf(vcf1))
    ##     v2 <- expand(readVcf(vcf2))
    ## }else if (is(vcf1, "VCF")){
    ##     v1 <- expand(vcf1)
    ##     v2 <- expand(vcf2)
    ## }

    if(is.character(vcf1)){
        vcf1 <- readVcf(vcf1)
    }
    if(is.character(vcf2)){
        vcf2 <- readVcf(vcf2)
    }

    if(pass_only){
        vcf1 <- vcf1[fixed(vcf1)$FILTER == "PASS",]
        vcf2 <- vcf2[fixed(vcf2)$FILTER == "PASS",]
    }

    if(nrow(vcf1) == 0 & nrow(vcf2) == 0){
        return(vcf1)
    }else if(nrow(vcf1) == 0){
        return(vcf2)
    }else if(nrow(vcf2) == 0){
        return(vcf1)
    }
    
    ## v1 <- unique(expand(vcf1))
    ## v2 <- unique(expand(vcf2))
    v1 <- expand(vcf1)
    v2 <- expand(vcf2)

    ## fix ids
    pid1 <- paste0(seqnames(v1), ":", start(v1), "_", ref(v1), "/", alt(v1))
    pid2 <- paste0(seqnames(v2), ":", start(v2), "_", ref(v2), "/", alt(v2))
    rownames(v1) <- pid1
    rownames(v2) <- pid2

    ## remove duplicates
    v1 <- v1[!duplicated(pid1)]
    v2 <- v2[!duplicated(pid2)]
    pid1 <- paste0(seqnames(v1), ":", start(v1), "_", ref(v1), "/", alt(v1))
    pid2 <- paste0(seqnames(v2), ":", start(v2), "_", ref(v2), "/", alt(v2))
    
    vars <- list(v1, v2)
    ids <- c(id_t, id_n)
    names(ids) <- c("TUMOR", "NORMAL")
    vars <- lapply(vars, function(x){
        if(ncol(x)==2){
            if(any(is.na(match(ids, colnames(x))))){
                colnames(x) <- ids[match(names(ids), colnames(x))]
            }
            x[, match(ids, colnames(x))]
        }else{
            colnames(x) <- id_t
            x
        }
    })
    vn <- unlist(lapply(vars, nrow))
    if(any(vn == 0)){
        return(vars[[which(vn!=0)]])
    }

    ## shared variants
    ## pid1 <- paste0(seqnames(v1), ":", start(v1), "_", ref(v1), "/", alt(v1))
    ## pid2 <- paste0(seqnames(v2), ":", start(v2), "_", ref(v2), "/", alt(v2))
    vcom <- intersect(pid1, pid2)
    message("Unique variants in vcf1: ", length(setdiff(pid1, vcom)))
    message("Unique variants in vcf2: ", length(setdiff(pid2, vcom)))
    message("Shared variants: ", length(vcom))
    
    ## rowRanges
    idx2 <- match(setdiff(pid2, pid1), pid2)
    rR <- c(rowRanges(vars[[1]])[,1], rowRanges(vars[[2]])[idx2, 1])
    names(rR) <- c(pid1, setdiff(pid2, pid1))
    rR <- sortSeqlevels(rR)
    rR <- sort(rR)
    vid_u <- names(rR)

    ## colData
    cData <- DataFrame(Samples = seq(ids), row.names = ids)

    ## fixed
    f1 <- fixed(vars[[1]])
    f2 <- fixed(vars[[2]])
    f1$ALT <- as.character(f1$ALT)
    f2$ALT <- as.character(f2$ALT)
    if(sources[1]!="" & !is.null(sources[1]) & !grepl(sources[1], f1$FILTER[1])){
        f1$FILTER <- paste(sources[1], f1$FILTER, sep = "_")
    }
    if(sources[2]!="" & !is.null(sources[2]) & !grepl(sources[2], f2$FILTER[1])){
        f2$FILTER <- paste(sources[2], f2$FILTER, sep = "_")
    }
    f1u <- f1[pid1 %in% setdiff(pid1, vcom),]
    f2u <- f2[pid2 %in% setdiff(pid2, vcom),]
    f12 <- f1[match(vcom, pid1),]
    f21 <- f2[match(vcom, pid2),]
    f12$FILTER <- paste(f12$FILTER, f21$FILTER, sep = "|")
    f_m <- rbind(f1u, f12, f2u)
    rownames(f_m) <- c(setdiff(pid1, vcom), vcom, setdiff(pid2, vcom))
    f_m <- f_m[match(vid_u, rownames(f_m)),]
    f_m <- DataFrame(f_m)
    
    ## INFO
    info1 <- info(vars[[1]])
    if(sources[1]!="" & !is.null(sources[1]) & !grepl(sources[1], colnames(info1)[1])){
        colnames(info1) <- paste(sources[1], colnames(info1), sep = "_")
    }
    info1 <- DataFrame(varId = pid1, info1)

    info2 <- info(vars[[2]])
    if(sources[1]!="" & !is.null(sources[1]) & !grepl(sources[2], colnames(info2)[1])){
        colnames(info2) <- paste(sources[2], colnames(info2), sep = "_")
    }
    info2 <- DataFrame(varId = pid2, info2)

    info_m <- merge(info1, info2, by = "varId", all = TRUE)
    ## rownames(info_m) <- info_m$varId
    info_m <- info_m[match(vid_u, info_m$varId),]
    info_m <- info_m[, -match("varId", colnames(info_m))]
    rownames(info_m) <- vid_u
    info_m <- DataFrame(info_m)

    lists <- sapply(info_m, function(x){
        is.list(x) || is(x, "List")
    })
    infoL <- lapply(info_m[lists], function(x)as(x, "CompressedList"))
    info_m[lists] <- DataFrame(infoL)
    
    
    ## GENO
    ## shared
    var_s <- list(vars[[1]][match(vcom, pid1),], vars[[2]][match(vcom, pid2),])

    geno_list <- lapply(names(GENO), function(gid){
        geno(var_s[[GENO[[gid]]]])[[gid]]
    })
    names(geno_list) <- names(GENO)

    ## uniq v1/v2
    geno_var_u <- list()
    pid_u <- c(vcom, setdiff(pid1, vcom), setdiff(pid2, vcom))
    for(i in seq(names(GENO))){
        gid <- names(GENO)[i]
        pids <- list(pid1, pid2)
        geno1 <- geno_list[[i]]
        gdim <- dim(geno_list[[gid]])
        for(j in seq(2)){
            var1 <- vars[[j]][!pids[[j]] %in% vcom,]
            if(names(GENO)[i] %in% names(geno(var1))){
                g1 <- geno(var1)[[gid]]
                if(ncol(var1) == 1){
                    gdim[1] <- nrow(var1)
                    g <- array(".", dim = gdim,
                               dimnames=list(NULL, c(ids)))
                    if(length(gdim) == 2){
                        g[,id_t] <- g1[,id_t]
                        if(nrow(geno1) == 0){
                            geno1 <- g
                        }else{
                            geno1 <- rbind(geno1, g)
                        }
                    }else if(length(gdim) == 3){
                        g[,id_t, ] <- g1[, id_t, ]
                        if(nrow(geno1) == 0){
                            geno1 <- g
                        }else{
                            geno1 <- abind(geno1, g, along = 1)
                        }
                    }
                }else{
                    if(nrow(geno1) == 0){
                        geno1 <- g1
                    }else{
                        if(length(gdim) == 2){
                            geno1 <- rbind(geno1, g1)
                        }else{
                            geno1 <- abind(geno1, g1, along = 1)
                        }
                    }
                }
            }else{
                gdim[1] <- nrow(var1)
                g <- array(".", dim = gdim,
                           dimnames=list(NULL, c(ids)))
                if(nrow(geno1) == 0){
                    geno1 <- g
                }else{
                    if(length(gdim) == 2){
                        geno1 <- rbind(geno1, g)
                    }else{
                        geno1 <- abind(geno1, g, along = 1)
                    }
                }
            }
        }
        if(length(gdim) == 2){
            geno1 <- geno1[match(vid_u, pid_u),,drop=F]
        }else if(length(gdim) == 3){
            geno1 <- geno1[match(vid_u, pid_u),,,drop=F]
        }
        geno_var_u[[i]] <- geno1
    }
    names(geno_var_u) <- names(GENO)

    ## header
    ref <- seqlevels(rR)
    df1 <- header(header(vars[[1]]))
    df2 <- header(header(vars[[2]]))
    df_n <- intersect(names(df1), names(df2))
    df1_u <- df1[!names(df1) %in% df_n]
    df2_u <- df2[!names(df2) %in% df_n]
    df_com <- list()
    for(i in seq(df_n)){
        if(df_n[i] == "INFO"){
            if(sources[1]!="" & !is.null(sources[1]) & !grepl(sources[1], rownames(df1[["INFO"]])[1])){
                rownames(df1[["INFO"]]) <- paste(sources[1], rownames(df1[["INFO"]]), sep = "_")
            }
            if(sources[2]!="" & !is.null(sources[2]) & !grepl(sources[2], rownames(df2[["INFO"]])[1])){
                rownames(df2[["INFO"]]) <- paste(sources[2], rownames(df2[["INFO"]]), sep = "_")
            }
        }
        if(df_n[i] == "FORMAT"){
            df1[["FORMAT"]] <- df1[["FORMAT"]][rownames(df1[["FORMAT"]]) %in% names(GENO)[GENO == 1],]
            df2[["FORMAT"]] <- df2[["FORMAT"]][rownames(df2[["FORMAT"]]) %in% names(GENO)[GENO == 2],]
        }
        if(df_n[i] == "contig"){
            ## only use the first contig
            df_com[[i]] <- df1[[df_n[i]]]
        }else{
            df_com[[i]] <- unique(rbind(df1[[df_n[i]]], df2[[df_n[i]]]))
        }
    }
    names(df_com) <- df_n
    df_h <- c(df_com, as.list(df1_u), as.list(df2_u))
    df_h$fileformat <- sort(df_h$fileformat)[1,,drop = FALSE]
    
    VH <- VCFHeader(reference = ref, samples = ids, header = DataFrameList(df_h))
    vcf <- VCF(rowRanges = rR, colData = cData, exptData = list(header = VH),
               fixed = f_m, info = info_m, geno = geno_var_u, collapsed = FALSE)
    return(vcf)
}

