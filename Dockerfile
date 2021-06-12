FROM bioconductor/bioconductor_docker:RELEASE_3_13

RUN R -e 'BiocManager::install("R.utils")'
RUN R -e 'BiocManager::install("hubentu/VariantCombiner", dependencies = TRUE)'
