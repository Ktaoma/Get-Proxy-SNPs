devtools::install_github("CBIIT/LDlinkR")
library(LDlinkR)
library(dplyr)
library(data.table)
library(splitstackshape)
library(ggplot2)

#CHeck if there any SNP on GSA array
df_query <- readxl::read_xlsx("input/Nature-2017-supp-PRS_313SNPs.xlsx")  %>% as.data.frame()

snp_mart = useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_snp")

SNPs <- readxl::read_xlsx("input/Nature-2017-supp-PRS_313SNPs.xlsx")[,c(2,3)]
SNPs$END <- SNPs$Positionb  
SNPs_id <- paste0(SNPs$Chromosome,":",SNPs$Positionb,":",SNPs$END)

get_rsid <- function(x){
  print(x)
  res <- getBM(attributes = c('refsnp_id', 'chr_name', 'chrom_start', 'chrom_end', 'allele'),
               filters = c('chromosomal_region'), 
               values = x, 
               mart = snp_mart,
               useCache = FALSE)  
  return(res)
}

res_tmp <- lapply(SNPs_id, get_rsid)

#############################################

get_proxy_snp <- function(x,ref){
  print(x)
  res <- LDproxy(x,pop=ref,r2d="r2",token="2692a80eeef4")
  res$query <- x
  res$rsid <- res[1,1]
  return(res)
}

pop <- c("EAS","SAS")
get_proxy_EAS_SAS <- lapply(df_query_list, get_proxy_snp,ref=pop) 
pop <- c("SAS")
get_proxy_SAS <- lapply(df_query_list, get_proxy_snp,ref=pop)

############################################

#find which SNPs are error
#If number of coluumn is less than 12, it implies those SNPs are monoallelic and not SNPs 
define_error <- function(i){
  error_snp <- which(ncol(as.data.frame(get_proxy_EAS_SAS[i])) < 12)
  return(error_snp)
}

detected_error_snp <- lapply(c(1:length(get_proxy_EAS_SAS)),define_error)
detected_error_snp_02 <- which(detected_error_snp == 1)

#detected_snp
detectable_snp <- do.call(rbind,get_proxy_EAS_SAS[-detected_error_snp_02])
write.csv(detectable_snp,"final_proxy_EAS_SAS.csv")

#not_detected_snp
not_detected_snp <- do.call(rbind,get_proxy_EAS_SAS[detected_error_snp_02])
write.csv(not_detected_snp,"SNPs_not_detected_EAS_SAS.csv")

