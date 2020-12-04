devtools::install_github("CBIIT/LDlinkR")
library(LDlinkR)
library(dplyr)
library(data.table)


#CHeck if there any SNP on GSA array
GSA_df <- fread("GSA-24v3-0_A1_b151_rsids.txt")
df_query <- readxl::read_xlsx("Nature-2017-supp-PRS_313SNPs.xlsx")  %>% as.data.frame()
df_query$id <- paste0(df_query$Chromosome,":",df_query$Positionb) 
df_query %>%  filter(id %in% GSA_df$Name) #No SNPs in summary statistic available in GSA array at all

#input list of rsID (example rs1234) and chr:position (example chr1:1234567)
df_query_list <- paste0("chr",df_query$Chromosome,":",df_query$Positionb) 
get_proxy_snp <- function(x){
  print(x)
  res <- LDproxy(x,pop="EAS",r2d="r2",token="2692a80eeef4")
  res$query <- x
  res$rsid <- res[1,1]
  return(res)
}

system.time(get_proxy <- lapply(df_query_list, get_proxy_snp)) ##313 SNPs takes around 40 minutes to finish


#find which SNPs are error
#If number of coluumn is less than 12, it implies those SNPs are monoallelic and not SNPs 
define_error <- function(i){
  error_snp <- which(ncol(as.data.frame(get_proxy[i])) < 12)
  return(error_snp)
}

detected_error_snp <- lapply(c(1:length(get_proxy)),define_error)
detected_error_snp_02 <- which(detected_error_snp == 1)

#detected_snp
detectable_snp <- do.call(rbind,df_tmp[-detected_error_snp_02])
write.csv(detectable_snp,"final_proxy.csv")

#not_detected_snp
not_detected_snp <- do.call(rbind,df_tmp[detected_error_snp_02])
write.csv(not_detected_snp,"SNPs_not_detected.csv")


