library(data.table)
library(dplyr)
library(splitstackshape)
library(tidyverse)


############## Check proxy available in GSA array #############################################################

#remove variant in GSA array without rsID ((RsID != ".")
GSA_array <- fread("input/GSA-24v3-0_A1_b151_rsids.txt") %>% filter(RsID != ".") %>% select(RsID) %>%  unlist()

#64 SNPs in summary statistics avaialbe in GSA array, No need proxy !!! 
SNP_avialable_in_GSA <- fread("result/final_proxy_EAS_SAS.csv") %>% filter(rsid %in% GSA_array) %>% select(query) %>% unique() %>% unlist()

#filter the optimal proxy
#regulome score index
rank_regulome_db <- data.frame(RegulomeDB = c("1a","1b","1c","1d","1e","1f","2a","2b","2c","3a","3b","4","5","6","7","."),
                               rank = c(1:16))

#filter proxy SNPs only avaiable in GSA array
# And select proxy avaialble in GSA array with rsID only 

SNP_not_avialable_in_GSA <- fread("result/final_proxy_EAS_SAS.csv") %>% 
  filter(!query %in% SNP_avialable_in_GSA) %>%
  filter(RS_Number %in% GSA_array)


# remove proxy without rsnumber (RS_Number != ".")
# filter Proxy with Dprime at least 0.9 (Dprime > 0.95)
# filter regulomeDB score at lowest value (more information annotated)
# filter proxy SNPs with lowest MAF  (lowest MAF, more effect )
# Dprime at 0.7 is set because it too strict some proxy SNPs will be removed

df_rank <- inner_join(rank_regulome_db,SNP_not_avialable_in_GSA,by="RegulomeDB") %>% 
  filter(Dprime > 0.70)  %>% 
  group_by(query) %>% 
  slice_min(order_by = rank, n = 1) %>%
  slice_min(order_by = MAF, n = 1) %>%  
  slice_max(order_by = R2, n = 1) 

#check which SNPs still have more than one candidate 
SNPs_two_proxy <- table(df_rank$query) %>% as.data.frame() %>% filter(Freq != 1) %>% select(Var1) %>% unlist()
df_rank$result <- ifelse(df_rank$query %in% SNPs_two_proxy,"More than one proxy SNPs","One proxy SNPs")
df_rank_final <- df_rank[c(1,4,6,7,9,13,14,15)]

####### Combine result between SNP available and not avaiable #############
SNP_avialable_in_GSA_02 <- fread("result/final_proxy_EAS_SAS.csv") %>% filter(Coord %in% SNP_avialable_in_GSA,Dprime == 1,R2==1)
SNP_avialable_in_GSA_02$result <- "Available in GSA array"
SNP_avialable_in_GSA_final <- SNP_avialable_in_GSA_02 %>% select(colnames(df_rank_final))

all_result <- rbind(df_rank_final %>% as.data.frame(),SNP_avialable_in_GSA_final%>% as.data.frame())
names(all_result) <- c("RegulomeDB","Proxy SNPs","Ref/Alt","Dprime","query variant","query","rsid","result")

write.csv(all_result,"proxy_SNPs_EAS_SAS.csv")



#plot

res_plot <- all_result %>% select(query,result) %>% unique() %>% select(result) %>% table() %>% as.data.frame() 
names(res_plot) <- c("Cat","Freq")
ggplot(res_plot) <- 
  
  
  