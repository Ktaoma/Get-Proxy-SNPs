library(data.table)
library(dplyr)
library(splitstackshape)
library(tidyverse)


#Get rsID from proxy_SNPs_EAS_SAS in not detected SNPs file by LDlink
not_detected_snp_EAS_SAS <- fread("result/SNPs_not_detected_EAS_SAS.csv",header = T)
rs_ID_not_detected_snp_EAS_SAS <- cSplit(not_detected_snp_EAS_SAS,"rsid"," ") %>% dplyr::select(query,rsid_2)
names(rs_ID_not_detected_snp_EAS_SAS) <- c("query","rsid")

#Get rsID from proxy_SNPs_EAS_SAS in detected SNPs file by LDlink
detected_snp_EAS_SAS <- fread("result/final_proxy_EAS_SAS.csv",header = T) %>% dplyr::select(query,rsid) %>% unique()

#combine two files
all_rs_id <- rbind(detected_snp_EAS_SAS,rs_ID_not_detected_snp_EAS_SAS)
all_rs_id$query <- gsub('chr', '',all_rs_id$query)
all_rs_id_02 <- all_rs_id %>% filter(rsid != ".")

#Check which SNPs avaialbe on GSA array
SNP_available_in_GSA_array_01 <- fread("input/GSA-24v3-0_A1_b151_rsids.txt") %>% filter(RsID %in% all_rs_id_02$rsid) %>% 
  select(RsID) %>% setnames("RsID","Name")
SNP_available_in_GSA_array_02 <- fread("input/GSA-24v3-0_A1_b151_rsids.txt") %>% filter(Name %in% all_rs_id_02$rsid) %>% 
  select(Name) 

SNP_available_in_GSA_array <- rbind(SNP_available_in_GSA_array_01,SNP_available_in_GSA_array_02) %>% unique()


#Annotate from LDlink result by adding column which proxy SNPs available in GSA array
detected_snp_EAS_SAS_df <- fread("result/final_proxy_EAS_SAS.csv",header = T)
detected_snp_EAS_SAS_df$type <- "SNPs"
not_detected_snp_EAS_SAS <- fread("result/SNPs_not_detected_EAS_SAS.csv",header = T)
not_detected_snp_EAS_SAS_02 <- cSplit(not_detected_snp_EAS_SAS,"rsid"," ") %>% dplyr::select(query,rsid_4,rsid_2)
names(not_detected_snp_EAS_SAS_02) <- c("query","type","rsid")

cdf <- bind_rows(detected_snp_EAS_SAS_df,not_detected_snp_EAS_SAS_02)%>% replace(is.na(.), 0)
cdf$Coord <- ifelse(cdf$Coord == 0, cdf$query,cdf$Coord)
cdf$RS_Number <- ifelse(cdf$RS_Number == 0,cdf$rsid,cdf$RS_Number)
cdf$cat <- ifelse(cdf$rsid %in% SNP_available_in_GSA_array$Name,"Available on GSA array",
                  "Not available on GSA array")

cdf[cdf ==  "not"] <- "Deletion mutation"
cdf[cdf ==  "SNPs"] <- "Proxy SNPs available"
cdf[cdf ==  "monoallelic"] <- "Monoallelic SNPs"

#Summarize data first before get the final set of proxy SNPs
cdf_summarize <- cdf %>% dplyr::select(query,cat,type) %>% unique() 
cdf_summarize$type_02 <- ifelse(cdf_summarize$cat == "Available on GSA array",
                                "Available on GSA array",
                                cdf_summarize$type)

table(cdf_summarize$type_02) %>% as.data.frame() %>% 
  ggplot(aes(reorder(Var1,Freq),Freq)) +
  geom_bar(stat="identity",position="dodge") +
  coord_flip() +
  geom_text(aes(label=Freq),
            position = position_dodge(width = 1),
            hjust = -0.5, size = 4,color="black") +
  ylim(c(0,260)) +
  theme_bw() +
  labs(color="",fill="") +
  ylab("Frequency") +
  xlab("Categories") 


#filter the optimal proxy
# select proxy SNPs avaialble in GSA array only

all_SNPs_in_GSA_1 <- fread("input/GSA-24v3-0_A1_b151_rsids.txt")[,1] %>% filter(Name != ".") 
all_SNPs_in_GSA_2 <- fread("input/GSA-24v3-0_A1_b151_rsids.txt")[,2] %>% filter(RsID != ".")

Proxy_SNP_avialable_in_GSA <- cdf %>% 
  filter(cat != "Available on GSA array",type == "Proxy SNPs available") %>%
  filter(RS_Number %in% c(all_SNPs_in_GSA_1$Name,all_SNPs_in_GSA_2$RsID)) %>%
  filter(RS_Number != ".") %>%
  dplyr::select(-type) 



Proxy_SNP_avialable_in_GSA$Distance <- ifelse(Proxy_SNP_avialable_in_GSA$Distance >= 0,Proxy_SNP_avialable_in_GSA$Distance,-1*Proxy_SNP_avialable_in_GSA$Distance)

#create regulome score index dataframe
rank_regulome_db <- data.frame(RegulomeDB = c("1a","1b","1c","1d","1e","1f","2a","2b","2c","3a","3b","4","5","6","7","."),
                               rank = c(1:16))

# annotate regulomedb score to ldlink result in order to easy manipulation
# remove proxy SNPs without rsnumber (RS_Number != ".")
# filter Proxy with Dprime at least 0.9 
# filter regulomeDB score at lowest value (more information annotated)
# filter proxy SNPs with lowest MAF  (lowest MAF, more effect )
index_df <- inner_join(rank_regulome_db,Proxy_SNP_avialable_in_GSA,by="RegulomeDB")
biallelic_index <- index_df$Alleles %>% unique()

df_rank <- inner_join(rank_regulome_db,Proxy_SNP_avialable_in_GSA,by="RegulomeDB") %>%
  filter(RS_Number != ".",Alleles %in% biallelic_index[1:12])  %>% 
  group_by(query) %>% 
  slice_max(order_by = R2, n = 1) %>%
  slice_min(order_by = Distance, n = 1) %>%
  slice_max(order_by = MAF, n = 1) %>% 
  slice_min(order_by = rank, n = 1) %>%
  as.data.frame()

#check which SNPs still have more than one candidate 
SNPs_two_proxy <- table(df_rank$query) %>% as.data.frame() %>% filter(Freq != 1) %>% dplyr::select(Var1) %>% unlist()
df_rank$result <- ifelse(df_rank$query %in% SNPs_two_proxy,"More than one proxy SNPs","One proxy SNPs")
df_rank_final <- df_rank[c(1,4,5,6,7,9,10,13,14,15,16)]

# Combine result
cdf$result <- "Available on GSA array"
cdf_02 <- cdf %>% filter(cat =="Available on GSA array") %>% dplyr::select(query) %>% unique()
SNP_avialable_in_GSA_final <- cdf %>% filter(cat == "Available on GSA array")%>% 
  dplyr::select(colnames(df_rank_final)) %>% filter(Coord %in% cdf_02$query,
                                                    Dprime %in% c(1,0))

all_result <- rbind(df_rank_final %>% as.data.frame(),SNP_avialable_in_GSA_final%>% as.data.frame())
names(all_result) <- c("RegulomeDB","Proxy SNPs","Coord","Ref/Alt","MAF",
                       "Dprime","R2","query variant","rsid","cat","result")

write.csv(all_result,"proxy_SNPs_EAS_SAS.csv")


