library(data.table)
library(dplyr)
library(ggplot2)
library(splitstackshape)


###explore detect SNPs across reference
detected_EAS <- fread("result/final_proxy_EAS.csv",header = T)[,13] %>% unique()
detected_EAS$type <- "EAS"
detected_SAS <- fread("result/final_proxy_SAS.csv",header = T)[,13] %>% unique()
detected_SAS$type <- "SAS"
detected_EAS_SAS <- fread("result/final_proxy_EAS_SAS.csv",header = T)[,13] %>% unique()
detected_EAS_SAS$type <- "EAS_SAS"
all_detect <- rbind(detected_EAS,detected_EAS_SAS,detected_SAS) %>% select(type) %>% table() %>% melt()
all_detect$Var1 <- "Detectable SNPs"

###explore not detect SNPs across reference
not_detected_EAS <- fread("result/SNPs_not_detected_EAS.csv",header = T)
not_detected_EAS_02 <- cSplit(not_detected_EAS,"data_out[1, 1]"," ") %>% 
  select(`data_out[1, 1]_4`) %>% setnames("data_out[1, 1]_4","type")
not_detected_EAS_02$ref <- "EAS"
not_detected_SAS <- fread("result/SNPs_not_detected_SAS.csv",header = T)
not_detected_SAS_02 <- cSplit(not_detected_SAS,"data_out[1, 1]"," ") %>% 
  select(`data_out[1, 1]_4`) %>% setnames("data_out[1, 1]_4","type")
not_detected_SAS_02$ref <- "SAS"
not_detected_EAS_SAS <- fread("result/SNPs_not_detected_EAS_SAS.csv",header = T)
not_detected_EAS_SAS_02 <- cSplit(not_detected_EAS_SAS,"data_out[1, 1]"," ") %>% 
  select(`data_out[1, 1]_4`) %>% setnames("data_out[1, 1]_4","type")
not_detected_EAS_SAS_02$ref <- "EAS_SAS"

#combine
not_detect_snp_df <- rbind(table(not_detected_EAS_02$type,not_detected_EAS_02$ref) %>% as.data.frame(),
                           table(not_detected_SAS_02$type,not_detected_SAS_02$ref) %>% as.data.frame(),
                           table(not_detected_EAS_SAS_02$type,not_detected_EAS_SAS_02$ref) %>% as.data.frame())

names(all_detect) <- c("reference","Freq","type")
names(not_detect_snp_df) <- c("type","reference","Freq")

plot_df <- rbind(all_detect,not_detect_snp_df)
plot_df[plot_df == "not"] <- "Not SNPs"
plot_df[plot_df == "monoallelic"] <- "Monoallelic SNPs"

#plot
ggplot(plot_df,aes(reference,Freq,color=type,fill=type)) +
  geom_bar(stat="identity",position="dodge") +
  labs(color="",fill="")+
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25,color="black") +
  theme_bw()+
  ylab("Frequency")

