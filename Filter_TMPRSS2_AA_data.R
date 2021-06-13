install.packages("tidyverse")
library(tidyverse)
library(readr)
library(scales)


#import datasets
Input_TMPRSS2 = read_tsv('Total_Input_Tmpress.txt') %>% 
  select(!X1) %>% 
  filter(!grepl('-', Variant), !grepl('S0', Variant), !grepl('B', Variant)) %>%
  na.omit() %>% rename(type = Phenotype) %>% 
  mutate(score = if_else(padj<0.05, "pass", "fail"), 
         log2FoldChange = -log2FoldChange)

Input_uPA = read_tsv('Total_Input_uPA.txt') %>% 
  select(!X1) %>% 
  filter(!grepl('-', Variant), !grepl('S0', Variant), !grepl('B', Variant)) %>%
  na.omit() %>% rename(type = Phenotype) %>% 
  mutate(score = if_else(padj<0.05, "pass", "fail"), 
         log2FoldChange = -log2FoldChange)

uPA_TMPRSS2 = read_tsv('Total_uPA_Tmpress.txt') %>% 
  select(!X1) %>% 
  filter(!grepl('-', Variant), !grepl('S0', Variant), !grepl('B', Variant)) %>%
  na.omit() %>% rename(type = Phenotype) %>% 
  mutate(score = if_else(padj<0.05, "pass", "fail"), 
         log2FoldChange = -log2FoldChange)

#Make MA Plot
ggplot()+
  scale_x_log10()+
  geom_hline(yintercept = 0, color = "black") +
  ggtitle('TMPRSS2 MA Plot')+
  theme_classic()+
  scale_shape_manual(values = c(1,19))+
  geom_jitter(data = subset(Input_TMPRSS2, type == 'Missense'),
              aes(x = baseMean, y = log2FoldChange, 
                  color = type, shape = score),
              stroke = 0.5, size = 1.2)+
  geom_jitter(data = subset(Input_TMPRSS2, type == 'Nonsense'),
              aes(x = baseMean, y = log2FoldChange, 
                  color = type, shape = score),
              stroke = 0.5, size = 1.2)+
  geom_jitter(data = subset(Input_TMPRSS2, type == 'wt'),
              aes(x = baseMean, y = log2FoldChange, 
                  color = type, shape = score),
              stroke = 0.5, size = 1.2)+
  geom_vline(xintercept = 100, color = "black", linetype = "dashed")


#Filter data, basemean = 100, padj < = 0.5
Input_TMPRSS2_filtered = Input_TMPRSS2 %>% filter(baseMean > 100, padj <= 0.05)
Input_uPA_filtered = Input_uPA %>% filter(baseMean > 100, padj <= 0.05)
uPA_TMPRSS2_filtered = uPA_TMPRSS2 %>% filter(baseMean > 100, padj <= 0.05)

write_tsv(Input_TMPRSS2_filtered,"Input_TMPRSS2_filtered.xls")
write_tsv(Input_uPA_filtered, "Input_uPA_filtered.xls")
write_tsv(uPA_TMPRSS2_filtered,"uPA_TMPRSS2_filtered.xls")

#Combine TMPRSS2 and uPA
T2_uPA_joined = full_join(Input_TMPRSS2_filtered, Input_uPA_filtered, 
                          by = "Variant") %>% na.omit() %>% 
  rename(TMPRSS2_log2 = log2FoldChange.x,
         uPA_log2 = log2FoldChange.y,
         TMPRSS2_padj = padj.x, 
         uPA_padj = padj.y,
         baseMean = baseMean.x,
         type = type.x,
         Amplicon_TMPRSS2 = Amplicon.x, Amplicon_uPA = Amplicon.y) %>% 
  select(Variant, Amplicon_TMPRSS2, Amplicon_uPA, baseMean, TMPRSS2_log2, 
         TMPRSS2_padj, uPA_log2, uPA_padj, type)
  
T2_preferred = T2_uPA_joined %>% filter(uPA_log2<=0, TMPRSS2_log2>0)
write_tsv(T2_preferred,"T2_preferred.xls")

T2_sigvsNOT = left_join(Input_TMPRSS2_filtered, Input_uPA, by = 'Variant') %>% 
  na.omit() %>% 
  rename(TMPRSS2_log2 = log2FoldChange.x,
         uPA_log2 = log2FoldChange.y,
         TMPRSS2_padj = padj.x, 
         uPA_padj = padj.y,
         baseMean = baseMean.x,
         type = type.x,
         Amplicon_TMPRSS2 = Amplicon.x, Amplicon_uPA = Amplicon.y) %>% 
  select(Variant, Amplicon_TMPRSS2, Amplicon_uPA, baseMean, TMPRSS2_log2, 
         TMPRSS2_padj, uPA_log2, uPA_padj, type) %>% 
  filter(TMPRSS2_log2>0,uPA_padj>0.05)
  
all_T2_preferred = rbind(T2_preferred,T2_sigvsNOT)

  
  