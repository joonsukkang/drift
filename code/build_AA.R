library(tidyverse)
library(here)

all_phase3_ns_maf0_05_pvar <- read_delim("~/Documents/data/1kgp/all_phase3_ns_maf0.05.pvar", 
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE, skip = 189)
INFO <- all_phase3_ns_maf0_05_pvar$INFO
loc <- gregexpr(patter='AA=', INFO)
AA <- substr(INFO, unlist(loc)+3, unlist(loc)+3)
all_phase3_ns_maf0_05_pvar %>%
  select(-QUAL, -FILTER, -INFO) %>%
  cbind(data.frame(AA=AA)) -> df.out

colnames(df.out)[1] <-'CHROM'
df.out %>%
  mutate(ID2 = paste0(CHROM, ':',POS)) -> df.out


# data from the 1000 Genomes Project paper
df2 <- read_delim("~/Documents/data/1kgp/ALL.wgs.phase3_shapeit2_filtered.20141217.maf0.05.pvar", 
                  delim = "\t", escape_double = FALSE, 
                  col_types = cols(ID = col_character()), 
                  trim_ws = TRUE)
df2 %>%
  select(ID, REF, ALT) %>%
  rename(ID2 = ID, REF2 = REF, ALT2 = ALT) -> df2

df2 %>%left_join(df.out, by='ID2') -> df.out2


saveRDS(df.out2, here('data', 'AA.rds'))
