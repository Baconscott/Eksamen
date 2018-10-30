library(tidyverse)
library(dplyr)

data <- readr::read_tsv("gene_expression_TPM.tsv")

data <- as_tibble(data)
View(data)

data %>% 
  gather(Behandling, TPM, -1) %>% 
  ggplot(aes(x = TPM)) +
  geom_histogram()

data %>% 
  ggplot(aes(x = `G0-1`, y = `G0-2`)) + geom_point()

plot(1:20000, log2(1:20000))

data <- bind_cols(geneID = data$geneID, log2(data[,-1]+1))
data

data.t <- data %>% 
  gather(Behandling, TPM, -geneID) %>% 
  spread(geneID, TPM) #Flipper datasette 

data.t %>% 
  separate(Behandling, into = c("Dose", "Replikat"), sep = "-") %>% 
  ggplot(aes(MA_1000049g0010, MA_10001337g0010, color = Dose)) + geom_point()


pca.res <- prcomp(data.t[,-1])

as_tibble(pca.res$x)
data.pca <- bind_cols(Behandling = data.t$Behandling, as_tibble(pca.res$x))

data.pca %>% 
  separate(Behandling, into = c("Dose", "Replikat"), sep = "-") %>%
  ggplot(aes(PC1, PC2, color = Dose)) + geom_point()

loadings <- bind_cols(geneID = row.names(pca.res$rotation), as.tibble(pca.res$rotation))

loadings %>%  
  select(geneID, PC1) %>% 
  arrange(desc(abs(PC1))) #Dette er genet som har størst innvirkinng for PC1

loadings %>%  
  select(geneID, PC2) %>% 
  arrange(desc(abs(PC2)))

data.t %>% 
  separate(Behandling, into = c("Dose", "Replikat"), sep = "-") %>% 
  ggplot(aes(MA_9438866g0010, MA_9293g0010, color = Dose)) + geom_point() #Nå har jeg plottet de to genene som er viktiskt for PC1 og PC2.


                  