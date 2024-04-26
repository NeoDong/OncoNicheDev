library(ggplot2)
library(cowplot)
library(ggpubr)
library(fitdistrplus)

obs <- read.csv("observed_info_change.csv")
rnd <- read.csv("random_info_change.csv")

norm_fit_p <- function(observed_value, random_dist){
  fd <- fitdist(random_dist, "norm") # fit random distribution with normal distribution
  # denscomp(list(fd), legendtext = c("normal"))
  e_mean <- fd$estimate["mean"] # estimated mean
  e_sd <- fd$estimate["sd"]# estimated sd
  
  # calculate 
  if(observed_value <= 0){
    e_p <- pnorm(observed_value, mean = e_mean, sd = e_sd)
  }else{
    e_p <- pnorm(observed_value, mean = e_mean, sd = e_sd, lower.tail=FALSE)
  }
  
  return(e_p)
}

targets <- colnames(rnd)
estimated_p_values <- c()
for (gene in obs$Target){
  obs_value <- obs[obs$Target==gene, ]$Delta_S
  gene <- gsub("-", ".", gene) # some "-"  in transcript factor name are translated into "." by R
  idx <- which(targets==gene)
  e_p  <- norm_fit_p(obs_value, rnd[, idx])
  estimated_p_values <- c(estimated_p_values, e_p)
}

obs$Estimated_P_Value <- estimated_p_values
obs$padj<- p.adjust(obs$Estimated_P_Value, method = "BH")

write.csv(obs, file="fitted_observed_info_change.csv", row.names = FALSE, quote = FALSE)

