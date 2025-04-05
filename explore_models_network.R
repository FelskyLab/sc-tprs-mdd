load("sc-tprs-mdd/og_white_all_tprs_immune.Rda")
scale_tprs <- og_white_all_tprs
for (event in events){
  scale_tprs[(scale_tprs$event_name == event),][, -c(1, 2, 3, 4, 6, 7, 8, 
                                                     (tprs_start_idx+2*n_ct), (tprs_start_idx+2*n_ct+1), 
                                                     (tprs_start_idx+2*n_ct+2), (tprs_start_idx+2*n_ct+3))] <- 
    scale(og_white_all_tprs[(og_white_all_tprs$event_name == event),][, -c(1, 2, 3, 4, 6, 7, 8, 
                                                                           (tprs_start_idx+2*n_ct), (tprs_start_idx+2*n_ct+1), 
                                                                           (tprs_start_idx+2*n_ct+2), (tprs_start_idx+2*n_ct+3))])
  
}
n_ct <- 22
tprs_start_idx <- 19
# PCA
pca_tprs_uw <- prcomp(t(scale_tprs[
  which(scale_tprs$event_name=="baseline_year_1_arm_1"), 
  tprs_start_idx:(n_ct + tprs_start_idx - 1)]))
pca_tprs_w <- prcomp(t(scale_tprs[
  which(scale_tprs$event_name=="baseline_year_1_arm_1"), 
  (n_ct + tprs_start_idx):(2*n_ct + tprs_start_idx - 1)]))
plot(pca_tprs_uw)
summary(pca_tprs_uw) #6 pcs to get 81% of var

# Plot projections along the components into a scatterplot.
# Axes for points are scaled as values, for vectors as variance
# Default for biplot() is the first and second component.

biplot(pca_tprs_uw)

plot(pca_tprs_w)
summary(pca_tprs_w)  #2 pcs to get 88% of var
biplot(pca_tprs_w)

# Examine the actual principal components in a parallel-coordinates

N <- 4
matplot(pca_tprs_uw$rotation[1:N, ],
        type="b", lwd=1,
        xlab = "ABCD participants", ylab="PCs")

matplot(pca_tprs_w$rotation[1:N, ],
type="b", lwd=1,
xlab = "Cell types", ylab="PCs")
pca_tprs_w$rotation <- t(pca_tprs_w$rotation)
cor_pc_w <- matrix(numeric(n_ct * 4), ncol = 4)
PC1 <- pca_tprs_w$rotation[ , 1]
PC2 <- pca_tprs_w$rotation[ , 2]
PC3 <- pca_tprs_w$rotation[ , 3]
PC4 <- pca_tprs_w$rotation[ , 4]

for (i in 1:n_ct) {
  cor_pc_w[i, 1] <- cor(scale_tprs[
      which(scale_tprs$event_name=="baseline_year_1_arm_1"), i + tprs_start_idx - 1], PC1)
  cor_pc_w[i, 2] <- cor(scale_tprs[
    which(scale_tprs$event_name=="baseline_year_1_arm_1"), i + tprs_start_idx - 1], PC2)
  cor_pc_w[i, 3] <- cor(scale_tprs[
    which(scale_tprs$event_name=="baseline_year_1_arm_1"), i + tprs_start_idx - 1], PC3)
  cor_pc_w[i, 4] <- cor(scale_tprs[
    which(scale_tprs$event_name=="baseline_year_1_arm_1"), i + tprs_start_idx - 1], PC4)
}
summary(cor_pc_w)

hist(cor_pc_w[ , 1])

# A function to find the best correlation:

bestCor <- function(V, dat = scale_tprs[
  which(scale_tprs$event_name=="baseline_year_1_arm_1"),], tprs_start = tprs_start_idx,
  tprs_end = tprs_start_idx + n_ct - 1) {
  # Identify that column in dataframe dat that has the highest correlation
  # or anticorrelation with vector V. Print some information about it.
  # Value: The index of the column (invisibly)
  #
  # ToDo: identify and return n-best correlations
  #
  cors <- numeric(tprs_end - tprs_start + 1)
  for (i in tprs_start:tprs_end) {
    cors[i-tprs_start+1] <- cor(V, dat[ , i]) # correlation between input vector and each
  }                              # column of the feature set in turn
  
  myBest <- max(abs(cors), na.rm = TRUE)
  sel <- abs(cors) == myBest   # Is TRUE for the highest absolute value in the
  # vector of correlation coefficients.
  idx <- which(sel)[1]       # which() finds the index of the TRUE value(s).
  # Pick only the first in case there are ties.
  print(sprintf("\n Highest correlation for index %d (%f): %s ",
              idx,
              cors[idx],
              colnames(dat)[idx+tprs_start-1]))
  return(idx)     # return the index.
}

bestCor(pca_tprs_uw$rotation[ , 1])  # This is PC1 ...
bestCor(pca_tprs_uw$rotation[ , 2])  # PC2
bestCor(pca_tprs_uw$rotation[ , 3])  # etc.
bestCor(pca_tprs_uw$rotation[ , 4])
bestCor(pca_tprs_uw$rotation[ , 5])  
bestCor(pca_tprs_uw$rotation[ , 6]) 
bestCor(pca_tprs_w$rotation[ , 1], tprs_start = tprs_start_idx + n_ct,
        tprs_end = tprs_start_idx + 2*n_ct - 1)  # etc.
bestCor(pca_tprs_w$rotation[ , 2], tprs_start = tprs_start_idx + n_ct,
        tprs_end = tprs_start_idx + 2*n_ct - 1)
bestCor(pca_tprs_w$rotation[ , 3], tprs_start = tprs_start_idx + n_ct,
        tprs_end = tprs_start_idx + 2*n_ct - 1)  # This is PC1 ...
bestCor(pca_tprs_w$rotation[ , 4], tprs_start = tprs_start_idx + n_ct,
        tprs_end = tprs_start_idx + 2*n_ct - 1)  # PC2
bestCor(pca_tprs_w$rotation[ , 5], tprs_start = tprs_start_idx + n_ct,
        tprs_end = tprs_start_idx + 2*n_ct - 1)  # etc.
bestCor(pca_tprs_w$rotation[ , 6], tprs_start = tprs_start_idx + n_ct,
        tprs_end = tprs_start_idx + 2*n_ct - 1)
pc_uw <- (pca_tprs_uw$rotation[, 1:6])
pc_w <- (pca_tprs_w$rotation[ ,1:6])
colnames(pc_uw) <- paste0(colnames(pca_tprs_uw$rotation[ ,1:6]), "_uw")
colnames(pc_w) <- paste0(colnames(pca_tprs_w$rotation[, 1:6]), "_w")
pc_df <- cbind(scale_tprs[
  which(scale_tprs$event_name=="baseline_year_1_arm_1"), c(1,
  tprs_start_idx:(2*n_ct + tprs_start_idx - 1))], pc_uw) 
pc_df <- cbind(pc_df,pc_w)
pc_df$anxdep <- base_bin[(base_bin$event_name == "baseline_year_1_arm_1"),]$anxdep_t
pc_df$withdep <- base_bin[(base_bin$event_name == "baseline_year_1_arm_1"),]$withdep_t
pc_df$anx_bin <- as.factor(as.integer(pc_df$anxdep > 65))
pc_df$with_bin <- as.factor(as.integer(pc_df$withdep > 65))
pc_df$anx_sqrt <- sqrt(pc_df$anxdep)
pc_df$with_sqrt <- sqrt(pc_df$anxdep)
pc_df$sex <- base_bin[(base_bin$event_name == "baseline_year_1_arm_1"),]$sex
colMeans(pc_df[which(pc_df$PC1_w > -0.02) , c(2:n_ct*2)])
colMeans(pc_df[which(pc_df$PC1_w < -0.02) , c(2:n_ct*2)])
colMeans((pc_df[ , c( "anxdep", "withdep")][which(pc_df$PC1_w> -0.02),]))
colMeans((pc_df[ , c("anxdep", "withdep")][which(pc_df$PC1_w < -0.02),]))
tabulate((pc_df[ , c( "sex")][which(pc_df$PC1_w> -0.02)]))/sum(pc_df$PC1_w> -0.02)
tabulate((pc_df[ , c("sex")][which(pc_df$PC1_w < -0.02)]))/sum(pc_df$PC1_w < -0.02)


plot(PC1, PC2)  # uninformative

ggplot(pc_df, aes(x = PC1_uw, y = PC2_uw)) +
  geom_point(aes(color = anxdep), alpha = 0.3)+
  scale_color_viridis_c()  
ggplot(pc_df[(which(pc_df$PC1_uw > -0.02)),], aes(x = PC1_uw, y = PC2_uw)) +
  geom_point(aes(color = naive_B_cell_uw), alpha = 0.3)+
  scale_color_viridis_c()  
ggplot(pc_df, aes(x = PC1_w, y = PC2_w, color = sex)) +
  geom_point( alpha = 0.3)
ggplot(pc_df[(which(pc_df$sex == 1)),], aes(x = PC1_uw, y = PC2_uw)) +
  geom_point(aes(color = anxdep), alpha = 0.3)+
  scale_color_viridis_c()  
ggplot(pc_df[(which(pc_df$sex == 2)),], aes(x = PC1_uw, y = PC2_uw)) +
  geom_point(aes(color = anxdep), alpha = 0.3)+
  scale_color_viridis_c()  

ggplot(pc_df, aes(x = PC1_uw, y = PC2_uw, color = with_bin)) +
  geom_point(alpha = 0.3)
ggplot(pc_df, aes(x = PC1_w, y = PC2_w, color = anx_bin)) +
  geom_point(alpha = 0.3) + 
  labs(title = "Weighted tPRS PC2 vs. PC1", x = "PC1", y = "PC2", color = "Borderline/clinically significant\nanxious/depressed score?") +
  scale_color_manual(labels = c("no", "yes"), values = c("blue", "red"))
ggplot(pc_df[(which(pc_df$sex == 1)),], aes(x = PC1_w, y = PC2_w, color = anx_bin)) +
  geom_point(alpha = 0.3) + 
  labs(title = "Weighted tPRS PC2 vs. PC1 in males", x = "PC1", y = "PC2", color = "Borderline/clinically significant\nanxious/depressed score?") +
  scale_color_manual(labels = c("no", "yes"), values = c("blue", "red"))
ggplot(pc_df[(which(pc_df$sex == 2)),], aes(x = PC1_w, y = PC2_w, color = anx_bin)) +
  geom_point(alpha = 0.3) + 
  labs(title = "Weighted tPRS PC2 vs. PC1 in females", x = "PC1", y = "PC2", color = "Borderline/clinically significant\nanxious/depressed score?") +
  scale_color_manual(labels = c("no", "yes"), values = c("blue", "red"))
ggplot(pc_df, aes(x = PC1_w, y = PC2_w)) +
  geom_point(aes(color = `CD4-positive_alpha-beta_T_cell_w`), alpha = 0.3)+
  scale_color_viridis_c() 
ggplot(pc_df, aes(x = PC1_w, y = PC2_w)) +
  geom_point(aes(color = double_negative_thymocyte_w), alpha = 0.3)+
  scale_color_viridis_c() 
ggplot(pc_df, aes(x = PC1_uw, y = PC2_uw)) +
  geom_point(aes(color = naive_B_cell_uw), alpha = 0.3)+
  scale_color_viridis_c() 
ggplot(pc_df, aes(x = PC1_uw, y = PC2_uw)) +
  geom_point(aes(color = double_negative_thymocyte_uw), alpha = 0.3)+
  scale_color_viridis_c() 


ggplot(pc_df, aes(x = PC2_uw, y = PC3_uw)) +
  geom_point(aes(color = anxdep), alpha = 0.3)+
  scale_color_viridis_c()  
ggplot(pc_df, aes(x = PC2_w, y = PC3_w)) +
  geom_point(aes(color = anxdep), alpha = 0.3)+
  scale_color_viridis_c()  
ggplot(pc_df, aes(x = PC2_w, y = PC3_w, color = sex)) +
  geom_point( alpha = 0.3)
ggplot(pc_df[(which(pc_df$sex == 1)),], aes(x = PC2_w, y = PC3_w)) +
  geom_point(aes(color = anxdep), alpha = 0.3)+
  scale_color_viridis_c()  
ggplot(pc_df[(which(pc_df$sex == 2)),], aes(x = PC2_w, y = PC3_w)) +
  geom_point(aes(color = anxdep), alpha = 0.3)+
  scale_color_viridis_c()  
ggplot(pc_df, aes(x = PC2_w, y = PC3_w, color = with_bin)) +
  geom_point(alpha = 0.3)
ggplot(pc_df, aes(x = PC2_w, y = PC3_w, color = anx_bin)) +
  geom_point(alpha = 0.3)

ggplot(pc_df, aes(x = PC3_w, y = PC4_w, color = anx_bin)) +
  geom_point(alpha = 0.3)
ggplot(pc_df, aes(x = PC3_w, y = PC4_w, color = sex)) +
  geom_point( alpha = 0.3)


# network to show overlap of sc-tPRSs
# Create a weighted graph

edges <- data.frame(
  from = sample(1:10, 20, replace = TRUE),
  to = sample(1:10, 20, replace = TRUE),
  weight = runif(20)
)
g <- graph_from_data_frame(edges, directed = FALSE)

# Plot the weighted graph
plot(g,
     edge.width = E(g)$weight * 8, # Edge width based on weight
     vertex.size = 20,
     vertex.color = "lightblue",
     vertex.label.color = "black",
     edge.color = "gray50",
     layout = layout_in_circle(g)
)



# histograms of every model's weight distr


load("sc-tprs-mdd/sc_tPRS_betas_immune.Rda")

sort_betas <- cbind(tprs_beta[, 2:(n_ct + 1)][,rev(order(tprs_beta[3,2:(n_ct + 1)]))],
                    tprs_beta[, 2:(n_ct + 1)+n_ct][,rev(order(tprs_beta[3,2:(n_ct + 1)+n_ct]))])

# by beta > 0

ggplot(stack( (scale_tprs[(colnames(sort_betas)[2:8 - 1 ])]), varying = colnames(sort_betas)[2:8 -1]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (scale_tprs[(colnames(sort_betas)[9:15 -1])]), varying = colnames(sort_betas)[9:15-1]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (scale_tprs[(colnames(sort_betas)[16:(n_ct+1)-1])]), varying = colnames(sort_betas)[16:(n_ct+1)-1]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))

normal_range <- c(-1, 1)
low_range <- c(-10,-2)
range_tprs<- scale_tprs[, c(19:(18+2*n_ct))]
range_tprs[,1:ncol(range_tprs)] <- NA
range_tprs[scale_tprs[, c(19:(18+2*n_ct))] > normal_range[1] & scale_tprs[, c(19:(18+2*n_ct))] < normal_range[2]] <- "normal"
range_tprs[scale_tprs[, c(19:(18+2*n_ct))] > low_range[1] & scale_tprs[, c(19:(18+2*n_ct))] < low_range[2]] <- "low"
if (!requireNamespace("RSQLite", quietly = TRUE)){
  install.packages("RSQLite")
}
library(RSQLite)
if (!requireNamespace("DBI", quietly = TRUE)){
  install.packages("DBI")
}
library(DBI)


model_files <- list.files(path="abcd/workspace/mding/immune_cell_models", pattern="*.db", full.names=TRUE, recursive=FALSE)
models <- list()
read_model <- function(filename){
  fn <- sub('\\.db$', '', basename(filename))
  sqlite.driver <- dbDriver("SQLite")
  db <- dbConnect(sqlite.driver, dbname = filename)
  
  model <- dbReadTable(db,"weights")
  dbDisconnect(db)
  #format chr and position
  #model_chr <- as.numeric(str_remove(str_extract(model$varID, "chr[0-9]+"), "chr"))
  #model_pos <- as.numeric(str_remove(str_extract(model$varID, "_[0-9]+"), "_"))
  
  models[[fn]] <<- data.frame(SNP = model$varID, 
                              ZSCORE = scale(model$weight), gene = model$gene)
}

sapply(model_files, read_model)
if (!requireNamespace("qqman", quietly = TRUE)){
  install.packages("qqman")
}
library(qqman)
normal_idx <- c(1,6,8,17)
mid_idx <- c(3, 13,16)
low_idx <- c(15,19,21)

for (i in 1:n_ct){
  print(paste(names(models)[i], length(unique(models[[i]]$SNP)), length(unique(models[[i]]$gene))))
}


manhattan(models[[normal_idx[1]]], ylim = c(-0.015, 0.015), logp = FALSE, p = "P", ylab = "Weight", genomewideline = FALSE,
                   suggestiveline = FALSE, main = paste0(names(models)[normal_idx[1]]," model weights"))
 manhattan(models[[normal_idx[2]]], ylim = c(-0.015, 0.015), logp = FALSE, p = "P", ylab = "Weight", genomewideline = FALSE,
                       suggestiveline = FALSE, main = paste0(names(models)[normal_idx[2]]," model weights"))
 manhattan(models[[normal_idx[3]]], ylim = c(-0.015, 0.015), logp = FALSE, p = "P", ylab = "Weight", genomewideline = FALSE,
                       suggestiveline = FALSE, main = paste0(names(models)[normal_idx[3]]," model weights"))
 manhattan(models[[normal_idx[4]]], ylim = c(-0.015, 0.015), logp = FALSE, p = "P", ylab = "Weight", genomewideline = FALSE,
                       suggestiveline = FALSE, main = paste0(names(models)[normal_idx[4]]," model weights"))

 ggplot(stack( (og_white_all_tprs[(colnames(tprs_beta)[1+c(normal_idx[1], normal_idx[1]+n_ct)])]),
                               varying = colnames(tprs_beta)[c(1+normal_idx[1], normal_idx[1]+n_ct)]), 
                aes(values)) +
   geom_histogram(bins = 20) + 
   facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
   theme(text=element_text(size=9))
 ggplot(stack( (og_white_all_tprs[(colnames(tprs_beta)[1+c(normal_idx[2], normal_idx[2]+n_ct)])]),
                             varying = colnames(tprs_beta)[1+c(normal_idx[2], normal_idx[2]+n_ct)]), 
               aes(values)) + 
   geom_histogram(bins = 20) + 
   facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
   theme(text=element_text(size=9))
 ggplot(stack( (og_white_all_tprs[(colnames(tprs_beta)[1+c(normal_idx[3], normal_idx[3]+n_ct)])]),
                             varying = colnames(tprs_beta)[1+c(normal_idx[3], normal_idx[3]+n_ct)]), 
               aes(values)) + 
   geom_histogram(bins = 20) + 
   facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
   theme(text=element_text(size=9))
 ggplot(stack( (og_white_all_tprs[(colnames(tprs_beta)[1+c(normal_idx[4], normal_idx[4]+n_ct)])]),
                             varying = colnames(tprs_beta)[1+c(normal_idx[4], normal_idx[4]+n_ct)]), 
               aes(values)) + 
   geom_histogram(bins = 20) + 
   facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
   theme(text=element_text(size=9))
 manhattan(models[[mid_idx[1]]], ylim = c(-0.015, 0.015), logp = FALSE, p = "P", ylab = "Weight", genomewideline = FALSE,
                     suggestiveline = FALSE, main = paste0(names(models)[mid_idx[1]]," model weights"))
 manhattan(models[[mid_idx[2]]], ylim = c(-0.015, 0.015), logp = FALSE, p = "P", ylab = "Weight", genomewideline = FALSE,
                     suggestiveline = FALSE, main = paste0(names(models)[mid_idx[2]]," model weights"))
 manhattan(models[[mid_idx[3]]], ylim = c(-0.015, 0.015), logp = FALSE, p = "P", ylab = "Weight", genomewideline = FALSE,
                     suggestiveline = FALSE, main = paste0(names(models)[mid_idx[3]]," model weights"))
 
   ggplot(stack( (og_white_all_tprs[(colnames(tprs_beta)[1+c(mid_idx[1], mid_idx[1]+n_ct)])]),
                               varying = colnames(tprs_beta)[c(1+mid_idx[1], mid_idx[1]+n_ct)]), 
                 aes(values)) + 
   geom_histogram(bins = 20) + 
   facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
   theme(text=element_text(size=9))
 ggplot(stack( (og_white_all_tprs[(colnames(tprs_beta)[1+c(mid_idx[2], mid_idx[2]+n_ct)])]),
                             varying = colnames(tprs_beta)[1+c(mid_idx[2], mid_idx[2]+n_ct)]), 
               aes(values)) + 
   geom_histogram(bins = 20) + 
   facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
   theme(text=element_text(size=9))
 ggplot(stack( (og_white_all_tprs[(colnames(tprs_beta)[1+c(mid_idx[3], mid_idx[3]+n_ct)])]),
                             varying = colnames(tprs_beta)[1+c(mid_idx[3], mid_idx[3]+n_ct)]), 
               aes(values)) + 
   geom_histogram(bins = 20) + 
   facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
   theme(text=element_text(size=9))
 manhattan(models[[low_idx[1]]], ylim = c(-0.015, 0.015), logp = FALSE, p = "P", ylab = "Weight", genomewideline = FALSE,
                     suggestiveline = FALSE, main = paste0(names(models)[low_idx[1]]," model weights"))
 manhattan(models[[low_idx[2]]], ylim = c(-0.015, 0.015), logp = FALSE, p = "P", ylab = "Weight", genomewideline = FALSE,
                     suggestiveline = FALSE, main = paste0(names(models)[low_idx[2]]," model weights"))
 manhattan(models[[low_idx[3]]], ylim = c(-0.015, 0.015), logp = FALSE, p = "P", ylab = "Weight", genomewideline = FALSE,
                     suggestiveline = FALSE, main = paste0(names(models)[low_idx[3]]," model weights"))
 
   ggplot(stack( (og_white_all_tprs[(colnames(tprs_beta)[1+c(low_idx[1], low_idx[1]+n_ct)])]),
                               varying = colnames(tprs_beta)[c(1+low_idx[1], low_idx[1]+n_ct)]), 
                 aes(values)) + 
   geom_histogram(bins = 20) + 
   facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
   theme(text=element_text(size=9))
 ggplot(stack( (og_white_all_tprs[(colnames(tprs_beta)[1+c(low_idx[2], low_idx[2]+n_ct)])]),
                             varying = colnames(tprs_beta)[1+c(low_idx[2], low_idx[2]+n_ct)]), 
               aes(values)) + 
   geom_histogram(bins = 20) + 
   facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
   theme(text=element_text(size=9))
 ggplot(stack( (og_white_all_tprs[(colnames(tprs_beta)[1+c(low_idx[3], low_idx[3]+n_ct)])]),
                             varying = colnames(tprs_beta)[1+c(low_idx[3], low_idx[3]+n_ct)]), 
               aes(values)) + 
   geom_histogram(bins = 20) + 
   facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
   theme(text=element_text(size=9))
 ggplot(stack( (scale_tprs[(colnames(sort_betas)[2:8 - 1 ])]), varying = colnames(sort_betas)[2:8 -1]), 
               aes(values)) + 
   geom_histogram(bins = 20) + 
   facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
   theme(text=element_text(size=9))
 ggplot(stack( (scale_tprs[(colnames(sort_betas)[9:15 -1])]), varying = colnames(sort_betas)[9:15-1]), 
               aes(values)) + 
   geom_histogram(bins = 20) + 
   facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
   theme(text=element_text(size=9))
 ggplot(stack( (scale_tprs[(colnames(sort_betas)[16:(n_ct+1)-1])]), varying = colnames(sort_betas)[16:(n_ct+1)-1]), 
               aes(values)) + 
   geom_histogram(bins = 20) + 
   facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
   theme(text=element_text(size=9))

 
tPRS_gene_diff <- function(filename, tprss){
   #filename <- "/nethome/kcni/mding/abcd/workspace/mding/abcd_imputed/with_sc_blood38/transitional_stage_B_cell"
   df <- fread(paste0(filename,"/abcd_pred_", 1, ".txt"), header = TRUE)
   df$IID <- mapply(sub, ".{9}_?", "", df$IID)
   overlap_genes <- intersect(witt_deg$ensembl_gene_id, colnames(df))
   df_sub <- df[(df$IID %in%  unique(tprss$IID)), ]
   df_sub <- df_sub[, c("FID", "IID", ..overlap_genes)]
   df_sub <- merge(df_sub, 
                        tprss[,c("IID", paste0(basename(filename), "_", "uw"), "PC1_uw")], 
                        by = c("IID"), sort = FALSE)
   g_diff <- vector()
   for (k in 3:(ncol(df_sub)-2)){
     a <- ncol(df_sub) - 1
     df_sub[,k] <- scale(df_sub[,..k])
     if (!is.na(df_sub[,..k][1,])){
       tt <- t.test(df_sub[(which(df_sub$PC1_uw > -0.02)),..k],
                    df_sub[(which(df_sub$PC1_uw < -0.02)),..k])
       if (!is.na(tt$p.value) & tt$p.value <= 0.05 ){
         g_diff[colnames(df_sub)[k]] <- tt$estimate[1] - tt$estimate[2]
       }
     }
   }
   return(g_diff)
}
 
ct_gene_diff_uw <- list()
for(c in cell_files){
  ct_gene_diff_uw[[paste0(basename(c), "_uw")]] <- tPRS_gene_diff(c, pc_df)
  
}

for (c in names(ct_gene_diff_uw)){
  ct_gene_diff_uw[[c]] <- ct_gene_diff_uw[[c]][order(abs(ct_gene_diff_uw[[c]]), decreasing = TRUE)]
}
lapply(ct_gene_diff_uw, head)
# install.packages("vroom")
# library(vroom)
vcf_id <- fread("/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_data/abcd_snps.txt", 
                header = TRUE, skip = 30)

wgt_snp_df <- data.frame(ct = c(), snp = c(), wgt = c(), CHR = c(), BP = c())
for(mod in names(models)){
  genes <- names(ct_gene_diff_uw[[paste0(mod,"_uw")]])
  idx <- which(models[[mod]]$gene %in% genes)
  idx2 <- which(models[[mod]]$SNP[idx] %in% vcf_id$ID)
  #print(idx2)
  df <- data.frame(ct = rep((paste0(mod,"_uw")), length(idx2)), 
                   SNP = models[[mod]]$SNP[idx][idx2], wgt = models[[mod]]$P[idx][idx2],
                   quant = ecdf(abs(models[[mod]]$P))(abs(models[[mod]]$P[idx][idx2])))
  model_chr <- as.numeric(str_remove(str_extract(df$SNP, "chr[0-9]+"), "chr"))
  model_pos <- as.numeric(str_remove(str_extract(df$SNP, "_[0-9]+"), "_"))
  df$CHR <- model_chr
  df$BP <- model_pos
  wgt_snp_df <- rbind(wgt_snp_df, df)
}

manhattan(wgt_snp_df, ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = "Aggregate model weights for SNPs that control expression\nlevel of influential genes")


balance_gene_expr <- unlist(lapply(ct_gene_diff_uw, sum))
balance_gene_expr <- balance_gene_expr[order(abs(balance_gene_expr), decreasing = TRUE)]
sort(table(wgt_snp_df$SNP),decreasing=TRUE)[1:52]
sort(table(names(unlist(unname(ct_gene_diff_uw)))),decreasing=TRUE)[1:59]
 fileConn<-file("infl_genes.txt")
 writeLines((names(unlist(unname(ct_gene_diff_uw)))), fileConn)
 close(fileConn)
ct_labs <- unique(wgt_snp_df$ct)
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[1]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[1], " model weights for SNPs that control expression\nlevel of influential genes"),
          highlight = c("chr12_47902703_A_G_b38"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[2]),], ylim = c(-0.001, 0.001), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, 
          main = paste0(ct_labs[2], " model weights for SNPs that control expression\nlevel of influential genes"),
          highlight = c("chr20_35291403_C_T_b38"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[3]),], ylim = c(-0.0005, 0.0005), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[3], " model weights for SNPs that control expression\nlevel of influential genes"),
          highlight = c("chr3_14978842_A_G_b38"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[4]),], ylim = c(-0.001, 0.001), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[4], " model weights for SNPs that control expression\nlevel of influential genes"),
          highlight = c("chr15_89091211_A_G_b38", "chr20_35291675_A_G_b38", "chr11_75398824_T_C_b38"), xlim = c(98000000, 150000000))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[5]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[5], " model weights for SNPs that control expression\nlevel of influential genes"),
          highlight = c("chr7_90656015_G_A_b38") , xlim = c(140000000, 190000000))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[6]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[6], " model weights for SNPs that control expression\nlevel of influential genes"),
          highlight =c("chr12_47902703_A_G_b38"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[7]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[7], " model weights for SNPs that control expression\nlevel of influential genes"),
          highlight = c("chr12_47902703_A_G_b38", "chr7_90656015_G_A_b38"), xlim = c(194000000, 400000000))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[8]),], ylim = c(-0.001, 0.001), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[8], " model weights for SNPs that control expression\nlevel of influential genes"),
          xlim = c(220000000, 270000000))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[9]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[9], " model weights for SNPs that control expression\nlevel of influential genes"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[10]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[10], " model weights for SNPs that control expression\nlevel of influential genes"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[11]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[11], " model weights for SNPs that control expression\nlevel of influential genes"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[12]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[12], " model weights for SNPs that control expression\nlevel of influential genes"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[13]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[13], " model weights for SNPs that control expression\nlevel of influential genes"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[14]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[14], " model weights for SNPs that control expression\nlevel of influential genes"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[15]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[15], " model weights for SNPs that control expression\nlevel of influential genes"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[16]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[16], " model weights for SNPs that control expression\nlevel of influential genes"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[17]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[17], " model weights for SNPs that control expression\nlevel of influential genes"),
          highlight = c("chr12_47902703_A_G_b38"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[18]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[18], " model weights for SNPs that control expression\nlevel of influential genes"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[19]),], ylim = c(-0.001, 0.001), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[19], " model weights for SNPs that control expression\nlevel of influential genes"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[20]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[20], " model weights for SNPs that control expression\nlevel of influential genes"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[21]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[21], " model weights for SNPs that control expression\nlevel of influential genes"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[22]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[22], " model weights for SNPs that control expression\nlevel of influential genes"),
          highlight = c("chr7_90656015_G_A_b38"))
iid <- data.frame(IID = scale_tprs[,2])
colnames(scale_tprs) <- gsub("-", "_", colnames(scale_tprs))
write.table(cbind(iid, scale_tprs)[(scale_tprs$event_name == events[1]),][, c(1,20:41)], 
            "/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_data/tprs_pheno.txt",
            sep="\t",row.names=FALSE, quote = F)

write.table(cbind(iid, scale_tprs)[(scale_tprs$event_name == events[1]),][, c(1,10:19)], 
            "/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_data/tprs_pheno_genpc.txt",
            sep="\t",row.names=FALSE, quote = F)

save(wgt_snp_df,file="sc-tprs-mdd/model_wgt_snp_df.Rda")
save(ct_gene_diff_uw,file="sc-tprs-mdd/ct_gene_diff_uw.Rda")

for (ct in names(ct_gene_diff_uw)){
  write.table(wgt_snp_df[which(wgt_snp_df$ct == ct), "SNP"], 
              paste0("/external/rprshnas01/kcni/mding/abcd/workspace/mding/immune_cell_models/snps/", ct,"_from_gene.txt"),
              sep="\t",row.names=FALSE, quote = F)
  
}