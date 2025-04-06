# required: load_abcd_data.R
# goal: PCA on Wittenberg tPRSs
# next: TODO interpret PCs

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



