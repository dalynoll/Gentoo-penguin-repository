setwd("path/to/dir")

## 1. Read data
morf <- read.table("morpho.tsv", 
                   header = TRUE, 
                   sep = "\t", 
                   stringsAsFactors = FALSE)


# 2. Exploration data
dim(morf)        
head(morf)       
str(morf)        

morf$Locality <- factor(morf$Locality)
morf$Lineage  <- factor(morf$Lineage)

vars_morf <- c("CL", "BWB", "BH", "BWG", "FW", "RL", "ML", "TML", "MTL")
# culmen length (CL; taken along the medial line)
# bill width at the base (BWB)
# bill height at gonys angle (BH) 
# bill width at gonys angle (BWG)
# flipper width (FW; shortest distance from anterior surface of flipper above the radiale to the posterior side of the flipper), 
# radius length (RL), 
# manus length (ML; indent at radiale/radius/ulna to distal flipper tip), 
# tarsus length (TML; anterior surface) and 
# middle toe length (MTL; digit I11 excluding nail)


# Check NA values
for (v in vars_morf) {
  if (any(morf[[v]] <= 0, na.rm = TRUE)) {
    stop(paste("Variable", v, "has values <= 0. Check before log-transformation."))
  }
}


# 2.1. Create log-transformed columns 
for (v in vars_morf) {
  morf[[paste0("log_", v)]] <- log(morf[[v]])  
}


# 2.2. Extract matrix with log-transformed variables
morf_log <- morf[, paste0("log_", vars_morf)]
head(morf_log)



# 3.1. PCA with prcomp 
pca_res <- prcomp(morf_log, center = TRUE, scale. = FALSE)

# 3.2. Summary
summary(pca_res)

## Importance of components:
##                          PC1    PC2     PC3     PC4     PC5     PC6     PC7     PC8
## Standard deviation     0.2292 0.1092 0.09355 0.07596 0.05176 0.04735 0.03428 0.03345
## Proportion of Variance 0.6041 0.1371 0.10062 0.06633 0.03079 0.02578 0.01351 0.01287
## Cumulative Proportion  0.6041 0.7412 0.84181 0.90814 0.93893 0.96471 0.97822 0.99109


pca_log <- prcomp(morf_log, center = TRUE, scale. = FALSE)

library(ggplot2)
library(dplyr)

pca_data <- data.frame(
  PC1      = pca_log$x[, 1],
  PC2      = pca_log$x[, 2],
  Locality = morf$Locality,
  Lineage  = morf$Lineage
)

pca_data_hull <- pca_data %>%
  group_by(Lineage) %>%
  slice(chull(PC1, PC2))  

# 2) Define colors
my_lineage_colors <- c(
  "Southern"       = "yellow",
  "South_Georgia"  = "orange",
  "Northern"       = "red",
  "Eastern"        = "blue",
  "Macquarie"      = "#7C217F",
  "Southeastern"   = "green"
)

my_locality_cols <- c(
  "Falklands" = "#DB153D",  
  "Kerguelen" = "#42AF71",
  "Crozet" = "#3267c0",
  "Macquarie" = "#b53bb9", 
  "Antarctic_Peninsula" = "#FFE05F",  
  "South_Georgia" = "#DD580E",  
  "South_Shetland" = "#EEE58A",  
  "South_Orkney" = "#E8B646"   
)

my_locality_shapes <- c(
  "Falklands" = 8,
  "Kerguelen" = 4,
  "Crozet" = 19,
  "Macquarie" = 19, 
  "Antarctic_Peninsula" = 0, 
  "South_Georgia" = 2,
  "South_Shetland" = 0, 
  "South_Orkney" = 0
)


# 3) Plot
mi_plot <- ggplot() +
  geom_polygon(
    data = pca_data_hull,
    aes(x = PC1, y = PC2, fill = Lineage, group = Lineage),
    alpha = 0.3
  ) +

geom_point(
    data = pca_data,
    aes(x = PC1, y = PC2, color = Locality, shape = Locality),
    size = 2
  ) +
  scale_fill_manual(values = my_lineage_colors) +     
  scale_color_manual(values = my_locality_cols) +     
  scale_shape_manual(values = my_locality_shapes) +   
  theme_minimal()

print(mi_plot)

ggsave("PCA_polygons.pdf", 
       plot = mi_plot, 
       device = "pdf", 
       width = 7, 
       height = 5, 
       units = "in")



########################################################### 
################## MULTIVARIATE ANALYSIS ################## 
########################################################### 

### MANOVA ###
mod_manova <- manova(as.matrix(morf_log) ~ Locality + Lineage, data = morf)
summary(mod_manova, test = "Pillai")

##           Df Pillai approx F num Df den Df    Pr(>F)    
## Locality   7 2.2762   2.1416     63    280 1.365e-05 ***
## Residuals 42                                            



library(vegan)

# Euclidian distance from log-transformed matrix
dist_morf <- dist(morf_log, method = "euclidean")


### PERMANOVA ### 
adonis_res <- adonis2(dist_morf ~ Locality + Lineage, data = morf, permutations = 999)
adonis_res

## Permutation test for adonis under reduced model
## Permutation: free
## Number of permutations: 999

## adonis2(formula = dist_morf ~ Locality + Lineage, data = morf, permutations = 999)
##          Df SumOfSqs      R2      F Pr(>F)    
## Model     7   1.4383 0.33745 3.0559  0.001 ***
## Residual 42   2.8240 0.66255                  
## Total    49   4.2623 1.00000                  




### LDA ###
library(MASS)

lda_df <- data.frame(Lineage = morf$Lineage, morf_log)
lda_res <- lda(Lineage ~ ., data = lda_df)

lda_res
## lda(Lineage ~ ., data = lda_df)
## Prior probabilities of groups:
##      Eastern     Macquarie      Northern South_Georgia  Southeastern      Southern 
##       0.02          0.18          0.18          0.26          0.16          0.20 


## Proportion of trace:
##   LD1    LD2    LD3    LD4    LD5 
## 0.5711 0.2999 0.0929 0.0321 0.0040 



predict_lda <- predict(lda_res)
table(Real = lda_df$Lineage, Predicted = predict_lda$class)

##                  Predicted
## Real            Eastern Macquarie Northern South_Georgia Southeastern Southern
##  Eastern             1         0        0             0            0        0
##  Macquarie           0         5        0             2            2        0
##  Northern            0         0        8             0            1        0
##  South_Georgia       0         1        1            11            0        0
##  Southeastern        0         0        0             0            7        1
##  Southern            0         1        0             0            0        9


############################################################################
#####################  ANTARCTIC POLAR FRONT INFLUENCE  #####################
############################################################################
morf$Region <- ifelse(morf$Locality %in% c("Crozet", "Macquarie", "Falklands"), "Norte",
                      ifelse(morf$Locality %in% c("Antarctic_Peninsula", "South_Georgia", "South_Orkney", "South_Shetland"), "Sur", "Kerguelen"))

# MANOVA
mod_manova <- manova(as.matrix(morf_log) ~ Region, data = morf)
summary(mod_manova, test = "Pillai")

##           Df  Pillai approx F num Df den Df   Pr(>F)    
## Region     2 0.81417   3.0515     18     80 0.000324 ***
## Residuals 47 


# PERMANOVA
library(vegan)
adonis_res <- adonis2(dist_morf ~ Region, data = morf, permutations = 999)
adonis_res

## Permutation test for adonis under reduced model
## Permutation: free
## Number of permutations: 999

## adonis2(formula = dist_morf ~ Region, data = morf, permutations = 999)
##          Df SumOfSqs      R2      F Pr(>F)   
## Model     2   0.6690 0.15695 4.3751  0.002 **
## Residual 47   3.5933 0.84305                 
## 
Total    49   4.2623 1.00000     

library(ggplot2)
library(dplyr)

morf$Region <- ifelse(morf$Locality %in% c("Crozet", "Macquarie", "Falklands"), "Norte",
                      ifelse(morf$Locality %in% c("Antarctic_Peninsula", "South_Georgia", "South_Orkney", "South_Shetland"), "Sur", "Kerguelen"))


morf$Region <- factor(morf$Region, levels = c("Norte", "Sur", "Kerguelen"))

pca_scores <- as.data.frame(pca_res$x)
pca_scores$Region <- morf$Region
pca_scores$Locality <- morf$Locality

pca_data_hull <- pca_scores %>%
  group_by(Region) %>%
  slice(chull(PC1, PC2))  

mi_plot <- ggplot() +
  geom_polygon(
    data = pca_data_hull,
    aes(x = PC1, y = PC2, fill = Region, group = Region),
    alpha = 0.3
  ) +
  geom_point(
    data = pca_scores,
    aes(x = PC1, y = PC2, color = Locality, shape = Locality),
    size = 2
  ) +
  scale_fill_manual(values = c("NorthAPF" = "blue", "SouthAPF" = "red", "Kerguelen" = "purple")) +
  scale_color_manual(values = my_locality_cols) +     
  scale_shape_manual(values = my_locality_shapes) +   
  theme_minimal() +
  labs(title = "PCA North vs. South APF")

print(mi_plot)
ggsave("PCA_polygon_N_S_APF.pdf", 
       plot = mi_plot, 
       device = "pdf", 
       width = 7, 
       height = 5, 
       units = "in")

