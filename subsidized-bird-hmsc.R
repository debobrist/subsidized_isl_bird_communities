#### Set up: ####
rm(list = ls())

library(Hmsc)
library(tidyverse)
library(corrplot)
library(parallel)
library(PNWColors)
library(greekLetters)

set.seed(999)

# Both raw data and cleaned data for this analysis are available in the Avian and paired Vegetation data from 100 Islands Project (Central Coast) in the  Hakai Data Catalogue:  
# https://catalogue.hakai.org/dataset/ca-cioos_12f951d4-4155-4c05-969d-a7158412f579

# For this analysis, read in data from hmsc-analysis-2022 sub-folder.
# Metadata can also be found in this sub-folder.
Y <- read.csv("hmsc-site-by-species.csv")
X <- read.csv("hmsc-env-params.csv")
Tr <- read.csv("hmsc-bird-traits.csv")
plots <- read.csv("hmsc-xycoords.csv")

# To run HMSC, need row names in X and Y to be point count ID (pcid). These are in the first column. Remove column.
row.names(Y) <- Y$pcid
Y$pcid <- NULL

row.names(X) <- X$pcid
X$pcid <- NULL

# Row name for traits data frame needs to be the species name. Remove column.
row.names(Tr) <- Tr$species
Tr$species <- NULL

# To make model spatially explicit, need a df of XY coordinates with pcid as the row name. Remove pcid and island ID from this df.
xycoords <- plots
row.names(xycoords) <- xycoords$pcid
xycoords$pcid <- NULL
xycoords$island <- NULL

# Set up random effect to represent study design:
studyDesign <- plots %>%
  dplyr::select(pcid, island)

# Not sure why this doesn't work without converting to a character and then back to a factor after but otherwise I get an error saying: Error in checkForRemoteErrors(val) : 
# 2 nodes produced errors; first error: subscript out of bounds
studyDesign$pcid <- as.character(studyDesign$pcid)
studyDesign$island <- as.character(studyDesign$island)
rownames(studyDesign) <- NULL

studyDesign <- studyDesign %>% 
  dplyr::select(island, pcid) 

studyDesign$pcid <- as.factor(studyDesign$pcid)
studyDesign$island <- as.factor(studyDesign$island)
studyDesign$spatial <- as.factor(studyDesign$pcid)

rL1 = HmscRandomLevel(units = unique(as.factor(plots$pcid)))
rL2 = HmscRandomLevel(units = unique(as.factor(plots$island)))
rL.spatial = HmscRandomLevel(sData = xycoords)

#### Step 1: Fit HMSC: ####
# Define formula for environmental parameters:
XFormula <- ~l.area.std * soild15N.std + l.isolation.90.std + rocky.prop.std + hab.het.std + sqrt.wrack.std + shoredist.std + year.std 

# Define formula for species traits: 
TrFormula <- ~mass.std + guild + feeding.height + nesting.height

# Set up the model: 
mod <- Hmsc(Y = Y,
            XData = X, 
            XFormula = XFormula,
            TrData = Tr, 
            TrFormula = TrFormula,
            studyDesign = studyDesign,
            ranLevels = list(pcid = rL1, island = rL2, spatial = rL.spatial),
            distr = "poisson")


# Number of chains to run: 
nChains <- 2

# Is this a test run or actual run? Use the test run to see if the model will actually run. The results will be completely meaningless. If test.run = FALSE, the 
# actual analysis will be replicated. Note - this took approximately a week to run on Compute Canada's Cedar computer at Simon Fraser University.
test.run = TRUE

if (test.run){
  thin = 1
  samples = 10
  transient = 1
  verbose = 1
} else {
  thin = 10
  samples = 40000
  transient = 1000
  verbose = 100
}

# Run the model: 
mod <- sampleMcmc(mod,
                  thin = thin,
                  samples = samples,
                  transient = transient,
                  nChains = nChains,
                  verbose = verbose,
                  updater = list(GammaEta = FALSE))

# Save the model output if running the full version. 
saveRDS(mod, "bird-hmsc.rds")

#### Step 2: Evaluate MCMC convergence: ####
mpost <- convertToCodaObject(mod)

# Summary of model outputs:
summary(mpost$Beta)

# From the HMSC 3.0 Univariate vignette: 
# Here, we want the two chains to yield approximately the same results. Second, they should mix well (go up and down without apparent autocorrelation). Finally, look to see that they seem to have reached a stationary distribution, as e.g. the first half of the recorded iterations looks statistically identical to the second half.
plot(mpost$Beta)

# More quantitatively, we can look at the effective sample sizes (ESS) and the potential scale reduction factors (PSRF). If the effective sample sizes are close to theoretical value of the number of samples, it means there is little autocorrelation among consecutive samples. If potential scale reductor fators give numbers close to 1, it means the two chains give consistent results.

# ESS and PSRF of Beta estimates (effects of environment):
ess.beta = effectiveSize(mpost$Beta)
psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf

# ESS and PSRF of Gamma estimates (biological traits):
ess.gamma = effectiveSize(mpost$Gamma)
psrf.gamma = gelman.diag(mpost$Gamma,multivariate=FALSE)$psrf

# ESS and PSRF of Omega estimates (biotic interactions):
sppairs = matrix(sample(x = 1:16^2, size = 100))
tmp = mpost$Omega[[1]]
for (chain in 1:length(tmp)){
  tmp[[chain]] = tmp[[chain]][,sppairs]
}
ess.omega = effectiveSize(tmp)
psrf.omega = gelman.diag(tmp, multivariate = FALSE)$psrf

par(mfrow = c(3,2))

hist(ess.beta, breaks = 30, cex.lab = 1.5, cex.axis = 1.5, xlab = "ESS (beta)", main = NULL)
hist(psrf.beta, breaks = 30, cex.lab = 1.5, cex.axis = 1.5, xlab = "PSRF (beta)", main = NULL)

hist(ess.gamma, breaks = 30, cex.lab = 1.5, cex.axis = 1.5, xlab = "ESS (gamma)", main = NULL)
hist(psrf.gamma, breaks = 30, cex.lab = 1.5, cex.axis = 1.5, xlab = "PSRF (gamma)", main = NULL)

hist(ess.omega, breaks = 30, cex.lab = 1.5, cex.axis = 1.5, xlab = "ESS (omega)", main = NULL)
hist(psrf.omega, breaks = 30, cex.lab = 1.5, cex.axis = 1.5, xlab = "PSRF (omega)", main = NULL)

dev.off()

#### Step 3: Evaluate model fit (explanatory power): ####
# Evaluate the explanatory power of the model: 
# Computes predicted values by taking the posterior medians.
preds_poisson <- computePredictedValues(mod, expected = FALSE) 

MF <- evaluateModelFit(hM = mod, predY = preds_poisson) # Look at model fit.

# From the vignette: 
# RMSE = Root-mean-square-error between the predicted and observed values times the sign of the correlation. 
# SR2 = pseudo-R2 for poisson models. This is computed as the squared Spearman correlation between observed and predicted values times the sign of the correlation.
# O.RMSE, O.AUC, O.Tjur2 = observed and predicted data are also truncated to occurrences (presence/absences). 
# C.RMSE and C.SR2 are conditional on presence.

par(mfrow = c(4,2))

hist(MF$RMSE, 
     main = paste0("Mean = ", round(mean(MF$RMSE), 2)),
     breaks = 10, 
     cex.lab = 1.5, 
     cex.axis = 1.5)
hist(MF$SR2, main = paste0("Mean = ", round(mean(MF$SR2), 2)),
     breaks = 10, 
     cex.lab = 1.5, 
     cex.axis = 1.5)
hist(MF$O.AUC, main = paste0("Mean = ", round(mean(MF$O.AUC), 2)),
     breaks = 10, 
     cex.lab = 1.5, 
     cex.axis = 1.5)
hist(MF$O.TjurR2, main = paste0("Mean = ", round(mean(MF$O.TjurR2), 2)),
     breaks = 10, 
     cex.lab = 1.5, 
     cex.axis = 1.5)
hist(MF$O.RMSE, main = paste0("Mean = ", round(mean(MF$O.RMSE), 2)),
     breaks = 10, 
     cex.lab = 1.5, 
     cex.axis = 1.5)
hist(MF$C.SR2, main = paste0("Mean = ", round(mean(MF$C.SR2), 2)),
     breaks = 10, 
     cex.lab = 1.5, 
     cex.axis = 1.5)
hist(MF$C.RMSE, main = paste0("Mean = ", round(mean(MF$C.RMSE), 2)),
     breaks = 10, 
     cex.lab = 1.5, 
     cex.axis = 1.5)

dev.off()


#### Step 4: Explore parameter estimates (results): ####
#### 4a: Variance partitioning: ####
# Group variables of interest. Here, I included IBT parameters Group 1, subsidy parameters in Group 2, and included a 3rd group to represent the interaction between these two. Year is considered separately. VP will also output variation attributed to plot-level, island-level, and spatiail random effects.
VP = computeVariancePartitioning(mod, 
                                 group = c(1,1,2,1,2,1,2,2,4,3), 
                                 groupnames = c("Island characteristics","Marine influence", "Interaction", "Year"))

# Save VP as RDS because it takes a long time to refit.
saveRDS(VP, "variance-partitioning.rds")

# Define colour palette from PNWColors package.
pnw_colors <- pnw_palette(name = "Lake", n = 7)

# Reorder values according to prop explained by island biogeography: 
VP$vals <- VP$vals[,order(VP$vals[1,]), drop = F]

tiff("variance-partitioning.tiff", units = "mm", height = 150, width = 168, compression = "lzw", res = 800)
plotVariancePartitioning(mod, 
                         VP = VP,
                         main = NULL,
                         cols = pnw_colors, 
                         cex.names = 0.75, 
                         las = 2,
                         args.legend = list(x = 19, 
                                            y = 0.265, 
                                            cex = 0.64))
dev.off()

# How much of the variation in Beta estimates is driven by traits?
knitr::kable(VP$R2T$Beta) # 0.21 - 0.59

# How much of the variation in species abundances is driven by traits?
VP$R2T$Y # 0.51

#### 4b: Environmental parameters (beta estimates): ####

# I used corrplot to make beta estimate plots because it seemed more customizable than the plotBeta function built into Hmsc. This involves pulling out means and support from the model.
postBeta = getPostEstimate(mod, parName = "Beta")

mbeta = postBeta$mean
betaP = postBeta$support
betaN = postBeta$supportNeg

supportLevel = 0.95

toPlot = mbeta
toPlot = toPlot * ((betaP > supportLevel) + (betaP < (1 - supportLevel)) > 0)
betaMat = matrix(toPlot, nrow = mod$nc, ncol = ncol(mod$Y))

colnames(betaMat) <- colnames(Y)
rownames(betaMat) <- c("Intercept", "Island area", "δ15N", "Isolation", "Prop. rocky shore", "Habitat heterogeneity", "Wrack biomass", "Distance to shore", "Year", "Island area * δ15N")

# Make delta symbols:
rownames(betaMat)[rownames(betaMat) == "d15N"] <- paste(greeks("delta"), "15N", sep = "")
rownames(betaMat)[rownames(betaMat) == "Island area * d15N"] <- paste("Island area:", greeks("delta"), "15N", sep = "")

# Sort the matrix by row of interest:
betaMat2 <- betaMat[,order(betaMat[3,]),drop=F] # Re-order by column of interest. In this case, by soil d15N. 

# Drop the intercept:
betaMat3 <- betaMat2[-1, ] # Drop the intercept:

# Plot:
corrplot(betaMat3, 
         method = "color",
         mar = c(0, 0, 0, 0), 
         tl.col = "black",
         tl.cex = 0.75,
         tl.srt = 45,
         cl.pos = "r",
         cl.ratio = 0.2,
         cl.length = 5,
         cl.lim = c(-1.2, 1.2), 
         bg = "white",
         addgrid.col = "grey",
         is.corr = FALSE) 

#### 4c: Species-to-species associations: ####
# Species-to-species interactions are evaluated through a latent variable approach.
# Plot-level:
OmegaCor = computeAssociations(mod) # This converts co-variances to the scale of correlations (-1 to 1)
supportLevel = 0 # Convention seems to be to show all co-occurrences, not just those at 95% cut-off.

# Plot-level:
toPlot = ((OmegaCor[[1]]$support > supportLevel)
          + (OmegaCor[[1]]$support < (1 - supportLevel)) > 0) * OmegaCor[[1]]$mean

# Island-level: 
toPlot2 = ((OmegaCor[[2]]$support > supportLevel)
           + (OmegaCor[[2]]$support < (1 - supportLevel)) > 0) * OmegaCor[[2]]$mean

# XY coordinate-level:
toPlot3 = ((OmegaCor[[3]]$support > supportLevel)
           + (OmegaCor[[3]]$support < (1 - supportLevel)) > 0) * OmegaCor[[3]]$mean


# Pull out median and sd for each level of the co-occurrences:
median(toPlot) 
sd(toPlot)

median(toPlot2)
sd(toPlot2)

median(toPlot3)
sd(toPlot3)

# Plotting:
tiff("omega-estimates-multipanel.tiff", res = 800, width = 168, height = 100, units = "mm", compression = "lzw")

par(mfrow = c(1, 2))

corrplot(toPlot,
         method = "color",
         title = "Plot-level co-occurences",
         mar = c(0, 0, 4, 0),
         order = "FPC",
         type = "lower",
         bg = "white",
         tl.col = "black",
         tl.cex = 0.6,
         tl.srt = 45,
         cl.cex = 0.6,
         cl.length = 3)
mtext("(a)", side = 3, at = -1)

corrplot(toPlot2,
         method = "color",
         title = "Island-level co-occurences",
         mar = c(0, 0, 4, 0),
         order = "FPC",
         type = "lower",
         bg = "white",
         tl.col = "black",
         tl.cex = 0.6,
         tl.srt = 45,
         cl.cex = 0.6,
         cl.length = 3)
mtext("(b)", side = 3, at = -1)
dev.off()


#### 4d: Traits: ####
postGamma = getPostEstimate(mod, parName = "Gamma")

mGamma = postGamma$mean
GammaP = postGamma$support
GammaN = postGamma$supportNeg

supportLevel = 0.95

# Means: 
# Mean coefficients: 
toPlot = mGamma
toPlot = toPlot * ((GammaP > supportLevel) + (GammaP < (1 - supportLevel)) > 0)
GammaMat = matrix(toPlot, nrow = mod$nc, ncol = ncol(mod$Tr))

# Note: Since many of these traits are categorical, the intercept in the columns represents ground-nester that feeds on the ground, and is an insectivore.
colnames(GammaMat) <- c("Nesting height - off ground", "Feeding height - upper canopy", "Feeding height - lower canopy", "Feeding height - ground/lower canopy", "Guild - omnivore", "mass.std", "Intercept")

rownames(GammaMat) <- c("Intercept", "Island area", "δ15N", "Isolation", "Prop. rocky shore", "Habitat heterogeneity", "Wrack biomass", "Distance to shore", "Year", "Island area * d15N")

corrplot(GammaMat, 
         method = "color",
         mar = c(0, 0, 0, 0), 
         tl.col = "black",
         tl.cex = 0.75,
         tl.srt = 45,
         cl.pos = "r",
         cl.ratio = 0.2,
         cl.length = 5,
         bg = "white",
         addgrid.col = "grey",
         is.corr = FALSE) 

dev.off()

# Built-in code (easier to see effects of changing supportLevel):
plotGamma(mod, 
          post = postGamma,
          supportLevel = 0.95,
          param = "Mean",
          colors = colorRampPalette(c("#6E0220", "#E6876B", "white", "#6BACD0", "#053061")),
          trNamesNumbers = c(TRUE, FALSE), 
          covNamesNumbers = c(TRUE, FALSE),
          cex = c(0.5, 1, 1))


#### Step 5: Evaluate predictive power ####
# Note: This takes 2 weeks to run on SFU's Cedar computer.
partition = createPartition(mod, nfolds = 2) 

preds = computePredictedValues(mod, partition = partition)

evaluateModelFit(hM = mod, predY = preds)