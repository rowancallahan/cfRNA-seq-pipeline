# Create Test Matrix and other metadata
require(RColorBrewer)
set.seed(seed = 1)
# Create gene names
sampleGeneNames_v <-  sprintf("Gene%d", 1:200)
sampleSampleNames_v <- sprintf("Sample%d", 1:150)

# Create 200 X 150 matrix of log-normal distribution
  # mean = 10
  # sd = 2.5
  # TO DO: what are appropriate values for these parameters?

sampleData_mat <- NULL
for (i in 1:150){
  currCol_v <- rlnorm(200, log(10), log(2.5))
  sampleData_mat <- cbind(sampleData_mat, currCol_v)
} # for
colnames(sampleData_mat) <- sampleSampleNames_v
rownames(sampleData_mat) <- sampleGeneNames_v

# Create sample grouping attributes
# TO DO: add different things to this list.

sampleAttribs_ls <- list("Treatment" = c(unlist(sapply(sprintf("Treat%s", LETTERS[6:10]), function(x) rep(x, times = 10), simplify = F, USE.NAMES = F)), # Out of order to help test ordering
                                         unlist(sapply(sprintf("Treat%s", LETTERS[1:5]), function(x) rep(x, times = 20), simplify = F, USE.NAMES = F))),
                         "Gender" = sample(c("M", "F"), size = 150, replace = T))

# Specify attribute for plotting
sampleOneClass_v <- "Treatment"

# Colors
sampleColors_v <- c(brewer.pal(8, "Dark2"), brewer.pal(3, "Paired")[2:3])

# Extra values for plot scaling
sampleVarPoints_v = sample(1:1000, size = 150, replace = T)
sampleVarPoints2_v = runif(150)

# Plot Stuff
sampleTitle_v <- "This is a Test"
sampleSubTitle_v <- "This should be below."

sampleNgenes_v <- 20

sampleLegendPos_v <- "topleft"
sampleLegendPos2_v <- "auto"

# Extra Params
samplePch_v = c(15:24)
sampleXLab = "Test.X"
sampleYLab = "Test.Y"
sampleXLab2 = "Test2.X"
sampleYLab2 = "Test2.Y"

sampleExtraParams_ls <- list(pch = samplePch_v, xlab = sampleXLab, ylab = sampleYLab)
sampleExtraParams2_ls <- list(xlab = sampleXLab2, ylab = sampleYLab2)

# Heatmap-specific
sampleColOrder_v <- unique(sampleAttribs_ls$Treatment)[sample(1:10, size = 10, replace = F)]
