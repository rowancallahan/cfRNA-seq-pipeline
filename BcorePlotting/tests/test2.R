library(testthat)
require(data.table)
arguments <- commandArgs(trailingOnly = T)
baseDir_v <- arguments[1]
source(paste0(baseDir_v, "BcorePlotting/tests/test_data.R"))
source(paste0(baseDir_v, "BcorePlotting/ClusteringPlots.R"))

# Test test - uncomment this to ensure tests are working
#expect_equal(2,3)

#############################################
### Check Column Assignment/Rearrangement ###
#############################################
  # Important input/output
      # all vectors related to columns must be same length and values in same order
      # ratio mat order must correspond
      # if no attribs supplied, setCol_v must be NULL

### Create object to test
colSpecs_lsv <- setColSpecs(ratiomat = sampleData_mat, 
                            attribs = sampleAttribs_ls, 
                            setCol_v = sampleOneClass_v, 
                            colOrder_v = sampleColOrder_v)

colSpecs2_lsv <- setColSpecs(ratiomat = sampleData_mat, 
                            attribs = NA, 
                            setCol_v = sampleOneClass_v, 
                            colOrder_v = sampleColOrder_v)

### Run Tests
test_that("Column Specifications are Correct", {
  # Same number of columns in and out
  expect_equal(ncol(colSpecs_lsv$ratiomat), ncol(sampleData_mat))
  # Same number of columns as labCol_v
  expect_equal(ncol(colSpecs_lsv$ratiomat), length(colSpecs_lsv$labCol_v))
  # Same column order as labCol_v
  expect_equal(colnames(colSpecs_lsv$ratiomat), colSpecs_lsv$labCol_v)
  # Mapping of attribs is samples is maintained
  test_dt <- data.table(colnames(colSpecs_lsv$ratiomat), colSpecs_lsv$attribs[[sampleOneClass_v]])
  base_dt <- data.table(colnames(sampleData_mat), sampleAttribs_ls[[sampleOneClass_v]])
  setkey(base_dt, `V1`); setkey(test_dt, `V1`) # If mapping between these two variables was conserved, should both be equal now
  expect_identical(test_dt, base_dt)
  # NULL column cluster (default) if no attribs
  expect_null(colSpecs2_lsv$setCol_v)
})






#######################################################################################################################

# Stuff for later tests
# Get Attributes

# Custom order
# makeHeatmap(sampleData_mat,
#             list("Treatment" = sampleAttribs_ls$Treatment),
#             plottitle = "This is a Heatmap",
#             subtitle = "Bottom of the Heatmap",
#             setCol_v = "Treatment",
#             colOrder_v = unique(sampleAttribs_ls$Treatment))
# 
# # Cluster
# makeHeatmap(sampleData_mat,
#             list("Treatment" = sampleAttribs_ls$Treatment),
#             plottitle = "This is a Heatmap",
#             subtitle = "Bottom of the Heatmap")
# 
# makeHeatmap(sampleData_mat,
#             attribs = NA,
#             plottitle = "This is a Heatmap",
#             subtitle = "Bottom of the Heatmap")