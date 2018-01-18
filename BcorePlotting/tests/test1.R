library(testthat)
require(limma)
require(data.table)
arguments <- commandArgs(trailingOnly = T)
baseDir_v <- arguments[1]
source(paste0(baseDir_v, "BcorePlotting/tests/test_data.R"))
source(paste0(baseDir_v, "BcorePlotting/ClusteringPlots.R"))

# Test test - uncomment this to ensure tests are working.
#expect_equal(2,3)

########################################################
### Check Mapping of plotting factors, samples, etc. ###
########################################################
  # Important output from this section:
      # There are the same number of grouping labels (sampClasses_v) as there are samples
      # There should be as many plotting colors as there are samples
      # There should be as many different plotting colors as there are different labels
      # Each label should only have one color assigned to it

### Create objects to test
mappingTest_lsv <- assignMappingSpecs(attribs = sampleAttribs_ls, 
                                      oneclass = sampleOneClass_v,
                                      colorspec = sampleColors_v,
                                      varPoints_v = sampleVarPoints_v)
mappingTest2_lsv <- assignMappingSpecs(attribs = sampleAttribs_ls, 
                                      oneclass = sampleOneClass_v,
                                      colorspec = sampleColors_v,
                                      varPoints_v = sampleVarPoints2_v)
mappingTest3_lsv <- assignMappingSpecs(attribs = sampleAttribs_ls, 
                                      oneclass = sampleOneClass_v,
                                      colorspec = sampleColors_v,
                                      varPoints_v = NULL)
compareMapping_dt <- data.table(mappingTest_lsv$sampClasses_v, mappingTest_lsv$plotColors_v)

### Test them!
test_that("Mapping Specs are correct", {
  # Total lengths are equal
  expect_equal(ncol(sampleData_mat), length(mappingTest_lsv$sampClasses_v))
  expect_equal(ncol(sampleData_mat), length(mappingTest_lsv$plotColors_v))
  # Groupings worked
  expect_equal(length(mappingTest_lsv$sampClasses_v), length(mappingTest_lsv$plotColors_v))
  expect_equal(length(mappingTest_lsv$uSampClasses_v), length(mappingTest_lsv$uPlotColors_v))
  # Mapping was appropriately done
  expect_equal(nrow(unique(compareMapping_dt)), length(mappingTest_lsv$uSampClasses_v))
  # Point scale is appropriate
  expect_true(unique(unique(mappingTest_lsv$plotCex_v) <= 1.9))
  expect_true(unique(unique(mappingTest_lsv$plotCex_v) >= 0.5))
  expect_true(unique(as.character(unique(mappingTest_lsv$plotCex_v)) %in% as.character(seq(.5,1.9, by = .1)))) # have to do char here
  expect_true(unique(as.character(unique(mappingTest2_lsv$plotCex_v)) %in% as.character(seq(.5,1.9, by = .1)))) 
  expect_equal(length(mappingTest_lsv$plotCex_v), length(mappingTest_lsv$plotColors_v))
  expect_equal(mappingTest3_lsv$plotCex_v, 0.6)
})

####################################
### Check Mapping of Label Stuff ###
####################################
  # Important output from this section:
      # Special plotting
          # There are as many point-type (pch) assignments as there are samples
          # plotLabels_v is null (to allow special pch's)
          # The point-types used for the legend are the same that are used for the samples (in terms of existence)
          # The point-types used for the legend are the same that are used for hte samples (in terms of mapping)
      # Default plotting
          # $pch is null
          # There are as many plot labels as there are samples
          # The plot labels are the same as the names of the samples

### Create objects to test
# Special plotting labels
labelTest1_lsv <- assignLabelSpecs(extraParams_ls = sampleExtraParams_ls,
                                   sampClasses_v = mappingTest_lsv$sampClasses_v,
                                   uSampClasses_v = mappingTest_lsv$uSampClasses_v,
                                   normmat = sampleData_mat)
labelTest1Pch_v <- sample(sampleExtraParams_ls$pch, size = 5)
# Default legend specs
labelTest2_lsv <- assignLabelSpecs(extraParams_ls = sampleExtraParams2_ls,
                                   sampClasses_v = mappingTest_lsv$sampClasses_v,
                                   uSampClasses_v = mappingTest_lsv$uSampClasses_v,
                                   normmat = sampleData_mat)

### Test them!
test_that("Special Label Specs are correct:", {
  # plotLabels_v is NULL
  expect_null(labelTest1_lsv$plotLabels_v)
  # Appropriate length
  expect_equal(ncol(sampleData_mat), length(labelTest1_lsv$pch))
  # Same values
  expect_equal(unique(labelTest1_lsv$pch), unname(labelTest1_lsv$pchLegend_v))
  # Same mapping
  for (i in 1:length(labelTest1Pch_v)){
    expect_equal(unique(names(labelTest1_lsv$pch[labelTest1_lsv$pch == labelTest1Pch_v[i]])),
                 names(labelTest1_lsv$pchLegend_v[labelTest1_lsv$pchLegend_v == labelTest1Pch_v[i]]))
  } # for
})

test_that("Default Label Specs are correct:", {
  # plotLabels_v is NULL
  expect_null(labelTest2_lsv$pch)
  # Same length
  expect_equal(ncol(sampleData_mat), length(labelTest2_lsv$plotLabels_v))
  # Same values
  expect_equal(labelTest2_lsv$plotLabels_v, colnames(sampleData_mat))
})

#############################
### Check Legend Creation ###
#############################
    # Important output from this section:
        # Always
            # Legend objects are lists of two lists. First specifies rectangle (4 elements), second specifies text (2 elements)
            # X and Y coordinates (second list) for text should be same length as the labels they're positioning
        # Default Location
            # Should be upper right corner
                # Upper coordinate should equal the upper limit of plotting region
                # Left coordinate should be less than right-most limit of plotting region
            
### Create objects to test
# Need a plot call for this to work...Is there a way to get the legend tests to work without this plot?
plot(x=1,y=1)
# Change position and use multiple pch's
legendTest1_lsv <- createLegend(legendPos_v = sampleLegendPos_v,
                               uSampClasses_v = mappingTest_lsv$uSampClasses_v,
                               uPlotColors_v = sampleColors_v,
                               pchLegend_v = labelTest1_lsv$pchLegend_v)

# Auto position and one pch
legendTest2_lsv <- createLegend(legendPos_v = sampleLegendPos2_v,
                                uSampClasses_v = mappingTest_lsv$uSampClasses_v,
                                uPlotColors_v = sampleColors_v,
                                pchLegend_v = labelTest2_lsv$pchLegend_v)

### Test them!
test_that("Legend Specifications are Appropriate", {
  # Legend object is list of 2 lists of length 4 and 2, respectively
  expect_length(legendTest1_lsv, 2)
  expect_length(legendTest2_lsv$rect, 4)
  expect_length(legendTest1_lsv$text, 2)

  # X and Y text positions are same length as attributes
  expect_equal(length(legendTest1_lsv$text$x), length(mappingTest_lsv$uSampClasses_v))
  expect_equal(length(legendTest2_lsv$text$y), length(mappingTest_lsv$uSampClasses_v))

  # Default plot location is in upper right corner
  expect_lt(legendTest2_lsv$rect$left, par("usr")[2])
  expect_equal(legendTest2_lsv$rect$top, par("usr")[4])
})

##########################
### Check MDS Creation ###
##########################

# Probably redundant to create tests here. If the three components pass tests...highly unlikely overall function is broken.
# Should still test that object is appropriate though
  # MDS type
  # 8 attributes
  # attributes 2, 3, 6, and 7 are same length as samples

plotMDSTest1 <- makeMDSplot(normmat = sampleData_mat,
                            attribs = sampleAttribs_ls,
                            oneclass = sampleOneClass_v,
                            colorspec = sampleColors_v,
                            plottitle = sampleTitle_v,
                            subtitle = sampleSubTitle_v,
                            ngenes = sampleNgenes_v,
                            legendPos_v = sampleLegendPos_v,
                            pch = samplePch_v,
                            xlab = sampleXLab,
                            varPoints_v = sampleVarPoints_v)

test_that("MDS object is appropriate", {
  # Should be MDS object
  expect_s4_class(plotMDSTest1, "MDS")
})

dev.off()