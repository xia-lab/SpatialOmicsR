testthat::test_that("regression: test-tissue_mask", {
checking_installed_package<-nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"));if(dir.exists(".lib")) .libPaths(c(normalizePath(".lib"),.libPaths()))
if(checking_installed_package)library(SpatialOmicsMSI)else{source(file.path(getwd(),"R","msi_pipeline.R"));source(file.path(getwd(),"R","real_data_adapters.R"))}
g<-expand.grid(x=1:9,y=1:9);inside<-g$x%in%3:7&g$y%in%3:7;tic<-ifelse(inside,1000,10);pc<-ifelse(inside,100,2)
r<-build_msi_tissue_mask(g,tic,pc,"kmeans_log_tic_peak_count",seed=17)
stopifnot(sum(r$mask$tissue)==25L,all(r$mask$tissue==inside),r$diagnostics$background_pixels==56L,!any(unlist(r$diagnostics[c("touch_left","touch_right","touch_bottom","touch_top")])))
s<-build_msi_tissue_mask(g,tic,pc,"score_quantile",.6,17);stopifnot(sum(s$mask$tissue)>0,all(s$mask$tissue<=s$mask$candidate))
cat("TISSUE_MASK_TEST_OK=TRUE\n")

  testthat::succeed()
})
