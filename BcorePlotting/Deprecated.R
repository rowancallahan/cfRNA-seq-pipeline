#####################################################################################################
# Set of functions useful for cell specificity analysis using the Fantom5 phase1 database.
#
# make_MA_plot_rows_cols             --- Replaced by maplot?
# make_heatmaps                      --- Replaced by makeHeatmap
# make_heatmap_versatile             --- Replaced by makeHeatmap
# hex.maplot                         --- Mothballed for now, but worth revisiting
# MA_plot_two_samples_same_treatment --- Replaced by maplot?
# make_heatmaps_exacloud             --- Replaced by makeHeatmap
# make_MA_plot                       --- Replaced by maplot
# chart.Corr                         --- Jessica has a better way of doing this
# panel.smoothcust                   --- Component of chart.Corr
#
# Author: Mark Fisher and Theresa Lusardi
# Started - October, 2016
#####################################################################################################


make_MA_plot_rows_cols <-
function(Expression_list, vec_of_treatments, cntrl_str){
  #Not really useful for knitr reports, but there's probably some salvageable code in here
  #A version of make_MA_plot that makes a matrix of plots, with nrow= #non-control treatments and ncol = #normalized matrices passed to the function
  ##Function to plot a panel of MA plots, with ncol = number of normalization algorithms and nrow=non-control treatments in array set.
  #Inputs are the list of expression matrices (each element represents a particular normalization algorithm), a vector with a treatment assigned to each array. E.g., if .CEL files were called 001.CEL, 002.CEL, 003.CEL, and 004.CEL, and 001 and 002 are treatment x and 003 and 004 are treatment y, vector would be c("x", "x", "y", "y"), and the string of the control treatment's name (e.g., "Mock" or "Control", or perhaps "X" in our example).
  #This could probably be split into two functions-one that generates log ratios to control for all non-control treatments, and one that plots MA plots.
  Treatment_wise_mean_expression_list = Generate_mean_by_treatment_matrix(Expression_list, vec_of_treatments)
  log_ratio_list = list()
  for(i in 1:length(Treatment_wise_mean_expression_list)){
    cntrl_col = which(colnames(Treatment_wise_mean_expression_list[[i]])==paste0(cntrl_str,"_Mn"))
    log_ratio_mat = NULL
    treatmnt_nm_vec = NULL
    for (j in 1:ncol(Treatment_wise_mean_expression_list[[i]])){
      if (j != cntrl_col){
        log_ratio_mat = cbind(log_ratio_mat, Treatment_wise_mean_expression_list[[i]][,j]-Treatment_wise_mean_expression_list[[i]][,cntrl_col])
        colnames(log_ratio_mat)[ncol(log_ratio_mat)] = paste0(colnames(Treatment_wise_mean_expression_list[[i]])[j], "_vs_", colnames(Treatment_wise_mean_expression_list[[i]])[cntrl_col])
        treatmnt_nm_vec = c(treatmnt_nm_vec,colnames(Treatment_wise_mean_expression_list[[i]])[j])
      }
    }
    rownames(log_ratio_mat) = rownames(Expression_list[[i]])
    log_ratio_list[[i]] = log_ratio_mat
    names(log_ratio_list)[i] = names(Expression_list)[i]
    rm(log_ratio_mat)
  }
  rs = ncol(Treatment_wise_mean_expression_list[[i]])-1 #number of non-control treatments
  cs = length(Treatment_wise_mean_expression_list) # number of normalization types
  par(mfcol = c(rs,cs), cex = 1.2, mar=c(4,4,0.5,0.5), oma=c(1,1,1,1), mgp=c(2,.6,0))
  for(e in 1:length(Treatment_wise_mean_expression_list)){
    for(i in 1:length(treatmnt_nm_vec)){
      xe = rowMeans(cbind(Treatment_wise_mean_expression_list[[e]][,paste0(cntrl_str,"_Mn")], Treatment_wise_mean_expression_list[[e]][,treatmnt_nm_vec[i]]))
      ye = log_ratio_list[[e]][,i]
      plot(x=xe, y=ye, xlab=paste0('Mn lg2 int, cntrl + ', treatmnt_nm_vec[i], ": ", names(Treatment_wise_mean_expression_list)[e]), ylab=colnames(log_ratio_list[[e]])[i], cex.lab=1.2, pch='.', col='blue', main=names(Expression_list)[i])
      abline(h=0,col=grey(.5),lty="22")
    }
  }
  #dev.off()
}
make_heatmaps <-
function(matrix_with_ID_symbol_and_betas_for_sig_IDs, beta_string_vec, ID_colname, Symbol_colname, gsub_str=''){
  #Assumes that matrix_with_ID_symbol_and_betas_for_sig_IDs already has rownames of ID (e.g., transcript cluster ID)
  subset_of_mat_of_interest = matrix_with_ID_symbol_and_betas_for_sig_IDs[,beta_string_vec]
  class(subset_of_mat_of_interest) = "numeric"
  c.lim = quantile(subset_of_mat_of_interest, probs=.95,na.rm=T)
  if(class(subset_of_mat_of_interest) != "matrix"){
    subset_of_mat_of_interest = t(as.matrix(subset_of_mat_of_interest))
  }
  if(ncol(subset_of_mat_of_interest)>2 & nrow(subset_of_mat_of_interest)>1){
    tmp=aheatmap(subset_of_mat_of_interest,color=paste0("-PiYG:", nrow(matrix_with_ID_symbol_and_betas_for_sig_IDs)), breaks=c(min(subset_of_mat_of_interest),seq(from=-c.lim,to=c.lim, by=c.lim/((nrow(matrix_with_ID_symbol_and_betas_for_sig_IDs)-2)/2)), max(subset_of_mat_of_interest)),annCol=gsub(gsub_str, '',beta_string_vec),labCol=gsub(gsub_str, '',beta_string_vec),labRow = rep(" ", nrow(subset_of_mat_of_interest)))
    print(tmp)
    indices_in_heat_map_order = rev(tmp$rowInd)
    ordered_probeset_ids = rownames(subset_of_mat_of_interest) [indices_in_heat_map_order]
    indices_in_original_matrix = matrix_with_ID_symbol_and_betas_for_sig_IDs[,ID_colname] %in% ordered_probeset_ids
    ordered_gene_symbols = matrix_with_ID_symbol_and_betas_for_sig_IDs[indices_in_original_matrix,Symbol_colname]
    return(list(ordered_probeset_ids, ordered_gene_symbols))
  }else{
    return("")
  }
}
