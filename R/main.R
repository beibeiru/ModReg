#' @title Modulation mode inference
#' @description Calculate the modulation relation between a modulator and a regulator.
#' @param exprMatrix Gene expression matrix with gene symbol (row) x sample (column).
#' @param modulator A gene symbol for the modulator.
#' @param regulator A gene symbol for the regulator
#' @param regulon A vector for the target genes of the regulator.
#' @param cutoff.modulator A number (1/5 <= cutoff <= 1/2) defining the samples with highly and lowly expressed modulator.
#' @param permut_num Permutation number for differential mutual information analysis.
#' @param mc.cores Number of CPU cores.
#' @return An ModReg object
#' @details
#' ModReg is comprised of two steps. At first, all tripartite relations of Modulator-Regulator-Target were estimated based on differential connectivity analysis. Specifically, all tumor samples were separated into the top (M+) and bottom (M−) of samples in which the modulator was most and least expressed. In these two sets (M+ and M−), the significance of differential correlation relationships between one regulator and its targets were quantified using the difference in z-scores by Fisher z-transformation. A two-sided p value can be calculated using the standard normal distribution, and Benjamini–Hochberg method was used to adjust P-values for multiple hypothesis tests. Significant tripartite relations were kept for the next step (adjusted p<0.05).
#'
#' Secondly, ModReg calculated the differential average mutual information (DAMI) of the significant regulated targets from step 1, which equals to the difference of average mutual information (i.e., dependence relationship) of all significant regulated targets between M+ and M− sample sets. Thus, DAMI quantified how much the modulator was able to affect the fucntion of the regulator (DAMI > 0 for activation and DAMI < 0 for inhibition, respectively). Permutation test is used to estimate the empirical distribution of DAMI by randomly selecting a target set of equal size with the significant regulated targets from the complete regulon of the regulator, and the associated p value was estimated as the percentage of random trials with an absolute value of DAMI greater than the absolute value of the measured DAMI. Based on the DAMI and p value here, users can determine the modulation mode of the modulator and the regulator.
#'
#' @examples
#' load(file.path(system.file(package = "ModReg"), "extdata/exampleData.rda"))
#' ModReg_obj <- ModReg(exprMatrix, modulator="SIRT7", regulator="NFE2L2", regulon)
#'
#' @rdname ModReg
#' @export
#'
ModReg <- function(exprMatrix, modulator, regulator, regulon, cutoff.modulator=1/4, permut_num=1000, mc.cores=1)
{
  exprMatrix <- as.matrix(exprMatrix)

  if(!modulator%in%rownames(exprMatrix)) stop("The modulator gene is not included in the expression matrix.")
  if(!regulator%in%rownames(exprMatrix)) stop("The regulator gene is not included in the expression matrix.")

  regulon_missing <- regulon[!regulon%in%rownames(exprMatrix)]
  if(length(regulon_missing)>0)
  {
    print("The following genes are not included in the expression matrix, and have been excluded for the current analysis.")
    print(paste0(regulon_missing,collapse=", "))
  }
  regulon <- setdiff(regulon, regulon_missing)

  if(cutoff.modulator > 1/2 | cutoff.modulator < 1/5) stop("The cutoff should be 1/5<=cutoff.modulator<=1/2.")


  sample_no <- round(dim(exprMatrix)[2] * cutoff.modulator)

  value_vector <- exprMatrix[modulator,]
  sorted_vector <- sort(value_vector)

  CRlow <- as.numeric( names(value_vector)%in%names(sorted_vector)[1:sample_no] )
  CRhigh <- as.numeric( names(value_vector)%in%names(sorted_vector)[(length(value_vector)-sample_no+1):length(value_vector)] )

  exprMatrix_CRlow <- exprMatrix[c(regulator,regulon),CRlow==1,drop=F]
  exprMatrix_CRhigh <- exprMatrix[c(regulator,regulon),CRhigh==1,drop=F]

  ###### 1. differential correlation
  res <- matrix(0,ncol=7,nrow=length(regulon))
  colnames(res) <- c("modulator_low_r","modulator_low_p","modulator_high_r","modulator_high_p","diff_z","diff_p","diff_padj")
  rownames(res) <- regulon

  for(target in regulon)
  {
    CRlow_cor <- cor.test(exprMatrix_CRlow[regulator,],exprMatrix_CRlow[target,])
    CRlow_r <- CRlow_cor$estimate
    CRlow_p <- CRlow_cor$p.value

    CRhigh_cor <- cor.test(exprMatrix_CRhigh[regulator,],exprMatrix_CRhigh[target,])
    CRhigh_r <- CRhigh_cor$estimate
    CRhigh_p <- CRhigh_cor$p.value

    CRlow_z <- atanh(CRlow_r)
    CRhigh_z <- atanh(CRhigh_r)

    dz <- (CRhigh_z-CRlow_z)/sqrt(2*(1/(sample_no-3)))
    dz_p <- 2*pnorm(-abs(dz))

    res[target,"modulator_low_r"] <- CRlow_r
    res[target,"modulator_low_p"] <- CRlow_p
    res[target,"modulator_high_r"] <- CRhigh_r
    res[target,"modulator_high_p"] <- CRhigh_p
    res[target,"diff_z"] <- dz
    res[target,"diff_p"] <- dz_p
  }

  res[,"diff_padj"] <- p.adjust(res[,"diff_p"],"BH")
  res <- res[order(res[,"diff_p"]),]


  ####### 2. differential average mutual information (DAMI) of the regulated TGs
  ####### 2.1. find the regulated TGs
  res_signif <- res[res[,"diff_padj"]<0.05,]
  regulon_signif <- rownames(res_signif)

  targets_num <- dim(res_signif)[1]
  signif_ratio <- dim(res_signif)[1]/dim(res)[1]

  act_targets <- c()
  act <- 0
  inh_targets <- c()
  inh <- 0

  for(target in regulon_signif)
  {
    if(cor(exprMatrix[regulator,],exprMatrix[target,])>0)
    {
      if( res_signif[rownames(res_signif)==target,"modulator_high_r"]>0 &
          res_signif[rownames(res_signif)==target,"modulator_high_r"]>res_signif[rownames(res_signif)==target,"modulator_low_r"])
      {
        act <- act+1
        act_targets <- c(act_targets,target)
      }
      if( res_signif[rownames(res_signif)==target,"modulator_low_r"]>0 &
          res_signif[rownames(res_signif)==target,"modulator_high_r"]<res_signif[rownames(res_signif)==target,"modulator_low_r"])
      {
        inh <- inh+1
        inh_targets <- c(inh_targets,target)
      }

    }else{
      if( res_signif[rownames(res_signif)==target,"modulator_high_r"]<0 &
          res_signif[rownames(res_signif)==target,"modulator_high_r"]<res_signif[rownames(res_signif)==target,"modulator_low_r"])
      {
        act <- act+1
        act_targets <- c(act_targets,target)
      }
      if( res_signif[rownames(res_signif)==target,"modulator_low_r"]<0 &
          res_signif[rownames(res_signif)==target,"modulator_high_r"]>res_signif[rownames(res_signif)==target,"modulator_low_r"])
      {
        inh <- inh+1
        inh_targets <- c(inh_targets,target)
      }

    }
  }

  act_ratio <- act/dim(res)[1]
  inh_ratio <- inh/dim(res)[1]
  uni_ratio <- 1-(act_ratio+inh_ratio)

  if(act_ratio > 2*inh_ratio)
  {
    mode <- "activate"
  }else if(inh_ratio > 2*act_ratio){
    mode <- "inhibit"
  }else{
    mode <- "undetermined"
  }


  ####### 2.2. dami and its p value of the regulated TGs
  mi_mat_low  <- parmigene::knnmi.all(as.matrix(exprMatrix_CRlow[regulon,]))
  mi_mat_high <- parmigene::knnmi.all(as.matrix(exprMatrix_CRhigh[regulon,]))
  dmmi_all <- mi_mat_high - mi_mat_low

  set.seed(123)
  dami_bg <- parallel::mclapply(
    rep(targets_num, permut_num),
    function(targets_num){
      targets_sub <- sample(regulon,targets_num);
      sum(dmmi_all[targets_sub,targets_sub])/(targets_num^2-targets_num)},
    mc.cores=mc.cores)

  dami <- sum(dmmi_all[regulon_signif,regulon_signif])/(targets_num^2-targets_num)
  dami_bg <- unlist(dami_bg)
  dami_p <- sum(abs(dami_bg)>=abs(dami))/permut_num

  res_summary <- data.frame()
  res_summary["Value", "modulator"] <- modulator
  res_summary["Value", "mode"] <- mode
  res_summary["Value", "regulator"] <- regulator
  res_summary["Value", "signif_ratio"] <- signif_ratio
  res_summary["Value", "activate_ratio"] <- act_ratio
  res_summary["Value", "inhibit_ratio"] <- inh_ratio
  res_summary["Value", "dami"] <- dami
  res_summary["Value", "dami_p"] <- dami_p

  list(
    modulator.regulator = res_summary,
    modulator.regulator.target = res,
    regulon = regulon,
    cutoff.modulator = cutoff.modulator,
    permut_num = permut_num
  )

}


#' @title Modulator-Regulator visualization
#' @description Visualize the modulation relation between a modulator and a regulator.
#' @param exprMatrix Gene expression matrix with gene symbol (row) x sample (column).
#' @param ModReg_obj An ModReg object
#' @return A ggplot2 object
#' @details
#' The left and right panels present the patient group with low (M−) and high (M+) level of SIRT7 expression, in which rows are targets and columns are patient samples, respectively. In each patient group (M− or M+), the samples are sorted based on NFE2L2 expression level.
#'
#' @examples
#' load(file.path(system.file(package = "ModReg"), "extdata/exampleData.rda"))
#' ModReg_obj <- ModReg(exprMatrix, modulator="SIRT7", regulator="NFE2L2", regulon)
#' plotModReg(exprMatrix, ModReg_obj)
#'
#' @rdname plotModReg
#' @export
#'
plotModReg <- function(exprMatrix, ModReg_obj)
{
  exprMatrix <- as.matrix(exprMatrix)

  modulator <- ModReg_obj$modulator.regulator$modulator
  regulator <- ModReg_obj$modulator.regulator$regulator
  cutoff.modulator <- ModReg_obj$cutoff.modulator

  res <- ModReg_obj$modulator.regulator.target
  res_signif <- res[res[,"diff_padj"]<0.05,]

  regulon_signif <- rownames(res_signif)
  cor_obj <- cor(t(exprMatrix[regulator,,drop=F]),t(exprMatrix[regulon_signif,]))
  regulon_signif <- names(sort(cor_obj[1,],decreasing=T))


  sample_no <- round(dim(exprMatrix)[2] * cutoff.modulator)

  value_vector <- exprMatrix[modulator,]
  sorted_vector <- sort(value_vector)

  CRlow <- as.numeric( names(value_vector)%in%names(sorted_vector)[1:sample_no] )
  CRhigh <- as.numeric( names(value_vector)%in%names(sorted_vector)[(length(value_vector)-sample_no+1):length(value_vector)] )

  exprMatrix_CRlow <- exprMatrix[c(regulator,regulon),CRlow==1,drop=F]
  exprMatrix_CRhigh <- exprMatrix[c(regulator,regulon),CRhigh==1,drop=F]

  exprMatrix_CRlow <- as.matrix(exprMatrix_CRlow[regulon_signif,order(exprMatrix_CRlow[regulator,])])
  exprMatrix_CRhigh <- as.matrix(exprMatrix_CRhigh[regulon_signif,order(exprMatrix_CRhigh[regulator,])])

  ###### scale start
  scale.matrix <- function(expr_mat)
  {
    for(i in 1:dim(expr_mat)[1])
    {
      temp <- expr_mat[i,]
      temp <- rank(temp)
      temp <- (temp-median(temp))*2/(max(temp)-min(temp))
      expr_mat[i,] <- temp
    }
    expr_mat
  }

  exprMatrix_CRlow <- scale.matrix(exprMatrix_CRlow)
  exprMatrix_CRhigh <- scale.matrix(exprMatrix_CRhigh)
  ###### scale end

  exprMatrix_CRlow.m <- reshape2::melt(exprMatrix_CRlow)
  exprMatrix_CRhigh.m <- reshape2::melt(exprMatrix_CRhigh)

  exprMatrix_CRlow.m[,1] <- factor(exprMatrix_CRlow.m[,1], levels=rev(rownames(exprMatrix_CRlow)) )
  exprMatrix_CRhigh.m[,1] <- factor(exprMatrix_CRhigh.m[,1], levels=rev(rownames(exprMatrix_CRhigh)) )

  library(ggplot2)
  library(patchwork)
  p1 <- ggplot(exprMatrix_CRlow.m,aes(x = Var2, y = Var1, fill= value))+
    geom_tile()+
    scale_fill_gradientn(colours=c("#209fff","white","#ff4d00"))+
    theme_void() +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(0, 6, 0, 0),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.position = "none"
    )
  p2 <- ggplot(exprMatrix_CRhigh.m,aes(x = Var2, y = Var1, fill= value))+
    geom_tile()+
    scale_fill_gradientn(colours=c("#209fff","white","#ff4d00"))+
    theme_void() +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(0, 0, 0, 6),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.position = "none"
    )

  p1|p2

}


#' @title Modulator-Regulator-Target visualization
#' @description ModReg
#' @param exprMatrix Gene expression matrix with gene symbol (row) x sample (column).
#' @param modulator A gene symbol for the modulator.
#' @param regulator A gene symbol for the regulator.
#' @param target A target gene of the regulator.
#' @param cutoff.modulator A number (1/5 <= cutoff <= 1/2) defining the samples with highly and lowly expressed modulator.
#' @return A ggplot2 object
#' @examples
#' load(file.path(system.file(package = "ModReg"), "extdata/exampleData.rda"))
#' plotModRegTar(exprMatrix, modulator="SIRT7", regulator="NFE2L2", target="IRF2")
#'
#' @rdname plotModRegTar
#' @export
#'
plotModRegTar <- function(exprMatrix, modulator, regulator, target, cutoff.modulator=1/4)
{
  exprMatrix <- as.matrix(exprMatrix)

  if(!modulator%in%rownames(exprMatrix)) stop("The modulator gene is not included in the expression matrix.")
  if(!regulator%in%rownames(exprMatrix)) stop("The regulator gene is not included in the expression matrix.")
  if(!target%in%rownames(exprMatrix)) stop("The target gene is not included in the expression matrix.")

  if(cutoff.modulator > 1/2 | cutoff.modulator < 1/5) stop("The cutoff should be 1/5<=cutoff.modulator<=1/2.")


  sample_no <- round(dim(exprMatrix)[2] * cutoff.modulator)

  value_vector <- exprMatrix[modulator,]
  sorted_vector <- sort(value_vector)

  CRlow <- as.numeric( names(value_vector)%in%names(sorted_vector)[1:sample_no] )
  CRhigh <- as.numeric( names(value_vector)%in%names(sorted_vector)[(length(value_vector)-sample_no+1):length(value_vector)] )

  exprMatrix_CRlow <- exprMatrix[c(regulator,regulon),CRlow==1,drop=F]
  exprMatrix_CRhigh <- exprMatrix[c(regulator,regulon),CRhigh==1,drop=F]

  CRlow_cor <- cor.test(exprMatrix_CRlow[regulator,],exprMatrix_CRlow[target,])
  CRlow_r <- round(CRlow_cor$estimate,3)
  CRlow_p <- signif(CRlow_cor$p.value,3)

  CRhigh_cor <- cor.test(exprMatrix_CRhigh[regulator,],exprMatrix_CRhigh[target,])
  CRhigh_r <- round(CRhigh_cor$estimate,3)
  CRhigh_p <- signif(CRhigh_cor$p.value,3)


  fg.df <- data.frame(
    regulators=c(exprMatrix_CRlow[regulator,],exprMatrix_CRhigh[regulator,]),
    targets=c(exprMatrix_CRlow[target,],exprMatrix_CRhigh[target,]),
    group=c(
      rep(paste0("\n", modulator," Low expr\n", "r=", CRlow_r, ", p=",CRlow_p , "\n"),sample_no),
      rep(paste0("\n", modulator," High expr\n", "r=", CRhigh_r, ", p=",CRhigh_p, "\n"),sample_no)
    )
  )

  library(ggplot2)
  ggplot(fg.df,aes(regulators, targets, color=factor(group)))+
    geom_point(alpha=0.8)+
    xlab(paste0(regulator," expression"))+
    ylab(paste0(target," expression"))+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 18, colour = "black"),
      axis.ticks = element_blank(),
      axis.text = element_text(size = 16,colour = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 12, colour = "black"),
      legend.position = "right"
    )+geom_smooth(method = "lm", se = FALSE)

}

