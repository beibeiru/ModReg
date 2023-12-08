#' @title Calculate the modulation relation between a modulator and a regulator
#' @description ModReg is a function that evaluate the differential connectivity of a regulator and its targets between high and low expression level of a modulator to identify the modulation relation
#' @param exprMatrix a numeric data frame containing data from expression profile where columns contain samples and rows contain genes
#' @param modulator a character variable which define the modulator gene
#' @param regulator a character variable which define the regulator gene
#' @param regulon a character vector with target genes of regulator
#' @param mc.cores number of cpu cores using in parallel computation
#' @usage ModReg(exprMatrix, modulator, regulator, regulon, mc.cores=1)
#' @return Typically, the returned list contains two objects. The first one is a data frame of the table of differential correlations between low (CS) and high (TS) expression of modulator. The second one is matirx showing relation between the modulator and regulator, percentage of influential targets of the regulator.
#' @examples
#' load(file.path(system.file(package = "ModReg"), "extdata/exampleData.rda"))
#' res <- ModReg(exprMatrix, modulator, regulator, regulon)
#' @rdname ModReg
#' @export
ModReg <- function(exprMatrix, modulator, regulator, regulon, mc.cores=1)
{
  exprMatrix <- as.matrix(exprMatrix)

  sample_no <- round(dim(exprMatrix)[2]/4)

  value_vector <- unlist(exprMatrix[modulator,])
  sorted_vector <- sort(value_vector)

  CS <- as.numeric( names(value_vector)%in%names(sorted_vector)[1:sample_no] )
  TS <- as.numeric( names(value_vector)%in%names(sorted_vector)[(length(value_vector)-sample_no+1):length(value_vector)] )

  exprMatrix_sub <- exprMatrix[c(regulator,regulon),,drop=F]
  exprMatrix_sub_CS <- exprMatrix_sub[,CS==1]
  exprMatrix_sub_TS <- exprMatrix_sub[,TS==1]


  ###### differential correlation
  res <- matrix(0,ncol=7,nrow=length(regulon))
  colnames(res) <- c("CS_r","CS_p","TS_r","TS_p","dz","dz_p","dz_padj")
  rownames(res) <- regulon

  for(target in regulon)
  {
    CS_cor <- cor.test(
      exprMatrix_sub_CS[regulator,],
      exprMatrix_sub_CS[target,],
      method="pearson")
    CS_r <- CS_cor$estimate
    CS_p <- CS_cor$p.value

    TS_cor <- cor.test(
      exprMatrix_sub_TS[regulator,],
      exprMatrix_sub_TS[target,],
      method="pearson")
    TS_r <- TS_cor$estimate
    TS_p <- TS_cor$p.value

    CS_z <- atanh(CS_r)
    TS_z <- atanh(TS_r)

    dz <- (TS_z-CS_z)/sqrt(2*(1/(sample_no-3)))
    dz_p <- 2*pnorm(-abs(dz))

    res[target,"CS_r"] <- CS_r
    res[target,"CS_p"] <- CS_p
    res[target,"TS_r"] <- TS_r
    res[target,"TS_p"] <- TS_p
    res[target,"dz"] <- dz
    res[target,"dz_p"] <- dz_p
  }

  res[,"dz_padj"] <- p.adjust(res[,"dz_p"],"BH")
  res <- res[order(res[,"dz_p"]),]


  ####### differential multi mutual information
  run.dmmi <- function(targets)
  {
    targets_num <- length(targets)

    mi_mat_low  <- knnmi.all(as.matrix(exprMatrix_sub_CS[targets,]))
    mi_mat_high <- knnmi.all(as.matrix(exprMatrix_sub_TS[targets,]))

    dmmi <- mi_mat_high-mi_mat_low
    dmmi
  }


  res_sig <- res[res[,"dz_padj"]<0.05,]
  regulon_sig <- rownames(res_sig)

  targets_num <- dim(res_sig)[1]
  sig_ratio <- dim(res_sig)[1]/dim(res)[1]

  act_targets <- c()
  act <- 0
  inh_targets <- c()
  inh <- 0

  for(target in regulon_sig)
  {
    if(cor(t(exprMatrix[regulator,,drop=F]),t(exprMatrix[target,,drop=F]))>0)
    {
      if( res_sig[rownames(res_sig)==target,"TS_r"]>0 &
          res_sig[rownames(res_sig)==target,"TS_r"]>res_sig[rownames(res_sig)==target,"CS_r"])
      {
        act <- act+1
        act_targets <- c(act_targets,target)
      }
      if( res_sig[rownames(res_sig)==target,"CS_r"]>0 &
          res_sig[rownames(res_sig)==target,"TS_r"]<res_sig[rownames(res_sig)==target,"CS_r"])
      {
        inh <- inh+1
        inh_targets <- c(inh_targets,target)
      }

    }else{
      if( res_sig[rownames(res_sig)==target,"TS_r"]<0 &
          res_sig[rownames(res_sig)==target,"TS_r"]<res_sig[rownames(res_sig)==target,"CS_r"])
      {
        act <- act+1
        act_targets <- c(act_targets,target)
      }
      if( res_sig[rownames(res_sig)==target,"CS_r"]<0 &
          res_sig[rownames(res_sig)==target,"TS_r"]>res_sig[rownames(res_sig)==target,"CS_r"])
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


  dami <- sum(run.dmmi(regulon_sig))/(targets_num^2-targets_num)

  dmmi_all <- run.dmmi(regulon)

  permut_num <- 10000
  bg_dami <- mclapply(rep(targets_num,permut_num),
                      function(targets_num){
                        targets_sub <- sample(regulon,targets_num);
                        sum(dmmi_all[targets_sub,targets_sub])/(targets_num^2-targets_num)},
                      mc.cores=mc.cores)
  bg_dami <- unlist(bg_dami)

  dami_p <- sum(abs(bg_dami)>=abs(dami))/permut_num
  dami_a_p <- sum(bg_dami>=dami)/permut_num
  dami_i_p <- sum(bg_dami<=dami)/permut_num

  res_summary <- matrix(0,ncol=10,nrow=1)
  colnames(res_summary) <- c("Modulator","mode","regulator","sig_ratio","act_ratio","inh_ratio",
                             "dami","dami_p","dami_a_p","dami_i_p")
  rownames(res_summary) <- "Value"

  res_summary[1,] <- c(modulator,mode,regulator,sig_ratio,act_ratio,inh_ratio,dami,dami_p,dami_a_p,dami_i_p)

  list(res,res_summary)
}


