#' Read data in MR format
#'
#' @param data input dataset
#' @param type whether exposure or outcome
#' @importFrom TwoSampleMR format_data clump_data
#' @export
readMRData = function(data, type, clump = FALSE){
  if(!type %in% c('outcome', 'exposure')) {
    stop('Please indicate type as "exposure" or "outcome"')
  }
  MR.data = TwoSampleMR::format_data(dat=data,
                                     type=type,
                                     phenotype_col = 'Phenotype',
                                     snp_col = 'SNP',
                                     beta_col = 'Beta',
                                     se_col = 'SE',
                                     eaf_col = 'EAF',
                                     effect_allele_col = 'EA',
                                     other_allele_col = 'OA',
                                     pval_col = 'P',
                                     chr_col = 'chromosome',
                                     pos_col = 'position')
  if (clump) {
    MR.data = clump_data(MR.data)
  }

  makeChrPos = function(data, type) {
    chr = paste0('chr.', type)
    pos = paste0('pos.', type)
    data$chr.pos = paste0(data[[chr]], ":", data[[pos]])
    return(data)
  }

  MR.data = makeChrPos(MR.data, type=type)
  return(MR.data)
}


#' Perform multivariate MR
#' @importFrom dplyr filter select
#' @importFrom TwoSampleMR mv_harmonise_data mv_multiple
#' @export
multivariate.MR = function(data, exposure, mediator, outcome) {

  exposure.data = dplyr::filter(data, Phenotype %in% c(exposure$name, mediator$name))
  outcome.data = dplyr::filter(data, Phenotype %in% outcome$name)

  # Read data in MR format
  MR.exposure.data = readMRData(exposure.data, 'exposure', clump = F)
  MR.outcome.data = readMRData(outcome.data, 'outcome')

  # Filter out by exposure chr.pos
  MR.outcome.data = filterByChrPos(MR.exposure.data, MR.outcome.data)

  # Harmonize for mvmr
  mv.data = TwoSampleMR::mv_harmonise_data(MR.exposure.data, MR.outcome.data)
  result = TwoSampleMR::mv_multiple(mv.data)
  result = result$result
  result = dplyr::select(result, c('exposure', 'outcome', 'nsnp', 'b', 'se', 'pval'))

  if (outcome$type=='binary') {
    result$mvmr.effect = paste0(
      round(exp(result$b),3), " (",
      round(exp(result$b - 1.96*result$se),2), ", ",
      round(exp(result$b + 1.96*result$se),2), ")"
    )
  } else{
    result$mvmr.effect = paste0(
      round((result$b),3), " (",
      round((result$b - 1.96*result$se),2), ", ",
      round((result$b + 1.96*result$se),2), ")"
    )
  }

  names(result) = c("exposure","outcome","mvmr.nsnp", "direct.effect",
                    "direct.se","direct.pval", "mvmr.direct.effect")

  return(result)
}


#' @importFrom dplyr filter
filterByChrPos = function(dt1, dt2) {
  dt1.chr.pos = dt1[['chr.pos']]
  return(dplyr::filter(dt2, chr.pos %in% dt1.chr.pos))
}

#' Perform univariate MR
#' @importFrom dplyr filter
#' @importFrom TwoSampleMR harmonise_data mr
#' @export
univariate.MR = function(data, exposure, outcome){

  extractSigNif = function(data, p=5e-8){
    return(dplyr::filter(data,P <= p))
  }

  exposure.data = extractSigNif(dplyr::filter(data, Phenotype == exposure$name))
  outcome.data = dplyr::filter(data, Phenotype %in% outcome$name)

  # Read data in MR format
  MR.exposure.data = readMRData(exposure.data, 'exposure',clump = F)
  MR.outcome.data = readMRData(outcome.data, 'outcome')
  # Filter out by exposure chr.pos
  MR.outcome.data = filterByChrPos(MR.exposure.data, MR.outcome.data)

  uni.mr.data <- TwoSampleMR::harmonise_data(MR.exposure.data, MR.outcome.data)
  uni.mr.res <- TwoSampleMR::mr(uni.mr.data, method_list = "mr_ivw")

  if(outcome$type == 'binary') {
    uni.mr.res$unimr.effect = paste0(
      round(exp(uni.mr.res$b),3), " (",
      round(exp(uni.mr.res$b - 1.96*uni.mr.res$se),2), ", ",
      round(exp(uni.mr.res$b + 1.96*uni.mr.res$se),2), ")"
    )
  } else{
    uni.mr.res$unimr.effect = paste0(
      round((uni.mr.res$b),3), " (",
      round((uni.mr.res$b - 1.96*uni.mr.res$se),2), ", ",
      round((uni.mr.res$b + 1.96*uni.mr.res$se),2), ")"
    )
  }

  uni.mr.res = uni.mr.res[, c('exposure', 'outcome', 'method', 'nsnp', 'b', 'se','pval', 'unimr.effect')]
  names(uni.mr.res) = c("exposure","outcome","method","unimr.nsnp","effect","se","pval","unimr.effect")
  return(uni.mr.res)
}

#' Perform mediator MR
#' @param data input data
#' @param exposure list of name of exposure and type
#' @param mediator list of name of mediator and type
#' @param outcome list of name of outcome and type
#' @export
#' @examples
#' fit = mrMediator(data=DTC,exposure = list(name='TSH', type = 'continuous'),mediator = list(name ='SHBG', type = 'continuous'),outcome = list(name='DTC', type ='binary'))
mrMediator = function(data, exposure, mediator, outcome) {
  expo.to.med = univariate.MR(data = data, exposure = exposure, outcome = mediator)
  expo.to.med$path = 'expo.to.med'
  expo.to.out = univariate.MR(data = data, exposure = exposure, outcome = outcome)
  expo.to.out$path = 'expo.to.out'
  med.to.out = univariate.MR(data = data, exposure = mediator, outcome = outcome)
  med.to.out$path = 'med.to.out'

  #combine the result
  unimr.res = rbind(expo.to.med, expo.to.out, med.to.out)
  mvmr.res = multivariate.MR(data=data, exposure = exposure, mediator = mediator, outcome = outcome)

  # diff
  diff.uni = unimr.res[unimr.res$exposure==exposure$name & unimr.res$outcome==outcome$name,]
  diff.mvmr = mvmr.res[mvmr.res$exposure==exposure$name & mvmr.res$outcome==outcome$name,]

  diff.res = merge(diff.uni, diff.mvmr, by=c('exposure', 'outcome'), all=T)
  diff.res$mediator = mediator$name
  diff.res$indirect.effect = diff.res$effect - diff.res$direct.effect
  diff.res$indirect.se = sqrt(diff.res$se^2 + diff.res$direct.se^2)
  diff.res$z.indirect = diff.res$indirect.effect / diff.res$indirect.se
  diff.res$indirect.pval = 2*pnorm(abs(diff.res$z.indirect), lower.tail = F)
  if(outcome$type == 'binary') {
    diff.res$diff.indirect.effect = paste0(
      round(exp(diff.res$indirect.effect),3), " (",
      round(exp(diff.res$indirect.effect - 1.96*diff.res$indirect.se),2), ", ",
      round(exp(diff.res$indirect.effect + 1.96*diff.res$indirect.se),2), ")"
    )
  } else{
    diff.res$diff.indirect.effect = paste0(
      round((diff.res$indirect.effect),3), " (",
      round((diff.res$indirect.effect - 1.96*diff.res$indirect.se),2), ", ",
      round((diff.res$indirect.effect + 1.96*diff.res$indirect.se),2), ")"
    )
  }
  diff.res = diff.res[, c('exposure', 'mediator', 'outcome',
                          # total effect
                          'unimr.nsnp', 'effect', 'se', 'pval', 'unimr.effect',
                          # direct effect
                          'mvmr.nsnp', 'direct.effect', 'direct.se', 'direct.pval', 'mvmr.direct.effect',
                          # indirect effect
                          'indirect.effect', 'indirect.se', 'indirect.pval', 'diff.indirect.effect')]

  # prod
  prod.res = data.frame(exposure = exposure$name)
  prod.res$mediator = mediator$name
  prod.res$outcome = outcome$name
  # indirect
  prod.res$indirect.effect = expo.to.med$effect * med.to.out$effect
  prod.res$indirect.se = sqrt(expo.to.med$effect^2 * med.to.out$se^2 + med.to.out$effect^2 + expo.to.med$se^2)
  prod.res$z.indirect = prod.res$indirect.effect / prod.res$indirect.se
  prod.res$indirect.pval = 2*pnorm(abs(prod.res$z.indirect), lower.tail = F)
  #direct
  prod.res$direct.effect = expo.to.out$effect - prod.res$indirect.effect
  prod.res$direct.se = sqrt(expo.to.out$se^2 + prod.res$indirect.se^2)
  prod.res$z.direct = prod.res$direct.effect / prod.res$direct.se
  prod.res$direct.pval = 2*pnorm(abs(prod.res$z.direct), lower.tail = F)

  if(outcome$type == 'binary') {
    prod.res$prod.direct.effect = paste0(
      round(exp(prod.res$direct.effect),3), " (",
      round(exp(prod.res$direct.effect - 1.96*prod.res$direct.se),2), ", ",
      round(exp(prod.res$direct.effect + 1.96*prod.res$direct.se),2), ")"
    )
    prod.res$prod.indirect.effect = paste0(
      round(exp(prod.res$indirect.effect),3), " (",
      round(exp(prod.res$indirect.effect - 1.96*prod.res$indirect.se),2), ", ",
      round(exp(prod.res$indirect.effect + 1.96*prod.res$indirect.se),2), ")"
    )
  } else{
    prod.res$prod.indirect.effect = paste0(
      round((prod.res$indirect.effect),3), " (",
      round((prod.res$indirect.effect - 1.96*prod.res$indirect.se),2), ", ",
      round((prod.res$indirect.effect + 1.96*prod.res$indirect.se),2), ")"
    )
    prod.res$prod.direct.effect = paste0(
      round((prod.res$direct.effect),3), " (",
      round((prod.res$direct.effect - 1.96*prod.res$direct.se),2), ", ",
      round((prod.res$direct.effect + 1.96*prod.res$direct.se),2), ")"
    )
  }


  result.list = list(diff =diff.res,
                     prod = prod.res)
  class(result.list) = 'mrMediator'
  return(result.list)

}

