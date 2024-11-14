#' WtCoxG method in GRAB package
#' 
#' WtCoxG method is an empirical approach to analyzing complex traits (including but not limited to time-to-event trait) by leveraging external MAFs. 
#'
#' @param GenoFile a character of genotype file. See Details section for more details.
#' @param GenoFileIndex additional index file(s) corresponding to GenoFile. See Details section for more details.
#' @param obj.WtCoxG a object with a class of "QCforBatchEffect"
#' @param control a list of parameters to decide which markers to extract.See \code{Details} section for more details
#'
#' @return an dataframe. including \code{WtCoxG.ext}, which utilizes external MAF, and \code{WtCoxG.noext} without external MAFs
#' 
#' @details
#' ## The details of \code{control} can be seen in \code{?GRAB.ReadGeno}
#' 
#' @export
#' @import dplyr, data.table
#' @examples
#' setwd(system.file("WtSPAG", package = "GRAB"))
#' PhenoData = read.table(system.file("WtSPAG", "simuPHENO_WtSPAG.txt", package = "GRAB"), header = T)
#' RefPrevalence = 0.1
#' #step0&1: fit a null model and estimate parameters according to batch effect p values
#' obj.WtCoxG = QCforBatchEffect(GenoFile = "simuBGEN1.bgen",
#'                              GenoFileIndex = c("simuBGEN1.bgen.bgi",
#'                                                 "simuBGEN1.sample"),
#'                              OutputFile = "qcBGEN1.txt",
#'                              control=list(AlleleOrder = "ref-first",
#'                                           AllMarkers = T,
#'                                           IndicatorColumn = "SurvEvent", SampleIDColumn = "IID"),
#'                              PhenoData=PhenoData,
#'                              RefAfFile = "RefMAFs.txt",
#'                              RefPrevalence = RefPrevalence,
#'                              SNPnum=1e4)
#' names(obj.WtCoxG)
#' #step2: conduct association testing
#' GRAB.WtCoxG(GenoFile = "simuBGEN1.bgen",
#'             GenoFileIndex = c("simuBGEN1.bgen.bgi", "simuBGEN1.sample"),
#'             obj.WtCoxG = obj.WtCoxG,
#'             OutputFile = "simuBGEN1.txt",
#'             control = list(AlleleOrder = "ref-first", AllMarkers=T))
GRAB.WtCoxG = function(GenoFile
                  , GenoFileIndex = NULL     # additional index file(s) corresponding to GenoFile
                  , obj.WtCoxG               # output list of QCforBatchEffect 
                  , GRM = NULL               # genetic relatedness matrix
                  , OutputFile               # output file path
                  , control=list(AlleleOrder = "ref-first", AllMarkers=T) 
                  , cutoff=0.1){
  
  PhenoData = obj.WtCoxG$PhenoData
  mergeGenoInfo = obj.WtCoxG$mergeGenoInfo
  RefPrevalence = obj.WtCoxG$RefPrevalence
  
  G = GRAB.ReadGeno(GenoFile = GenoFile
                    ,GenoFileIndex = GenoFileIndex
                    ,SampleIDs = PhenoData$SampleID
                    ,control = control)$GenoMat
  
  
  ## GWAS analysis ----------------------------------------------------------------------
  cat("Start GWAS analysis ########################################################## \n")
 
  t1=Sys.time()
  
  GWAS = lapply(1:ncol(G),function(i){
    if(i%%1000==0)cat("Complete ",i,"/",ncol(G),"\n")
    
    g = G[,i]
    R = PhenoData$R
    w = PhenoData$weight
    
    mu.ext = mergeGenoInfo$AF_ref[i]
    n.ext = mergeGenoInfo$AN_ref[i]/2
    TPR = mergeGenoInfo$TPR[i]
    sigma2 = mergeGenoInfo$sigma2[i]
    p_bat = mergeGenoInfo$pvalue_bat[i]
    w.ext = mergeGenoInfo$w.ext[i]
    var.ratio.w0 = mergeGenoInfo$var.ratio.w0[i]
    var.ratio.int = mergeGenoInfo$var.ratio.int[i]

    if(!is.null(GRM)){
      R_tilde = R - mean(R) * w.ext
      var.ratio0= (t(R_tilde) %*% GRM.ds %*% R_tilde  + w.ext^2 * sum(R)^2/n.ext)/(sum(R_tilde^2) + w.ext^2 * sum(R)^2/n.ext) 
      
    }else{var.ratio0=1}
    
      WtCoxG.ext = WtCoxG.test(g = g,
                          R = R,
                          w = w,
                          TPR=TPR,
                          sigma2 = sigma2,
                          b = w.ext,
                          var.ratio.w0=var.ratio.w0,
                          var.ratio.w1=var.ratio.w0,
                          var.ratio0 = var.ratio0,
                          var.ratio1 = var.ratio0,
                          mu.ext = mu.ext,
                          n.ext = n.ext,
                          p_bat = p_bat,
                          p_cut = cutoff)
      
      WtCoxG.noext = WtCoxG.test(g = g,
                               R = R,
                               w = w,
                               var.ratio.int = var.ratio.int,
                               p_bat = p_bat,
                               p_cut = cutoff)
      
    
    
    return( cbind(WtCoxG.ext, WtCoxG.noext) )
  }) %>%
    do.call("rbind",.) %>%
    cbind(.,mergeGenoInfo)
  
  t2=Sys.time()
  print(t2-t1)
  
  
  fwrite(GWAS, file = OutputFile)
  
}           


###SPA method----------------------------------------------------------------------
#####calculate p-value --------------------------------------------------

WtCoxG.test = function(g,
                       R,
                       w,
                       p_bat,
                       TPR=NA,
                       sigma2=NA,
                       b=0,
                       var.ratio.int = 1,              ### variance ratio of S.int
                       var.ratio.w0 = 1,           ### variance ratio of S.bat when bathceffect = 0
                       var.ratio.w1 = 1,           ### variance ratio of S.bat when bathceffect = 1
                       var.ratio0 = 1,             ### variance ratio of Score when bathceffect = 0
                       var.ratio1 = 1,             ### variance ratio of Score when bathceffect = 1
                       mu.ext=NA,
                       n.ext=NA,
                       p_cut=0.1                   ### batcheffect p value cut off
){

  ##imputation missing SNP
  missing.rate = mean(is.na(g))
  pos.na = which(is.na(g))  
  if(missing.rate != 0){
    g[pos.na] = mean(na.omit(g))
  }
  
  
  ### if external MAF is unavailable, then S.int = sum(R*(g-mu.int))
  if(is.na(mu.ext)){
    
    mu.int = mean(g)/2
    p.con = SPA_G.one.SNP_homo(g=g, R=R, mu.ext=NA, n.ext= 0,sigma2=0, var.ratio = var.ratio.int)[1]
    p_deno =NA
    return(p.con)
    
  }
  
  if(p_bat<p_cut|is.na(p_bat)|sum(g)<10|sum(2-g)<10 ){
    
    p.con=NA
    return(p.con)
    
  }else{
    
    meanR= mean(R)
    sumR = sum(R)
    mu.int = mean(g)/2
    N = length(g)
    
    
    mu = (1-b) * mu.int + b * mu.ext
    S= sum(R*(g-2*mu))
    # S=sum((R-(1-b)*meanR)*g)-sumR*2*b*mu.ext
    
    w1 = w/(2*sum(w))
    
    var_mu_ext = mu * (1 - mu)/(2 * n.ext)
    var_Sbat = sum(w1^2) * 2 * mu * (1 - mu) + var_mu_ext
    
    
    lb = -qnorm(1-p_cut/2) * sqrt(var_Sbat) * sqrt(var.ratio.w0)
    ub = qnorm(1-p_cut/2) *sqrt(var_Sbat) * sqrt(var.ratio.w0)
    c = pnorm(ub/sqrt(var.ratio.w1), 0, sqrt(var_Sbat+sigma2), log.p = T)
    d = pnorm(lb/sqrt(var.ratio.w1), 0, sqrt(var_Sbat+sigma2), log.p = T)
    p_deno = TPR*(  exp(d) * (exp(c-d) - 1) )+(1-TPR)*(1-p_cut)
    
    ##sigma2=0

    p_spa_s0 = SPA_G.one.SNP_homo(g=g, R=R, b=b ,mu.ext=mu.ext, n.ext= n.ext, sigma2=0,
                                  var.ratio = var.ratio0)[1]
    var_S = S^2 / var.ratio0 / qchisq(p_spa_s0, 1, ncp = 0, lower.tail = F)

    
    var.int = sum((R-(1-b) * meanR)^2) * 2 * mu * (1-mu)
    # var_S = var.int +  4*b^2 *sumR^2* var_mu_ext
    cov_Sbat_S = sum(w1 * (R-(1-b) * meanR)) * 2*mu * (1-mu) + 2 * b * sumR * var_mu_ext
    cov_Sbat_S = cov_Sbat_S*sqrt(var_S/(var.int +  4 * b^2 * sumR^2 * var_mu_ext))
    VAR = matrix(c(var_S, cov_Sbat_S, cov_Sbat_S, var_Sbat),nrow=2)
    p0 = max(0,min(1, pmvnorm(lower=c( -Inf, lb/sqrt(var.ratio.w0)),
                              upper=c( -abs(S/sqrt(var.ratio0)), ub/sqrt(var.ratio.w0)), mean=c(0,0), sigma=VAR)))
    
    ##sigma2!=0
    p_spa_s1 = SPA_G.one.SNP_homo(g=g, R=R, b=b, mu.ext=mu.ext, n.ext= n.ext, sigma2=sigma2,
                                  var.ratio = var.ratio1)[1]

    var_S1 = S^2/var.ratio1/qchisq(p_spa_s1, 1, ncp = 0, lower.tail = F)
    
    #var_S1 = var.int +  4*b^2 *sumR^2* (var_mu_ext+sigma2)
    cov_Sbat_S1 = sum(w1*(R-(1-b)*meanR))*2*mu*(1-mu)+2*b*sumR*(var_mu_ext+sigma2)
    cov_Sbat_S1 = cov_Sbat_S1*sqrt(var_S1/(var.int +  4*b^2 *sumR^2* (var_mu_ext+sigma2)))
    var_Sbat1 = var_Sbat+sigma2
    VAR1 = matrix(c(var_S1, cov_Sbat_S1, cov_Sbat_S1, var_Sbat1),nrow=2)
    p1= max(0,min(1,pmvnorm(lower=c( -Inf, lb/sqrt(var.ratio.w1)),
                            upper=c( -abs(S/sqrt(var.ratio1)), ub/sqrt(var.ratio.w1)), mean=c(0,0), sigma=VAR1)))
    
    p.con = 2*(TPR*p1+(1-TPR)*p0)/p_deno
    
    
    
    return(p.con) 
    
    
  }
  
  
}

#############################################################################
library(mvtnorm)
M_G0 = function(t, MAF){
  re = (1 - MAF + MAF * exp(t))^2
  return(re)
}
# The first derivative of the MGF of G (genotype)
M_G1 = function(t, MAF){
  re = 2*(MAF * exp(t))*(1 - MAF + MAF * exp(t))
  return(re)                           
}

# The second derivative of the MGF of G (genotype)
M_G2 = function(t, MAF){
  re = 2*(MAF * exp(t))^2 + 2*(MAF * exp(t))*(1 - MAF + MAF * exp(t))
  return(re)
}

# The CGF of G (genotype)
K_G0 = function(t, MAF){
  re = log(M_G0(t, MAF))
  return(re)
}

K_G1 = function(t, MAF){
  re = M_G1(t, MAF)/M_G0(t, MAF)
  return(re)
}

# The second derivative of the CGF of G (genotype)
K_G2 = function(t, MAF){
  re = M_G0(t, MAF)/M_G0(t, MAF) * M_G2(t, MAF)/M_G0(t, MAF) -(M_G1(t, MAF)/M_G0(t, MAF))^2
  return(re)
}

# The CGF of score test statistic 
H_org = function(t, R, MAF, n.ext, N.all, sumR, var_mu_ext, g.var.est, meanR, b){
  
  n.t = length(t)
  out = rep(0,n.t)
  R_hat  = sumR/N.all
  
  mu.adj = -2 * b * sumR * MAF
  var.adj = 4 * b^2 * sumR^2 * var_mu_ext
  for(i in 1:n.t){
    t1 = t[i]
    out[i] = sum(K_G0(t1 * (R - (1 - b) * meanR), MAF)) + mu.adj * t1 + var.adj/2 * t1^2  
    
  }
  
  return(out)
}

# The first derivative of the CGF of score test statistic 
H1_adj = function(t, R, s, MAF, n.ext, N.all, sumR,var_mu_ext, g.var.est,meanR,b)
{
  n.t = length(t)
  out = rep(0,n.t)
  R_hat  = sumR/N.all
  
  mu.adj = -2*b*sumR*MAF
  var.adj = 4*b^2 *sumR^2* var_mu_ext
  
  for(i in 1:n.t){
    t1 = t[i]
    out[i] = sum(( R - (1 - b) * meanR) *K_G1(t1 * ((R - (1 - b) * meanR)), MAF)) + mu.adj + var.adj * t1 - s
  }
  return(out)
}

# The second derivative of the CGF of score test statistic 
H2 = function(t, R, MAF, n.ext, N.all, sumR, var_mu_ext, g.var.est, meanR, b)
{
  n.t = length(t)
  out = rep(0,n.t)
  R_hat  = sumR/N.all
  
  mu.adj = -n.ext*R_hat*2*MAF
  var.adj = n.ext*R_hat^2*2*MAF*(1-MAF)
  
  
  for(i in 1:n.t){
    t1 = t[i]
    out[i] = sum((R - (1 - b) * meanR)^2 * K_G2(t1 * (R - (1 - b) * meanR) , MAF)) + var.adj
  }
  return(out)
}

GetProb_SPA_G = function(MAF, R, s, n.ext, N.all, sumR,var_mu_ext, g.var.est,meanR,b, lower.tail){
  
  out = uniroot(H1_adj, c(-1,1), extendInt = "yes",
                R=R, s=s, MAF=MAF, n.ext=n.ext, N.all = N.all, sumR = sumR,
                var_mu_ext = var_mu_ext, g.var.est=g.var.est, meanR =meanR, b=b)
  zeta = out$root
  
  k1 = H_org(zeta, R=R, MAF=MAF, n.ext=n.ext, N.all = N.all, sumR = sumR,
             var_mu_ext = var_mu_ext, g.var.est=g.var.est, meanR =meanR,b=b)
  k2 = H2(zeta, R=R, MAF=MAF, n.ext=n.ext, N.all = N.all, sumR = sumR,
          var_mu_ext = var_mu_ext, g.var.est=g.var.est, meanR =meanR, b=b)
  
  temp1 = zeta * s - k1
  
  w = sign(zeta) * (2 *temp1)^{1/2}
  v = zeta * (k2)^{1/2}
  
  pval = pnorm(w + 1/w * log(v/w), lower.tail = lower.tail)
  # pval.norm = pnorm(q2, lower.tail = lower.tail)
  re = pval
  return(re)
}


###saddlepoint approximation (SPA) to calicrate the p value in the case that batcheffect = 0 or 1
SPA_G.one.SNP_homo = function(g,                ### genotype vector 
                              R,                ### null model residual vector
                              mu.ext = NA,      ### external MAF
                              n.ext = NA,       ### external sample size
                              b = 0,            ### weight of external MAF
                              sigma2 = NA,      ###
                              var.ratio = 1,
                              Cutoff = 2,
                              impute.method = "fixed",
                              missing.cutoff = 0.15,
                              min.mac = 10,          # update on 2022-08-16 : replace 0.0001 by 0.000001
                              G.model = "Add")
{
  ## calculate MAF and update genotype vector
  
  ##imputation missing SNP
  missing.rate = mean(is.na(g))
  pos.na = which(is.na(g))
  if(missing.rate != 0){
    g[pos.na] = mean(na.omit(g))
  }
  
  
  if(is.na(mu.ext)){
    mu.ext =0
    n.ext=0
    
  }
  
  
  if(sum(g)<min.mac|sum(2-g)<min.mac|missing.rate>missing.cutoff){
    
    MAF= mean(na.omit(g))/2
    
    return(c(NA, NA))      
  }
  
  
  ######################################################################
  
  ## Score statistic
  N=length(g)
  mu.int = mean(g)/2
  MAF = (1-b) * mu.int + b * mu.ext 
  sumR = sum(R)
  N.all = N + n.ext
  S = sum(R *(g-2*MAF))
  S = S / var.ratio
  
  ## estimated variance without adjusting for covariates
  g.var.est = 2 * MAF * (1 - MAF)
  var_mu_ext = ifelse(n.ext==0,0, MAF*(1-MAF)/(2*n.ext)+sigma2) 
  
  
  #  S.var = sum(R^2 * g.var.est + sum(R)^2 * g.var.est/n.ext)
  meanR  = mean(R)
  S.var = sum((R-(1-b)*meanR)^2)*g.var.est +  4*b^2 *sumR^2* var_mu_ext
  
  z = S/sqrt(S.var)
  
  if(abs(z) < Cutoff){
    pval.norm = pnorm(abs(z), lower.tail = FALSE)*2
    return(c(pval.norm, pval.norm))  # update on 2022-10-05 : MAF.est.negative.num 
  }else{
    pval1 = GetProb_SPA_G(MAF, R = R, abs(S), n.ext=n.ext, N.all=N.all,
                          var_mu_ext = var_mu_ext, g.var.est=g.var.est,meanR =meanR,
                          sumR=sumR, b=b, lower.tail = FALSE) # EmpSPA-G p value 
    pval2 = GetProb_SPA_G(MAF, R = R, -abs(S), n.ext=n.ext, N.all=N.all,
                          var_mu_ext = var_mu_ext, g.var.est=g.var.est,meanR =meanR,
                          sumR=sumR, b=b, lower.tail = TRUE) # EmpSPA-G p value 
    
    pval3 = pnorm(abs(z), lower.tail = FALSE) # Normal
    pval4 = pnorm(-abs(z), lower.tail = TRUE) # Normal
    
    pval.spa.G = pval1 + pval2
    pval.norm = pval3 + pval4
    
    # if(abs(z) < Cutoff){
    #   pval.spa.G = pval.norm
    # }
    
    return(c(pval.spa.G, pval.norm))
    
    
  }
  
}


