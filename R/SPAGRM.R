
#' SPAGRM method in GRAB package
#' 
#' SPAGRM method is to analysis longitudinal phenotype for related samples in a large-scale biobank. 
#' 
#' @details 
#' Please check \code{?GRAB.control} for the generic list of \code{control} in \code{GRAB.NullModel()} and \code{GRAB.Marker()}.
#' 
#' @examples 
#' # Step 1: fit a null model
#' PhenoData = read.table(system.file("extdata", "example.pheno", package = "GRAB"), header = T)
#' obj.SPAGRM = GRAB.NullModel(survival::Surv(time,event) ~ Cova1 + Cova2, 
#'                             data = PhenoData, subjData = PhenoData$IID, 
#'                             method = "SPAGRM", traitType = "time-to-event")
#' 
#' # Step 2: perform score test
#' GenoFile = system.file("extdata", "nSNPs-10000-nsubj-1000-ext.bed", package = "GRAB")
#' OutputDir = system.file("results", package = "GRAB")
#' OutputFile = paste0(OutputDir, "/SPAGRMMarkers.txt")
#' GRAB.Marker(obj.SPAGRM, GenoFile = GenoFile,
#'             OutputFile = OutputFile)
#' head(read.table(OutputFile, header=T))
#' @export
GRAB.SPAGRM = function(){
  print("Check ?GRAB.SPAGRM for more details about 'SPAGRM' method.")
}

################### This file includes the following functions

# ------------ used in 'GRAB_Marker.R' -----------
# 1. checkControl.Marker.SPAGRM(control)
# 2. setMarker.SPAGRM(objNull, control)
# 3. mainMarker.SPAGRM()

# check the control list in marker-level testing
checkControl.Marker.SPAGRM = function(control)
{
  default.control = list(SPA_Cutoff = 2,
                         zeta = 0,
                         tol = 1e-5)
  
  control = updateControl(control, default.control)  # This file is in 'control.R'
  
  return(control)
}

checkControl.SPAGRM.NullModel = function(control,
                                         ResidMat,
                                         SparseGRM,
                                         PairwiseIBD)
{
  if(control$MaxQuantile < control$MinQuantile)
    stop("MaxQuantile(default is 0.75) should be larger than MinQuantile(default is 0.25).")
  
  if(control$OutlierRatio < 0)
    stop("OutlierRatio should be larger than or equal 0 (default is 1.5).")
  
  SubjID.In.Resid = ResidMat$SubjID
  SubjID.In.GRM = unique(c(SparseGRM$ID1, SparseGRM$ID2))
  SubjID.In.IBD = unique(c(PairwiseIBD$ID1, PairwiseIBD$ID2))
  
  if(any(!SubjID.In.Resid %in% SubjID.In.GRM))
    stop("At least one subject in residual matrix does not have GRM information.")
  
  if(any(!SubjID.In.IBD %in% SubjID.In.GRM))
    stop("At least one subject has IBD information but does not have GRM information.")
  
  return(control)
}

setMarker.SPAGRM = function(objNull, control)
{
  # the below function is in 'Main.cpp'
  setSPAGRMobjInCPP(objNull$Resid,
                    objNull$Resid.unrelated.outliers,
                    objNull$sum_R_nonOutlier,
                    objNull$R_GRM_R_nonOutlier,
                    objNull$R_GRM_R_TwoSubjOutlier,
                    objNull$R_GRM_R,
                    objNull$MAF_interval,
                    objNull$TwoSubj_list,
                    objNull$ThreeSubj_list,
                    control$SPA_Cutoff,
                    control$zeta,
                    control$tol)
  
  print(paste0("The current control$nMarkersEachChunk is ", control$nMarkersEachChunk,".")) # This file is in 'control.R'
}

mainMarker.SPAGRM = function(genoType, genoIndex, outputColumns)
{
  OutList = mainMarkerInCPP("SPAGRM", genoType, genoIndex);
  
  obj.mainMarker = data.frame(Marker = OutList$markerVec,           # marker IDs
                              Info = OutList$infoVec,               # marker information: CHR:POS:REF:ALT
                              AltFreq = OutList$altFreqVec,         # alternative allele frequencies
                              AltCounts = OutList$altCountsVec,     # alternative allele counts
                              MissingRate = OutList$missingRateVec, # alternative allele counts
                              zScore = OutList$zScore,              # standardized score statistics
                              Pvalue = OutList$pvalVec,             # marker-level p-value
                              hwepval = OutList$hwepvalVec)
                              
  return(obj.mainMarker)
}

SPAGRM.NullModel = function(ResidMatFile,    # two columns: column 1 is subjID, column 2 is Resid
                            SparseGRMFile,   # a path of SparseGRMFile get from getSparseGRM() function.
                            PairwiseIBDFile, # a path of PairwiseIBDFile get from getPairwiseIBD() function.
                            control = list(MaxQuantile = 0.75,
                                           MinQuantile = 0.25,
                                           OutlierRatio = 1.5,
                                           MaxNuminFam = 5,
                                           MAF_interval = c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)))
{
  ResidMat = data.table::fread(ResidMatFile)
  SparseGRM = data.table::fread(SparseGRMFile)
  PairwiseIBD = data.table::fread(PairwiseIBDFile)
  
  ResidMat$SubjID = as.character(ResidMat$SubjID)
  SparseGRM$ID1 = as.character(SparseGRM$ID1); SparseGRM$ID2 = as.character(SparseGRM$ID2)
  PairwiseIBD$ID1 = as.character(PairwiseIBD$ID1); PairwiseIBD$ID2 = as.character(PairwiseIBD$ID2)
  
  control = checkControl.SPAGRM.NullModel(control, ResidMat, SparseGRM, PairwiseIBD)
  
  MaxQuantile = control$MaxQuantile;
  MinQuantile = control$MinQuantile;
  OutlierRatio = control$OutlierRatio;
  MaxNuminFam = control$MaxNuminFam;
  MAF_interval = control$MAF_interval;
  
  SubjID = ResidMat$SubjID
  SparseGRM = SparseGRM %>% filter(ID1 %in% SubjID & ID2 %in% SubjID)
  PairwiseIBD = PairwiseIBD %>% filter(ID1 %in% SubjID & ID2 %in% SubjID)
  
  # Use residual information to define outliers / non-outliers
  Resid = ResidMat$Resid
  Quant = quantile(Resid, probs = c(MinQuantile, MaxQuantile))
  Range = max(Quant) - min(Quant)
  cutoffVec = c(min(Quant) - OutlierRatio * Range, max(Quant) + OutlierRatio * Range)
  
  cat("cutoffVec:\t",cutoffVec,"\n")
  ResidMat$Outlier = ifelse(Resid < cutoffVec[1] | Resid > cutoffVec[2],
                            TRUE, FALSE)
  
  cat("Outliers information is as below\n")
  print(ResidMat %>% filter(Outlier == TRUE) %>% dplyr::select(SubjID, Resid, Outlier) %>% arrange(Resid))
  
  # Decompose the subjects based on family structure and use a greedy algorithm to reduce family size if needed
  SparseGRM1 = SparseGRM
  SparseGRM1$pos1 = ResidMat$Resid[match(SparseGRM$ID1, ResidMat$SubjID)]
  SparseGRM1$pos2 = ResidMat$Resid[match(SparseGRM$ID2, ResidMat$SubjID)]
  SparseGRM1 = SparseGRM1 %>% mutate(Cov = abs(Value * pos1 * pos2))
  
  edges = t(SparseGRM1[, c("ID1", "ID2")])
  graph_GRM = make_graph(edges, directed = F)
  graph_list_all = graph_GRM %>% decompose()
  graph_length = lapply(graph_list_all, length)
  
  graph_list_1 = graph_list_all[graph_length == 1]
  SubjID.unrelated = lapply(graph_list_1, get.vertex.attribute) %>% unlist(use.names = FALSE)
  ResidMat.unrelated = ResidMat %>% filter(SubjID %in% SubjID.unrelated)
  SubjID.unrelated.nonOutlier = ResidMat.unrelated %>% filter(Outlier == FALSE) %>% select(SubjID) %>% unlist(use.names = F)
  
  # Values used in association analysys
  R_GRM_R = SparseGRM1 %>% filter(ID1 %in% SubjID.unrelated) %>% select(Cov) %>% sum
  sum_R_nonOutlier = ResidMat.unrelated %>% filter(Outlier == FALSE) %>% select(Resid) %>% sum
  R_GRM_R_nonOutlier = SparseGRM1 %>% filter(ID1 %in% SubjID.unrelated.nonOutlier) %>% select(Cov) %>% sum
  Resid.unrelated.outliers = ResidMat.unrelated %>% filter(Outlier == TRUE) %>% select(Resid) %>% unlist(use.names = F)
  R_GRM_R_TwoSubjOutlier = 0; TwoSubj_list = ThreeSubj_list = list();
  
  # initialize parameters
  graph_list_updated = list()
  graph_list = graph_list_all[graph_length > 1]
  nGraph = length(graph_list)
  index.outlier = 1
  
  if(nGraph != 0)
  {
    cat("Start process the related residual information.\n")
    
    for(i in 1:nGraph)
    {
      if(i %% 1000 == 0)
        cat("Processing the related residual information:\t", i,"/",nGraph,"\n")
      
      comp1 = graph_list[[i]]
      comp3 = V(comp1)$name
      
      # Step 0: calculate variance for the family
      pos1 = match(comp3, SubjID)
      outlierInFam = any(ResidMat$Outlier[pos1])
      
      block_GRM = make.block.GRM(comp1, SparseGRM)
      
      R_GRM_R.temp = as.numeric(t(Resid[pos1]) %*% block_GRM %*% Resid[pos1])
      R_GRM_R = R_GRM_R + R_GRM_R.temp
      
      if(!outlierInFam)
      {
        sum_R_nonOutlier = sum_R_nonOutlier + sum(ResidMat$Resid[pos1])
        R_GRM_R_nonOutlier = R_GRM_R_nonOutlier + R_GRM_R.temp
        next
      }
      
      # cat("Family ", i, " (with outliers) includes ", length(comp3), " subjects:", comp3, "\n")
      
      vcount = vcount(comp1)   # number of vertices 
      
      if(vcount <= MaxNuminFam)
      {
        graph_list_updated[[index.outlier]] = comp1
        index.outlier = index.outlier + 1
        next
      }
      
      # Step 1: remove the edges until the largest family size is <= MaxNuminFam, default is 5.
      
      comp1.temp = comp1
      tempGRM1 = SparseGRM1 %>% filter(ID1 %in% comp3 | ID2 %in% comp3) %>% arrange(Cov)
      for(j in 1:nrow(tempGRM1))
      {
        # cat("j:\t",j,"\n")
        edgesToRemove = paste0(tempGRM1$ID1[j],"|",tempGRM1$ID2[j])
        comp1.temp = delete.edges(comp1.temp, edgesToRemove)
        vcount = decompose(comp1.temp) %>% sapply(vcount)  # vertices count for the new graph after edge removal
        # cat("vcount:\t",vcount,"\n")
        if(max(vcount) <= MaxNuminFam)
          break;
      }
      
      # cat("Edge removal complete. Counts of vertices:\t", vcount,"\n")
      
      # Step 2: add the (removed) edges while keeping the largest family size <= MaxNuminFam, default is 5.
      
      tempGRM1 = tempGRM1[1:j,] %>% arrange(desc(Cov))
      comp1 = comp1.temp
      for(k in 1:nrow(tempGRM1))
      {
        # cat("k:\t",k,"\n")
        edgesToAdd = c(tempGRM1$ID1[k], tempGRM1$ID2[k])
        comp1.temp = add.edges(comp1, edgesToAdd)
        
        vcount = decompose(comp1.temp) %>% sapply(vcount)  # vertices count for the new graph after edge removal
        # cat("vcount:\t",vcount,"\n")
        
        if(max(vcount) <= MaxNuminFam)
          comp1 = comp1.temp
      }
      
      comp1 = decompose(comp1)
      
      # cat("Edge add complete. Counts of vertices:\t", comp1 %>% sapply(vcount),"\n")
      
      for(k in 1:length(comp1))
      {
        comp11 = comp1[[k]]
        comp13 = V(comp11)$name
        
        pos2 = match(comp13, SubjID)
        outlierInFam = any(ResidMat$Outlier[pos2])
        
        block_GRM = make.block.GRM(comp11, SparseGRM)
        
        R_GRM_R.temp = as.numeric(t(Resid[pos2]) %*% block_GRM %*% Resid[pos2])
        
        if(!outlierInFam){
          sum_R_nonOutlier = sum_R_nonOutlier + sum(ResidMat$Resid[pos2])
          R_GRM_R_nonOutlier = R_GRM_R_nonOutlier + R_GRM_R.temp
        }else{
          graph_list_updated[[index.outlier]] = comp11
          index.outlier = index.outlier + 1;
        }
      }
    }
    
    cat("Start process the Chow-Liu tree.\n")
    
    # Make a list of array index.
    arr.index = list()
    for(n in 1:MaxNuminFam)
    {
      temp = c()
      for(i in 1:n)
      {
        indexString = rep("c(1, 1, 1)", n)
        indexString[i] = "0:2"
        indexString = paste0(indexString, collapse = "%o%")
        cmd = paste0("temp = c(temp, list(arr.index", i, "=", indexString, "))")
        eval(parse(text = cmd))
      }
      arr.index[[n]] = temp
    }
    
    # build chou-liu-tree.
    n.outliers = length(graph_list_updated)
    if(n.outliers != 0)
    {
      ## The below values are only used in chou.liu.tree
      TwofamID.index = ThreefamID.index = 0
      for(index.outlier in 1:n.outliers)
      {
        if(index.outlier %% 1000 == 0)
          cat("Processing the CLT for families with outliers:\t", TwofamID.index, ", ", ThreefamID.index, "/", nGraph, "\n")
        
        comp1 = graph_list_updated[[index.outlier]]
        comp3 = V(comp1)$name
        n1 = length(comp3)
        pos3 = match(comp3, SubjID)
        
        Resid.temp = ResidMat$Resid[pos3]
        
        if(n1 == 1)
        {
          Resid.unrelated.outliers = c(Resid.unrelated.outliers, Resid.temp)
          next;
        }
        
        block_GRM = make.block.GRM(comp1, SparseGRM)
        
        tempIBD = PairwiseIBD %>% filter(ID1 %in% comp3 & ID2 %in% comp3)
        
        if(n1 == 2)
        {
          TwofamID.index = TwofamID.index + 1;
          
          R_GRM_R_TwoSubjOutlier.temp = as.numeric(t(Resid.temp) %*% block_GRM %*% Resid.temp)
          R_GRM_R_TwoSubjOutlier = R_GRM_R_TwoSubjOutlier + R_GRM_R_TwoSubjOutlier.temp
          
          Rho.temp = tempIBD$pa + 0.5*tempIBD$pb
          midterm = sqrt(Rho.temp^2 - tempIBD$pa)
          
          TwoSubj_list[[TwofamID.index]] = list(Resid = Resid.temp,
                                                Rho = c(Rho.temp + midterm, Rho.temp - midterm))
          next;
        }
        
        ThreefamID.index = ThreefamID.index + 1;
        
        CLT = chow.liu.tree(N = n1,
                            IBD = tempIBD,
                            IDs = comp3,
                            MAF_interval = MAF_interval)
        
        stand.S.temp = array(rowSums(mapply(function(x, y) x*y, arr.index[[n1]], Resid.temp)), rep(3, n1))
        
        ThreeSubj_list[[ThreefamID.index]] = list(CLT = CLT,
                                                  stand.S = c(stand.S.temp))
      }
      cat("Completed processing the CLT for families with outliers:\t", TwofamID.index, ", ", ThreefamID.index, "/", nGraph, "\n")
    }
  }
  
  obj = list(Resid = Resid, subjData = SubjID, N = length(SubjID), Resid.unrelated.outliers = Resid.unrelated.outliers,
             R_GRM_R = R_GRM_R, R_GRM_R_TwoSubjOutlier = R_GRM_R_TwoSubjOutlier,
             sum_R_nonOutlier = sum_R_nonOutlier, R_GRM_R_nonOutlier = R_GRM_R_nonOutlier,
             TwoSubj_list = TwoSubj_list, ThreeSubj_list = ThreeSubj_list, 
             MAF_interval = MAF_interval)
  
  class(obj) = "SPAGRM_NULL_Model"
  
  return(obj)
}

make.block.GRM = function(graph, 
                          GRM)    # three columns: "ID1", "ID2", and "Value"
{
  comp2 = get.data.frame(graph)
  
  # igraph gives an unexpected additional loop, which may change the block GRM
  # the below is to remove the additional loop
  comp2 = comp2[!duplicated(comp2),]
  
  comp3 = V(graph)$name
  
  colnames(GRM) = c("to", "from", "Value")
  
  n1 = nrow(comp2)
  comp2 = merge(comp2, GRM)
  n2 = nrow(comp2)
  
  if(n1 != n2)
    stop("Ask Wenjian Bi (wenjianb@pku.edu.cn) to check why 'n1 != n2'.")
  
  block_GRM = sparseMatrix(i = match(comp2$from, comp3),
                           j = match(comp2$to, comp3),
                           x = comp2$Value,
                           symmetric = T)
  return(block_GRM)
}

chow.liu.tree = function(N,
                         IBD,
                         IDs,
                         MAF_interval)
{
  CLT = c()
  
  # MAF_interval = c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
  for(index in 1:length(MAF_interval))
  {
    mu = MAF_interval[index]
    
    # p = c(G0, G1, G2)
    p0 = c((1-mu)^2, 2*mu*(1-mu), mu^2)
    
    # p = c(G00, G10, G20, G01, G11, G21, G02, G12, G22)
    pa.allele2 = c((1-mu)^2, 0, 0, 0, 2*mu*(1-mu), 0, 0, 0, mu^2)
    
    pb.allele1 = c((1-mu)^3, mu*(1-mu)^2, 0, mu*(1-mu)^2, mu*(1-mu), mu^2*(1-mu), 0, mu^2*(1-mu), mu^3)
    
    pc.allele0 = c((1-mu)^4, 2*mu*(1-mu)^3, mu^2*(1-mu)^2, 2*mu*(1-mu)^3, 4*mu^2*(1-mu)^2, 
                   2*mu^3*(1-mu), mu^2*(1-mu)^2, 2*mu^3*(1-mu), mu^4)
    
    # calculate entropy I(Gi, Gj). Noting that entropy of unrelated pairs is zero.
    for(j in 1:nrow(IBD))
    {
      pro = IBD$pa[j] * pa.allele2 + IBD$pb[j] * pb.allele1 + IBD$pc[j] * pc.allele0
      
      entropy = sum(pro * log(pro/pc.allele0), na.rm = T)
      IBD$entropy[j] = entropy
    }
    
    # use the "prim" lgorithm to bulid a maximum spanning tree.
    Max_span_tree = IBD %>% graph_from_data_frame(directed = T) %>% 
      mst(weights = - IBD$entropy, algorithm = "prim") %>% get.edgelist() %>% 
      data.table::as.data.table() %>% rename(ID1 = V1, ID2 = V2)
    
    mst.IBD = merge(Max_span_tree, IBD, all.x = T) %>%
      mutate(idxID1 = match(ID1, IDs), idxID2 = match(ID2, IDs))
    
    arr.prob = array(1, dim = rep(3, N))
    for(i in 1:N)
      dimnames(arr.prob)[[i]] = paste0("ID",i,":",0:2) 
    
    vec = c(mst.IBD$idxID1, mst.IBD$idxID2); vec = vec[duplicated(vec)]
    
    for(k in 1:(N - 1))
    {
      pro = mst.IBD$pa[k] * pa.allele2 + mst.IBD$pb[k] * pb.allele1 + mst.IBD$pc[k] * pc.allele0
      
      matrix.prob = matrix(pro, 3, 3)
      matrix.index1 = mst.IBD$idxID1[k]; matrix.index2 = mst.IBD$idxID2[k]
      for(i in 1:3){
        for(j in 1:3){
          indexString = rep("", N)
          indexString[matrix.index1] = i
          indexString[matrix.index2] = j
          indexString = paste0(indexString, collapse = ",")
          cmd = paste0("arr.prob[",indexString,"] = arr.prob[", indexString, "] * matrix.prob[",i,",",j,"]")
          # "arr.prob[1,1,] = arr.prob[1,1,] * matrix.prob[1,1]"
          eval(parse(text = cmd))
        }
      }
    }
    
    for(k in 1:(N - 2))
    {
      vector.prob = p0
      vector.index = vec[k]
      for(i in 1:3){
        indexString = rep("", N)
        indexString[vector.index] = i
        indexString = paste0(indexString, collapse = ",")
        cmd = paste0("arr.prob[",indexString,"] = arr.prob[", indexString, "] / vector.prob[",i,"]")
        # arr.prob[,,1] = arr.prob[,,1] / vector.prob[1]"
        eval(parse(text = cmd))
      }
    }
    
    CLT = cbind(CLT, c(arr.prob))
  }
  
  return(CLT)
}

