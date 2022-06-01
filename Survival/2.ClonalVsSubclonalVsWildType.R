
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/OS.byDriverGene.ClonalSubclonal.RData");
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/DFS.byDriverGene.ClonalSubclonal.RData");
source("0.CodeLibrary/Cox.function.R")
#处理清洗OS COX生存数据
{
  OS.significant <- c("ANK1","ARID1A","CASP8","GRIN2A","SMAD2")
  pos <- match(OS.significant,names(OS.byDriverGene.ClonalSubclonal))
  COX.OS.byDriverGene.ClonalSubclonal <- vector(mode="list",length=length(OS.significant))
  for(i in 1:length(pos)){
    COX.OS.byDriverGene.ClonalSubclonal[[i]] <- OS.byDriverGene.ClonalSubclonal[[pos[i]]]
  }
  names(COX.OS.byDriverGene.ClonalSubclonal) <- OS.significant
  head(COX.OS.byDriverGene.ClonalSubclonal[[1]])
  levels(COX.OS.byDriverGene.ClonalSubclonal[[1]]$sample.label)
  ## 对临床字段进行处理
  COX.OS.byDriverGene.ClonalSubclonal <- lapply(COX.OS.byDriverGene.ClonalSubclonal,function(x){
    #temp <- x[,c(1,10,11,2,3,4,5,6,7,8,9,12,13,14)] #所有临床变量
    temp <- x[,c(1,10,11,2,3,4,5,6,7,9,13,14)] #剔除亚型信息
    #colnames(temp) <- c("Patient_ID","event","time","sex","MSS_state","T_stage","N_stage","M_stage","AJCC_stage","Residual","age","Subtype","Tumour_site","Mutation_status")
    #colnames(temp) <- c("Patient_ID","event","time","sex","MSS_state","T_stage","N_stage","M_stage","AJCC_stage","Residual","age","Tumour_site","Mutation_status")
    #对协变量进行factor转换，控制levels的第一项为参考项,目前默认第一个level为参考项 
    temp$sex <- factor(temp$sex,levels=c("Female","Male"))
    temp$MSS_state[grep("MSI",temp$MSS_state)] <- "MSI"
    temp$MSS_state <- factor(temp$MSS_state,levels=c("MSI","MSS"))
    temp$T_stage[grep("T1|T2",temp$T_stage)] <- "T1/T2"
    temp$T_stage[grep("T3|T4",temp$T_stage)] <- "T3/T4"
    temp$T_stage <- factor(temp$T_stage,levels=c("T1/T2","T3/T4"))
    temp$N_stage[grep("N1|N2",temp$N_stage)] <- "N1/N2"
    temp$N_stage <- factor(temp$N_stage,levels=c("N0","N1/N2"))
    temp$M_stage <- factor(temp$M_stage,levels=c("M0","M1"))
    temp$TNM_stage[grep("Stage 1|Stage 2",temp$TNM_stage)] <- "Stage 1/Stage 2"
    temp$TNM_stage[grep("Stage 3|Stage 4",temp$TNM_stage)] <- "Stage 3/Stage 4"
    temp$TNM_stage <- factor(temp$TNM_stage,levels=c("Stage 1/Stage 2","Stage 3/Stage 4"))
    #temp$Residual[grep("R1|R2",temp$Residual)] <- "R1/R2"
    #temp$Residual <- factor(temp$Residual,levels=c("R0","R1/R2"))
    temp$age <- factor(temp$age,levels=c("<67.5",">=67.5"))
    #temp$Subtype <- factor(temp$Subtype,levels=c("CIN","GS","HM-SNV","MSI"))
    temp$Tumour_site <- factor(temp$Tumour_site,levels=c("right-sided","left-sided"))
    temp
  })
  head(COX.OS.byDriverGene.ClonalSubclonal[[1]])
  levels(COX.OS.byDriverGene.ClonalSubclonal[[1]]$TNM_stage)
  names(COX.OS.byDriverGene.ClonalSubclonal) <- OS.significant
  head(COX.OS.byDriverGene.ClonalSubclonal[[1]])
  library(mice)
  md.pattern(COX.OS.byDriverGene.ClonalSubclonal[[1]])
  library("VIM") #观察NA分布，发现Residual tumour非常多,所以删除Residual tumour
  dev.new();aggr(COX.OS.byDriverGene.ClonalSubclonal[[1]],prop=FALSE,numbers=TRUE)
  dim(na.omit(COX.OS.byDriverGene.ClonalSubclonal[[1]]))
}
# 进行Test.proportional.hazards分析
{
  ## 只对基因突变克隆状态进行等比例风险分析
  uni.proportional.hazards <- list()
  for(i in names(COX.OS.byDriverGene.ClonalSubclonal)){
    uni.proportional.hazards[[i]] <- Test.proportional.hazards(COX.OS.byDriverGene.ClonalSubclonal[[i]][,c(1,2,3,12)],paste0("F:/",i,".pdf"))
  }
  ## 对所有临床字段进行等比例风险假设分析
  proportional.hazards <- list()
  for(i in names(COX.OS.byDriverGene.ClonalSubclonal)){
    proportional.hazards[[i]] <- Test.proportional.hazards(COX.OS.byDriverGene.ClonalSubclonal[[i]],paste0("F:/",i,".pdf"))
  }
}
# 对每个基因单独进行OS单多因素COX分析 得到结果COX.OS.ClonalSubclonal.result
{
  COX.OS.ClonalSubclonal.result <- lapply(COX.OS.byDriverGene.ClonalSubclonal,function(x){
    Cox.function(x)
  })
  names(COX.OS.ClonalSubclonal.result) <- OS.significant
  save(COX.OS.ClonalSubclonal.result,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/COX.OS.ClonalSubclonal.result.RData")
}
