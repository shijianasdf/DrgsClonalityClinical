#Assessing the relationship between clonality of driver genes and clinical variables
load("D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/CRC.clonal.subclonal.data.list.RData");
##对CRC.clonal.subclonal.data.list进行处理，挑选要用的临床变量来分析
Correlation.table.list <- lapply(CRC.clonal.subclonal.data.list,function(x){
	temp <- x[,c(2,3,4,5,6,7,9,17,18)]
	colnames(temp) <- c("sex","MSS_state","T_stage","N_stage","M_stage","AJCC_stage","age","Tumour_site",colnames(x)[18])
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
	temp$AJCC_stage[grep("Stage 1|Stage 2",temp$AJCC_stage)] <- "Stage 1/Stage 2"
	temp$AJCC_stage[grep("Stage 3|Stage 4",temp$AJCC_stage)] <- "Stage 3/Stage 4"
	temp$AJCC_stage <- factor(temp$AJCC_stage,levels=c("Stage 1/Stage 2","Stage 3/Stage 4"))
	#temp$Subtype <- factor(temp$Subtype,levels=c("CIN","GS","HM-SNV","MSI"))
	temp$Tumour_site <- factor(temp$Tumour_site,levels=c("right-sided","left-sided"))
	temp
})
# make table data need by fisher.test or chisq.test
Correlation.table.list <- lapply(Correlation.table.list,function(x){
	tablename <- colnames(x)
	chisq.fisher.list <- list()
	for(i in tablename[1:8]){
	   chisq.fisher.list[[paste(tablename[9],i,sep="&")]] <- unclass(table(x[,i],x[,tablename[9]]))
	}
	chisq.fisher.list
})
Filter.Correlation.table.list <- unlist(Filter.Correlation.table.list,recursive=F,use.names=T)
noncorrect.chisq.result <- lapply(Filter.Correlation.table.list, function(x){chisq.test(x, correct = F)})
correct.chisq.result <- lapply(Filter.Correlation.table.list, function(x){chisq.test(x, correct = T)})
fisher.result <- lapply(Filter.Correlation.table.list, function(x){fisher.test(x, workspace=100000000)})
fisher.p.result <-  lapply(fisher.result, function(x){x$p.value})
correct.p.result <- lapply(correct.chisq.result, function(x){x$p.value})
noncorrect.p.result <- lapply(noncorrect.chisq.result, function(x){x$p.value})
result <- cbind(names(noncorrect.chisq.result), do.call(rbind.data.frame, correct.p.result), do.call(rbind.data.frame, noncorrect.p.result), do.call(rbind.data.frame, fisher.p.result));
colnames(result) <- c("name", "correct.p", "noncorrect.p", "fisher.p");
save(result,file="D:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/CloneAndVariables.RData");
