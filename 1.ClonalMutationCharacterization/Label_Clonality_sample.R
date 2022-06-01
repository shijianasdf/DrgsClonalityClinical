#基于驱动基因克隆状态对样本分成3类(Clone,Subclone,WT)
grouping.byclonality.oneMarker <- function(marker.symbol, clonality.data, clinical.data)
{    
    #STEP 1: 处理样本，保留在突变数据和临床数据中共同出现的样本
    common.sample <- intersect(unique(clonality.data$Patient_ID), clinical.data$Patient_ID);
    clinical.data <- clinical.data[clinical.data$Patient_ID %in% common.sample, ];
    clonality.data <- clonality.data[clonality.data$Patient_ID %in% common.sample, ];
    
	  #STEP 2:识别marker基因是否出现在克隆性数据中
    markers <- intersect(marker.symbol,unique(clonality.data$Gene_Symbol));
    if(length(markers) == 0)
    {
  		sample.label <- rep(paste(marker.symbol,"WT"),length(common.sample));
  		return(cbind.data.frame(clinical.data,sample.label));
    }  
	
    #STEP 3: 得到克隆、亚克隆样本
    clonality.cellmarker.data <- unique(clonality.data[clonality.data$Gene_Symbol %in% marker.symbol,])   
    clone.sample <- unique(clonality.cellmarker.data[which(clonality.cellmarker.data$"CI95.timing" == "Clonal"),1])
    subclone.sample <- unique(clonality.cellmarker.data[which(clonality.cellmarker.data$"CI95.timing" == "Subclonal"),1])
    #####找到对于某个基因即是克隆又是亚克隆的样本,在我们分析中将这样的克隆亚克隆样本去掉,并且给样本打上标签
    temp.sample <- intersect(clone.sample,subclone.sample)
    if(length(temp.sample)==0){
  	 label1 <- c(paste(markers, "WT",sep = " "),paste(markers, "Clonal",sep = " "));
  	 label2 <- c(paste(markers, "WT",sep = " "),paste(markers, "Subclonal",sep = " "));
  	 level <- c(paste(markers, "WT",sep = " "),paste(markers, "Clonal",sep = " "),paste(markers, "Subclonal",sep = " ")); #设置因子顺序，对应画图颜色
  	 sample.label1 <- label1[as.numeric(clinical.data$Patient_ID %in% clone.sample)+1];
  	 sample.label2 <- label2[as.numeric(clinical.data$Patient_ID %in% subclone.sample)+1];
  	 pos2 <- grep("Subclonal",sample.label2);
  	 sample.label1[pos2] <- paste(markers,"Subclonal",sep=" ");
  	 clinical.data$'sample.label' <- factor(sample.label1,levels=level); #数据框新建一列	 
    }else{
     label1 <- c(paste(markers, "WT",sep = " "),paste(markers, "Clonal",sep = " "));
	   label2 <- c(paste(markers, "WT",sep = " "),paste(markers, "Subclonal",sep = " "));
     level <- c(paste(markers, "WT",sep = " "),paste(markers, "Clonal",sep = " "),paste(markers, "Subclonal",sep = " ")); #设置因子顺序，对应画图颜色	 
     sample.label1 <- label1[as.numeric(clinical.data$Patient_ID %in% clone.sample)+1];
	   sample.label2 <- label2[as.numeric(clinical.data$Patient_ID %in% subclone.sample)+1];
	   pos1 <- grep("Clonal",sample.label1);
	   pos2 <- grep("Subclonal",sample.label2);
	   Clonal.Subclonal.pos <- intersect(pos1,pos2);
	   sample.label1[pos2] <- paste(markers,"Subclonal",sep=" ");
	   sample.label1 <- sample.label1[-Clonal.Subclonal.pos] #删除对于该基因即有克隆突变又有亚克隆突变的标签
	   clinical.data <- clinical.data[!(clinical.data$Patient_ID %in% temp.sample),]; #删除对于该基因即有克隆突变又有亚克隆突变的样本
	   #sample.label1[Clonal.Subclonal.pos] <- paste(markers,"WT",sep = " "); 
	   clinical.data$'sample.label' <- factor(sample.label1,levels=level); #数据框新建一列	 
	}              
   return(clinical.data)
}
