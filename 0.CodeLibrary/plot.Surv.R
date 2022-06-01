#函数三 plot.surv
#' @description 绘制KM曲线并且进行logrank检验
#' @param clinical.data 具有分类标签的临床数据信息，至少包含四列：Patient_ID表示临床样本ID；event表示样本事件结局(数值型或逻辑型)；time表示样本事件的时间（数值型）；最后一列是通过上面的分类函数产生的分类标签sample.label，必须是factor
#' @param upper.time: 数值型，生存时间上限，必须与clinical.data中的单位一致，默认为NULL；如果使用，超过年限的样本将被去掉
#' @param xscale: 字符型，允许的选项包括"d_m"，"d_y"，"m_d"，"m_y"，"y_d"和"y_m"，其中d =天，m =月和y =年。 例如，xscale =“d_m”会将x轴单位从几天转换为几个月
#' @param unit.xlabel: x轴label展示的时间单位，四种选择：c("year", "month", "week", "day"),分别代表年，月，周，天，必须与最终图中展示使用的时间单位一致
#' @param surv.median.line: 在中位生存期时绘制水平/垂直线的字符向量，允许的值包括c（"none"，"hv"，"h"，"v"）中的一个，v：垂直，h：水平
#' @param risk.table: 逻辑向量，是否绘制risk.table，默认为TRUE
#' @param pval: 逻辑向量，是否给出log rank p值，默认为TRUE
#' @param survival.event: 结局事件的类型，允许的值包括c(“Overall Survival”, “Progress Free Survival”)等
#' @param main: 主标题的名字
#' @param inputFilePath: KMplot存储路径
#' @param picture.name: 存储的KMplot名字
#' @return 中位生存期，1,3,5年生存率，KM曲线图片

plot.surv <- function(clinical.data,
                      time = "time",
					  event = "event",
                      sample.label = "sample.label",
                      upper.time=NULL,
                      xscale = c( "d_m", "d_y", "m_d", "m_y", "y_d", "y_m"),#用于转换x轴坐标刻度
                      unit.xlabel = c("year", "month", "week", "day"),#所用生存时间的单位,
                      surv.median.line = c("none", "hv", "h", "v"),#是否画出中位生存时间，默认不给出
                      risk.table = c(TRUE, FALSE),#是否显示risk table
                      pval = c(TRUE, FALSE),#是否给出log-rank的p值
                      conf.int = c(FALSE, TRUE),#是否画出置信区间
                      main=NULL,#主标题名字
                      survival.event=c("Overall Survival","Progress Free Survival"),#事件类型
					            inputFilePath = NULL, #KMplot存储路径
					            picture.name = NULL  #存储的KMplot名字
)
{
    options(stringsAsFactors = FALSE);
    #载入生存分析包
    require(survival);
	  # ggplot2 survival plot ggsurvplot and pairwise_survdiff
    require(survminer);
    require(RColorBrewer);
	  # Provides a number of user-level functions to work with "grid" graphics, notably to arrange multiple grid-based plots on a page, and draw tables.
    require(gridExtra);
    
    survival.event <- survival.event[1];
    unit.xlabel <- unit.xlabel[1];
    
	# 重构生存数据
	clinical.data <- clinical.data[,c(time,event,sample.label)]## data1.FastSurvplot一定要为全局变量
    colnames(clinical.data)[1:3] <- c("time","event","sample.label")
    clinical.data$time <- as.numeric(as.character(clinical.data$time))
    clinical.data$event <- as.numeric(as.character(clinical.data$event))

    #去除生存时间超过upper.time的样本
    if(!is.null(upper.time))
    {
        clinical.data <- clinical.data[clinical.data$time <= upper.time,]
    }
    
    #选择日期格式 
    xSL <- data.frame(xScale=c(1,7,30,365.25),xLab=c("Days","Weeks","Months","Years"), stringsAsFactors=FALSE)
    switch(unit.xlabel, year={xScale <- 365.25;}, month={xScale <- 30;}, week={xScale <- 7;}, day={xScale <- 1})
    xLab <- xSL[which(xSL[,1]==xScale),2];
    
	
	
    #构造颜色
    t.name <- levels(clinical.data$sample.label);
    colors <- c("#808080","#E7B800","#2E9FDF","#FBBC05","#34A853","#000000");#顺序：灰，红，蓝，黄，绿，黑
    t.col <- colors[1:length(t.name)]; #palette = c("#E7B800", "#2E9FDF")
    
	
    print("fbhfth")
    #构造生存对象Surv(time, event)并且依据sample.label分层绘制KM生存曲线
    km.curves <- survfit(Surv(time, event)~sample.label, data=clinical.data)
    #km.curves <- survfit(Surv(time, event)~sample.label, data=OS.byDriverGene.ClonalSubclonal[[1]]);
	  
    #生成log-rank检验的p值,保留两位有效数字
	pval <- round(surv_pvalue(km.curves, clinical.data)$pval,digits=4)
    
	#计算中位生存期以及CI
    median.os <- surv_median(fit=km.curves)
	print("中位生存期")
	
	#计算5年生存率
    if(unit.xlabel == "month")
	  survival_rate <- summary(km.curves,times=c(12,36,60)) 
    if(unit.xlabel == "year")
      survival_rate <- summary(km.curves,times=c(1,3,5))
    if(unit.xlabel == "week")
      survival_rate <- summary(km.curves,times=c(104.35,156.53,260.89))
    if(unit.xlabel == "day")
      survival_rate <- summary(km.curves,times=c(365.25,1095.75,1826.25))
    print("五年生存率")
    #返回中位生存期和5年生存率
    if(pval <= 0.05){
      result <- list(median.os=median.os,survival_rate=survival_rate)
    }else{
      result <- NULL;
    }
      
	#构造生存图像中图例显示文字
    legend.content <- substr(names(km.curves$strata),start = 14,stop = 1000)
    print("legend.content")
  	#设置生存图片地址及名字
  	if(is.null(inputFilePath)){
  	   inputFilePath <- getwd()
  	}
  	if(!file.exists(inputFilePath)){
  	   dir.create(inputFilePath,recursive=T)
  	}
  	if(is.null(picture.name)){
  		picture.name <- "default"
  	}
  	
	# 存储地址
	pdf(paste0(inputFilePath,picture.name,"_",pval,".pdf"))
    cat(paste0(inputFilePath,picture.name,"_",pval,".pdf"))
	#ggsurvplot绘制生存图像
    ggsurv <- ggsurvplot(                         #注意此时并没有画图，只是存在ggsurv变量中
                         km.curves,               # survfit object with calculated statistics.
                         data = clinical.data,             # data used to fit survival curves.
                         palette = t.col,
                         
                         #图的主题构架
                         risk.table = risk.table[1],       # 是否展示风险table.
                         pval = pval[1],             # 是否展示log-rank检验的p值.
                         surv.median.line = surv.median.line[1],  # 是否展示中位生存期.surv_median
						             conf.int = conf.int[1], #是否画出置信区间
                         title = main,     #主标题名字
                         font.main = 15,       #主标题字体大小              
                         xlab = paste("Time","in",xLab,sep = " "),   # customize X axis label.
                         ylab = survival.event,   # customize Y axis label.
                         
                         #图例设置
                         legend.title = "", #图例标题，一般不用，设为空
                         legend.labs = legend.content, #图例文字描述
                         #legend = c(0.8,0.9), #图例的位置，取值在【0,1】之间
                         font.legend = 9,     #图例字体大小
                         
                         #risk table设置
                         tables.theme = theme_cleantable(),#table主题
                         risk.table.title = "No. at risk:",#table标题
                         risk.table.y.text.col = T, # 使用颜色代替Y轴文字
                         risk.table.y.text = FALSE, # Y轴不使用文字注释
                         tables.height = 0.15,      # table的高度
                         risk.table.fontsize = 3,    # risk table内文字的大小
						 ggtheme = theme_bw() # Change ggplot2 theme
                        );
						print("画图")
    ggsurv$plot <- ggsurv$plot + theme(plot.title = element_text(hjust = 0.5), plot.margin = unit(c(5.5, 5.5, 5.5, 50), "points"));#KM曲线图的标题居中
    ggsurv$table <- ggsurv$table + theme(plot.title = element_text(hjust = -0.04), plot.margin = unit(c(5.5, 5.5, 5.5, 50), "points"));#风险table图的标题居中
    print("画图")
    #判断分类的类数，如果只有两类，就不必计算两两之间的log rank p值
    if(length(t.name) > 2)
    {
    
    #计算pairwise的log rank的p值
    res <- pairwise_survdiff(Surv(time, event)~sample.label, data=clinical.data);
    pairwise.pvalue <- round(res$p.value, digits = 4);
    pairwise.pvalue[which(pairwise.pvalue < 0.0001)] <- "<0.0001";
    pairwise.pvalue[is.na(pairwise.pvalue)] <- "-";
    
    #添加表格
    tt <- ttheme_minimal(
                         core = list(fg_params = list(col = "black"),bg_params = list(fill = NA, col = "black")),
                         colhead = list(fg_params = list(col = NA),bg_params = list(fill = t.col, col = "black")),
                         rowhead = list(fg_params = list(col = NA, hjust = 1),bg_params = list(fill = c("white",t.col[-1]), col = "black"))
                        );
    pairwise.table <- tableGrob(pairwise.pvalue, theme = tt);
     print(ggarrange(ggarrange(ggsurv$plot, ggsurv$table,nrow=2,heights=c(2,0.5)), pairwise.table, nrow=2,heights = c(2,0.5),labels = c("","p from pairwise comparisons"), hjust = 0, font.label = list(size = 15, face = "plain"))); #最终的绘图函数ggarange(Arrange multiple ggplots on the same page)，包括3部分ggsurv$plot(KM生存曲线)，ggsurv$table(风险table)，成对比较table，print是为了让循环画图可以进行
    }
    else
    {
      print(ggsurv)
    }
	dev.off();
	#返回中位生存期和5年生存率
	return(result);
	#ggsave("shijian.pdf", device="pdf",dpi=300, plot=a, path="D:/");也可以用ggsave保存图片 
}
