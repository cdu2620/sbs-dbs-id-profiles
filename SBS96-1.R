
load(file = 'data.RData')
load(file="data1.RData")

sbs1 <- read.csv("sigProfiler_SBS_signatures_2019_05_22.csv")
dbs1 <- read.csv("sigProfiler_DBS_signatures.csv")

plot_sbs_96_profile <- function(data,samplename=NULL,totalmut=NULL,samplename_plot=TRUE, totalmut_plot=TRUE,filename=NULL){
  
  
  require(tidyverse)
  require(hrbrthemes)
  require(scales)
  require(ggpubr)  
  
  #### data is a frame with the following columns: Type, SubType, MutationType, Value ###
  
  if(is.null(samplename)){
    samplename <- colnames(data)[4]
  }
  
  colnames(data)[4] <- "Value"
  
  data <- data %>% mutate(Seq=seq_along(MutationType)) %>% mutate(Type=fct_inorder(Type),SubType=fct_inorder(SubType),MutationType=fct_inorder(MutationType))
  
  stype <- unique(data$Type)
  
  if(!is.null(totalmut)){
    totalmut <- comma_format()(totalmut)
  }
  
  COLORS6 = c(
    "#2EBAED", "#000000", "#DE1C14",
    "#D4D2D2", "#ADCC54", "#F0D0CE")
  names(COLORS6) <- stype
  
  ymax <- max(data$Value)/0.9
  ymax <- 0.04*ceiling(ymax/0.04)
  
  data0 <- bind_cols(
    data %>% group_by(Type) %>% filter(row_number() == 1) %>% ungroup() %>% mutate(Value1=Seq+0.25) %>% select(Value1,Type),
    data %>% group_by(Type) %>% filter(row_number() == 16) %>% ungroup() %>% mutate(Value2=Seq+0.25) %>% select(Value2)
  )
  ylabp <- percent(pretty_breaks(n = 5)(data$Value))
  ylabp <- ylabp[length(ylabp)]
  
  data00 <- data %>%  mutate(B1=str_sub(SubType,1,1),B2=str_sub(SubType,2,2),B3=str_sub(SubType,3,3))
  
  p1 <- data0 %>% ggplot(aes(xmin=Value1,xmax=Value2,ymin=0,ymax=1,fill=Type))+
    geom_rect()+
    scale_fill_manual(values=COLORS6)+
    scale_x_continuous(expand = c(0,0),sec.axis = dup_axis(labels = NULL,name = NULL),limits=c(1,96.5))+
    scale_y_continuous(expand = c(0,0),labels =ylabp,breaks = 1,limits = c(0,2.2),sec.axis = dup_axis(labels = NULL,name = NULL))+
    labs(x="",y=" ")+
    theme(axis.title.y=element_text(size=12),text = element_text(family = "Arial"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=10,color = "white"),
          axis.text.x=element_blank(),
          strip.text=element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid = element_blank(),
          #axis.line = element_blank(),
          axis.line.x = element_line(colour = 'white',size=0.6),
          axis.line.y = element_line(colour = 'white',size=0.6),
          plot.margin=margin(b=-1.2,l = 0.2,unit="cm")
    )+
    annotate("text", x = 8.75+16*0, y = 1.5, label = stype[1],size=5,fontface =2)+
    annotate("text", x = 8.75+16*1, y = 1.5, label = stype[2],size=5,fontface =2)+
    annotate("text", x = 8.75+16*2, y = 1.5, label = stype[3],size=5,fontface =2)+
    annotate("text", x = 8.75+16*3, y = 1.5, label = stype[4],size=5,fontface =2)+
    annotate("text", x = 8.75+16*4, y = 1.5, label = stype[5],size=5,fontface =2)+
    annotate("text", x = 8.75+16*5, y = 1.5, label = stype[6],size=5,fontface =2)
   p2 <- data %>% 
    ggplot(aes(Seq,Value,fill=Type))+
    geom_col(color="white",width = 0.5,size=0)+
    scale_fill_manual(values=COLORS6)+
    scale_x_continuous(expand = c(0,0),labels = data$SubType,breaks = data$Seq,sec.axis = dup_axis(labels = NULL,name = NULL))+
    scale_y_continuous(expand = c(0,0),labels = percent,breaks = seq(0,ymax,length.out = 5),limits = c(0,ymax),sec.axis = dup_axis(labels = NULL,name = NULL))+
    labs(x="",y="Percentage of Single Base Subsitutions")+
    theme(axis.title.y=element_text(size=12),text = element_text(family = "Arial"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=10),
          #axis.text.x=element_text(angle = 90,hjust = 1,vjust = 0.5),
          axis.text.x=element_blank(),
          strip.text=element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid.major.y = element_line(colour = "gray90",size = 0.3,linetype="solid"),
          panel.grid.minor.y = element_line(colour = "gray96",size = 0.3,linetype="solid"),
          panel.grid.major.x = element_blank(),
          axis.line.x = element_line(colour = 'gray80',size=0.6),
          axis.line.y = element_line(colour = 'gray80',size=0.6),
          plot.margin=margin(t=-1.2,b=-1,l=0.2,unit="cm")
          
    )
  
  if(samplename_plot & !is.null(samplename)){
    p2 <- p2+annotate("text", x = 1.5, y = ymax*0.9, label = samplename,size=7,fontface =2,hjust = 0)
  }
  
  if(totalmut_plot & !is.null(totalmut)){
    p2 <- p2+annotate("text", x = 78, y = ymax*0.88, label = paste0("Total Mutations: ",totalmut),size=6.5,fontface =2,hjust = 0)
  }
  
  invisible(capture.output(p2 <- flush_ticks(p2)+theme(axis.text.x = element_blank())))
  
  
  p3 <- data00 %>% 
    ggplot(aes(Seq))+
    geom_text(aes(y=0.4,label=B1),col="gray60",angle=90,hjust = 0,vjust = 0.5,size=3.5,family = "Arial")+
    geom_text(aes(y=1,label=B2,col=Type),angle=90,hjust = 0,vjust = 0.5,size=3.5,family = "Arial")+
    geom_text(aes(y=1.6,label=B3),col="gray60",angle=90,hjust = 0,vjust = 0.5,size=3.5,family = "Arial")+
    scale_color_manual(values=COLORS6)+
    scale_x_continuous(expand = c(0,0),sec.axis = dup_axis(labels = NULL,name = NULL),limits = c(0.5,96.5))+
    scale_y_continuous(expand = c(0,0),labels =ylabp,breaks = 1,limits = c(0,2.2),sec.axis = dup_axis(labels = NULL,name = NULL))+
    labs(x="",y=" ")+
    theme(axis.title.y=element_text(size=12),text = element_text(family = "Arial"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=10,colour = "white"),
          #axis.text.x=element_text(angle = 90,hjust = 1,vjust = 0.5),
          axis.text.x=element_blank(),
          strip.text=element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid.major.y = element_line(colour = "white",size = 0.4,linetype="solid"),
          panel.grid.minor.y = element_line(colour = "white",size = 0.4,linetype="solid"),
          panel.grid.major.x = element_blank(),
          axis.line.x = element_line(colour = 'white',size=0.6),
          axis.line.y = element_line(colour = 'white',size=0.6),
          plot.margin=margin(t=-1,l=0.2,unit="cm")
    )
  
  
  #cairo_pdf(file = 'tmp.pdf',width = 16,height = 4)
  pcomb <- ggarrange(p1,p2,p3,nrow = 3,align = "h",heights = c(1.8,10,1.5))
  
  if(is.null(filename)){
    return(pcomb)
  } else {
    ggsave(filename,width = 18,height = 4,device = cairo_pdf)
  }
  
}

plot_sbs_96_profile(data,filename = "tmpx.pdf",totalmut = 43145)

plot_dbs_78_profile <- function(data,samplename=NULL,totalmut=NULL,samplename_plot=TRUE, totalmut_plot=TRUE,filename=NULL, col){
  
  dataRep = data[1]
  dataRep$Type = paste0(substr(as.character(data$MutationType), 1, 3), "NN")
  dataRep$SubType = substr(as.character(data$MutationType), 4, 5)
  dataRep = dataRep[-1]
  dataRep$MutationType = data$MutationType
  temp = data[col]
  dataRep <- cbind(dataRep, temp)
  
  if(is.null(samplename)){
    samplename <- colnames(dataRep)[4]
  }

  colnames(dataRep)[4] <- "Value"
  data <- dataRep %>% mutate(Seq=seq_along(MutationType)) %>% mutate(Type=fct_inorder(Type),SubType=fct_inorder(SubType),MutationType=fct_inorder(MutationType))
  stype <- unique(data$Type)
  
  if(!is.null(totalmut)){
    totalmut <- comma_format()(totalmut)
  }
  
  COLORS10 = c(
    "#2EBAED", "#000000", "#DE1C14",
    "#D4D2D2", "#ADCC54", "#F0D0CE", "#1D2451", "#145E0F", "#FF9900",
    "#3D2384")
  names(COLORS10) <- stype
  
  ymax <- max(data$Value)
  
  data0 <- bind_cols(
    data %>% group_by(Type) %>% filter(row_number() == 1) %>% ungroup() %>% mutate(Value1=Seq+0.25) %>% select(Value1,Type),
    data %>% group_by(Type) %>% "["(.,c(which(data$Type != lag(data$Type))-1, 78),) %>% ungroup() %>% mutate(Value2=Seq+0.25) %>% select(Value2)
  )
  
  ylabp <- pretty_breaks(n = 5)(data$Value)
  ylabp <- ylabp[length(ylabp)]
  data00 <- data %>%  mutate(B1=str_sub(SubType,1,1),B2=str_sub(SubType,2,2),B3=str_sub(SubType,3,3))
  p1 <- data0 %>% ggplot(aes(xmin=Value1,xmax=Value2,ymin=0,ymax=1,fill=Type))+
    geom_rect()+
    scale_fill_manual(values=COLORS10)+
    #scale_x_continuous(expand = c(0,0),sec.axis = dup_axis(labels = NULL,name = NULL),limits=c(1,96.5))+
    scale_y_continuous(expand = c(0,0),breaks = 1,limits = c(0,2.2),sec.axis = dup_axis(labels = NULL,name = NULL))+
    labs(x="",y=" ")+
    theme(axis.title.y=element_text(size=12),text = element_text(family = "Arial"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=10,color = "white"),
          axis.text.x=element_blank(),
          strip.text=element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid = element_blank(),
          #axis.line = element_blank(),
          axis.line.x = element_line(colour = 'white',size=0.6),
          axis.line.y = element_line(colour = 'white',size=0.6),
          plot.margin=margin(b=-1.2,l = 0.2,unit="cm")
    )+
    annotate("text", x = 5+8*0, y = 1.5, label = stype[1],size=5,fontface =2)+
    annotate("text", x = 5+8*1, y = 1.5, label = stype[2],size=5,fontface =2)+
    annotate("text", x = 5+8*1.7, y = 1.5, label = stype[3],size=5,fontface =2)+
    annotate("text", x = 5+8*2.8, y = 1.5, label = stype[4],size=5,fontface =2)+
    annotate("text", x = 5+8*3.8, y = 1.5, label = stype[5],size=5,fontface =2)+
    annotate("text", x = 5+8*4.8, y = 1.5, label = stype[6],size=5,fontface =2)+
  annotate("text", x = 5+8*5.5, y = 1.5, label = stype[7],size=5,fontface =2)+
    annotate("text", x = 5+8*6.5, y = 1.5, label = stype[8],size=5,fontface =2)+
    annotate("text", x = 5+8*7.5, y = 1.5, label = stype[9],size=5,fontface =2)+
    annotate("text", x = 5+8*8.5, y = 1.5, label = stype[10],size=5,fontface =2)

  p2 <- data %>% 
    ggplot(aes(Seq,Value,fill=Type))+
    geom_col(color="white",width = 0.4,size=0)+
    scale_fill_manual(values=COLORS10)+
    scale_x_continuous(labels = data$SubType,breaks = data$Seq)+
    scale_y_continuous(expand = c(0,0),breaks = seq(0,ymax,length.out = 5),limits = c(0,ymax),sec.axis = dup_axis(labels = NULL,name = NULL))+
    labs(x="",y="Number of Double Base Subsitutions")+
    theme(axis.title.y=element_text(size=12),text = element_text(family = "Arial"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=10),
          #axis.text.x=element_text(angle = 90,hjust = 1,vjust = 0.5),
          axis.text.x=element_blank(),
          strip.text=element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid.major.y = element_line(colour = "gray90",size = 0.3,linetype="solid"),
          panel.grid.minor.y = element_line(colour = "gray96",size = 0.3,linetype="solid"),
          panel.grid.major.x = element_blank(),
          axis.line.x = element_line(colour = 'gray80',size=0.6),
          axis.line.y = element_line(colour = 'gray80',size=0.6),
          plot.margin=margin(t=-1.2,b=-1,l=0.2,unit="cm")
          
    )
  
  if(samplename_plot & !is.null(samplename)){
    p2 <- p2+annotate("text", x = 1.5, y = ymax*0.9, label = samplename,size=7,fontface =2,hjust = 0)
  }
  
  if(totalmut_plot & !is.null(totalmut)){
    p2 <- p2+annotate("text", x = 78, y = ymax*0.88, label = paste0("Total Mutations: ",totalmut),size=6.5,fontface =2,hjust = 0)
  }
  
  invisible(capture.output(p2 <- flush_ticks(p2)+theme(axis.text.x = element_blank())))
  
  
  p3 <- data00 %>% 
    ggplot(aes(Seq))+
    geom_text(aes(y=0.4,label=B1),col="gray60",angle=90,hjust = 0,vjust = 0.5,size=3.5,family = "Arial")+
    geom_text(aes(y=1,label=B2,col=Type),angle=90,hjust = 0,vjust = 0.5,size=3.5,family = "Arial")+
    geom_text(aes(y=1.6,label=B3),col="gray60",angle=90,hjust = 0,vjust = 0.5,size=3.5,family = "Arial")+
    scale_color_manual(values=COLORS10)+
    #scale_x_continuous(expand = c(0,0),sec.axis = dup_axis(labels = NULL,name = NULL),limits = c(0.5,96.5))+
    scale_y_continuous(expand = c(0,0),labels =ylabp,breaks = 1,limits = c(0,2.2),sec.axis = dup_axis(labels = NULL,name = NULL))+
    labs(x="",y=" ")+
    theme(axis.title.y=element_text(size=12),text = element_text(family = "Arial"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=10,colour = "white"),
          #axis.text.x=element_text(angle = 90,hjust = 1,vjust = 0.5),
          axis.text.x=element_blank(),
          strip.text=element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid.major.y = element_line(colour = "white",size = 0.4,linetype="solid"),
          panel.grid.minor.y = element_line(colour = "white",size = 0.4,linetype="solid"),
          panel.grid.major.x = element_blank(),
          axis.line.x = element_line(colour = 'white',size=0.6),
          axis.line.y = element_line(colour = 'white',size=0.6)
          #plot.margin=margin(t=-1,l=2.4,unit="cm")
    )
  
  
  #cairo_pdf(file = 'tmp.pdf',width = 16,height = 4)
  pcomb <- ggarrange(p1,p2,p3,nrow = 3,align = "h",heights = c(1.8,10,1.5))
  
  if(is.null(filename)){
    return(pcomb)
  } else {
    ggsave(filename,width = 18,height = 4,device = cairo_pdf)
  }
  
}
plot_dbs_78_profile(data_dbs,filename = "tmp2.pdf", col=2)

plot_id_profile <- function(data,samplename=NULL,totalmut=NULL,samplename_plot=TRUE, totalmut_plot=TRUE,filename=NULL, col){
  
  dataRep = data[1]
  for (i in 1:83) {
    dataRep$Type[[i]] = substr(as.character(data$MutationType[[i]]), 1, 7)
  }
  
  for (i in 1:83) {
    if (substr(as.character(data$MutationType[[i]]), 3, 5) == "Del") {
      if (substr(as.character(data$MutationType[[i]]), 7, 7) == "M") {
        if (substr(as.character(data$MutationType[[i]]), 9, 9) == "5") {
          dataRep$SubType[[i]] = "5+"
        } else {
          dataRep$SubType[[i]] = substr(as.character(data$MutationType[[i]]), 9, 9)
        }
      } else {if (substr(as.character(data$MutationType[[i]]), 9, 9) == "5") {
        dataRep$SubType[[i]] = "6+"
      } else {
        dataRep$SubType[[i]] = as.character(as.numeric(substr(as.character(data$MutationType[[i]]), 9, 9)) +1)
      }}
      
    } else {
      if (substr(as.character(data$MutationType[[i]]), 9, 9) == "5") {
        dataRep$SubType[[i]] = "5+"
      } else {
        dataRep$SubType[[i]] = substr(as.character(data$MutationType[[i]]), 9, 9)
      }
    }
  }
  
  dataRep = dataRep[-1]
  dataRep$MutationType = data$MutationType
  temp = data[col]
  dataRep <- cbind(dataRep, temp)
  
  if(is.null(samplename)){
    samplename <- colnames(dataRep)[4]
  }
  
  colnames(dataRep)[4] <- "Value"
  
  stuff = c("C", "T", "C", "T", "2", "3", "4", "5+", "2", 
            "3", "4", "5+", "2", "3", "4", "5+")
  
  coordinates = c(3, 10, 16, 22, 27, 34, 40, 45, 52, 57, 63, 
                  70, 73, 75, 77, 81)
  
  captions = c("Homopolymer length", "Homopolymer Length", 
               "Number of Repeat Units", "Number of Repeat Units",
               "Microhomology Length")
  
  data <- dataRep %>% mutate(Seq=seq_along(MutationType)) %>% mutate(Type=fct_inorder(Type),SubType=fct_inorder(SubType),MutationType=fct_inorder(MutationType))
  stype <- unique(data$Type)
   stype1 <- c("1bp Deletion", "1bp Insertion", ">1bp Deletion at Repeats (Deletion Length)",
           ">1bp Insertions at Repeats (Insertion Length)", "Microhomology (Deletion Length")
  
  if(!is.null(totalmut)){
    totalmut <- comma_format()(totalmut)
  }
  
  COLORS16 = c("#EB911B", "#FF9100", "#6DCE91", "#007C2D", "#FCEEF6", "#EC84A0",
               "#FF0404", "#8A0A0A", "#BAE2F6", "#72ACC9", 
               "#2384B5", "#054389", "#BE9BE6", "#A784C1",
               "#704492", "#380473")
  names(COLORS16) <- stype
  
  ymax <- max(data$Value)
  
  data0 <- bind_cols(
    data %>% group_by(Type) %>% filter(row_number() == 1) %>% ungroup() %>% mutate(Value1=Seq+0.25) %>% select(Value1,Type),
    data %>% group_by(Type) %>% "["(.,c(which(data$Type != lag(data$Type))-1, 83),) %>% ungroup() %>% mutate(Value2=Seq+0.25) %>% select(Value2)
  )
  ylabp <- pretty_breaks(n = 5)(data$Value)
  ylabp <- ylabp[length(ylabp)]
  
  p1 <- data0 %>% ggplot(aes(xmin=Value1,xmax=Value2,ymin=0,ymax=1,fill=Type))+
    geom_rect()+
    geom_text(aes(x= coordinates,y=0.4,label=stuff))+
    geom_label(aes(x= coordinates,y=0.4,label=stuff), fill="white")+
    scale_fill_manual(values=COLORS16)+
    scale_x_continuous(expand = c(0,0),sec.axis = dup_axis(labels = NULL,name = NULL),limits=c(0,83.5))+
    scale_y_continuous(expand = c(0,0),breaks = 1,limits = c(0,2.2),sec.axis = dup_axis(labels = NULL,name = NULL))+
    labs(x="",y=" ")+
    theme(axis.title.y=element_text(size=12),text = element_text(family = "Arial"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=10,color = "white"),
          axis.text.x=element_blank(),
          strip.text=element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid = element_blank(),
          #axis.line = element_blank(),
          axis.line.x = element_line(colour = 'white',size=0.6),
          axis.line.y = element_line(colour = 'white',size=0.6),
          plot.margin=margin(b=-1.2,l = 0.2,unit="cm")
    )+
    annotate("text", x = 5+8*0, y = 1.5, label = stype1[1],size=3,fontface =2)+
    annotate("text", x = 5+8*1.5, y = 1.5, label = stype1[2],size=3,fontface =2)+
    annotate("text", x = 5+8*4, y = 1.5, label = stype1[3],size=3,fontface =2)+
    annotate("text", x = 5+8*7, y = 1.5, label = stype1[4],size=3,fontface =2)+
    annotate("text", x = 5+8*9, y = 1.5, label = stype1[5],size=3,fontface =2)
  
  p1Copy <- data0 %>% ggplot(aes(xmin=Value1,xmax=Value2,ymin=0,ymax=1,fill=Type))+
    geom_rect()+
    scale_fill_manual(values=COLORS16)+
    scale_x_continuous(expand = c(0,0),sec.axis = dup_axis(labels = NULL,name = NULL),limits=c(0,83.5))+
    scale_y_continuous(expand = c(0,0),breaks = 1,limits = c(0,2.2),sec.axis = dup_axis(labels = NULL,name = NULL))+
    labs(x="",y=" ")+
    theme(axis.title.y=element_text(size=12),text = element_text(family = "Arial"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=10,color = "white"),
          axis.text.x=element_blank(),
          strip.text=element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid = element_blank(),
          #axis.line = element_blank(),
          axis.line.x = element_line(colour = 'white',size=0.6),
          axis.line.y = element_line(colour = 'white',size=0.6),
          plot.margin=margin(t=0,l = 0.2,unit="cm")
    )
  
  p2 <- data %>% 
    ggplot(aes(Seq,Value,fill=Type))+
    geom_col(color="white",width = 0.4,size=0)+
    scale_fill_manual(values=COLORS16)+
    scale_x_continuous(expand = c(0,0),labels = data$SubType,breaks = data$Seq,sec.axis = dup_axis(labels = NULL,name = NULL))+
    scale_y_continuous(expand = c(0,0),breaks = seq(0,ymax,length.out = 5),limits = c(0,ymax),sec.axis = dup_axis(labels = NULL,name = NULL))+
    labs(x="",y="Number of Indels")+
    theme(axis.title.y=element_text(size=12),text = element_text(family = "Arial"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=10),
          #axis.text.x=element_text(angle = 90,hjust = 1,vjust = 0.5),
          axis.text.x=element_blank(),
          strip.text=element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid.major.y = element_line(colour = "gray90",size = 0.3,linetype="solid"),
          panel.grid.minor.y = element_line(colour = "gray96",size = 0.3,linetype="solid"),
          panel.grid.major.x = element_blank(),
          axis.line.x = element_line(colour = 'gray80',size=0.6),
          axis.line.y = element_line(colour = 'gray80',size=0.6),
          plot.margin=margin(t=-1.2,b=0,l=0.2,unit="cm")
          
    )
  
  if(samplename_plot & !is.null(samplename)){
    p2 <- p2+annotate("text", x = 1.5, y = ymax*0.9, label = samplename,size=7,fontface =2,hjust = 0)
  }
  
  if(totalmut_plot & !is.null(totalmut)){
    p2 <- p2+annotate("text", x = 78, y = ymax*0.88, label = paste0("Total Mutations: ",totalmut),size=6.5,fontface =2,hjust = 0)
  }
  
  invisible(capture.output(p2 <- flush_ticks(p2)+theme(axis.text.x = element_blank())))
  p3 <- data %>% 
    ggplot(aes(Seq))+
    geom_text(aes(y=1.5,label=data$SubType),col="gray60",size=5,family = "Arial")+
    scale_color_manual(values=COLORS16)+
    scale_x_continuous(expand = c(0,0),sec.axis = dup_axis(labels = NULL,name = NULL),limits=c(0,83.5))+
    scale_y_continuous(expand = c(0,0),labels =ylabp,breaks = 1,limits = c(0,2.2),sec.axis = dup_axis(labels = NULL,name = NULL))+
    labs(x="",y=" ")+
    theme(axis.title.y=element_text(size=12),text = element_text(family = "Arial"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=10,colour = "white"),
          #axis.text.x=element_text(angle = 90,hjust = 1,vjust = 0.5),
          axis.text.x=element_blank(),
          strip.text=element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid.major.y = element_line(colour = "white",size = 0.4,linetype="solid"),
          panel.grid.minor.y = element_line(colour = "white",size = 0.4,linetype="solid"),
          panel.grid.major.x = element_blank(),
          axis.line.x = element_line(colour = 'white',size=0.6),
          axis.line.y = element_line(colour = 'white',size=0.6),
          plot.margin=margin(t=-1,l=.1,unit="cm")
    )+ annotate("text", x = 5+8*0, y = .5, label = captions[1],size=3,fontface =2)+
    annotate("text", x = 5+8*1.5, y = .5, label = captions[2],size=3,fontface =2)+
    annotate("text", x = 5+8*4, y = .5, label = captions[3],size=3,fontface =2)+
    annotate("text", x = 5+8*7, y = .5, label = captions[4],size=3,fontface =2)+
    annotate("text", x = 5+8*9, y = .5, label = captions[5],size=3,fontface =2)
  
  
  #cairo_pdf(file = 'tmp.pdf',width = 16,height = 4)
  pcomb <- ggarrange(p1,p2,p1Copy,p3,nrow = 4,align = "h",heights = c(3,10,3,4))
  
  if(is.null(filename)){
    return(pcomb)
  } else {
    ggsave(filename,width = 18,height = 4,device = cairo_pdf)
  }
  
}
plot_id_profile(data_indels,filename = "tmp3.pdf", col=2)

