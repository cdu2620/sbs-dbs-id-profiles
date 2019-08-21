setwd('"/Users/cdu2620"')

library(tidyverse)
library(hrbrthemes)
library(scales)
library(ggpubr)
library(ggsci)
library(BSgenome.Hsapiens.UCSC.hg19)
library(janitor)
library(factoextra)
library(ggplot2)
library(dplyr)
library(rlist)
library(entropy)
library(MutationalPatterns)
library(devtools)
library(ggnomics)

mut_mat2 <- read.delim("breast_cancer_samples_example.SBS96.all", header= TRUE, sep="\t")
mut_mat2 <- mut_mat2[-6]
mut_matTest <- read.delim("breast_cancer_samples_example.SBS192.all", header=TRUE, sep="\t")
mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
                               package="MutationalPatterns"))
mut_mat3 <- read.delim("breast_cancer_samples_example.DBS78.all", header=TRUE, sep="\t")
mut_mat3 <- mut_mat3[-6]
mut_mat <- read.delim("breast_cancer_samples_example.ID94.all", header=TRUE, sep="\t")

COLORS6 = c(
  "#2EBAED", "#000000", "#DE1C14",
  "#D4D2D2", "#ADCC54", "#F0D0CE")

COLORS7 = c(
  "#2EBAED", "#000000", "#DE1C14",
  "#E98C7B", "#D4D2D2", "#ADCC54",
  "#F0D0CE")

COLORS10 = c(
  "#2EBAED", "#000000", "#DE1C14",
  "#D4D2D2", "#ADCC54", "#F0D0CE", "#1D2451", "#145E0F", "#FF9900",
  "#3D2384")

COLORS16 = c("#EB911B", "#FF9100", "#6DCE91", "#007C2D", "#FCEEF6", "#EC84A0",
             "#FF0404", "#8A0A0A", "#BAE2F6", "#72ACC9", 
             "#2384B5", "#054389", "#D8C1E9", "#A784C1",
             "#704492", "#380473")

SUBSTITUTIONS = c('C>A','C>G','C>T','T>A','T>C','T>G')
SUBSTITUTIONS_96 = rep(SUBSTITUTIONS, each=16)
SUBSTITUTIONS_192 = rep(SUBSTITUTIONS, each=32)

C_TRIPLETS = c(
  "ACA", "ACC", "ACG", "ACT",
  "CCA", "CCC", "CCG", "CCT",
  "GCA", "GCC", "GCG", "GCT",
  "TCA", "TCC", "TCG", "TCT")

T_TRIPLETS = c(
  "ATA", "ATC", "ATG", "ATT",
  "CTA", "CTC", "CTG", "CTT",
  "GTA", "GTC", "GTG", "GTT",
  "TTA", "TTC", "TTG", "TTT")

CONTEXTS_96 = c(rep(C_TRIPLETS, 3), rep(T_TRIPLETS, 3))

# combine substitutions and context in one 
TRIPLETS_96 = paste(substr(CONTEXTS_96,1,1), "[", SUBSTITUTIONS_96, "]", substr(CONTEXTS_96,3,3), sep = "")

STRAND = rep(c("U","T"), 96)
DNA_BASES = c("A", "C", "G", "T")

plot_96_profile2 = function(mut_matrix, colors, ymax = 100)
{
  # Relative contribution
  norm_mut_matrix = apply(mut_matrix, 2, function(x) x )
  
  # Check color vector length
  # Colors for plotting
  if(missing(colors)){colors=COLORS6}
  if(length(colors) != 6){stop("Provide colors vector with length 6")}
  context = CONTEXTS_96
  substitution = rep(SUBSTITUTIONS, each=16)
  
  # Replace mutated base with dot to get context
  substring(context, 2, 2) = "C"
  
  # Construct dataframe
  df = data.frame(substitution = substitution, context = context)
  rownames(norm_mut_matrix) = NULL
  df2 = cbind(df, as.data.frame(norm_mut_matrix))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  
  # These variables will be available at run-time, but not at compile-time.
  # To avoid compiling trouble, we initialize them to NULL.
  value = NULL
  names = as.character(unique(df3$variable))
  for (i in 1:length(names)) {
    names[[i]] = paste(names[[i]], ":")
  }
  df3$annotations = c(paste(names[[1]], sum(df3$value[1:96]), "subs"), rep("",95), paste(names[[2]], sum(df3$value[97:192]), "subs"), 
                      rep("",95), paste(names[[3]], sum(df3$value[193:288]), "subs"),
                      rep("",95), paste(names[[4]], sum(df3$value[289:384]), "subs"),
                      rep("",95))
    plot = ggplot(data=df3, aes(x=context,
                                y=value,
                                fill=substitution,
                                width=0.6)) +
      geom_bar(stat="identity", colour="black", size=.2) + 
      scale_fill_manual(values=colors) + 
      geom_text(aes(label=annotations), size=2.2, vjust=-7, hjust =-.2) +
      facet_grid(variable ~ substitution, scales="free") + 
      ylab("Number of Single Base Substitutions") +
      # no legend
      guides(fill=FALSE) +
      # white background
      theme_bw() +
      # format text
      theme(axis.title.y=element_text(size=12,vjust=1),
            axis.text.y=element_text(size=8),
            axis.title.x=element_text(size=12),
            axis.text.x=element_text(size=5,angle=90,vjust=.3),
            strip.text.x=element_text(size=9),
            strip.text.y=element_text(size=9),
            panel.grid.major.x = element_blank())
  
  ggarrange(plot)
  ggsave("test.pdf",height = 10,width = 20,device = cairo_pdf)
}
plot_96_profile2(mut_mat2[-1])

plot_dbs_78_profile2 = function(mut_matrix, colors)
{
  dfTest = mut_matrix[1]
  context = c()
  substitution = c()
  for (i in 1:nrow(dfTest)) {
    answer = substr(as.character(mut_matrix[[1]][[i]]), 1, 2)
    test = paste(answer, "> NN", sep=" ")
    substitution = append(substitution, test)
    dfTest$substitution[[i]] = test
    test2 = substr(as.character(mut_matrix[[1]][[i]]), 4, 5)
    context = append(context, test2)
    dfTest$context[[i]] = test2
  }
  dfTest = dfTest[-1]
  norm_mut_matrix = apply(mut_matrix[-1], 2, function(x) x )
  
  # Check color vector length
  # Colors for plotting
  if(missing(colors)){colors=COLORS10}
  if(length(colors) != 10){stop("Provide colors vector with length 10")}
  #context = CONTEXTS_96
  #substitution = rep(SUBSTITUTIONS, each=16)
  
  # Replace mutated base with dot to get context
  #substring(context, 2, 2) = "C"
  
  # Construct dataframe
  df2 = cbind(dfTest, as.data.frame(norm_mut_matrix))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  # These variables will be available at run-time, but not at compile-time.
  # To avoid compiling trouble, we initialize them to NULL.
  value = NULL
  names = as.character(unique(df3$variable))
  for (i in 1:length(names)) {
    names[[i]] = paste(names[[i]], ":")
  }
  df3$annotations = c(paste(names[[1]], sum(df3$value[1:78]), "double subs"), rep("",77), paste(names[[2]], sum(df3$value[79:156]), "double subs"), 
                           rep("",77), paste(names[[3]], sum(df3$value[157:234]), "double subs"),
                           rep("",77), paste(names[[4]], sum(df3$value[235:312]), "double subs"),
                           rep("",77))
    plot = ggplot(data=df3, aes(x=context,
                                y=value,
                                fill=(substitution),
                                width=0.6)) +
      geom_bar(stat="identity", colour="black", size=.2) + 
      scale_fill_manual(values=colors) + 
      facet_grid(variable ~ substitution, scales="free") + 
      geom_text(aes(label=annotations), size=3, vjust=-15, hjust =0) +
      ylab("Number of double base substitutions") + 
      # no legend
      guides(fill=FALSE) +
      # white background
      theme_bw() +
      # format text
      theme(axis.title.y=element_text(size=12,vjust=1),
            axis.text.y=element_text(size=8),
            axis.title.x=element_text(size=12),
            axis.text.x=element_text(size=5,angle=90,vjust=0.4),
            strip.text.x=element_text(size=9),
            strip.text.y=element_text(size=9),
            panel.grid.major.x = element_blank())
  ggarrange(plot)
  ggsave("test2.pdf",height = 10,width = 20,device = cairo_pdf)
}

plot_dbs_78_profile2(mut_mat3)

plot_192_profile2 = function(mut_matrix, colors, ymax = 0.2, condensed = FALSE)
{
  # Relative contribution
  
  # Check color vector length
  # Colors for plotting
  if(missing(colors)){colors=COLORS6}
  if(length(colors) != 6){stop("Provide colors vector with length 6")}
  
  norm_mut_matrix = apply(mut_matTest[1:96, 2:5], 2, function(x) x )
  
  context = CONTEXTS_96
  substitution = rep(SUBSTITUTIONS, each=16)
  
  # Replace mutated base with dot to get context
  substring(context, 2, 2) = "C"
  
  # Construct dataframe
  df = data.frame(substitution = substitution, context = context)
  rownames(norm_mut_matrix) = NULL
  df2 = cbind(df, as.data.frame(norm_mut_matrix))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  
  norm_mut_matrix = apply(mut_matTest[97:192, 2:5], 2, function(x) x )
  
  context = CONTEXTS_96
  substitution = rep(SUBSTITUTIONS, each=16)
  
  # Replace mutated base with dot to get context
  substring(context, 2, 2) = "C"
  
  # Construct dataframe
  df = data.frame(substitution = substitution, context = context)
  rownames(norm_mut_matrix) = NULL
  df2 = cbind(df, as.data.frame(norm_mut_matrix))
  df4 = melt(df2, id.vars = c("substitution", "context"))
  
  df5 = rbind(df3, df4)
  
  for (i in 1:768) {
    if (i >= 1 & i <=384) {
      df5$transcribe[[i]] = "Transcribed strand"
    } else {
      df5$transcribe[[i]] = "Untranscribed strand"
    }
  }
  # These variables will be available at run-time, but not at compile-time.
  # To avoid compiling trouble, we initialize them to NULL.
  value = NULL
  names = as.character(unique(df5$variable))
  for (i in 1:length(names)) {
    names[[i]] = paste(names[[i]], ":")
  }
  df5$annotations = c(paste(names[[1]], sum(df5$value[1:96], df5$value[385:480]), "transcribed subs"), rep("",191), paste(names[[2]], sum(df5$value[97:192], df5$value[481:576]), "transcribed subs"), 
                      rep("",191), paste(names[[3]], sum(df5$value[193:288], df5$value[577:672]), "transcribed subs"),
                      rep("",191), paste(names[[4]], sum(df5$value[289:384], df5$value[673:768]), "transcribed subs"),
                      rep("",191))
    plot = ggplot(data=df5, aes(x=context,
                                y=value,
                                color=df5$transcribe,
                                fill=df5$transcribe,
                                width=0.6)) +
      geom_bar(position = position_dodge(.9), stat="identity",  size=.2) +
      #scale_color_manual(values=colors) + 
      facet_grid(variable ~ substitution, scales ="free") + 
      ylab("Number of single base substitutions") + 
      geom_text(aes(label=annotations), size=3, vjust=-15, hjust =0) +
      # no legend
      guides(fill=FALSE) + 
      # white background
      theme_bw() +
      # format text
      theme(axis.title.y=element_text(size=12,vjust=1),
            axis.text.y=element_text(size=8),
            axis.title.x=element_text(size=12),
            axis.text.x=element_text(size=5,angle=90,vjust=0.4),
            strip.text.x=element_text(size=9),
            strip.text.y=element_text(size=9),
            panel.grid.major.x = element_blank())
  
  ggarrange(plot)
  ggsave("test3.pdf", width=20, height = 20)
}
plot_192_profile2(mut_matTest)

plot_id_profile <- function(mut_matrix, colors) {
  norm_mut_matrix = apply(mut_matrix[-1], 2, function(x) x )
  if(missing(colors)){colors=COLORS16}
  if(length(colors) != 16){stop("Provide colors vector with length 6")}
  df = mut_matrix[1]
  for (i in 1:83) {
    if (i >=1 & i <=6) {
      df$substitution[[i]] = "C"
      df$context[[i]] = as.character(i)
    } else if (i >= 7 & i <= 12) {
      df$substitution[[i]] = "T"
      df$context[[i]] = as.character(i - 6)
    } else if (i >=13 & i <= 18) {
      df$substitution[[i]] = "C"
      df$context[[i]] = as.character(i - 13)
    } else if (i >=19 & i <= 24) {
      df$substitution[[i]] = "T"
      df$context[[i]] = as.character(i - 19)
    } else if (i >=25 & i <= 30) {
      df$substitution[[i]] = "2"
      df$context[[i]] = as.character(i-24)
    } else if (i >=31 & i <= 36) {
      df$substitution[[i]] = "3"
      df$context[[i]] = as.character(i-30)
    } else if (i >=37 & i <= 42) {
      df$substitution[[i]] = "4"
      df$context[[i]] = as.character(i-36)
    } else if (i >=43 & i <= 48) {
      df$substitution[[i]] = "5+"
      df$context[[i]] = as.character(i-42)
    } else if (i >=49 & i <= 54) {
      df$substitution[[i]] = "2"
      df$context[[i]] = as.character(i-49)
    } else if (i >=55 & i <= 60) {
      df$substitution[[i]] = "3"
      df$context[[i]] = as.character(i-55)
    } else if (i >=61 & i <= 66) {
      df$substitution[[i]] = "4"
      df$context[[i]] = as.character(i-61)
    } else if (i >=67 & i <= 72) {
      df$substitution[[i]] = "5+"
      df$context[[i]] = as.character(i-67)
    } else if (i == 73) {
      df$substitution[[i]] = "2"
      df$context[[i]] = as.character(1)
    } else if (i >= 74 & i <= 75) {
      df$substitution[[i]] = "3"
      df$context[[i]] = as.character(i-73)
    } else if (i >=76 & i <= 78) {
      df$substitution[[i]] = "4"
      df$context[[i]] = as.character(i-75)
    } else {
      df$substitution[[i]] = "5+"
      df$context[[i]] = as.character(i-78)
    }
  }
  df = df[-1]
  dfTest = cbind(df, as.data.frame(norm_mut_matrix))
  testdata = melt(dfTest, id.vars = c("substitution", "context"))
  names = as.character(unique(testdata$variable))
  for (i in 1:length(names)) {
    names[[i]] = paste(names[[i]], ":")
  }
  testdata$annotations = c(paste(names[[1]], sum(testdata$value[1:83]), "indels"), rep("",85), paste(names[[2]], sum(testdata$value[84:166]), "indels"), 
                           rep("",80), paste(names[[3]], sum(testdata$value[167:249]), "indels"),
                           rep("",85), paste(names[[4]], sum(testdata$value[250:332]), "indels"),
                           rep("",78))
   value = NULL
   testdata$f3 = c(rep("1bp Deletion",12), rep("1bp Insertion",12), rep(">1bp Deletion at Repeats (Deletion Length)", 24),
                   rep(">1bp Insertions at Repeats (Insertion Length)", 24),
                   rep("Microhomology (Deletion Length)", 11),rep("1bp Deletion",12), rep("1bp Insertion",12), rep(">1bp Deletion at Repeats (Deletion Length)", 24),
                   rep(">1bp Insertions at Repeats (Insertion Length)", 24),
                   rep("Microhomology (Deletion Length)", 11),rep("1bp Deletion",12), rep("1bp Insertion",12), rep(">1bp Deletion at Repeats (Deletion Length)", 24),
                   rep(">1bp Insertions at Repeats (Insertion Length)", 24),
                   rep("Microhomology (Deletion Length)", 11),rep("1bp Deletion",12), rep("1bp Insertion",12), rep(">1bp Deletion at Repeats (Deletion Length)", 24),
                   rep(">1bp Insertions at Repeats (Insertion Length)", 24),
                   rep("Microhomology (Deletion Length)", 11))
   testdata$f3_f <- factor(testdata$f3, levels=as.character(unique(testdata$f3)))
  plot = ggplot(data=testdata, aes(x=context,
                              y=value,
                              fill=(substitution),
                              width=0.6)) +
    geom_bar(stat="identity", colour="black", size=.2) + 
    scale_fill_manual(values=colors) +
    facet_nested(variable ~ f3_f + substitution, scales="free") +
    #facet_grid(variable ~ substitution_f, scales= "free") + 
    ylab("Number of indels") + 
    geom_text(aes(x=3.5,y=100,label=annotations), size=2.2) +
    #coord_cartesian(ylim=c(0,100)) +
    #scale_y_continuous(breaks=seq(0, 100, 10)) +
    # no legend
    guides(fill=FALSE) +
    # white background
    theme_bw() +
    # format text
    theme(axis.title.y=element_text(size=12,vjust=1),
          axis.text.y=element_text(size=8),
          axis.title.x=element_text(size=12),
          axis.text.x=element_text(size=10),
          strip.text.x=element_text(size=9),
          strip.text.y=element_text(size=9),
          panel.grid.major.x = element_blank(),
          panel.spacing = unit(1, "mm"),
          panel.border = element_rect(color = "white", fill = NA, size = 1))
  ggarrange(plot)
  ggsave("test4.pdf",height = 20,width = 20,device = cairo_pdf)
}
plot_id_profile(mut_mat)
