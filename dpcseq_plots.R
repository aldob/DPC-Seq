# R script
library(ggplot2)
library(ggrepel)
library(gplots)
library(dplyr)
library(stringr)
library(reshape2)
library(pbapply)
library(viridis)
'%notin%' <- function(x,y) !('%in%'(x,y))
ttesting <- function(x, nonpar=FALSE) {
    wt0 <- as.numeric(c(x["wt_0_2"], x["wt_0_3"], x["wt_0_4"]))
    wt6 <- as.numeric(c(x["wt_6_2"], x["wt_6_3"], x["wt_6_4"]))
    csb0 <- as.numeric(c(x["csb_0_2"], x["csb_0_3"], x["csb_0_4"]))
    csb6 <- as.numeric(c(x["csb_6_2"], x["csb_6_3"], x["csb_6_4"]))

    if (nonpar==TRUE) {
        wt0csb0 <- wilcox.test(wt0, csb0, paired=TRUE)
        wt6csb6 <- wilcox.test(wt6, csb6, paired=TRUE)
        wt06 <- wilcox.test(wt0, wt6, paired=TRUE)
        csb06 <- wilcox.test(csb0, csb6, paired=TRUE)    
    } else {
        wt0csb0 <- t.test(wt0, csb0, paired=TRUE, var.equal=TRUE)
        wt6csb6 <- t.test(wt6, csb6, paired=TRUE, var.equal=TRUE)
        wt06 <- t.test(wt0, wt6, paired=TRUE, var.equal=TRUE)
        csb06 <- t.test(csb0, csb6, paired=TRUE, var.equal=TRUE)        
    }

    x["avg0_pval"] <- wt0csb0$p.value
    x["avg0_stat"] <- wt0csb0$statistic
    x["avg6_pval"] <- wt6csb6$p.value
    x["avg6_stat"] <- wt6csb6$statistic
    x["wt06_pval"] <- wt06$p.value
    x["wt06_stat"] <- wt06$statistic
    x["csb06_pval"] <- csb06$p.value
    x["csb06_stat"] <- csb06$statistic
    return(x)
}
theme3 <- theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    element_line(size = 0.5),
    axis.line = element_line(colour = 'black', size = 0.5),
    plot.title = element_text(hjust = 0.5, size = 10),
    title = element_text(size = 8),
    axis.text = element_text(size=8, colour='black'),
    axis.ticks = element_line(colour='black', size=0.5),
    axis.ticks.length=unit(0.1, "cm")
)
theme4 <- theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    element_line(size = 0.5),
    axis.line = element_line(colour = 'black', size = 0.5),
    plot.title = element_text(hjust = 0.5, size = 10),
    title = element_text(size = 8),
    axis.text = element_text(size=8, colour='black'),
    axis.ticks = element_line(colour='black', size=0.5),
    axis.ticks.length=unit(0.1, "cm"),
    strip.background = element_blank(),
    strip.text.x = element_blank()
)
theme5 <- theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    element_line(size = 0.5),
    axis.line = element_line(colour = 'black', size = 0.5),
    plot.title = element_text(hjust = 0.5, size = 10),
    title = element_text(size = 8),
    axis.text = element_text(size=8, colour='black'),
    axis.text.x = element_text(angle=-45, vjust=-0.2),
    axis.ticks = element_line(colour='black', size=0.5),
    axis.ticks.length=unit(0.1, "cm"),
    strip.background = element_blank(),
    strip.text.x = element_blank()
)
theme6 <- theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    element_line(size = 0.5),
    axis.line = element_line(colour = 'black', size = 0.5),
    plot.title = element_text(hjust = 0.5, size = 10),
    title = element_text(size = 8),
    axis.text = element_text(size=8, colour='black'),
    axis.ticks.y = element_line(colour='black', size=0.5),
    axis.ticks.x = element_blank(),
    axis.ticks.length=unit(0.1, "cm"),
    strip.background = element_blank(),
    strip.text.x = element_blank()
)

########################################## Import ##########################################

setwd("C:/Users/Aldo Bader/My Drive/science/Experiments/DPCs/dpc_seq/csb/analysis")
df <- read.delim("combined_new.cov", sep="\t")

df2 <- as.data.frame(t(pbapply(df, MARGIN=1, FUN=ttesting)))
df2[,c(2,3,8:ncol(df2))] <- lapply(df2[,c(2,3,8:ncol(df2))], as.numeric)

# write.table(df2, "combined_new_proc.cov", sep="\t", quote=FALSE)
# df2 <- read.delim("combined_new_proc.cov", sep="\t")










########################################## Dotplots ##########################################

df2$group <- 'none'
df2$group[df2["avg6_stat"]>0 & df2["avg6_pval"]<0.05 & abs(df2["csb6_avg"]/df2["wt6_avg"]-1)>0.2] <- "wt6"
df2$group[df2["avg6_stat"]<0 & df2["avg6_pval"]<0.05 & abs(df2["csb6_avg"]/df2["wt6_avg"]-1)>0.2] <- "csb6"
df2$group[abs(df2["csb6_avg"]/df2["wt6_avg"]-1)<=0.2] <- "none"
df2$group[df2["avg6_pval"]>0.99 & abs(df2["csb6_avg"]/df2["wt6_avg"]-1)<0.1] <- "absnone"

df2$act <- "low"
df2$act[df2$activity_pkb > 11.92] <- "mid-low"
df2$act[df2$activity_pkb > 29.03] <- "mid-high"
df2$act[df2$activity_pkb > 144.64] <- "high"

df2$atc <- "low"
df2$atc[df2$atac > 1000] <- "medium"
df2$atc[df2$atac > 5000] <- "high"

##### colour by group #####
plot_temp <- ggplot(arrange(df2, group), aes(y=wt6_avg, x=csb6_avg, col=group)) +
    theme3+
    guides(fill=FALSE, colour=FALSE, alpha=FALSE)+

    scale_x_log10(limits=c(0.03,5), expand=c(0,0), name="CSBko 6H DPC-seq coverage (RPM)")+
    scale_y_log10(limits=c(0.03,5), expand=c(0,0), name="WT 6H DPC-seq coverage (RPM)")+

    geom_point(size=2, aplha=0.3) + 
    scale_color_manual(values = c('grey60', '#ff5e00', '#005ef6')) + 

    geom_abline(slope=1, intercept=0, color='black', alpha=0.3, linetype='dashed')+

    ggsave(filename=("plots/final/biplot_group_6hr.png"), width = 3.2, height = 3, dpi=500, bg='white')



##### colour by activity #####

plot_temp <- ggplot(df2) +
    theme3+
    guides(fill=FALSE, colour=FALSE, alpha=FALSE)+
    scale_x_log10(limits=c(0.01,10), expand=c(0,0))+
    scale_y_log10(limits=c(0.01,10), expand=c(0,0))+
 
    geom_point(alpha = 0.1, aes(x=wt0_avg, y=wt6_avg, col=(-1)*(log10(activity_pkb)))) + 
    scale_color_distiller(palette = "YlOrRd")+

    geom_abline(slope=1, intercept=0, color='black', alpha=0.4, linetype="dashed")+

    ggsave(filename=("plots/biplot_act_wt06.png"), width = 3.2, height = 3, dpi=500, bg='white')



##### correlate features plot #####

plot_temp <- ggplot(filter(arrange(df2, desc(group)), group %in% c("absnone", "csb6")), aes(x=atac, y=activity_pkb)) +
    theme3+
    guides(fill=FALSE, colour=FALSE, alpha=FALSE)+

    scale_x_log10(name="Relative DNA accessibility (ATAC-seq)")+
    scale_y_log10(name="Relative RNA-Pol II occupancy (ChIP-seq)")+
    coord_cartesian(ylim=c(0.1,10000))+

    geom_point(aes(col=group), size=2, alpha=0.1) + 
    scale_color_manual(values = c('grey60', '#005ef6')) + 

    geom_hline(yintercept = 40, linetype='dashed', col='grey40')+
    geom_vline(xintercept = 1400, linetype='dashed', col='grey40')+

    ggsave(filename=("plots/final/biplot_atac_vs_act_groups_quads.png"), width = 3.2, height = 3, dpi=500, bg='white')










########################################## Boxplot ##########################################

boxdf <- data.frame()
interdf <- df2[,c("wt0_avg", "act", "atc")]
colnames(interdf) <- c("val", "act", "atc")
interdf$Sample <- "wt0"
boxdf <- rbind(boxdf, interdf)
interdf <- df2[,c("wt6_avg", "act", "atc")]
colnames(interdf) <- c("val", "act", "atc")
interdf$Sample <- "wt6"
boxdf <- rbind(boxdf, interdf)
interdf <- df2[,c("csb0_avg", "act", "atc")]
colnames(interdf) <- c("val", "act", "atc")
interdf$Sample <- "csb0"
boxdf <- rbind(boxdf, interdf)
interdf <- df2[,c("csb6_avg", "act", "atc")]
colnames(interdf) <- c("val", "act", "atc")
interdf$Sample <- "csb6"
boxdf <- rbind(boxdf, interdf)

plot_temp <- ggplot(filter(boxdf, Sample %in% c("csb0", "csb6", "wt0", "wt6")), aes(y=val, x=factor(Sample, levels=c("wt0", "csb0", "wt6", "csb6")), fill=factor(act, levels=c("low", "mid-low", "mid-high", "high")))) +
    theme5+
    guides(fill=FALSE, colour=FALSE, alpha=FALSE)+
    stat_boxplot(geom='errorbar',coef=1.5, width=0.7)+
    geom_boxplot(colour='black', outlier.shape = NA)+

    scale_y_continuous(name="DPC-seq coverage (RPM)")+
    facet_wrap(~factor(atc, levels=c("low", "mid-low", "mid-high", "high")), nrow=1)+
    coord_cartesian(ylim=c(0,1.04))+
    scale_fill_manual(values=c("yellow", "#ffbb00", "#ff8000", "red"))+

    ggsave(filename=("plots/boxplot_act_all.png"), width = 6, height = 3, dpi=500, bg='white')










########################################## Genome Browser Plots ##########################################

library(ggplot2)
library(dplyr)
library(gplots)
library(zoo)
library(dunn.test)
library(data.table)
'%notin%' <- function(x,y)!('%in%'(x,y))

metatheme2 <- theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    element_line(size = 1),
    axis.line = element_line(colour = 'black', size = 1),
    plot.title = element_text(hjust = 0.5, size = 30),
    title = element_text(face='bold', size = 15),
    axis.text = element_text(size=15, face="bold", colour='black'),
    axis.ticks = element_line(colour='black', size=1),
    axis.ticks.length=unit(0.15, "cm")
    )



region0 <- fread("Q:/seq/dpcseq/processed/regions/wt0_comb_r3.bed", sep="\t")[,1:4]
colnames(region0) <- c("chr", "coord1", "coord2", "Cov")
region0$Sample <- "RNA-Pol II"

region6 <- fread("Q:/seq/dpcseq/processed/regions/wt6_comb_r3.bed", sep="\t")[,1:4]
colnames(region6) <- c("chr", "coord1", "coord2", "Cov")
region6$Sample <- "RNA-Pol II"

rnpii <- fread("Q:/seq/dpcseq/processed/regions/rnpii_comb_r3.bed", sep="\t")[,1:4]
colnames(rnpii) <- c("chr", "coord1", "coord2", "Cov")
rnpii$Sample <- "RNA-Pol II"

region_comb <- rbind(region0, region6, rnpii)
rolled <- region_comb %>%
    group_by(Sample) %>%
    mutate(Cov = rollapply(Cov, 15, FUN=mean, partial=TRUE)) %>%
    select(chr, coord1, Cov, Sample)

rolled <- rbind(as.data.frame(rolled), c("chr6", 5, 16, "WT0"))
rolled <- rbind(rolled, c("chr6", 5, 16, "WT6"))



plot_temp <- ggplot(rolled) +
  metatheme2+
  geom_histogram(position=position_dodge(), stat='identity', aes(x=as.numeric(coord1), y=as.numeric(Cov), col=factor(Sample, levels = c("WT0", "WT6", "RNA-Pol II")))) +
  facet_wrap(~factor(Sample, levels=c("WT0", "WT6", "RNA-Pol II")), ncol=1, scales='free')+ 
  guides(color=FALSE)+
  scale_x_continuous(name = "", expand=c(0,0), breaks=c(53063000, 53075000, 53090000, 53104000), labels=c()) +
  scale_y_continuous(name = "DPC-seq coverage (RPM)", expand=c(0,0)) +
  coord_cartesian(xlim=c(53063000,53104000)) +
  theme(strip.text.x = element_blank()) +
  scale_colour_manual(values = c('#00aaff', '#003cff', "grey")) +
  ggsave(filename = "plots/final/browser_region3_wt.png", width = 8, height = 4, dpi=500, bg="white")





