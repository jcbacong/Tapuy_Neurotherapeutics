## Import packages
library(survival)
library("survminer")
library(dplyr)
library("ggplot2")
library(ggsignif)
library(ggpubr)
library(rstatix)
library(cowplot)
library("lemon")
library(reprex)
library(AICcmodavg)
library(DescTools)
library("plotrix")
library(data.table)
library("rockchalk")
library(truncnorm)
library("BSDA")
library(ggtext)
library("scales")

## Set working directory to the file location using Rstudio API
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

######################### READ FILES #############################
day5.csv = read.csv("Day 5-output.csv", header = T)
day5.csv$control = "OP50 Day 5"
day6.csv = read.csv("Acute-output.csv", header = T)
chronic.csv = read.csv("Chronic-output.csv", header = T)
acute.csv = rbind(day5.csv,day6.csv)
acute.csv = subset(acute.csv, select = -c(X))
chronic.csv = subset(chronic.csv, select = -c(X))


######################### ACUTE FILE #############################
main.acute.df <- acute.csv %>% 
  group_by(control,day) %>%
  summarise(mean.rawint = mean(RawIntDen),
            sem.rawint = std.error(RawIntDen),
            mean.Area = mean(Area),
            sem.Area = std.error(Area),
            obs.rawint = n()) %>%
  as.data.frame()

acute.df = main.acute.df
acute.df = subset(acute.df, !(control %in% c('Field Lees','Lab Lees','Lab Wine')))
acute.df$day = factor(acute.df$day, levels = c(5, 6), labels=c(5,6))

acute.df$control = factor(acute.df$control, levels=c("OP50 Day 5", "OP50", "DMSO", "EGCG",
                                             "Itogon Low", "Itogon Mid", "Itogon High",
                                             "Kapangan Low", "Kapangan Mid", "Kapangan High",
                                             "La Trinidad Low", "La Trinidad Mid", "La Trinidad High",
                                             "Sablan Low", "Sablan Mid", "Sablan High"),
                          
                          labels=c("OP50 <i>E. coli</i><br>(Day 5)", "OP50 <i>E. coli</i><br>(Day 6)", "0.1% DMSO", "200 uM EGCG",
                                   "Itogon Low", "Itogon Mid", "Itogon High",
                                   "Kapangan Low", "Kapangan Mid", "Kapangan High",
                                   "La Trinidad Low", "La Trinidad Mid", "La Trinidad High",
                                   "Sablan Low", "Sablan Mid", "Sablan High"))
acute.df <- acute.df %>% arrange(control)

color_scheme = c("#84d94f", "#ffa600", "#f25c78",  "#b64cdf", "#fa5cf2")

acute.df <- acute.df %>%
  mutate(type = case_when(
    grepl("E. coli|DMSO|EGCG", control) ~ "Control",
    grepl("Itogon", control) ~ "Itogon",
    grepl("Kapangan", control) ~ "Kapangan",
    grepl("La Trinidad", control) ~ "La Trinidad",
    grepl("Sablan", control) ~ "Sablan"
  ),
  concentration = case_when(
    grepl("Day 5", control) ~ "OP50 <i>E. coli</i><br>(Day 5)",
    grepl("Day 6", control) ~ "OP50 <i>E. coli</i><br>(Day 6)",
    grepl("DMSO", control) ~ "0.1% DMSO",
    grepl("EGCG", control) ~ "200 uM EGCG",
    grepl("Low", control) ~ "Low",
    grepl("Mid", control) ~ "Mid",
    grepl("High", control) ~ "High"
  ))

acute.df$type = factor(acute.df$type, levels=c("Control", "Itogon", "Kapangan", "La Trinidad", "Sablan"),
                    labels=c("Control", "Itogon", "Kapangan", "La Trinidad", "Sablan"))

acute.df$concentration = factor(acute.df$concentration, levels=c("OP50 <i>E. coli</i><br>(Day 5)",
                                                                 "OP50 <i>E. coli</i><br>(Day 6)",
                                                                 "0.1% DMSO", "200 uM EGCG",
                                                                 "Low", "Mid", "High"),
                       labels=c("OP50<br><i>E. coli</i><br>(Day 5)",
                                "OP50<br><i>E. coli</i><br>(Day 6)",
                                "0.1% <br>DMSO", "200 uM <br>EGCG",
                                "Low", "Mid", "High"))

write.csv(acute.df, file="acute-values.csv")
## Itogon
# ascorbic.acid = subset(df.csv.summary, Treatment == "Ascorbic acid")
q1 = ( ggplot(acute.df, aes(x=concentration, y=mean.rawint, fill=type)) +
         theme_classic()  + labs(x = '', y = 'Mean RFP Intensity (A.U)') +
         geom_bar(stat="identity", color="black") +
         geom_errorbar(aes(ymin=mean.rawint - sem.rawint, ymax=mean.rawint + sem.rawint), width=.2,position=position_dodge(.9)) +
         facet_grid('. ~ type', scales="free", switch = "x", space='free') +
         # # geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
         scale_fill_manual(values=color_scheme) +
         theme(strip.placement = "outside") +
         theme(axis.title.x = element_blank(),
               axis.title.y =element_text(size=15),
               plot.title=element_text(hjust=0.5),
               axis.text.x = element_markdown(size=8),
               axis.text.y = element_text(size=11),
               legend.position="none",
               strip.text.x = element_text(size=13, face = "bold"),
               strip.background = element_blank()) 

         # scale_y_continuous(limits = c(0,100))
       
       # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
       #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
       #       strip_text_x = element_text(size=12), strip_background = element_blank(),
       #       axis_line_y = element_line(colour = 'black', linetype='solid'))
       
       # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
       
)
q1

## Statistical tests
## 1. T-test between day 5 vs day 6 (OP50) not sig )p = 0.6596
dat.test = subset(acute.csv, control %in% c("OP50 Day 5","OP50"))
result = t.test(dat.test[dat.test$control=="OP50", "RawIntDen"],
             dat.test[dat.test$control=="OP50 Day 5", "RawIntDen"],
             alternative = "greater")
result
summary_table <- data.frame(
  Statistic = result$statistic,
  p_value = result$p.value,
  CI_lower = result$conf.int[1],
  CI_upper = result$conf.int[2]
)
write.csv(summary_table, file="ttest-stat-results.csv")



dat.test = subset(acute.csv, !(control %in% c("OP50 Day 5", "Field Lees", "Lab Wine", "Lab Lees")))
dat.test$control = factor(dat.test$control, levels=c( "OP50", "DMSO", "EGCG",
                                                     "Itogon Low", "Itogon Mid", "Itogon High",
                                                     "Kapangan Low", "Kapangan Mid", "Kapangan High",
                                                     "La Trinidad Low", "La Trinidad Mid", "La Trinidad High",
                                                     "Sablan Low", "Sablan Mid", "Sablan High"),
                          
                          labels=c("OP50", "0.1% DMSO", "200 uM EGCG",
                                   "Itogon Low", "Itogon Mid", "Itogon High",
                                   "Kapangan Low", "Kapangan Mid", "Kapangan High",
                                   "La Trinidad Low", "La Trinidad Mid", "La Trinidad High",
                                   "Sablan Low", "Sablan Mid", "Sablan High"))

res2 = DunnettTest(x=dat.test$RawIntDen, g=dat.test$control)
write.csv(as.data.frame(res2$`OP50`), file="dunnet-stat-results.csv")

####### VISUALIZATION WITH STAT
color_scheme = c("#84d94f", "#ffa600", "#f25c78",  "#b64cdf", "#fa5cf2")
dat = subset(acute.df, type == "Control")
h1 = ( ggplot(dat, aes(x=concentration, y=mean.rawint, fill=type)) +
         theme_classic()  + labs(x = '', y = 'Mean RFP Intensity (A.U)', title = "Control") +
         geom_bar(stat="identity", color="black") +
         geom_errorbar(aes(ymin=mean.rawint - sem.rawint, ymax=mean.rawint + sem.rawint), width=.2,position=position_dodge(.9)) +
         # facet_grid('. ~ type', scales="free", switch = "x", space='free') +
         # # geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
         scale_fill_manual(values=color_scheme[1]) +
         theme(axis.title.x = element_blank(),
               axis.title.y =element_text(size=15),
               plot.margin = unit(c(0, 0, 0, 0), "cm"),
               plot.title=element_text(hjust=0.5),
               axis.text.x = element_markdown(size=8),
               axis.text.y = element_text(size=11),
               legend.position="none") +
       
       scale_y_continuous(limits = c(0,60000))       # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
       #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
       #       strip_text_x = element_text(size=12), strip_background = element_blank(),
       #       axis_line_y = element_line(colour = 'black', linetype='solid'))
       
       # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
       
)
h1

op50 <- tibble(
  x = c(1, 1, 2, 2),
  y = c(48000, 52000, 52000, 46000)
) ## p = 0

h1 = h1 +
  geom_line(data = op50, aes(x= x, y = y, group=1), inherit.aes = F) +
  annotate("text", x = 1.5, y = 55000, label = 'n.s', size = 5, color = "black") +
  annotate("text", x = 4, y = 25000, label = '**', size = 5, color = "#960101")
h1

dat = subset(acute.df, type == "Itogon")
h2 = ( ggplot(dat, aes(x=concentration, y=mean.rawint, fill=type)) +
         theme_classic()  + labs(x = '', y = '', title = "Itogon") +
         geom_bar(stat="identity", color="black") +
         geom_errorbar(aes(ymin=mean.rawint - sem.rawint, ymax=mean.rawint + sem.rawint), width=.2,position=position_dodge(.9)) +
         # facet_grid('. ~ type', scales="free", switch = "x", space='free') +
         # # geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
         scale_fill_manual(values=color_scheme[2]) +
         theme(axis.title.x = element_blank(),
               axis.title.y =element_text(size=15),
               plot.title=element_text(hjust=0.5),
               axis.text.x = element_markdown(size=8),
               axis.text.y = element_text(size=11),
               legend.position="none") +
         scale_x_discrete(labels = c("Low","Mid", "High")) +
         scale_y_continuous(breaks = c(0, 20000, 40000, 60000),
                            labels = c("", "", "", ""),
                            limits = c(0,60000))+
         theme(legend.position="none",
               plot.title = element_text(hjust = 0.5),
               plot.margin = unit(c(0, 0, 0, 0), "cm"),
               axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))      # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
       #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
       #       strip_text_x = element_text(size=12), strip_background = element_blank(),
       #       axis_line_y = element_line(colour = 'black', linetype='solid'))
       
       # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
       
)
h2

h2 = h2 +
  annotate("text", x = 1, y = 7000, label = '***', size = 5, color = "#960101") +
  annotate("text", x = 2, y = 5500, label = '***', size = 5, color = "#960101") +
  annotate("text", x = 3, y = 6000, label = '***', size = 5, color = "#960101")
h2


dat = subset(acute.df, type == "Kapangan")
h3 = ( ggplot(dat, aes(x=concentration, y=mean.rawint, fill=type)) +
         theme_classic()  + labs(x = '', y = '', title = "Kapangan") +
         geom_bar(stat="identity", color="black") +
         geom_errorbar(aes(ymin=mean.rawint - sem.rawint, ymax=mean.rawint + sem.rawint), width=.2,position=position_dodge(.9)) +
         # facet_grid('. ~ type', scales="free", switch = "x", space='free') +
         # # geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
         scale_fill_manual(values=color_scheme[3]) +
         theme(axis.title.x = element_blank(),
               axis.title.y =element_text(size=15),
               plot.title=element_text(hjust=0.5),
               axis.text.x = element_markdown(size=8),
               axis.text.y = element_text(size=11),
               legend.position="none") +
         scale_x_discrete(labels = c("Low","Mid", "High")) +
         scale_y_continuous(breaks = c(0, 20000, 40000, 60000),
                            labels = c("", "", "", ""),
                            limits = c(0,60000))+
         theme(legend.position="none",
               plot.title = element_text(hjust = 0.5),
               plot.margin = unit(c(0, 0, 0, 0), "cm"),
               axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))      # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
       #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
       #       strip_text_x = element_text(size=12), strip_background = element_blank(),
       #       axis_line_y = element_line(colour = 'black', linetype='solid'))
       
       # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
       
)
h3

h3 = h3 +
  annotate("text", x = 1, y = 6300, label = '***', size = 5, color = "#960101") +
  annotate("text", x = 2, y = 8500, label = '***', size = 5, color = "#960101") +
  annotate("text", x = 3, y = 6500, label = '***', size = 5, color = "#960101")
h3

dat = subset(acute.df, type == "La Trinidad")
h4 = ( ggplot(dat, aes(x=concentration, y=mean.rawint, fill=type)) +
         theme_classic()  + labs(x = '', y = '', title = "La Trinidad") +
         geom_bar(stat="identity", color="black") +
         geom_errorbar(aes(ymin=mean.rawint - sem.rawint, ymax=mean.rawint + sem.rawint), width=.2,position=position_dodge(.9)) +
         # facet_grid('. ~ type', scales="free", switch = "x", space='free') +
         # # geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
         scale_fill_manual(values=color_scheme[4]) +
         theme(axis.title.x = element_blank(),
               axis.title.y =element_text(size=15),
               plot.title=element_text(hjust=0.5),
               axis.text.x = element_markdown(size=8),
               axis.text.y = element_text(size=11),
               legend.position="none") +
         scale_x_discrete(labels = c("Low","Mid", "High")) +
         scale_y_continuous(breaks = c(0, 20000, 40000, 60000),
                            labels = c("", "", "", ""),
                            limits = c(0,60000))+
         theme(legend.position="none",
               plot.title = element_text(hjust = 0.5),
               plot.margin = unit(c(0, 0, 0, 0), "cm"),
               axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))      # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
       #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
       #       strip_text_x = element_text(size=12), strip_background = element_blank(),
       #       axis_line_y = element_line(colour = 'black', linetype='solid'))
       
       # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
       
)
h4

h4 = h4 +
  annotate("text", x = 1, y = 9500, label = '***', size = 5, color = "#960101") +
  annotate("text", x = 2, y = 5000, label = '***', size = 5, color = "#960101") +
  annotate("text", x = 3, y = 5500, label = '***', size = 5, color = "#960101")
h4

dat = subset(acute.df, type == "Sablan")
h5 = ( ggplot(dat, aes(x=concentration, y=mean.rawint, fill=type)) +
         theme_classic()  + labs(x = '', y = '', title = "Sablan") +
         geom_bar(stat="identity", color="black") +
         geom_errorbar(aes(ymin=mean.rawint - sem.rawint, ymax=mean.rawint + sem.rawint), width=.2,position=position_dodge(.9)) +
         # facet_grid('. ~ type', scales="free", switch = "x", space='free') +
         # # geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
         scale_fill_manual(values=color_scheme[5]) +
         theme(axis.title.x = element_blank(),
               axis.title.y =element_text(size=15),
               plot.title=element_text(hjust=0.5),
               axis.text.x = element_markdown(size=8),
               axis.text.y = element_text(size=11),
               legend.position="none") +
         scale_x_discrete(labels = c("Low","Mid", "High")) +
         scale_y_continuous(breaks = c(0, 20000, 40000, 60000),
                            labels = c("", "", "", ""),
                            limits = c(0,60000))+
         theme(legend.position="none",
               plot.title = element_text(hjust = 0.5),
               plot.margin = unit(c(0, 0, 0, 0), "cm"),
               axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))      # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
       #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
       #       strip_text_x = element_text(size=12), strip_background = element_blank(),
       #       axis_line_y = element_line(colour = 'black', linetype='solid'))
       
       # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
       
)
h5

h5 = h5 +
  annotate("text", x = 1, y = 7000, label = '***', size = 5, color = "#960101") +
  annotate("text", x = 2, y = 8300, label = '***', size = 5, color = "#960101") +
  annotate("text", x = 3, y = 7300, label = '***', size = 5, color = "#960101")
h5


# ## Saving the figure for short term memory
acute.plot = plot_grid(h1, h2, h3, h4, h5,  align = "h", ncol = 5, rel_widths = c(0.28, 0.18, 0.18,0.18,0.18)) +
  theme(plot.margin = unit(c(0.0,0,0,0), "cm")) 

png(filename = file.path(getwd(),"/Acute-All.png"), width = 9, height = 3.5, units = "in", res = 600)
plot(acute.plot)
dev.off()
print(paste("Saved at ", getwd()))
print(acute.plot)


######################### Chronic FILE #############################
main.chronic.df <- chronic.csv %>% 
  group_by(control,day) %>%
  summarise(mean.rawint = mean(RawIntDen),
            sem.rawint = std.error(RawIntDen),
            mean.Area = mean(Area),
            sem.Area = std.error(Area),
            obs.rawint = n()) %>%
  as.data.frame()

main.chronic.df <- main.chronic.df %>% add_row(control='DMSO', 
                                             day=6, 
                                             mean.rawint = 0,
                                            sem.rawint = 0)
main.chronic.df <- main.chronic.df %>% add_row(control='EGCG', 
                                     day=6, 
                                     mean.rawint = 0,
                                     sem.rawint = 0)
main.chronic.df <- main.chronic.df %>% add_row(control='Kapangan Low', 
                                     day=6, 
                                     mean.rawint = 0,
                                     sem.rawint = 0)

chronic.df = main.chronic.df
chronic.df = subset(chronic.df, !(control %in% c('Field Lees','Lab Lees','Lab Wine')))

chronic.df$control = factor(chronic.df$control, levels=c("OP50", "DMSO", "EGCG",
                                                     "Itogon Low", "Itogon Mid", "Itogon High",
                                                     "Kapangan Low", "Kapangan Mid", "Kapangan High",
                                                     "La Trinidad Low", "La Trinidad Mid", "La Trinidad High",
                                                     "Sablan Low", "Sablan Mid", "Sablan High"),
                          
                          labels=c( "OP50 <i>E. coli</i><br>(Day 6)", "0.1% DMSO", "200 uM EGCG",
                                   "Itogon Low", "Itogon Mid", "Itogon High",
                                   "Kapangan Low", "Kapangan Mid", "Kapangan High",
                                   "La Trinidad Low", "La Trinidad Mid", "La Trinidad High",
                                   "Sablan Low", "Sablan Mid", "Sablan High"))

chronic.df <- chronic.df %>% arrange(control)
write.csv(chronic.df, file="chronic-values.csv")

color_scheme = c("#84d94f", "#ffa600", "#f25c78",  "#b64cdf", "#fa5cf2")

chronic.df <- chronic.df %>%
  mutate(type = case_when(
    grepl("E. coli|DMSO|EGCG", control) ~ "Control",
    grepl("Itogon", control) ~ "Itogon",
    grepl("Kapangan", control) ~ "Kapangan",
    grepl("La Trinidad", control) ~ "La Trinidad",
    grepl("Sablan", control) ~ "Sablan"
  ),
  concentration = case_when(
    grepl("Day 5", control) ~ "OP50 <i>E. coli</i><br>(Day 5)",
    grepl("Day 6", control) ~ "OP50 <i>E. coli</i><br>(Day 6)",
    grepl("DMSO", control) ~ "0.1% DMSO",
    grepl("EGCG", control) ~ "200 uM EGCG",
    grepl("Low", control) ~ "Low",
    grepl("Mid", control) ~ "Mid",
    grepl("High", control) ~ "High"
  ))

chronic.df$type = factor(chronic.df$type, levels=c("Control", "Itogon", "Kapangan", "La Trinidad", "Sablan"),
                       labels=c("Control", "Itogon", "Kapangan", "La Trinidad", "Sablan"))

chronic.df$concentration = factor(chronic.df$concentration, levels=c("OP50 <i>E. coli</i><br>(Day 5)",
                                                                 "OP50 <i>E. coli</i><br>(Day 6)",
                                                                 "0.1% DMSO", "200 uM EGCG",
                                                                 "Low", "Mid", "High"),
                                labels=c("OP50<br><i>E. coli</i><br>(Day 5)",
                                         "OP50<br><i>E. coli</i><br>(Day 6)",
                                         "0.1% <br>DMSO", "200 uM <br>EGCG",
                                         "Low", "Mid", "High"))

## Itogon
q2 = ( ggplot(chronic.df, aes(x=concentration, y=mean.rawint, fill=type)) +
         theme_classic()  + labs(x = '', y = 'Mean RFP Intensity (A.U)') +
         geom_bar(stat="identity", color="black") +
         geom_errorbar(aes(ymin=mean.rawint - sem.rawint, ymax=mean.rawint + sem.rawint), width=.2,position=position_dodge(.9)) +
         facet_grid('. ~ type', scales="free", switch = "x", space='free') +
         # # geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
         scale_fill_manual(values=color_scheme) +
         theme(strip.placement = "outside") +
         theme(axis.title.x = element_blank(),
               axis.title.y =element_text(size=15),
               plot.title=element_text(hjust=0.5),
               axis.text.x = element_markdown(size=8),
               axis.text.y = element_text(size=11),
               legend.position="none",
               strip.text.x = element_text(size=13, face = "bold"),
               strip.background = element_blank()) 
       
       # scale_y_continuous(limits = c(0,100))
       
       # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
       #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
       #       strip_text_x = element_text(size=12), strip_background = element_blank(),
       #       axis_line_y = element_line(colour = 'black', linetype='solid'))
       
       # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
       
)
q2

#### STATISTICAL TEST
dat.test = subset(chronic.csv, !(control %in% c("Field Lees", "Lab Wine", "Lab Lees")))
dat.test$control = factor(dat.test$control, levels=c( "OP50", "DMSO", "EGCG",
                                                      "Itogon Low", "Itogon Mid", "Itogon High",
                                                      "Kapangan Low", "Kapangan Mid", "Kapangan High",
                                                      "La Trinidad Low", "La Trinidad Mid", "La Trinidad High",
                                                      "Sablan Low", "Sablan Mid", "Sablan High"),
                          
                          labels=c("OP50", "0.1% DMSO", "200 uM EGCG",
                                   "Itogon Low", "Itogon Mid", "Itogon High",
                                   "Kapangan Low", "Kapangan Mid", "Kapangan High",
                                   "La Trinidad Low", "La Trinidad Mid", "La Trinidad High",
                                   "Sablan Low", "Sablan Mid", "Sablan High"))

res2 = DunnettTest(x=dat.test$RawIntDen, g=dat.test$control)
res2
write.csv(as.data.frame(res2$`OP50`), file="dunnet-stat-results-chronic.csv")



####### VISUALIZATION WITH STAT
color_scheme = c("#84d94f", "#ffa600", "#f25c78",  "#b64cdf", "#fa5cf2")
dat = subset(chronic.df, type == "Control")
h1 = ( ggplot(dat, aes(x=concentration, y=mean.rawint, fill=type)) +
         theme_classic()  + labs(x = '', y = 'Mean RFP Intensity (A.U)', title = "Control") +
         geom_bar(stat="identity", color="black") +
         geom_errorbar(aes(ymin=mean.rawint - sem.rawint, ymax=mean.rawint + sem.rawint), width=.2,position=position_dodge(.9)) +
         # facet_grid('. ~ type', scales="free", switch = "x", space='free') +
         # # geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
         scale_fill_manual(values=color_scheme[1]) +
         theme(axis.title.x = element_blank(),
               axis.title.y =element_text(size=15),
               plot.margin = unit(c(0, 0, 0, 0), "cm"),
               plot.title=element_text(hjust=0.5),
               axis.text.x = element_markdown(size=8),
               axis.text.y = element_text(size=11),
               legend.position="none") +
         
         scale_y_continuous(limits = c(0,35000))       # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
       #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
       #       strip_text_x = element_text(size=12), strip_background = element_blank(),
       #       axis_line_y = element_line(colour = 'black', linetype='solid'))
       
       # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
       
)
h1

dat = subset(chronic.df, type == "Itogon")
h2 = ( ggplot(dat, aes(x=concentration, y=mean.rawint, fill=type)) +
         theme_classic()  + labs(x = '', y = '', title = "Itogon") +
         geom_bar(stat="identity", color="black") +
         geom_errorbar(aes(ymin=mean.rawint - sem.rawint, ymax=mean.rawint + sem.rawint), width=.2,position=position_dodge(.9)) +
         # facet_grid('. ~ type', scales="free", switch = "x", space='free') +
         # # geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
         scale_fill_manual(values=color_scheme[2]) +
         theme(axis.title.x = element_blank(),
               axis.title.y =element_text(size=15),
               plot.title=element_text(hjust=0.5),
               axis.text.x = element_markdown(size=8),
               axis.text.y = element_text(size=11),
               legend.position="none") +
         scale_x_discrete(labels = c("Low","Mid", "High")) +
         scale_y_continuous(breaks = c(0, 10000, 20000, 30000),
                            labels = c("", "", "", ""),
                            limits = c(0,35000))+
         theme(legend.position="none",
               plot.title = element_text(hjust = 0.5),
               plot.margin = unit(c(0, 0, 0, 0), "cm"),
               axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))      # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
       #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
       #       strip_text_x = element_text(size=12), strip_background = element_blank(),
       #       axis_line_y = element_line(colour = 'black', linetype='solid'))
       
       # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
       
)
h2


dat = subset(chronic.df, type == "Kapangan")
h3 = ( ggplot(dat, aes(x=concentration, y=mean.rawint, fill=type)) +
         theme_classic()  + labs(x = '', y = '', title = "Kapangan") +
         geom_bar(stat="identity", color="black") +
         geom_errorbar(aes(ymin=mean.rawint - sem.rawint, ymax=mean.rawint + sem.rawint), width=.2,position=position_dodge(.9)) +
         # facet_grid('. ~ type', scales="free", switch = "x", space='free') +
         # # geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
         scale_fill_manual(values=color_scheme[3]) +
         theme(axis.title.x = element_blank(),
               axis.title.y =element_text(size=15),
               plot.title=element_text(hjust=0.5),
               axis.text.x = element_markdown(size=8),
               axis.text.y = element_text(size=11),
               legend.position="none") +
         scale_x_discrete(labels = c("Low","Mid", "High")) +
         scale_y_continuous(breaks = c(0, 10000, 20000, 30000),
                            labels = c("", "", "", ""),
                            limits = c(0,35000))+
         theme(legend.position="none",
               plot.title = element_text(hjust = 0.5),
               plot.margin = unit(c(0, 0, 0, 0), "cm"),
               axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))      # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
       #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
       #       strip_text_x = element_text(size=12), strip_background = element_blank(),
       #       axis_line_y = element_line(colour = 'black', linetype='solid'))
       
       # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
       
)
h3

dat = subset(chronic.df, type == "La Trinidad")
h4 = ( ggplot(dat, aes(x=concentration, y=mean.rawint, fill=type)) +
         theme_classic()  + labs(x = '', y = '', title = "La Trinidad") +
         geom_bar(stat="identity", color="black") +
         geom_errorbar(aes(ymin=mean.rawint - sem.rawint, ymax=mean.rawint + sem.rawint), width=.2,position=position_dodge(.9)) +
         # facet_grid('. ~ type', scales="free", switch = "x", space='free') +
         # # geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
         scale_fill_manual(values=color_scheme[4]) +
         theme(axis.title.x = element_blank(),
               axis.title.y =element_text(size=15),
               plot.title=element_text(hjust=0.5),
               axis.text.x = element_markdown(size=8),
               axis.text.y = element_text(size=11),
               legend.position="none") +
         scale_x_discrete(labels = c("Low","Mid", "High")) +
         scale_y_continuous(breaks = c(0, 10000, 20000, 30000),
                            labels = c("", "", "", ""),
                            limits = c(0,35000))+
         theme(legend.position="none",
               plot.title = element_text(hjust = 0.5),
               plot.margin = unit(c(0, 0, 0, 0), "cm"),
               axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))      # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
       #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
       #       strip_text_x = element_text(size=12), strip_background = element_blank(),
       #       axis_line_y = element_line(colour = 'black', linetype='solid'))
       
       # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
       
)
h4


dat = subset(chronic.df, type == "Sablan")
h5 = ( ggplot(dat, aes(x=concentration, y=mean.rawint, fill=type)) +
         theme_classic()  + labs(x = '', y = '', title = "Sablan") +
         geom_bar(stat="identity", color="black") +
         geom_errorbar(aes(ymin=mean.rawint - sem.rawint, ymax=mean.rawint + sem.rawint), width=.2,position=position_dodge(.9)) +
         # facet_grid('. ~ type', scales="free", switch = "x", space='free') +
         # # geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
         scale_fill_manual(values=color_scheme[5]) +
         theme(axis.title.x = element_blank(),
               axis.title.y =element_text(size=15),
               plot.title=element_text(hjust=0.5),
               axis.text.x = element_markdown(size=8),
               axis.text.y = element_text(size=11),
               legend.position="none") +
         scale_x_discrete(labels = c("Low","Mid", "High")) +
         scale_y_continuous(breaks = c(0, 10000, 20000, 30000),
                            labels = c("", "", "", ""),
                            limits = c(0,35000))+
         theme(legend.position="none",
               plot.title = element_text(hjust = 0.5),
               plot.margin = unit(c(0, 0, 0, 0), "cm"),
               axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))      # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
       #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
       #       strip_text_x = element_text(size=12), strip_background = element_blank(),
       #       axis_line_y = element_line(colour = 'black', linetype='solid'))
       
       # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
       
)
h5


# ## Saving the figure for short term memory
chronic.plot = plot_grid(h1, h2, h3, h4, h5,  align = "h", ncol = 5, rel_widths = c(0.24, 0.19, 0.19,0.19,0.19)) +
  theme(plot.margin = unit(c(0.0,0,0,0), "cm")) 

png(filename = file.path(getwd(),"/Chronic-All.png"), width = 9, height = 3.5, units = "in", res = 600)
plot(chronic.plot)
dev.off()
print(paste("Saved at ", getwd()))
print(chronic.plot)

########################################################################### 
########################### BCH 205 Results  ##############################
######################### Junelle & Mark ACUTE ############################
###########################################################################  
jm.acute.df = main.acute.df
jm.acute.df = subset(jm.acute.df, (control %in% c('OP50','OP50 Day 5','DMSO','EGCG', 'Itogon High', 'Field Lees','Lab Lees','Lab Wine')))

jm.acute.df$control = factor(jm.acute.df$control, levels=c("OP50 Day 5", "OP50", "DMSO", "EGCG",
                                                    "Itogon High",
                                                    'Lab Wine',
                                                    'Field Lees','Lab Lees'),
                          
                          labels=c("OP50 <i>E. coli</i><br>(Day 5)", "OP50 <i>E. coli</i><br>(Day 6)", "0.1% DMSO", "200 uM EGCG",
                                   "Field Wine",
                                   "Lab Wine",
                                   "Field Lees",
                                   "Lab Lees"))
jm.acute.df <- jm.acute.df %>% arrange(control)

color_scheme = c("#84d94f", "#ffa600", "#f25c78",  "#b64cdf", "#fa5cf2")

jm.acute.df <- jm.acute.df %>%
  mutate(type = case_when(
    grepl("E. coli|DMSO|EGCG", control) ~ "Control",
    grepl("Lab", control) ~ "Lab",
    grepl("Field", control) ~ "Field"
  ),
  concentration = case_when(
    grepl("Day 5", control) ~ "OP50 <i>E. coli</i><br>(Day 5)",
    grepl("Day 6", control) ~ "OP50 <i>E. coli</i><br>(Day 6)",
    grepl("DMSO", control) ~ "0.1% DMSO",
    grepl("EGCG", control) ~ "200 uM EGCG",
    grepl("Wine", control) ~ "Wine",
    grepl("Lees", control) ~ "Lees"
  ))

jm.acute.df$type = factor(jm.acute.df$type, levels=c("Control", "Field", "Lab"),
                         labels=c("Control", "Field", "Lab"))

jm.acute.df$concentration = factor(jm.acute.df$concentration, levels=c("OP50 <i>E. coli</i><br>(Day 5)",
                                                                     "OP50 <i>E. coli</i><br>(Day 6)",
                                                                     "0.1% DMSO", "200 uM EGCG",
                                                                     "Wine", "Lees"),
                                  labels=c("OP50 <i>E. coli</i><br>(Day 5)",
                                           "OP50 <i>E. coli</i><br>(Day 6)",
                                           "0.1% <br>DMSO", "200uM <br>EGCG",
                                           "Wine", "Lees"))


## Statistics
dat.test = subset(acute.csv, (control %in% c("OP50","DMSO", "EGCG",
                                             "Itogon High","Lab Wine",
                                             "Field Lees",  "Lab Lees")))

dat.test$control = factor(dat.test$control, levels=c("OP50","DMSO", "EGCG",
                                                     "Itogon High","Lab Wine",
                                                     "Field Lees",  "Lab Lees"),
                          
                          labels=c("OP50","DMSO", "EGCG",
                                   "Field Wine","Lab Wine",
                                   "Field Lees",  "Lab Lees"))
dat.test = arrange(dat.test, control)
res2 = DunnettTest(x=dat.test$RawIntDen, g=dat.test$control)
write.csv(as.data.frame(res2$`OP50`), file="BCH205-dunnet-stat-results.csv")
write.csv(jm.acute.df, "BCH205-acute-values.csv")

### Facet plotss
jm.acute.df = subset(jm.acute.df, control != "0.1% DMSO")
# color_scheme = c("#6BC2F4","#fa5cf2", "#ffa600", "#84d94e","#f25c78", "#b64cdf")
color_scheme = c("#f25c78","#f25c78",  "#84d94e","#ffa600", "#b64cdf")
q3 = ( ggplot(jm.acute.df, aes(x=concentration, y=mean.rawint, fill=concentration)) +
         theme_classic()  + labs(x = '', y = 'Mean RFP Intensity (A.U)') +
         geom_bar(stat="identity", color="black") +
         geom_errorbar(aes(ymin=mean.rawint - sem.rawint, ymax=mean.rawint + sem.rawint), width=.2,position=position_dodge(.9)) +
         facet_grid('. ~ type', scales="free", switch = "x", space='free') +
         # # geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
         scale_fill_manual(values=color_scheme) +
         theme(strip.placement = "outside") +
         theme(axis.title.x = element_blank(),
               axis.title.y =element_text(size=15),
               plot.title=element_text(hjust=0.5),
               axis.text.x = element_markdown(size=11),
               axis.text.y = element_text(size=11),
               legend.position="none",
               strip.text.x = element_text(size=13, face = "bold"),
               strip.background = element_blank())  +
         # geom_signif(comparisons = list(c("OP50 <i>E. coli</i><br>(Day 5)", "OP50 <i>E. coli</i><br>(Day 6)"),
         #                                c("OP50 <i>E. coli</i><br>(Day 6)", "0.1% <br>DMSO"),
         #                                c("OP50 <i>E. coli</i><br>(Day 6)", "200 uM EGCG")),
         #             test = "wilcox.test", step_increase = 0.075,
         #             map_signif_level = TRUE, tip_length = 0)
       
       scale_y_continuous(limits = c(0,60000))
       
       # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
       #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
       #       strip_text_x = element_text(size=12), strip_background = element_blank(),
       #       axis_line_y = element_line(colour = 'black', linetype='solid'))
       
       # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
       
)
q3



op50 <- tibble(
  type = "Control", 
  x = c(1, 1, 2, 2),
  y = c(48000, 52000, 52000, 46000)
) ## p = 0

h1 = q3 +
  geom_line(data = op50,
            aes(x= x, y = y), inherit.aes = F) +
  # annotate("text", x = 1.5, y = 55000, label = 'n.s', size = 5, color = "black") 
  geom_text(data = subset(jm.acute.df, type=="Control"),
            aes(x = 1.5, y = 55000, label = "n.s"), size=5, color = "black") +
  geom_text(data = subset(jm.acute.df, type=="Control"),
            aes(x = 3, y = 24000, label = "*"), size=7, color = "#960101") +
  geom_text(data = subset(jm.acute.df, type=="Field"),
          aes(x = 1, y = 6000, label = "**"), size=7, color = "#960101") +
  geom_text(data = subset(jm.acute.df, type=="Field"),
            aes(x = 2, y = 4500, label = "***"), size=7, color = "#960101") +
  geom_text(data = subset(jm.acute.df, type=="Lab"),
            aes(x = 1, y = 10000, label = "*"), size=7, color = "#960101") +
  geom_text(data = subset(jm.acute.df, type=="Lab"),
            aes(x = 2, y = 6000, label = "***"), size=7, color = "#960101") 

h1
ggsave(plot=h1, filename="Final-acute-jkm-barplot.png", width=8, height=4, dpi=600)
write.csv(jm.acute.df, "acute-jkm-high.csv")

# 
# ####### VISUALIZATION WITH STAT
# color_scheme = c("#6BC2F4","#fa5cf2", "#ffa600", "#84d94e","#f25c78", "#b64cdf")
# dat = subset(jm.acute.df, type == "Control")
# h1 = ( ggplot(dat, aes(x=concentration, y=mean.rawint, fill=concentration)) +
#          theme_classic()  + labs(x = '', y = 'Mean RFP Intensity (A.U)', title = "Control") +
#          geom_bar(stat="identity", color="black") +
#          geom_errorbar(aes(ymin=mean.rawint - sem.rawint, ymax=mean.rawint + sem.rawint), width=.2,position=position_dodge(.9)) +
#          # facet_grid('. ~ type', scales="free", switch = "x", space='free') +
#          # # geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
#          scale_fill_manual(values=c("black", "black", "#6BC2F4", "#247447")) +
#          theme(axis.title.x = element_blank(),
#                axis.title.y =element_text(size=15),
#                plot.margin = unit(c(0, 0, 0, 0), "cm"),
#                plot.title=element_text(hjust=0.5),
#                axis.text.x = element_markdown(size=8),
#                axis.text.y = element_text(size=11),
#                legend.position="none") +
#          scale_y_continuous(breaks = c(0, 10000, 20000, 30000, 40000, 50000),
#                             
#                             limits = c(0,55000))
#          
#          # scale_y_continuous(limits = c(0,55000))       # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
#        #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
#        #       strip_text_x = element_text(size=12), strip_background = element_blank(),
#        #       axis_line_y = element_line(colour = 'black', linetype='solid'))
#        
#        # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
#        
# )
# h1
# 
# op50 <- tibble(
#   x = c(1, 1, 2, 2),
#   y = c(48000, 52000, 52000, 46000)
# ) ## p = 0
# 
# h1 = h1 +
#   geom_line(data = op50, aes(x= x, y = y, group=1), inherit.aes = F) +
#   annotate("text", x = 1.5, y = 55000, label = 'n.s', size = 5, color = "black") +
#   annotate("text", x = 4, y = 24000, label = '*', size = 5, color = "#960101")
# h1
# 
# dat = subset(jm.acute.df, type == "Wine")
# h2 = ( ggplot(dat, aes(x=concentration, y=mean.rawint, fill=concentration)) +
#          theme_classic()  + labs(x = '', y = '', title = "Wine") +
#          geom_bar(stat="identity", color="black") +
#          geom_errorbar(aes(ymin=mean.rawint - sem.rawint, ymax=mean.rawint + sem.rawint), width=.2,position=position_dodge(.9)) +
#          # facet_grid('. ~ type', scales="free", switch = "x", space='free') +
#          # # geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
#          scale_fill_manual(values=c("#8C363A", "#206999")) +
#          theme(axis.title.x = element_blank(),
#                axis.title.y =element_text(size=15),
#                plot.title=element_text(hjust=0.5),
#                axis.text.x = element_markdown(size=8),
#                axis.text.y = element_text(size=11),
#                legend.position="none") +
#          scale_x_discrete(labels = c("Field","Lab")) +
#          scale_y_continuous(breaks = c(0, 10000, 20000, 30000, 40000, 50000),
#                             labels = c("", "", "", "", "", ""),
#                             limits = c(0,60000))+
#          theme(legend.position="none",
#                plot.title = element_text(hjust = 0.5),
#                plot.margin = unit(c(0, 0, 0, 0), "cm"),
#                axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))      # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
#        #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
#        #       strip_text_x = element_text(size=12), strip_background = element_blank(),
#        #       axis_line_y = element_line(colour = 'black', linetype='solid'))
#        
#        # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
#        
# )
# h2
# 
# h2 = h2 +
#   annotate("text", x = 1, y = 6000, label = '**', size = 5, color = "#960101") +
#   annotate("text", x = 2, y = 10000, label = '***', size = 5, color = "#960101")
# h2
# 
# dat = subset(jm.acute.df, type == "Lees")
# h3 = ( ggplot(dat, aes(x=concentration, y=mean.rawint, fill=concentration)) +
#          theme_classic()  + labs(x = '', y = '', title = "Lees") +
#          geom_bar(stat="identity", color="black") +
#          geom_errorbar(aes(ymin=mean.rawint - sem.rawint, ymax=mean.rawint + sem.rawint), width=.2,position=position_dodge(.9)) +
#          # facet_grid('. ~ type', scales="free", switch = "x", space='free') +
#          # # geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
#          scale_fill_manual(values=c("#D1A979", "#BC858A")) +
#          theme(axis.title.x = element_blank(),
#                axis.title.y =element_text(size=15),
#                plot.title=element_text(hjust=0.5),
#                axis.text.x = element_markdown(size=8),
#                axis.text.y = element_text(size=11),
#                legend.position="none") +
#          scale_x_discrete(labels = c("Field","Lab")) +
#          scale_y_continuous(breaks = c(0, 10000, 20000, 30000, 40000, 50000),
#                             labels = c("", "", "", "", "", ""),
#                             limits = c(0,60000))+
#          theme(legend.position="none",
#                plot.title = element_text(hjust = 0.5),
#                plot.margin = unit(c(0, 0, 0, 0), "cm"),
#                axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))      # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
#        #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
#        #       strip_text_x = element_text(size=12), strip_background = element_blank(),
#        #       axis_line_y = element_line(colour = 'black', linetype='solid'))
#        
#        # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
#        
# )
# h3
# 
# h3 = h3 +
#   annotate("text", x = 1, y = 4500, label = '**', size = 5, color = "#960101") +
#   annotate("text", x = 2, y = 6000, label = '***', size = 5, color = "#960101")
# h3
# 
# # ## Saving the figure for short term memory
# jm.acute.plot = plot_grid(h1, h2, h3,  align = "h", ncol = 3, rel_widths = c(0.50, 0.25, 0.25)) +
#   theme(plot.margin = unit(c(0.0,0,0,0), "cm")) 
# 
# png(filename = file.path(getwd(),"/BCH205-Acute-All.png"), width = 9, height = 3.5, units = "in", res = 600)
# plot(jm.acute.plot)
# dev.off()
print(paste("Saved at ", getwd()))
# plot(jm.acute.plot)


######################### Junelle & Mark CHRONIC #############################
jm.chronic.df = main.chronic.df
jm.chronic.df = subset(jm.chronic.df, (control %in% c('OP50','DMSO','EGCG', 'Itogon High', 'Field Lees','Lab Lees','Lab Wine')))

jm.chronic.df$control = factor(jm.chronic.df$control, levels=c( "OP50", "DMSO", "EGCG",
                                                           "Itogon High",
                                                           'Lab Wine',
                                                           'Field Lees','Lab Lees'),
                             
                             labels=c("OP50 <i>E. coli</i><br>(Day 6)", "0.1% DMSO", "200uM EGCG",
                                      "Field Wine",
                                      "Lab Wine",
                                      "Field Lees",
                                      "Lab Lees"))
jm.chronic.df <- jm.chronic.df %>% arrange(control)

color_scheme = c("#84d94f", "#ffa600", "#f25c78",  "#b64cdf", "#fa5cf2")

jm.chronic.df <- jm.chronic.df %>%
  mutate(type = case_when(
    grepl("E. coli|DMSO|EGCG", control) ~ "Control",
    grepl("Lab", control) ~ "Lab",
    grepl("Field", control) ~ "Field"
  ),
  concentration = case_when(
    grepl("Day 6", control) ~ "OP50 <i>E. coli</i><br>(Day 6)",
    grepl("DMSO", control) ~ "0.1% DMSO",
    grepl("EGCG", control) ~ "200 uM EGCG",
    grepl("Wine", control) ~ "Wine",
    grepl("Lees", control) ~ "Lees"
  ))

jm.chronic.df$type = factor(jm.chronic.df$type, levels=c("Control", "Field", "Lab"),
                          labels=c("Control", "Field", "Lab"))

jm.chronic.df$concentration = factor(jm.chronic.df$concentration, levels=c("OP50 <i>E. coli</i><br>(Day 6)",
                                                                       "0.1% DMSO", "200 uM EGCG",
                                                                       "Wine", "Lees"),
                                   labels=c("OP50 <i>E. coli</i><br>(Day 6)",
                                            "0.1% <br>DMSO", "200uM <br>EGCG",
                                            "Wine", "Lees"))

jm.chronic.df = subset(jm.chronic.df, control != "0.1% DMSO")
color_scheme = c("#f25c78", "#84d94e","#ffa600", "#b64cdf")
q4 = ( ggplot(jm.chronic.df, aes(x=concentration, y=mean.rawint, fill=concentration)) +
         theme_classic()  + labs(x = '', y = 'Mean RFP Intensity (A.U)') +
         geom_bar(stat="identity", color="black") +
         geom_errorbar(aes(ymin=mean.rawint - sem.rawint, ymax=mean.rawint + sem.rawint), width=.2,position=position_dodge(.9)) +
         facet_grid('. ~ type', scales="free", switch = "x", space='free') +
         # # geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
         scale_fill_manual(values=color_scheme) +
         theme(strip.placement = "outside") +
         theme(axis.title.x = element_blank(),
               axis.title.y =element_text(size=15),
               plot.title=element_text(hjust=0.5),
               axis.text.x = element_markdown(size=11),
               axis.text.y = element_text(size=11),
               legend.position="none",
               strip.text.x = element_text(size=13, face = "bold"),
               strip.background = element_blank())  +
         # geom_signif(comparisons = list(c("OP50 <i>E. coli</i><br>(Day 5)", "OP50 <i>E. coli</i><br>(Day 6)"),
         #                                c("OP50 <i>E. coli</i><br>(Day 6)", "0.1% <br>DMSO"),
         #                                c("OP50 <i>E. coli</i><br>(Day 6)", "200 uM EGCG")),
         #             test = "wilcox.test", step_increase = 0.075,
         #             map_signif_level = TRUE, tip_length = 0)
         
         scale_y_continuous(limits = c(0,60000))
       
       # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
       #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
       #       strip_text_x = element_text(size=12), strip_background = element_blank(),
       #       axis_line_y = element_line(colour = 'black', linetype='solid'))
       
       # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
       
)
q4


ggsave(plot=q4, filename="Final-chronic-jkm-barplot.png", width=8, height=4, dpi=600)
write.csv(jm.chronic.df, "chronic-jkm-high.csv")



# 
# 
# 
# 
# 
# color_scheme = c("#fa5cf2", "#ffa600", "#84d94e","#f25c78", "#b64cdf")
# q4 = ( ggplot(jm.chronic.df, aes(x=concentration, y=mean.rawint, fill=concentration)) +
#          theme_classic()  + labs(x = '', y = 'Mean RFP Intensity (A.U)') +
#          geom_bar(stat="identity", color="black") +
#          geom_errorbar(aes(ymin=mean.rawint - sem.rawint, ymax=mean.rawint + sem.rawint), width=.2,position=position_dodge(.9)) +
#          facet_grid('. ~ type', scales="free", switch = "x", space='free') +
#          # # geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
#          scale_fill_manual(values=color_scheme) +
#          theme(strip.placement = "outside") +
#          theme(axis.title.x = element_blank(),
#                axis.title.y =element_text(size=15),
#                plot.title=element_text(hjust=0.5),
#                axis.text.x = element_markdown(size=8),
#                axis.text.y = element_text(size=11),
#                legend.position="none",
#                strip.text.x = element_text(size=13, face = "bold"),
#                strip.background = element_blank())
#        
#        # scale_y_continuous(limits = c(0,100))
#        
#        # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
#        #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
#        #       strip_text_x = element_text(size=12), strip_background = element_blank(),
#        #       axis_line_y = element_line(colour = 'black', linetype='solid'))
#        
#        # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
#        
# )
# q4
# 
# ## Statistics
# ## Statistics
# dat.test = subset(chronic.csv, (control %in% c("OP50","DMSO", "EGCG",
#                                              "Itogon High","Lab Wine",
#                                              "Field Lees",  "Lab Lees")))
# 
# dat.test$control = factor(dat.test$control, levels=c("OP50","DMSO", "EGCG",
#                                                      "Itogon High","Lab Wine",
#                                                      "Field Lees",  "Lab Lees"),
#                           
#                           labels=c("OP50","DMSO", "EGCG",
#                                    "Field Wine","Lab Wine",
#                                    "Field Lees",  "Lab Lees"))
# dat.test = arrange(dat.test, control)
# res2 = DunnettTest(x=dat.test$RawIntDen, g=dat.test$control)
# write.csv(as.data.frame(res2$`OP50`), file="BCH205-dunnet-stat-results-chronic.csv")
# write.csv(jm.chronic.df, "BCH205-chronic-values.csv")
# 
# ####### VISUALIZATION WITH STAT
# color_scheme = c("#6BC2F4","#fa5cf2", "#ffa600", "#84d94e","#f25c78", "#b64cdf")
# dat = subset(jm.chronic.df, type == "Control")
# h1 = ( ggplot(dat, aes(x=concentration, y=mean.rawint, fill=concentration)) +
#          theme_classic()  + labs(x = '', y = 'Mean RFP Intensity (A.U)', title = "Control") +
#          geom_bar(stat="identity", color="black") +
#          geom_errorbar(aes(ymin=mean.rawint - sem.rawint, ymax=mean.rawint + sem.rawint), width=.2,position=position_dodge(.9)) +
#          # facet_grid('. ~ type', scales="free", switch = "x", space='free') +
#          # # geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
#          scale_fill_manual(values=c("black", "#6BC2F4", "#247447")) +
#          theme(axis.title.x = element_blank(),
#                axis.title.y =element_text(size=15),
#                plot.margin = unit(c(0, 0, 0, 0), "cm"),
#                plot.title=element_text(hjust=0.5),
#                axis.text.x = element_markdown(size=8),
#                axis.text.y = element_text(size=11),
#                legend.position="none") +
#          scale_y_continuous(breaks = c(0, 10000, 20000, 30000),
#                             
#                             limits = c(0,30000))
#        
#        # scale_y_continuous(limits = c(0,55000))       # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
#        #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
#        #       strip_text_x = element_text(size=12), strip_background = element_blank(),
#        #       axis_line_y = element_line(colour = 'black', linetype='solid'))
#        
#        # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
#        
# )
# h1
# 
# 
# 
# dat = subset(jm.chronic.df, type == "Wine")
# h2 = ( ggplot(dat, aes(x=concentration, y=mean.rawint, fill=concentration)) +
#          theme_classic()  + labs(x = '', y = '', title = "Wine") +
#          geom_bar(stat="identity", color="black") +
#          geom_errorbar(aes(ymin=mean.rawint - sem.rawint, ymax=mean.rawint + sem.rawint), width=.2,position=position_dodge(.9)) +
#          # facet_grid('. ~ type', scales="free", switch = "x", space='free') +
#          # # geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
#          scale_fill_manual(values=c("#8C363A", "#206999")) +
#          theme(axis.title.x = element_blank(),
#                axis.title.y =element_text(size=15),
#                plot.title=element_text(hjust=0.5),
#                axis.text.x = element_markdown(size=8),
#                axis.text.y = element_text(size=11),
#                legend.position="none") +
#          scale_x_discrete(labels = c("Field","Lab")) +
#          scale_y_continuous(breaks = c(0, 10000, 20000, 30000),
#                             labels = c("", "", "", ""),
#                             limits = c(0,30000))+
#          theme(legend.position="none",
#                plot.title = element_text(hjust = 0.5),
#                plot.margin = unit(c(0, 0, 0, 0), "cm"),
#                axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))      # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
#        #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
#        #       strip_text_x = element_text(size=12), strip_background = element_blank(),
#        #       axis_line_y = element_line(colour = 'black', linetype='solid'))
#        
#        # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
#        
# )
# h2
# 
# 
# dat = subset(jm.chronic.df, type == "Lees")
# h3 = ( ggplot(dat, aes(x=concentration, y=mean.rawint, fill=concentration)) +
#          theme_classic()  + labs(x = '', y = '', title = "Lees") +
#          geom_bar(stat="identity", color="black") +
#          geom_errorbar(aes(ymin=mean.rawint - sem.rawint, ymax=mean.rawint + sem.rawint), width=.2,position=position_dodge(.9)) +
#          # facet_grid('. ~ type', scales="free", switch = "x", space='free') +
#          # # geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
#          scale_fill_manual(values=c("#D1A979", "#BC858A")) +
#          theme(axis.title.x = element_blank(),
#                axis.title.y =element_text(size=15),
#                plot.title=element_text(hjust=0.5),
#                axis.text.x = element_markdown(size=8),
#                axis.text.y = element_text(size=11),
#                legend.position="none") +
#          scale_x_discrete(labels = c("Field","Lab")) +
#          scale_y_continuous(breaks = c(0, 10000, 20000, 30000),
#                             labels = c("", "", "", ""),
#                             limits = c(0,30000))+
#          theme(legend.position="none",
#                plot.title = element_text(hjust = 0.5),
#                plot.margin = unit(c(0, 0, 0, 0), "cm"),
#                axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))      # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
#        #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
#        #       strip_text_x = element_text(size=12), strip_background = element_blank(),
#        #       axis_line_y = element_line(colour = 'black', linetype='solid'))
#        
#        # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
#        
# )
# h3
# 
# 
# # ## Saving the figure for short term memory
# jm.chronic.plot = plot_grid(h1, h2, h3,  align = "h", ncol = 3, rel_widths = c(0.50, 0.25, 0.25)) +
#   theme(plot.margin = unit(c(0.0,0,0,0), "cm")) 
# 
# png(filename = file.path(getwd(),"/BCH205-Chronic-All.png"), width = 9, height = 3.5, units = "in", res = 600)
# plot(jm.chronic.plot)
# dev.off()
print(paste("Saved at ", getwd()))
# plot(jm.chronic.plot)





