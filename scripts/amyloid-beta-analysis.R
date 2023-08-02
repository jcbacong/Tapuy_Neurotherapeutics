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


############################################################# 
#########################  ACUTE ############################
############################################################# 
jm.acute.df = main.acute.df
jm.acute.df = subset(jm.acute.df, (control %in% c('OP50','DMSO','EGCG', 'Itogon High', 'Field Lees','Lab Lees','Lab Wine')))

jm.acute.df$control = factor(jm.acute.df$control, levels=c("OP50", "DMSO", "EGCG",
                                                           "Itogon High",
                                                           'Lab Wine',
                                                           'Field Lees','Lab Lees'),
                             
                             labels=c("Untreated", "0.1% DMSO", "200 uM EGCG",
                                      "Field Wine",
                                      "Lab Wine",
                                      "Field Lees",
                                      "Lab Lees"))
jm.acute.df <- jm.acute.df %>% arrange(control)

color_scheme = c("#84d94f", "#ffa600", "#f25c78",  "#b64cdf", "#fa5cf2")

jm.acute.df <- jm.acute.df %>%
  mutate(type = case_when(
    grepl("Untreated|DMSO|EGCG", control) ~ "Control",
    grepl("Lab", control) ~ "Lab",
    grepl("Field", control) ~ "Traditional"
  ),
  concentration = case_when(
    grepl("Untreated", control) ~ control,
    grepl("DMSO", control) ~ "0.1% DMSO",
    grepl("EGCG", control) ~ "200 uM EGCG",
    grepl("Wine", control) ~ "Wine",
    grepl("Lees", control) ~ "Lees"
  ))

jm.acute.df$type = factor(jm.acute.df$type, levels=c("Control", "Traditional", "Lab"),
                          labels=c("Control", "Traditional", "Lab"))

jm.acute.df$concentration = factor(jm.acute.df$concentration, levels=c("Untreated",
                                                                       "0.1% DMSO", "200 uM EGCG",
                                                                       "Wine", "Lees"),
                                   labels=c("Untreated",
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
color_scheme = c("#f25c78","#f25c78",  "#84d94e","#ffa600", "#b64cdf")
q3 = ( ggplot(jm.acute.df, aes(x=concentration, y=mean.rawint, fill=concentration)) +
         theme_classic()  + labs(x = '', y = 'Mean RFP Intensity (A.U)') +
         geom_bar(stat="identity", color="black") +
         geom_errorbar(aes(ymin=mean.rawint - sem.rawint, ymax=mean.rawint + sem.rawint), width=.2,position=position_dodge(.9)) +
         facet_grid('. ~ type', scales="free", switch = "x", space='free') +
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
         scale_y_continuous(limits = c(0,50000))
)
q3


h1 = q3 +
  geom_text(data = subset(jm.acute.df, type=="Control"),
            aes(x = 2, y = 24000, label = "*"), size=7, color = "#960101") +
  geom_text(data = subset(jm.acute.df, type=="Traditional"),
            aes(x = 1, y = 6000, label = "**"), size=7, color = "#960101") +
  geom_text(data = subset(jm.acute.df, type=="Traditional"),
            aes(x = 2, y = 4500, label = "***"), size=7, color = "#960101") +
  geom_text(data = subset(jm.acute.df, type=="Lab"),
            aes(x = 1, y = 10000, label = "*"), size=7, color = "#960101") +
  geom_text(data = subset(jm.acute.df, type=="Lab"),
            aes(x = 2, y = 6000, label = "***"), size=7, color = "#960101") 

h1
ggsave(plot=h1, filename="Final-acute-jkm-barplot.png", width=6.5, height=4, dpi=600)
write.csv(jm.acute.df, "acute-jkm-high.csv")
print(paste("Saved at ", getwd()))



######################### Junelle & Mark CHRONIC #############################
jm.chronic.df = main.chronic.df
jm.chronic.df = subset(jm.chronic.df, (control %in% c('OP50','DMSO','EGCG', 'Itogon High', 'Field Lees','Lab Lees','Lab Wine')))

jm.chronic.df$control = factor(jm.chronic.df$control, levels=c( "OP50", "DMSO", "EGCG",
                                                                "Itogon High",
                                                                'Lab Wine',
                                                                'Field Lees','Lab Lees'),
                               
                               labels=c("Untreated", "0.1% DMSO", "200uM EGCG",
                                        "Field Wine",
                                        "Lab Wine",
                                        "Field Lees",
                                        "Lab Lees"))
jm.chronic.df <- jm.chronic.df %>% arrange(control)

jm.chronic.df <- jm.chronic.df %>%
  mutate(type = case_when(
    grepl("Untreated|DMSO|EGCG", control) ~ "Control",
    grepl("Lab", control) ~ "Lab",
    grepl("Field", control) ~ "Traditional"
  ),
  concentration = case_when(
    grepl("Untreated", control) ~ "Untreated",
    grepl("DMSO", control) ~ "0.1% DMSO",
    grepl("EGCG", control) ~ "200 uM EGCG",
    grepl("Wine", control) ~ "Wine",
    grepl("Lees", control) ~ "Lees"
  ))

jm.chronic.df$type = factor(jm.chronic.df$type, levels=c("Control", "Traditional", "Lab"),
                            labels=c("Control", "Traditional", "Lab"))

jm.chronic.df$concentration = factor(jm.chronic.df$concentration, levels=c("Untreated",
                                                                           "0.1% DMSO", "200 uM EGCG",
                                                                           "Wine", "Lees"),
                                     labels=c("Untreated",
                                              "0.1% <br>DMSO", "200uM <br>EGCG",
                                              "Wine", "Lees"))

jm.chronic.df = subset(jm.chronic.df, !control %in% c("0.1% DMSO", "200uM EGCG"))
color_scheme = c("#f25c78", "#ffa600", "#b64cdf")
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
         scale_y_continuous(limits = c(0,30000))
)
q4

ggsave(plot=q4, filename="Final-chronic-jkm-barplot.png", width=6.5, height=4, dpi=600)
write.csv(jm.chronic.df, "chronic-jkm-high.csv")
print(paste("Saved at ", getwd()))
