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

## Set working directory to the file location using Rstudio API
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

######################### READ FILES #############################
df.csv = read.csv("DPPH 205 data.csv", header = T)
df.csv$Control = factor(df.csv$Control, 
                                 levels = c("Ascorbic Acid", "Field Wine", "Lab Wine","Field Lees", "Lab Lees"))

df.csv.summary =  df.csv %>% 
  group_by(Control) %>%
  summarise(Mean= mean(DPPH_value),
            Std = std.error(DPPH_value),
            count = n()) %>%
  as.data.frame()
write.csv(df.csv.summary, "DPPH-values.csv")
color_scheme = c("#84d94f", "#ffa600", "#f25c78",  "#b64cdf", "#fa5cf2")

## Itogon
# ascorbic.acid = subset(df.csv.summary, Treatment == "Ascorbic acid")
q1 = ( ggplot(df.csv.summary, aes(x=Control, y=Mean, fill=Control)) +
         theme_classic()  + labs(x = '', y = '% DPPH') +
         geom_bar(stat="identity", color="black") +
         geom_errorbar(aes(ymin=Mean-Std, ymax=Mean+Std), width=.2,position=position_dodge(.9)) +
         # facet_grid('. ~ Sample', scales="free", switch = "x", space='free') +
         # geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
         scale_fill_manual(values=color_scheme) +
         theme(strip.placement = "outside") +
         theme(axis.title.x = element_blank(),
               axis.title.y =element_text(size=15),
               plot.title=element_text(hjust=0.5),
               axis.text.x = element_text(size=13, face="bold"),
               axis.text.y = element_text(size=11),
               legend.position="none",
               strip.text.x = element_text(size=13, face = "bold"),
               strip.background = element_blank()) + 

         scale_y_continuous(limits = c(0,100))
       
       # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
       #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
       #       strip_text_x = element_text(size=12), strip_background = element_blank(),
       #       axis_line_y = element_line(colour = 'black', linetype='solid'))
       
       # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
       
)
q1

## Statistics
# df.stat.summary = read.csv("DPPH Stat Data.csv", header = T)
# df.stat.summary$Control = factor(df.stat.summary$Control, 
#                                  levels = c("Ascorbic Acid", "Itogon", "Kapangan","La Trinidad", "Sablan"))
res = DunnettTest(x=df.csv$DPPH_value, g=df.csv$Control)
write.csv(res$`Ascorbic Acid`, "DPPH-stat-results.csv")
res


q1 = q1 +
  annotate("text", x = 2, y = 25, label = "***", size = 7, color = "#960101") +
  annotate("text", x = 3, y = 51, label = "***", size = 7, color = "#960101") +
  annotate("text", x = 4, y = 43.5, label = "***", size = 7, color = "#960101") +
  annotate("text", x = 5, y = 58, label = "***", size = 7, color = "#960101")
q1

png(filename = file.path(getwd(),"/DPPH-All.png"), width = 8, height = 3.5, units = "in", res = 600)
plot(q1)
dev.off()
print(paste("Saved at ", getwd()))


######################### READ FILES #############################
df.csv = read.csv("FRAP 205 data.csv", header = T)
df.csv$Control = factor(df.csv$Control, 
                        levels = c("Ascorbic Acid", "Field Wine", "Lab Wine","Field Lees", "Lab Lees"))

df.csv.summary =  df.csv %>% 
  group_by(Control) %>%
  summarise(Mean= mean(FRAP_reduce),
            Std = std.error(FRAP_reduce),
            count = n()) %>%
  as.data.frame()
write.csv(df.csv.summary, "FRAP-values.csv")
color_scheme = c("#84d94f", "#ffa600", "#f25c78",  "#b64cdf", "#fa5cf2")



## Itogon
# ascorbic.acid = subset(df.csv.summary, Treatment == "Ascorbic acid")
q2 = ( ggplot(df.csv.summary, aes(x=Control, y=Mean, fill=Control)) +
         theme_classic()  + labs(x = '', y = '% Iron reduced') +
         geom_bar(stat="identity", color="black") +
         geom_errorbar(aes(ymin=Mean-Std, ymax=Mean+Std), width=.2,position=position_dodge(.9)) +
         # facet_grid('. ~ Sample', scales="free", switch = "x", space='free') +
         # geom_hline(aes(yintercept = 0.9, color = "U"), linetype='dashed', color='#960101') +
         scale_fill_manual(values=color_scheme) +
         theme(strip.placement = "outside") +
         theme(axis.title.x = element_blank(),
               axis.title.y =element_text(size=15),
               plot.title=element_text(hjust=0.5),
               axis.text.x = element_text(size=13, face="bold"),
               axis.text.y = element_text(size=11),
               legend.position="none",
               strip.text.x = element_text(size=13, face = "bold"),
               strip.background = element_blank()) + 
         
         scale_y_continuous( limits = c(0,1))
       
       # theme(figure_size=c(16,4), legend_position = "bottom", legend_title = element_blank(), legend_entry_spacing_x = 20,
       #       axis_title_y = element_text(size = 14), axis_title_x = element_text(size = 14),
       #       strip_text_x = element_text(size=12), strip_background = element_blank(),
       #       axis_line_y = element_line(colour = 'black', linetype='solid'))
       
       # scale_fill_manual(name = 'Toxin', values=gradient_scheme, labels=['$t = 0$h', '$t = 3$h', '$t = 6$h', '$t = 9$h'])
       
)
q2

## Statistics
# df.stat.summary = read.csv("FRAP Stat Data.csv", header = T)
# df.stat.summary$Control = factor(df.stat.summary$Control, 
#                                  levels = c("Ascorbic Acid", "Itogon", "Kapangan","La Trinidad", "Sablan"))
res = DunnettTest(x=df.csv$FRAP_reduce, g=df.csv$Control)
write.csv(res$`Ascorbic Acid`, "FRAP-stat-results.csv")
res

q2 = q2 + 
  annotate("text", x = 2, y = 68.2, label = "***", size = 7, color = "#960101") +
  annotate("text", x = 3, y = 81.5, label = "***", size = 7, color = "#960101") +
  annotate("text", x = 4, y = 81.5, label = "***", size = 7, color = "#960101") +
  annotate("text", x = 5, y = 85.5, label = "***", size = 7, color = "#960101")
q2



png(filename = file.path(getwd(),"/FRAP-All.png"), width = 8, height = 3.5, units = "in", res = 600)

plot(q2)
dev.off()
print(paste("Saved at ", getwd()))