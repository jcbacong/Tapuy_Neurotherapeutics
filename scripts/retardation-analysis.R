
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
library(grid)


## Set working directory to the file location using Rstudio API
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


## Prepare group name
rename_wine_or_lees = function(x) {
  if( toupper(substr(x,start=1,stop=2)) == "FW"){
    return("Field Wine")
  } else if ( toupper(substr(x,start=1,stop=2)) == "FL") {
    return("Field Lees")
  } else if( toupper(substr(x,start=1,stop=2)) == "LW") {
    return("Lab Wine")
  } else if( toupper(substr(x,start=1,stop=2)) == "LL") {
    return("Lab Lees")
  } else {
    return("Control")
  }
}

rename_concentration = function(x) {
  if( toupper(substr(x,start=3,stop=3)) == "L"){
    return("Low")
  } else if ( toupper(substr(x,start=3,stop=3)) == "M") {
    return("Mid")
  } else if( toupper(substr(x,start=3,stop=3)) == "H") {
    return("High")
  } else {
    return("Control")
  }
}

## Prepare group name
rename_fill = function(x) {
  if( toupper(substr(x,start=1,stop=2)) == "FW"){
    return("Wine")
  } else if ( toupper(substr(x,start=1,stop=2)) == "FL") {
    return("Lees")
  } else if( toupper(substr(x,start=1,stop=2)) == "LW") {
    return("Wine")
  } else if( toupper(substr(x,start=1,stop=2)) == "LL") {
    return("Lees")
  } 
  else if( x == "Food") {
    return("OP50")
  } else if( x == "EGCG") {
    return("EGCG")
  }
  # } else if( x == "EGCG") {
  #   return("200uM EGCG")
  # }
  else {
    return("Control")
  }
}


rename_lab_or_field = function(x) {
  if(substr(x,start=1,stop=2) == "FW"){
    return("Field")
  } else if (substr(x,start=1,stop=2) == "FL") {
    return("Field")
  } else if(substr(x,start=1,stop=2) == "LW") {
    return("Lab")
  } else if(substr(x,start=1,stop=2) == "LL") {
    return("Lab")
  } else {
    return("Control")
  }
}


strain_comparison_t_test_2 = function(strain.summary.norm.9h, alternative = "greater") {
  strain.summary.norm.9h$sem.norm = strain.summary.norm.9h$mean.rawint * strain.summary.norm.9h$rel.error
  strain.summary.norm.9h$p.value = NaN
  tox.mean = strain.summary.norm.9h[strain.summary.norm.9h$control=="Food", mean.rawint]
  tox.sem = strain.summary.norm.9h[strain.summary.norm.9h$control=="Food", sem.norm]
  tox.obs = 15
  
  for(ctrl in strain.summary.norm.9h$control) {
    if(ctrl!="Food") {
      var.mean = strain.summary.norm.9h[strain.summary.norm.9h$control==ctrl, mean.rawint]
      var.sem = strain.summary.norm.9h[strain.summary.norm.9h$control==ctrl, sem.norm]
      var.obs = 15
      res = tsum.test(var.mean, var.sem, var.obs, tox.mean, tox.sem, tox.obs, alternative=alternative)
      strain.summary.norm.9h[strain.summary.norm.9h$control == ctrl, "p.value"] = res$p.value
    }
  }
  return(strain.summary.norm.9h)
}


######################### READ FILES #############################
### UA57
strain.ua57 = "UA57"
controls.ua57.csv = read.csv(file.path(getwd(),"controls.csv"), header = T)
treatment.ua57.csv = read.csv(file.path(getwd(),"treatment.csv"), header = T)
treatment.ua57.csv$control = lapply(treatment.ua57.csv$control, toupper)
ua57.df = rbind(controls.ua57.csv, treatment.ua57.csv)

ua57.df$concentration = lapply(ua57.df$control, rename_concentration)
ua57.df$group = lapply(ua57.df$control, rename_wine_or_lees)
ua57.df[ua57.df == "OP50"] = "Food"

## Factor categories
ua57.df$control <- factor(ua57.df$control, levels = c('Food', 'DMSO', 'EGCG', 'Metab',
                                                                'FWL', 'FWM', 'FWH',
                                                                'LWL', 'LWM', 'LWH',
                                                                'FLL', 'FLM', 'FLH',
                                                                'LLL', 'LLM', 'LLH'))

ua57.df$concentration <- factor(ua57.df$concentration, levels = c('Control','Low', 'Mid','High'))

ua57.df$group <- factor(ua57.df$group, levels = c('Control','Field Wine', 'Lab Wine',
                                                        'Field Lees', 'Lab Lees'))
ua57.df$neuron = 'ua57'

## Remove metab
ua57.df = subset(ua57.df, control != "Metab")
ua57.df


## Divide 45 total number of worms per concentration 
ua57.summary <- ua57.df %>% 
                group_by(control, time, group, neuron, concentration) %>%
                summarise(mean.rawint = mean(RawIntDen),
                          sem.rawint = std.error(RawIntDen),
                          mean.area = mean(Area),
                          sem.area = std.error(Area),
                          obs.rawint = n()/15) %>% 
                as.data.frame()

ua57.summary

##########################################################################################################
################################## COMPARISON OF NUMBER NEURONS  #########################################
##########################################################################################################

neuron.summary <- ua57.df %>% 
  group_by(control, time, group, neuron, concentration, replicate) %>%
  summarise(mean.area = mean(Area),
            obs.rawint = n()/15) %>% 
  as.data.frame()
neuron.summary

neuron.summary2 <- neuron.summary %>% 
  group_by(control, time, group, neuron, concentration) %>%
  summarise(mean.dops = mean(obs.rawint),
            sem.dops= std.error(obs.rawint)) %>% 
  as.data.frame()
neuron.summary2$rel.error = neuron.summary2$sem.dops/neuron.summary2$mean.dops
neuron.summary2$mean.rawint = neuron.summary2$mean.dops/8.0


neuron.summary2 = neuron.summary2[neuron.summary2$time=="t3",]
setDT(neuron.summary2)
set(neuron.summary2, j = "mean.rawint", value = as.numeric(neuron.summary2[["mean.rawint"]]))
set(neuron.summary2, j = "rel.error", value = as.numeric(neuron.summary2[["rel.error"]]))

neuron.summary2 = subset(neuron.summary2, control %in% c("Food", "EGCG", "FWH", "LWH", "FLH", "LLH"))

neuron.summary2$name = lapply(neuron.summary2$control, rename_fill)
neuron.summary2$name = factor(neuron.summary2$name, 
                      levels = c("OP50", "EGCG","Wine", "Lees"),
                      labels = c("OP50 <i>E. coli</i>", "200uM EGCG","Wine", "Lees"))

neuron.summary2$type = lapply(neuron.summary2$control, rename_lab_or_field)
neuron.summary2$type = factor(neuron.summary2$type, 
                      levels = c("Control","Field", "Lab"),
                      labels = c("Control","Field", "Lab"))
neuron.summary2

color_scheme = c("#f25c78",  "#84d94e","#ffa600", "#b64cdf")
p = ggplot(neuron.summary2, aes(name, mean.rawint, fill=name)) +
  theme_classic() + labs(x = '', y = '% Dopaminergic neurons') + 
  theme(axis.title.y = element_text(size = 14), axis.text.x=element_text(size=13), axis.text.y=element_text(size=12)) +
  geom_bar(stat="identity", color="black",position=position_dodge()) +
  geom_errorbar(aes(ymin=mean.rawint*(1-rel.error), ymax=mean.rawint*(1+rel.error)), width=.2,position=position_dodge(.9)) +
  facet_grid('. ~ type', scales="free", switch = "x", space='free') +
  scale_fill_manual(values=color_scheme) +
  theme(legend.position="none") + 
  theme(strip.placement = "outside") +
  theme(axis.title.x = element_blank(),
        axis.title.y =element_text(size=16),
        plot.title=element_text(hjust=0.5),
        axis.text.x = element_markdown(size=10),
        axis.text.y = element_text(size=11),
        legend.position="none",
        strip.text.x = element_text(size=13, face = "bold"),
        strip.background = element_blank()) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1))

res = strain_comparison_t_test_2(neuron.summary2)
res

q = p +
  geom_text(data = subset(neuron.summary2, type=="Control"),
            aes(x = 2, y = 0.74, label = "***"), size=7, color = "#960101") +
  geom_text(data = subset(neuron.summary2, type=="Field"),
            aes(x = 1, y = 0.66, label = "***"), size=7, color = "#960101") +
  geom_text(data = subset(neuron.summary2, type=="Field"),
            aes(x = 2, y = 0.675, label = "***"), size=7, color = "#960101") +
  geom_text(data = subset(neuron.summary2, type=="Lab"),
            aes(x = 1, y = 0.61, label = "***"), size=7, color = "#960101") +
  geom_text(data = subset(neuron.summary2, type=="Lab"),
            aes(x = 2, y = 0.56, label = "***"), size=7, color = "#960101")
q

ggsave(plot=p, filename="Final-neuronal-loss-barplot.png", width=8, height=4, dpi=600) 

print(paste("Saved at ", getwd()))

neuron.summary2$error = neuron.summary2$sem.rawint * neuron.summary2$rel.error
write.csv(select(neuron.summary2, c("control", "group", "mean.rawint", "error")), "retardation-neurons.csv")



print(paste("Saved at ", getwd()))

