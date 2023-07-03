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
library("data.table")
library("BSDA")
library(ggtext)

## Set working directory to the file location using Rstudio API
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

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
  else if( x == "CTRL") {
    return("OP50")
  } else if( x == "0.01% EGCG") {
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

######################################  GENTLE HEAD TOUCH  ######################################
## Import file name
## Change the title
df.csv = read.csv(file.path(getwd(),"Gentle Head Touch.csv"), header = T)
df.csv$Treatment <- factor(df.csv$Treatment, levels = unique(df.csv$Treatment))

## Summarize the data
df.csv.summary <- df.csv %>% 
  group_by(Treatment) %>%
  summarise(mean.val = mean(Mean)/100,
            sem.val = std.error(Mean)/100,
            obs.val = n()) %>%
  as.data.frame()



## Using Dunnett's Tests (Apply statistics)
## Run first & determine the significant value
## Only those with stars
comp = lm( Mean ~ Treatment, data = df.csv)
summary(aov(comp))
res = DunnettTest(x=df.csv$Mean, g=df.csv$Treatment)


### Bar plots
df.bar = subset(df.csv.summary, Treatment %in% c("CTRL", "0.01% EGCG", "FWH", 
                                                 "FLH", "LWH", "LLH"))

df.bar$name = lapply(df.bar$Treatment, rename_fill)
df.bar$name = factor(df.bar$name, 
                      levels = c("OP50", "EGCG","Wine", "Lees"),
                      labels = c("OP50<br><i>E. coli</i>", "200uM<br>EGCG","Wine", "Lees"))

df.bar$type = lapply(df.bar$Treatment, rename_lab_or_field)
df.bar$type = factor(df.bar$type, 
                      levels = c("Control","Field", "Lab"),
                      labels = c("Control","Field", "Lab"))
df.bar

color_scheme = c("#f25c78",  "#84d94e","#ffa600", "#b64cdf")
p = ggplot(df.bar, aes(name,mean.val,fill=name)) +
  theme_classic() + labs(x = '', y = 'Responsiveness') + 
  theme(axis.title.y = element_text(size = 14), axis.text.x=element_text(size=13), axis.text.y=element_text(size=12)) +
  geom_bar(stat="identity", color="black",position=position_dodge()) +
  geom_errorbar(aes(ymin=mean.val - sem.val, ymax=mean.val + sem.val), width=.2,position=position_dodge(.9)) +
  facet_grid('. ~ type', scales="free", switch = "x", space='free') +
  scale_fill_manual(values=color_scheme) +
  theme(legend.position="none") + 
  theme(strip.placement = "outside") +
  theme(axis.title.x = element_blank(),
        axis.title.y =element_text(size=16),
        plot.title=element_text(hjust=0.5),
        axis.text.x = element_markdown(size=12),
        axis.text.y = element_text(size=11),
        legend.position="none",
        strip.text.x = element_text(size=13, face = "bold"),
        strip.background = element_blank()) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1.02))

p

q = p +
  geom_text(data = subset(df.bar, type=="Control"),
            aes(x = 2, y = 1.0, label = "*"), size=7, color = "#960101") +
  geom_text(data = subset(df.bar, type=="Field"),
            aes(x = 2, y = 1.0, label = "*"), size=7, color = "#960101") +
  geom_text(data = subset(df.bar, type=="Lab"),
            aes(x = 2, y = 1.02, label = "*"), size=7, color = "#960101")
q

### Saving the plot
png(filename = file.path(getwd(),"/Gentle-Head-Touch.png"), width = 5, height = 3, units = "in", res = 600)
plot(q)
dev.off()

write.csv(as.data.frame(res$`CTRL`), file="Gentle-head-touch-stat-results.csv")
write.csv(df.csv.summary, file="Gentle-head-touch-values.csv")

print(paste("Saved at ", getwd()))

####################################################################################
######################################  GENTLE TAIL TOUCH  ######################################
## Import file name
## Change the title
df.csv = read.csv(file.path(getwd(),"Gentle Tail Touch.csv"), header = T)
df.csv$Treatment <- factor(df.csv$Treatment, levels = unique(df.csv$Treatment))

## Summarize the data
df.csv.summary <- df.csv %>% 
  group_by(Treatment) %>%
  summarise(mean.val = mean(Mean)/100,
            sem.val = std.error(Mean)/100,
            obs.val = n()) %>%
  as.data.frame()

## Using Dunnett's Tests (Apply statistics)
## Run first & determine the significant value
## Only those with stars
comp = lm( Mean ~ Treatment, data = df.csv)
summary(aov(comp))
res = DunnettTest(x=df.csv$Mean, g=df.csv$Treatment)
res

### Bar plots
df.bar = subset(df.csv.summary, Treatment %in% c("CTRL", "0.01% EGCG", "FWH", 
                                                 "FLH", "LWH", "LLH"))

df.bar$name = lapply(df.bar$Treatment, rename_fill)
df.bar$name = factor(df.bar$name, 
                     levels = c("OP50", "EGCG","Wine", "Lees"),
                     labels = c("OP50<br><i>E. coli</i>", "200uM<br>EGCG","Wine", "Lees"))

df.bar$type = lapply(df.bar$Treatment, rename_lab_or_field)
df.bar$type = factor(df.bar$type, 
                     levels = c("Control","Field", "Lab"),
                     labels = c("Control","Field", "Lab"))
df.bar

color_scheme = c("#f25c78",  "#84d94e","#ffa600", "#b64cdf")
p = ggplot(df.bar, aes(name,mean.val,fill=name)) +
  theme_classic() + labs(x = '', y = 'Responsiveness') + 
  theme(axis.title.y = element_text(size = 14), axis.text.x=element_text(size=13), axis.text.y=element_text(size=12)) +
  geom_bar(stat="identity", color="black",position=position_dodge()) +
  geom_errorbar(aes(ymin=mean.val - sem.val, ymax=mean.val + sem.val), width=.2,position=position_dodge(.9)) +
  facet_grid('. ~ type', scales="free", switch = "x", space='free') +
  scale_fill_manual(values=color_scheme) +
  theme(legend.position="none") + 
  theme(strip.placement = "outside") +
  theme(axis.title.x = element_blank(),
        axis.title.y =element_text(size=16),
        plot.title=element_text(hjust=0.5),
        axis.text.x = element_markdown(size=12),
        axis.text.y = element_text(size=11),
        legend.position="none",
        strip.text.x = element_text(size=13, face = "bold"),
        strip.background = element_blank()) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1.02))

p


### Saving the plot
png(filename = file.path(getwd(),"/Gentle-Tail-Touch.png"), width = 5, height = 3, units = "in", res = 600)
plot(p)
dev.off()

write.csv(as.data.frame(res$`CTRL`), file="Gentle-tail-touch-stat-results.csv")
write.csv(df.csv.summary, file="Gentle-tail-touch-values.csv")

print(paste("Saved at ", getwd()))

####################################################################################
######################################  NOSE TAP TOUCH  ######################################
## Import file name
## Change the title
df.csv = read.csv(file.path(getwd(),"Nose Touch Test.csv"), header = T)
df.csv$Treatment <- factor(df.csv$Treatment, levels = unique(df.csv$Treatment))

## Summarize the data
df.csv.summary <- df.csv %>% 
  group_by(Treatment) %>%
  summarise(mean.val = mean(Mean)/100,
            sem.val = std.error(Mean)/100,
            obs.val = n()) %>%
  as.data.frame()


## Using Dunnett's Tests (Apply statistics)
## Run first & determine the significant value
## Only those with stars
comp = lm( Mean ~ Treatment, data = df.csv)
summary(aov(comp))
res = DunnettTest(x=df.csv$Mean, g=df.csv$Treatment)
res

### Bar plots
df.bar = subset(df.csv.summary, Treatment %in% c("CTRL", "0.01% EGCG", "FWH", 
                                                 "FLH", "LWH", "LLH"))

df.bar$name = lapply(df.bar$Treatment, rename_fill)
df.bar$name = factor(df.bar$name, 
                     levels = c("OP50", "EGCG","Wine", "Lees"),
                     labels = c("OP50<br><i>E. coli</i>", "200uM<br>EGCG","Wine", "Lees"))

df.bar$type = lapply(df.bar$Treatment, rename_lab_or_field)
df.bar$type = factor(df.bar$type, 
                     levels = c("Control","Field", "Lab"),
                     labels = c("Control","Field", "Lab"))
df.bar

color_scheme = c("#f25c78",  "#84d94e","#ffa600", "#b64cdf")
p = ggplot(df.bar, aes(name,mean.val,fill=name)) +
  theme_classic() + labs(x = '', y = 'Responsiveness') + 
  theme(axis.title.y = element_text(size = 14), axis.text.x=element_text(size=13), axis.text.y=element_text(size=12)) +
  geom_bar(stat="identity", color="black",position=position_dodge()) +
  geom_errorbar(aes(ymin=mean.val - sem.val, ymax=mean.val + sem.val), width=.2,position=position_dodge(.9)) +
  facet_grid('. ~ type', scales="free", switch = "x", space='free') +
  scale_fill_manual(values=color_scheme) +
  theme(legend.position="none") + 
  theme(strip.placement = "outside") +
  theme(axis.title.x = element_blank(),
        axis.title.y =element_text(size=16),
        plot.title=element_text(hjust=0.5),
        axis.text.x = element_markdown(size=12),
        axis.text.y = element_text(size=11),
        legend.position="none",
        strip.text.x = element_text(size=13, face = "bold"),
        strip.background = element_blank()) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1.02))

p

### Saving the plot
png(filename = file.path(getwd(),"/Nose-Touch.png"), width = 5, height = 3, units = "in", res = 600)
plot(p)
dev.off()

write.csv(as.data.frame(res$`CTRL`), file="Nose-tap-touch-stat-results.csv")
write.csv(df.csv.summary, file="Nose-tap-touch-values.csv")

print(paste("Saved at ", getwd()))

####################################################################################

####################################################################################
