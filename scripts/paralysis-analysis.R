## Install packages for Kaplan-Meier plots
# install.packages("survivalAnalysis")
# install.packages("survminer")

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
library("plotrix")     
library(AICcmodavg)
library(DescTools)
library(ggpattern)
library("BSDA")


## Set working directory to the file location using Rstudio API
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

## Load survival csv file
data = read.csv("CL4176-data.csv", header = T)

controls = c("Neg", "OP50", "DMSO", "EGCG")
treatments = c("FWL", "FWM", "FWH",
               "LWL", "LWM", "LWH",
               "FLL", "FLM", "FLH",
               "LLL", "LLM", "LLH")

## Pattern the table from the previous KM plots.
## Drop the missing and consider only those alive and dead
## Arrange the tables as necessary
colnames(data) = c("group", "alive", "dead", "missing","day", "trial", "type")
## Important note here:
## Combined missing + dead data for consistency
data$comb = data$dead + data$missing
ndata = subset(data, select = -c(missing))


##########################################################
###################### CONTROL ###########################
##########################################################
## Define a function that will create a survival/event table
## Per row will be a biological sample that either survived (0) or dead (1).
## No censored will be applied here. 
newdat.controls = ndata[ndata$type=="Control",]
newdat.controls.transform = data.frame(newdat.controls[rep(seq_len(dim(newdat.controls)[1]), newdat.controls$comb), c('day','comb','group', 'trial', 'type'), drop = FALSE], row.names=NULL)
newdat.controls.transform$group = factor(newdat.controls.transform$group, levels=c("Neg", "OP50", "DMSO", "EGCG"), labels=c("Negative", "Food", "Vehicle (DMSO)", "Positive (EGCG)"))
newdat.controls.transform$status = 1

## Fit survival function
fit <- survfit(Surv(newdat.controls.transform$day, newdat.controls.transform$status) ~ group, data = newdat.controls.transform)
fit

## Graph the K-M plot.
color_scheme = c("#D99879", "#BA135D", "#5F9DF7","black")
g = ggsurvplot(fit, 
               # risk.table = TRUE, # Add risk table
               # risk.table.col = "group", # Change risk table color by groups
               # linetype = "strata", # Change line type by groups
               ggtheme = theme_classic(), # Change ggplot2 theme
               palette = color_scheme,
               xlab="Time in days",
               ylab="Paralysis probability",
               size = 0.8,
               legend = c(0.2, 0.73),
               legend.title = element_blank(),
               legend.labs = c("Food + 15\u00B0C", "Food + 25\u00B0C", "Vehicle (DMSO) + 25\u00B0C", "Antioxidant (EGCG) + 25\u00B0C"),
               xlim = c(0,5),
               font.legend = c(11, "plain", "black"),
               font.x = c(15, 'plain', 'black'),
               font.y = c(15, 'plain', 'black'),
               font.tickslab = c(11, 'plain', 'black'),
               surv.scale = 'percent',
               fun = "event",
) + guides(colour = guide_legend(nrow = 4))
g


#### Mean paralysis index
newdat.controls.transform
final.alive = newdat.controls[newdat.controls$day==5,]
final.alive.transform = data.frame(final.alive [rep(seq_len(dim(final.alive)[1]), final.alive$alive), c('day','alive','group', 'trial', 'type'), drop = FALSE], row.names=NULL)
final.alive.transform$group = factor(final.alive.transform$group, levels=c("Neg", "OP50", "DMSO", "EGCG"), labels=c("Negative", "Food", "Vehicle (DMSO)", "Positive (EGCG)"))
final.alive.transform$lifetime = final.alive.transform$day

newdat.controls.transform$lifetime = newdat.controls.transform$day -1

final.alive.transform = subset(final.alive.transform, select = -c(alive))
newdat.controls.transform = subset(newdat.controls.transform, selec = -c(comb, status))

ndf = rbind(newdat.controls.transform, final.alive.transform)

df.trial = ndf %>% group_by(group,trial) %>% summarise(lifetime = mean(lifetime), count = n())
df = df.trial %>% group_by(group) %>% summarise(mean = mean(lifetime), sem = std.error(lifetime), count = n())
df$upper = df$mean + df$sem
df$lower = df$mean - df$sem
df

## Visualize plot
## Using bar plots of means. Test for significance will be t-test
## 1. Checked anova. p-val = 1.34e-6 
anova.res = oneway.test(lifetime ~ group, data = ndf, var.equal = FALSE) # assuming equal variances
anova.res

## 2. Used pairwise independent two-tailed t-test with hochberg correction
stat.test = ndf %>% t_test(lifetime ~ group, p.adjust.method = "hochberg", var.equal = F)
stat.test


## 1. Manual addition of significance bars (preferred)
p = ggplot(df, aes(group,mean,fill=group)) + 
  theme_classic() + labs(x = '', y = 'Ave. time of paralysis onset') + 
  theme(axis.title.y = element_text(size = 14), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12)) +
  geom_bar(stat="identity", color="black",position=position_dodge()) +
  geom_errorbar(aes(ymin=mean, ymax=upper), width=.2,position=position_dodge(.9)) +
  scale_fill_manual(values=color_scheme) +
  theme(legend.position="none") + ylim(c(0,7)) +
  scale_x_discrete(labels = c("Food \n+ 15\u00B0C", "Food \n+ 25\u00B0C", "Vehicle (DMSO) \n+ 25\u00B0C", "Antioxidant (EGCG) \n+ 25\u00B0C"))
p

## The following results are given:
## Negative vs Food: p < 0.001
## Negative vs DMSO: p < 0.001
## Negative vs EGCG: p < 0.01
## Food vs DMSO: p = 0.221 ns
## Food vs EGCG: p = 0.048 
## DMSO vs EGCG: p < 0.01 **

## Significance bars
## Only with ***
neg.vs.food.p <- tibble(
  x = c(1, 1, 1.98, 1.98),
  y = c(5.1, 5.7, 5.7, 3.6)
) ## p < 0.001

neg.vs.egcg.p <- tibble(
  x = c(1, 1, 4, 4),
  y = c(6.2, 6.8, 6.8, 6.5)
) ## p = 0.000081

food.vs.dmso.p <- tibble(
  x = c(2.02, 2.02, 2.98, 2.98),
  y = c(3.6, 5.7, 5.7, 3.1)
) ## p = n.s

food.vs.egcg.p <- tibble(
  x = c(2, 2, 4, 4),
  y = c(6.0, 6.3, 6.3, 4.5)
) ## p = 0.02


p = p + 
  geom_line(data = neg.vs.food.p, aes(x= x, y = y, group=1), inherit.aes = F) +
  annotate("text", x = 1.5, y = 5.8, label = "***", size = 5, color = "#22292F") +
  
  geom_line(data = neg.vs.egcg.p, aes(x= x, y = y, group=1), inherit.aes = F) +
  annotate("text", x = 2.5, y = 6.9, label = "**", size = 5, color = "#22292F") +
  
  geom_line(data = food.vs.dmso.p, aes(x= x, y = y, group=1), inherit.aes = F) +
  annotate("text", x = 2.5, y = 5.9, label = "n.s", size = 4, color = "#22292F") +
  
  geom_line(data = food.vs.egcg.p, aes(x= x, y = y, group=1), inherit.aes = F) +
  annotate("text", x = 3.0, y = 6.4, label = "*", size = 5, color = "#22292F") 

p

## Combine

tiff(filename = file.path(getwd(),"/KM-Plot-Paralysis.tiff"), width = 7, height = 7, units = "in", res = 600)
kmbar = ggarrange(plotlist = list(g$plot, p), ncol = 1, nrow = 2,
                  labels = c("A", "B"))
plot(kmbar)
dev.off()

png(filename = file.path(getwd(),"/KM-Plot-Paralysis.png"), width = 7, height = 7, units = "in", res = 600)
kmbar = ggarrange(plotlist = list(g$plot, p), ncol = 1, nrow = 2,
                  labels = c("A", "B"))
plot(kmbar)
dev.off()
print(paste("Saved at ", getwd()))


##################################################################
###################### TREATMENTS ################################
##################################################################
newdat.treatment = ndata[!(ndata$type %in% c("Sablan", "Kapangan", "La Trinidad")),]
newdat.treatment = newdat.treatment[!(newdat.treatment$group %in% c("Neg", "DMSO")),]
newdat.treatment.transform = data.frame(newdat.treatment[rep(seq_len(dim(newdat.treatment)[1]),newdat.treatment$comb), c('day','comb','group', 'trial', 'type'), drop = FALSE], row.names=NULL)
newdat.treatment.transform$group = factor(newdat.treatment.transform$group,
                                          levels=c("OP50", "EGCG", "FWL", "FWM", "FWH",
                                                   "LWL", "LWM", "LWH",
                                                   "FLL", "FLM", "FLH",
                                                   "LLL", "LLM", "LLH"),
                                          labels=c("OP50", "EGCG", "FWL", "FWM", "FWH",
                                                   "LWL", "LWM", "LWH",
                                                   "FLL", "FLM", "FLH",
                                                   "LLL", "LLM", "LLH"))
newdat.treatment.transform$type= factor(newdat.treatment.transform$type,
                                          levels=c("Control", "Field Wine",
                                                   "Lab Wine",
                                                   "Field Lees",
                                                   "Lab Lees"),
                                          labels=c("Control", "Field Wine",
                                                   "Lab Wine",
                                                   "Field Lees",
                                                   "Lab Lees"))
newdat.treatment.transform$lifetime = newdat.treatment.transform$day -1


## Get the contribution of worms who were alive after 5 days
final.alive.treatment = newdat.treatment[newdat.treatment$day==5,]
final.alive.treatment.transform = data.frame(final.alive.treatment[rep(seq_len(dim(final.alive.treatment)[1]), final.alive.treatment$alive), c('day','alive','group', 'trial', 'type'), drop = FALSE], row.names=NULL)
final.alive.treatment.transform$group = factor(final.alive.treatment.transform$group,
                                          levels=c("OP50", "EGCG","FWL", "FWM", "FWH",
                                                   "LWL", "LWM", "LWH",
                                                   "FLL", "FLM", "FLH",
                                                   "LLL", "LLM", "LLH"),
                                          labels=c("OP50", "EGCG","FWL", "FWM", "FWH",
                                                   "LWL", "LWM", "LWH",
                                                   "FLL", "FLM", "FLH",
                                                   "LLL", "LLM", "LLH"))
final.alive.treatment.transform$type= factor(final.alive.treatment.transform$type,
                                        levels=c("Control", "Field Wine",
                                                 "Lab Wine",
                                                 "Field Lees",
                                                 "Lab Lees"),
                                        labels=c("Control", "Field Wine",
                                                 "Lab Wine",
                                                 "Field Lees",
                                                 "Lab Lees"))

final.alive.treatment.transform$lifetime = final.alive.treatment.transform$day
final.alive.treatment.transform = subset(final.alive.treatment.transform, select = -c(alive))
newdat.treatment.transform = subset(newdat.treatment.transform, selec = -c(comb))

## Combine
ndf.treatment = rbind(newdat.treatment.transform, final.alive.treatment.transform)
df.treatment.trial = ndf.treatment %>% group_by(group, type, trial) %>% summarise(lifetime = mean(lifetime), count = n())
df.treatment = df.treatment.trial %>% group_by(group, type) %>% summarise(mean = mean(lifetime), sem = std.error(lifetime), count = n())

df.treatment$upper = df.treatment$mean + df.treatment$sem
df.treatment$lower = df.treatment$mean - df.treatment$sem
df.treatment

## Visualize plot
##Use bar plots

## Controls (OP50 & EGCG)
color_scheme = c("#fa5cf2", "#ffa600", "#84d94e","#f25c78", "#b64cdf")

# # fill_scheme = c("#F8EDE3", "#DFD3C3", "#DFD3C3", "#FBF8F1")
# fill_scheme = c("#F8EDE3", "#F8EDE3", "#F8EDE3", "#F8EDE3")
# df.treatment.ctrl = df.treatment[df.treatment$type == "Control", ] 
# y1 = ggplot(df.treatment.ctrl, aes(group, mean, fill=group)) +
#   theme_classic() + labs(x = '', y = 'Average time of paralysis onset', title="Control") +
#   geom_bar(stat="identity", color="black",position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean, ymax=upper), width=.2,position=position_dodge(.9)) +
#   scale_x_discrete(labels = c("Food\n+ 25\u00B0C", "EGCG")) +
#   scale_fill_manual(values=c(color_scheme[2],"black")) +
#   scale_y_continuous(breaks = c(0, 2, 4, 6),
#                      labels = c(0, 2, 4, 6),
#                      limits = c(0,6))+
#   theme(legend.position="none",
#         plot.title = element_text(hjust = 0.5),
#         plot.margin = unit(c(1, 0, 0, 0), "cm"),
#         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
# 
# y1
# 
# ## Field Wine
# df.treatment.fw = df.treatment[df.treatment$type == "Field Wine", ] 
# y2 = ggplot(df.treatment.fw, aes(group, mean)) +
#   theme_classic() + labs(x = '', y = '', title="Field Wine") +
#   geom_bar_pattern(stat = "identity",
#                    fill = fill_scheme[1],
#                    colour="black",
#                    pattern_color = "black",
#                    pattern_fill = "black",
#                    aes(pattern = group)) +
#   scale_x_discrete(labels = c("Low","Mid", "High")) +
#   scale_y_continuous(breaks = c(0, 2, 4, 6),
#                      labels = c("", "", "", ""),
#                      limits = c(0,6))+
#   geom_errorbar(aes(ymin=mean, ymax=upper), width=.2,position=position_dodge(.9)) +
#   theme(legend.position="none",
#         plot.title = element_text(hjust = 0.5),
#         plot.margin = unit(c(0, 0, 0, 0), "cm"),
#         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
# 
# y2
# 
# 
# 
# ## Field Lees
# df.treatment.fl = df.treatment[df.treatment$type == "Field Lees", ] 
# y3 = ggplot(df.treatment.fl, aes(group, mean)) +
#   theme_classic() + labs(x = '', y = '', title="Field Lees") +
#   geom_bar_pattern(stat = "identity",
#                    fill = fill_scheme[2],
#                    colour="black",
#                    pattern_color = "black",
#                    pattern_fill = "black",
#                    aes(pattern = group)) +
#   scale_x_discrete(labels = c("Low","Mid", "High")) +
#   scale_y_continuous(breaks = c(0, 2, 4, 6),
#                      labels = c("", "", "", ""),
#                      limits = c(0,6))+
#   geom_errorbar(aes(ymin=mean, ymax=upper), width=.2,position=position_dodge(.9)) +
#   theme(legend.position="none",
#         plot.margin = unit(c(0, 0, 0, 0), "cm"),
#         plot.title = element_text(hjust = 0.5),
#         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
# 
# y3
# 
# ## Lab Wine
# df.treatment.lw = df.treatment[df.treatment$type == "Lab Wine", ] 
# y4 = ggplot(df.treatment.lw, aes(group, mean)) +
#   theme_classic() + labs(x = '', y = '', title="Lab Wine") +
#   geom_bar_pattern(stat = "identity",
#                    fill = fill_scheme[3],
#                    colour="black",
#                    pattern_color = "black",
#                    pattern_fill = "black",
#                    aes(pattern = group)) +
#   scale_x_discrete(labels = c("Low","Mid", "High")) +
#   scale_y_continuous(breaks = c(0, 2, 4, 6),
#                      labels = c("", "", "", ""),
#                      limits = c(0,6))+
#   geom_errorbar(aes(ymin=mean, ymax=upper), width=.2,position=position_dodge(.9)) +
#   theme(legend.position="none",
#         plot.margin = unit(c(0, 0, 0, 0), "cm"),
#         plot.title = element_text(hjust = 0.5),
#         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
# 
# y4
# 
# ## Lab Lees
# df.treatment.ll = df.treatment[df.treatment$type == "Lab Lees", ] 
# y5 = ggplot(df.treatment.ll, aes(group, mean)) +
#   theme_classic() + labs(x = '', y = '', title="Lab Lees") +
#   geom_bar_pattern(stat = "identity",
#                    fill = fill_scheme[4],
#                    colour="black",
#                    pattern_color = "black",
#                    pattern_fill = "black",
#                    aes(pattern = group)) +
#   scale_x_discrete(labels = c("Low","Mid", "High")) +
#   scale_y_continuous(breaks = c(0, 2, 4, 6),
#                      labels = c("", "", "", ""),
#                      limits = c(0,6))+
#   geom_errorbar(aes(ymin=mean, ymax=upper), width=.2,position=position_dodge(.9)) +
#   theme(legend.position="none",
#         plot.margin = unit(c(0, 0.5, 0, 0), "cm"),
#         plot.title = element_text(hjust = 0.5),
#         axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
# 
# y5



## Using bar plots of means. Test for significance will be t-test
## 1. Checked anova. 
comp.treatment.group = lm( lifetime ~ group, data = ndf.treatment)
comp.treatment.type = lm( lifetime ~ type, data = ndf.treatment)
comp.treatment.comb = lm( lifetime ~ group*type, data = ndf.treatment)
comp.treatement.notox.models <- list(
  comp.treatment.group,
  comp.treatment.type,
  comp.treatment.comb)
aictab(cand.set = comp.treatement.notox.models, modnames = c('group', 'type', 'comb'))

tx.mod = aov( lifetime ~ group, data = ndf.treatment)
summary(tx.mod)

op50.tukey.res = TukeyHSD(tx.mod, conf.level=.95)
op50.tukey.res

## Results for the group
# DunnettTest(x=df.treatment.trial$lifetime, g=df.treatment.trial$group)

mean.x = df.treatment[df.treatment$group=="OP50",]$mean
s.x = df.treatment[df.treatment$group=="OP50",]$sem

tx.vector = c( "EGCG","FWL", "FWM", "FWH",
                      "LWL", "LWM", "LWH",
                      "FLL", "FLM", "FLH",
                      "LLL", "LLM", "LLH")

res.pval = data.frame(matrix(ncol = 8, nrow = length(tx.vector)))
colnames(res.pval) <- c('Comparison','sig',  'mean.x', 'sem.x','mean.y', 'sem.y', 't-value', 'p.val' )
for (i in c(1:length(tx.vector))) {
  tx.name = tx.vector[i]
  mean.y = df.treatment[df.treatment$group==tx.name,]$mean
  s.y = df.treatment[df.treatment$group==tx.name,]$sem
  res = tsum.test(
    mean.y,s.y,n.y = 3,
    mean.x,s.x,n.x=3,
    alternative = "greater",
    var.equal = FALSE,
    conf.level = 0.95
  )
  res.pval[i,1] = tx.name
  res.pval[i,3] = mean.x
  res.pval[i,4] = s.x
  res.pval[i,5] = mean.y
  res.pval[i,6] = s.y
  res.pval[i,7] = res$statistic
  res.pval[i,8] = res$p.val
  if(res$p.val >= 0.05) {
    res.pval[i,2] = "n.s"
  }
  if(res$p.val < 0.05) {
    res.pval[i,2] = "*"
  }
  if(res$p.val <= 0.01) {
    res.pval[i,2] = "**"
  }
  if(res$p.val <= 0.001) {
    res.pval[i,2] = "***"
  }
}

## Anova has showed independence through different groups
## Using independent t-test, we tested whether each group would have a greater mean time of paralysis onset
## compared to the control (Food + 25C)
res.pval
## p < 0.001 : EGCG, FWL
## p < 0.01 : FWM, LWH, FLL, LLM, LLH

## EGCG
y1 = y1 +
  annotate("text", x = 2, y = 4.2, label = "***", size = 5, color = color_scheme[2])
y1

## FWL, FWM
y2 = y2 +
  annotate("text", x = 1, y = 4.0, label = "***", size = 5, color = color_scheme[2]) +
  annotate("text", x = 2, y = 3.8, label = "**", size = 5, color = color_scheme[2])
y2

## FLL
y3 = y3 + 
  annotate("text", x = 1, y = 3.8, label = "**", size = 5, color = color_scheme[2])
y3

## LWH
y4 = y4 + 
  annotate("text", x = 3, y = 4.1, label = "**", size = 5, color = color_scheme[2])
y4

## LLM, LLH
y5 = y5 +
  annotate("text", x = 2, y = 3.8, label = "**", size = 5, color = color_scheme[2]) +
  annotate("text", x = 3, y = 4.0, label = "**", size = 5, color = color_scheme[2])
y5


## Visualize part 2
df.treatment
df.high = subset(df.treatment, group %in% c("OP50", "EGCG", "FWH", "LWH", "FLH", "LLH"))

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
  else if( x == "OP50") {
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

df.high$name = lapply(df.high$group, rename_fill)
df.high$name = factor(df.high$name, 
                       levels = c("OP50", "EGCG","Wine", "Lees"),
                       labels = c("OP50 <i>E. coli</i>", "200uM EGCG","Wine", "Lees"))

df.high$type = lapply(df.high$group, rename_lab_or_field)
df.high$type = factor(df.high$type, 
                      levels = c("Control","Field", "Lab"),
                      labels = c("Control","Field", "Lab"))
df.high

color_scheme = c("#f25c78",  "#84d94e","#ffa600", "#b64cdf")
p = ggplot(df.high, aes(name,mean,fill=name)) +
  theme_classic() + labs(x = '', y = 'Ave. time of paralysis onset') + 
  theme(axis.title.y = element_text(size = 14), axis.text.x=element_text(size=13), axis.text.y=element_text(size=12)) +
  geom_bar(stat="identity", color="black",position=position_dodge()) +
  geom_errorbar(aes(ymin=mean - sem, ymax=mean + sem), width=.2,position=position_dodge(.9)) +
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
  ylim(0,5)

p

p = p +
  geom_text(data = subset(df.high, type=="Control"),
            aes(x = 2, y = 4.2, label = "**"), size=7, color = "#960101") +
  geom_text(data = subset(df.high, type=="Lab"),
          aes(x = 2, y = 3.9, label = "*"), size=7, color = "#960101") +
  geom_text(data = subset(df.high, type=="Lab"),
            aes(x = 1, y = 4.0, label = "**"), size=7, color = "#960101")
p
  
ggsave(plot=p, filename="Final-paralysis-barplot.png", width=8, height=4, dpi=600)
write.csv(df.treatment, "paralysis-high.csv")

# ###### 2, Check among those red dots were greater or equal to EGCG
# ###### 
# ###### 
# mean.x = df.treatment[df.treatment$group=="EGCG",]$mean
# s.x = df.treatment[df.treatment$group=="EGCG",]$sem
# 
# tx.vector2 = c("FWL", "FWM",
#                "LWH", 
#                "FLL",
#                "LLM", "LLH")
# 
# res.pval2 = data.frame(matrix(ncol = 8, nrow = length(tx.vector2)))
# 
# colnames(res.pval) <- c('Comparison','sig',  'mean.x', 'sem.x','mean.y', 'sem.y', 't-value', 'p.val' )
# for (i in c(1:length(tx.vector2))) {
#   tx.name = tx.vector2[i]
#   mean.y = df.treatment[df.treatment$group==tx.name,]$mean
#   s.y = df.treatment[df.treatment$group==tx.name,]$sem
#   res = tsum.test(
#     mean.y,s.y,n.y = 3,
#     mean.x,s.x,n.x=3,
#     alternative = "two.sided",
#     var.equal = FALSE,
#     conf.level = 0.95
#   )
#   res.pval2[i,1] = tx.name
#   res.pval2[i,3] = mean.x
#   res.pval2[i,4] = s.x
#   res.pval2[i,5] = mean.y
#   res.pval2[i,6] = s.y
#   res.pval2[i,7] = res$statistic
#   res.pval2[i,8] = res$p.val
#   if(res$p.val >= 0.05) {
#     res.pval2[i,2] = "n.s"
#   }
#   if(res$p.val < 0.05) {
#     res.pval2[i,2] = "*"
#   }
#   if(res$p.val <= 0.01) {
#     res.pval2[i,2] = "**"
#   }
#   if(res$p.val <= 0.001) {
#     res.pval2[i,2] = "***"
#   }
# }
# res.pval2
# 
# ## LWH
# y4 = y4 + 
#   annotate("text", x = 3, y = 4.7, label = "n.s", size = 4, color = "black")
# y4
# 
# ## LLM, LLH
# y5 = y5 +
#   annotate("text", x = 3, y = 4.6, label = "n.s", size = 4, color = "black")
# y5
# 
# 
# # ## Saving the figure for short term memory
# paralysis.plot = plot_grid(y1, y2, y3, y4, y5, align = "h", ncol = 5) +
#   theme(plot.margin = unit(c(0,0,0,0.15), "cm")) 
# 
# tiff(filename = file.path(getwd(),"/Treatments.tiff"), width = 8, height = 4, units = "in", res = 600)
# plot(paralysis.plot)
# dev.off()
# print(paste("Saved at ", getwd()))
# 
# png(filename = file.path(getwd(),"/Treatments.png"), width = 8, height = 4, units = "in", res = 600)
# plot(paralysis.plot)
# dev.off()
# print(paste("Saved at ", getwd()))
# paralysis.plot
# 
# 
# 
# fit<- survfit(Surv(time, status) ~ sex, data = lung)
# 
# # Basic survival curves
# # ggsurvplot(fit, data = lung,linetype = "strata")$plot
# 
# #your own choice of linetype
# ggsurvplot(fit, data = lung,
#            legend.labs = c("Male","Female"),
#            linetype = "strata")$plot +
#   scale_linetype_manual(values = c("dotted",
#                                    "dashed")) +
#   scale_colour_manual(values = c("black",
#                                  "red")) +
#   theme(legend.title = element_blank(),
#         legend.text = element_markdown(size=12))

