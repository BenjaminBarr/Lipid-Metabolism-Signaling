library(here) # here makes a project transportable
library(janitor) # clean_names
library(readxl) # read excel, duh!
library(data.table) # magical data frames
library(magrittr) # pipes
library(stringr) # string functions
library(forcats) # factor functions

# analysis packages
library(emmeans) # the workhorse for inference
library(nlme) # gls and some lmm
library(lme4) # linear mixed models
library(lmerTest) # linear mixed model inference
library(afex) # ANOVA linear models
library(glmmTMB) # generalized linear models
library(MASS) # negative binomial and some other functions
library(car) # model checking and ANOVA
library(DHARMa) # model checking

# graphing packages
library(ggsci) # color palettes
library(ggpubr) # publication quality plots
library(ggforce) # better jitter
#library(cowplot) # combine plots
library(knitr) # kable tables
library(kableExtra) # kable_styling tables

library(insight)
library(lazyWeave)

# Data Read In

library(ggplot2)
library(reshape2)
library(matrixStats)
library(dplyr)
library(AICcmodavg)
library(tidyverse)
library(rstatix)
library(ggdark)
library(multcompView)
library(RColorBrewer)
library(patchwork)
library(cowplot)
library(gridExtra)

getwd()

gluc <- read_xlsx("./Metabolites/Metabolites Data/Metabolites Neg.xlsx", sheet = "Glucose")
gluc <- data.frame(gluc)

lactate <- read_xlsx("./Metabolites/Metabolites Data/Metabolites Neg.xlsx", sheet = "Lactate")
lactate <- data.frame(lactate)

ala <- read_xlsx("./Metabolites/Metabolites Data/Metabolites Pos.xlsx", sheet = "Alanine Pos")
ala <- data.frame(ala)

gln <- read_xlsx("./Metabolites/Metabolites Data/Metabolites Pos.xlsx", sheet = "Glutamine Pos")
gln <- data.frame(gln)

carn <- read_xlsx("./Metabolites/Metabolites Data/Metabolites Pos.xlsx", sheet = "Carnitine Pos")
carn <- data.frame(carn)

gluc <- melt(gluc)
lactate <- melt(lactate)
ala <- melt(ala)
gln <- melt(gln)
carn <- melt(carn)


met <- cbind(gluc, lactate[3], ala[3], gln[3], carn[3])
names(met) <- c("label", "variable", "glucose", "lactate", 
                "alanine", "glutamine", "carnitine")

met.labs <- str_split_fixed(met$label, " ", 3)
met.labs <- data.frame(met.labs)
names(met.labs) <- c("age", "group", "sex")
met.labs

met <- cbind(met, met.labs)
met <- na.omit(met)
met <- met[, -2]
met

met$fat <- "HF"
met$fat[which(met$group %in% c("G1", "G2"))] <- "C"
met

met$protein <- "C"
met$protein[which(met$group %in% c("G3", "G4"))] <- "Beef"
met

met$trt <- "yes"
met$trt[which(met$group %in% c("G2", "G4", "G6"))] <- "no"
met

f.met <- met[which(met$sex == "F"),] # Females
hff.met <- f.met[which(f.met$group %in% c("G3", "G4", "G5", "G6")),] #HF only

hff.met$group[which(hff.met$group == "G3")] <- "group 3"
hff.met$group[which(hff.met$group == "G4")] <- "group 4"
hff.met$group[which(hff.met$group == "G5")] <- "group 5"
hff.met$group[which(hff.met$group == "G6")] <- "group 6"

hff.met

hff6 <- hff.met[which(hff.met$age == "H1"),]
hff12 <- hff.met[which(hff.met$age == "H2"),]
hff18 <- hff.met[which(hff.met$age == "H3"),]

#Glucose First#
ggqqplot(hff6, x = "glucose", color = "group")

hf.glucose.f6.lm <- lm(glucose ~ trt * protein, data = hff6)

anova(hf.glucose.f6.lm)
summary(hf.glucose.f6.lm)


f6.glucose.em <- emmeans(hf.glucose.f6.lm, specs = c("trt", "protein"))
f6.glucose.em

#f6.glucose.dt <- data.table(summary(f6.glucose.em))
#f6.glucose.dt

#emmeans(hf.glucose.f6.lm, specs = c("trt", "protein"))

## Main Effects Only ##

hf.glu.ahe.f6.lm <- lm(glucose ~ trt, data = hff6)

summary(hf.glu.ahe.f6.lm)
f6.glu.ahe.em <- emmeans(hf.glu.ahe.f6.lm, specs = "trt")
emmeans(f6.glu.ahe.em, pairwise ~ trt)
hff6 %>% group_by(trt) %>%
  summarise(emmean = mean(glucose),
            SD = sd(glucose),
            SE = sd(glucose)/sqrt(length(glucose)),
            n = length(glucose))


hf.glu.p.f6.lm <- lm(glucose ~ protein, data = hff6)

summary(hf.glu.p.f6.lm)
f6.glu.p.em <- emmeans(hf.glu.p.f6.lm, specs = "protein")
emmeans(f6.glu.p.em, pairwise ~ protein)
hff6 %>% group_by(protein) %>%
  summarise(emmean = mean(glucose),
            SD = sd(glucose),
            SE = sd(glucose)/sqrt(length(glucose)),
            n = length(glucose))
anova(hf.glu.ahe.f6.lm)
anova(hf.glu.p.f6.lm)

f6.glu.me <- rbind(anova(hf.glu.ahe.f6.lm),anova(hf.glu.p.f6.lm))

# Graphing labels ()

### Customize Here ###
f6g.lab <- paste("group", c(3:6))
f6g.lab <- data.frame(f6g.lab)
f6g.lab$lab <- c("", "B", "A", "A")
names(f6g.lab)[1] <- "group1"
f6g.lab$group2 <- NA
f6g.lab$group <- f6g.lab$group1


# Customize for each set of comp
glab <- c(`group 3` = "HFBN",
          `group 4` = "HFB",
          `group 5` = "HFCN",
          `group 6` = "HFC")

#Plot
f6.glu.sum <- hff6 %>% group_by(group) %>%
  summarise(emmean = mean(glucose),
            SD = sd(glucose),
            SE = sd(glucose)/sqrt(length(glucose)),
            n = length(glucose))

g3 <- data.frame(group = "group 3", emmean = 0, SD = 0, SE = 0, n = 0)
f6.glu.sum <- rbind(g3, f6.glu.sum)

f6g.lab$y <- f6.glu.sum$emmean + f6.glu.sum$SD + 5

f6.glu.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
#f6.glucose.dt <- f6.glucose.dt[order(f6.glucose.dt$group, decreasing = T),]
f6.glu.sum <- f6.glu.sum[order(f6.glu.sum$group, decreasing = T),]

names(f6.glu.sum)[c(1,6)] <- c("lab", "group")

f6.glu.plot <- ggplot(data = f6.glu.sum, 
                      aes(x = factor(group, level = c("group 6", "group 5", 
                                                      "group 4", "group 3")),
                          y = emmean, fill = group))+
geom_col(fill = c("brown1", "deepskyblue3", "chartreuse4","purple"), width = 0.5, 
           color = "black") +
  scale_x_discrete(label = glab) +
  scale_y_continuous(limits = c(0,40), expand = expansion(mult = c(0,.10))) +
  theme_pubr() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title.y = element_text(size = 16))+
  geom_errorbar(data = f6.glu.sum, aes(y = emmean, ymin = emmean, ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  geom_point(data = hff6, aes(x = group, y = glucose), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 4),
             size = 1, stroke = 1) +
  annotate("segment", x = 1, xend = 2, y = 26, yend = 26) +
  annotate("text", x = 1.5, y = 28, label = "A", size = 5) +
  annotate("text", x = 3, y = 20, label = "B", size = 5) +
  annotate("text", x = .6, y = 40, fontface = 3, label = "DPS (Casein > Beef)", hjust = 0, size = 5) +
  #coord_cartesian(clip = "off") +
  #annotate("text", x = .6, y = 33, fontface = 3, label = "Main Effect: DPS (Casein > Beef)", hjust = 0) +
  labs(title = "Glucose",
       subtitle = "Females 6-Months",
       y ="[Glucose] 
(mg/mL)",
       x = "Diet")
f6.glu.plot

#### Glucose 12 Months ####
ggqqplot(hff12, x = "glucose", color = "group")

hf.glucose.f12.lm <- lm(glucose ~ trt * protein, data = hff12)
anova(hf.glucose.f12.lm)
summary(hf.glucose.f12.lm)

f12.glu.res <- residuals(hf.glucose.f12.lm)
ggqqplot(f12.glu.res)

f12.glucose.em <- emmeans(hf.glucose.f12.lm, specs = c("trt", "protein"))
f12.glucose.em

#f12.glucose.dt <- data.table(summary(f12.glucose.em))
#f12.glucose.dt

# Proper/Simple Post-Hoc

tp.plot <- emmip(hf.glucose.f12.lm, trt ~ protein)
tp <- emmeans(hf.glucose.f12.lm, pairwise ~ protein)
tp

pairs(f12.glucose.em, simple = "protein")

# Graphing labels ()
#f12.glucose.dt$group <- c("group 4", "group 3", "group 6", "group 5")

### Customize Here ###
f12.glu.sum <- hff12 %>% group_by(group) %>%
  summarise(emmean = mean(glucose),
            SD = sd(glucose),
            SE = sd(glucose)/sqrt(length(glucose)),
            n = length(glucose))

f12.glu.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
#f12.glucose.dt <- f12.glucose.dt[order(f12.glucose.dt$group, decreasing = T),]
f12.glu.sum <- f12.glu.sum[order(f12.glu.sum$group, decreasing = T),]

names(f12.glu.sum)[c(1,6)] <- c("lab", "group")

f12.glu.plot <- ggplot(data = f12.glu.sum, 
                       aes(x = factor(group, level = c("group 6", "group 5", 
                                                       "group 4", "group 3")),
                           y = emmean, fill = group))+
geom_col(fill = c("brown1", "deepskyblue3", "chartreuse4","purple"), width = 0.5, 
           color = "black") +
  scale_y_continuous(limits = c(0,40), expand = expansion(mult = c(0,.10))) +
  scale_x_discrete(label = glab) +
  theme_pubr() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24),
        plot.subtitle = element_text(hjust = 0.5))+
  geom_errorbar(data = f12.glu.sum, aes(y = emmean, ymin = emmean, ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  geom_point(data = hff12, aes(x = group, y = glucose), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 4),
             size = 1, stroke = 1) +
  annotate("segment", x = 1, xend = 2, y = 24, yend = 24) +
  annotate("text", x = 1.5, y = 26, label = "A", size = 5) +
  annotate("segment", x = 3, xend = 4, y = 30, yend = 30) +
  annotate("text", x = 3.5, y = 32, label = "B", size = 5) +
  annotate("text", x = .6, y = 40, fontface = 3, label = "DPS (Beef > Casein)",
           hjust = 0, size = 5) +
  coord_cartesian(clip = "off") +
  labs(title = "Glucose",
       subtitle = "Females 12-Months",
       y ="Concentration (mg/mL)",
       x = "Diet")
f12.glu.plot

#### Glucose 18 Months ####
ggqqplot(hff18, x = "glucose", color = "group")

hf.glucose.f18.lm <- lm(glucose ~ trt * protein, data = hff18)
anova(hf.glucose.f18.lm)
summary(hf.glucose.f18.lm)

f18.glu.res <- residuals(hf.glucose.f18.lm)
ggqqplot(f18.glu.res)

f18.glucose.em <- emmeans(hf.glucose.f18.lm, specs = c("trt", "protein"))
f18.glucose.em

#f18.glucose.dt <- data.table(summary(f18.glucose.em))
#f18.glucose.dt

# Proper/Simple Post-Hoc

tp.plot <- emmip(hf.glucose.f18.lm, trt ~ protein)
tp <- emmeans(hf.glucose.f18.lm, pairwise ~ protein)

pairs(f18.glucose.em, simple = "protein")
pairs(f18.glucose.em, simple = c("trt", "protein"))

# Graphing labels ()
#f18.glucose.dt$group <- c("group 4", "group 3", "group 6", "group 5")

# Within Group Indications #

### Customize Here ###
f18.glu.sum <- hff18 %>% group_by(group) %>%
  summarise(emmean = mean(glucose),
            SD = sd(glucose),
            SE = sd(glucose)/sqrt(length(glucose)),
            n = length(glucose))

f18.glu.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
#f18.glucose.dt <- f18.glucose.dt[order(f18.glucose.dt$group, decreasing = T),]
f18.glu.sum <- f18.glu.sum[order(f18.glu.sum$group, decreasing = T),]

names(f18.glu.sum)[c(1,6)] <- c("lab", "group")

f18.glu.plot <- ggplot(data = f18.glu.sum, 
                       aes(x = factor(group, level = c("group 6", "group 5", 
                                                       "group 4", "group 3")),
                           y = emmean, fill = group))+
geom_col(fill = c("brown1", "deepskyblue3", "chartreuse4","purple"), width = 0.5, 
           color = "black") +
  scale_y_continuous(limits = c(0,40), expand = expansion(mult = c(0,.10))) +
  scale_x_discrete(label = glab) +
  theme_pubr() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24),
        plot.subtitle = element_text(hjust = 0.5))+
  geom_errorbar(data = f18.glu.sum, aes(y = emmean, ymin = emmean, ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  geom_point(data = hff18, aes(x = group, y = glucose), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 4),
             size = 1, stroke = 1) +
  coord_cartesian(clip = "off") +
  labs(title = "Glucose",
       subtitle = "Females 18-Months",
       y ="Concentration (mg/mL)",
       x = "Diet")
f18.glu.plot

# Create your plot grid first (with cleaned-up plots)
f.glu.plot <- plot_grid(
  f6.glu.plot + 
    theme(
      plot.subtitle = element_blank(),
      plot.title = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    ),
  f12.glu.plot +
    theme(
      plot.subtitle = element_blank(),
      plot.title = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    ),
  f18.glu.plot +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    ),
  nrow = 1,
  align = "h",
  axis = "tblr",
  rel_widths = c(1, 1, 1)
)

glu.label.plot <- ggplot() + 
  theme_void() +
  annotate("text", 
           x = 0, y = 0.5, 
           label = "[Glucose]\n(mg/mL)",  # mg/mL
           angle = 90, 
           size = 5, 
           fontface = "bold",
           hjust = 0.5)
# Step 3: Combine the label and plot using plot_grid again
f.glu.plot.labeled <- plot_grid(
  glu.label.plot, f.glu.plot,
  ncol = 2,
  rel_widths = c(0.07, 1))  # Adjust label space vs plot


f.glu.plot.labeled

#Lactate
hf.lactate.f6.lm <- lm(lactate ~ trt * protein, data = hff6)
anova(hf.lactate.f6.lm)
summary(hf.lactate.f6.lm)

ggqqplot(hff6, x = "lactate", color = "group")

f6.lac.res <- residuals(hf.lactate.f6.lm)
ggqqplot(f6.lac.res)

f6.lactate.em <- emmeans(hf.lactate.f6.lm, specs = c("trt", "protein"))
f6.lactate.em

#f6.lactate.dt <- data.table(summary(f6.lactate.em))
#f6.lactate.dt

# Proper/Simple Post-Hoc

tp.plot <- emmip(hf.lactate.f6.lm, trt ~ protein)
tp <- emmeans(hf.lactate.f6.lm, pairwise ~ trt)
tp
pt.plot <- emmip(hf.lactate.f6.lm, protein ~ trt)
pt <- emmeans(hf.lactate.f6.lm, pairwise ~ protein)
pt

pairs(f6.lactate.em, simple = "trt")

## Main Effects Only ##

hf.lac.ahe.f6.lm <- lm(lactate ~ trt, data = hff6)
anova(hf.lac.ahe.f6.lm)
summary(hf.lac.ahe.f6.lm)
f6.lac.ahe.em <- emmeans(hf.lac.ahe.f6.lm, specs = "trt")
emmeans(f6.lac.ahe.em, pairwise ~ trt)
hff6 %>% group_by(trt) %>%
  summarise(emmean = mean(lactate),
            SD = sd(lactate),
            SE = sd(lactate)/sqrt(length(lactate)),
            n = length(lactate))


hf.lac.p.f6.lm <- lm(lactate ~ protein, data = hff6)
anova(hf.lac.p.f6.lm)
summary(hf.lac.p.f6.lm)
f6.lac.p.em <- emmeans(hf.lac.p.f6.lm, specs = "protein")
emmeans(f6.lac.p.em, pairwise ~ protein)
hff6 %>% group_by(protein) %>%
  summarise(emmean = mean(lactate),
            SD = sd(lactate),
            SE = sd(lactate)/sqrt(length(lactate)),
            n = length(lactate))
f6.lac.me <- rbind(anova(hf.lac.ahe.f6.lm),anova(hf.lac.p.f6.lm))
### Customize Here ###
f6g.lab <- paste("group", c(3:6))
f6g.lab <- data.frame(f6g.lab)
f6g.lab$lab <- c("", "A", "B", "A")
names(f6g.lab)[1] <- "group1"
f6g.lab$group2 <- NA
f6g.lab$group <- f6g.lab$group1

#Plot
f6.lac.sum <- hff6 %>% group_by(group) %>%
  summarise(emmean = mean(lactate),
            SD = sd(lactate),
            SE = sd(lactate)/sqrt(length(lactate)),
            n = length(lactate))

g3 <- data.frame(group = "group 3", emmean = 0, SD = 0, SE = 0, n = 0)
f6.lac.sum <- rbind(g3, f6.lac.sum)

f6g.lab$y <- f6.lac.sum$emmean + f6.lac.sum$SD + 0.3

f6.lac.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
#f6.lactate.dt <- f6.lactate.dt[order(f6.lactate.dt$group, decreasing = T),]
f6.lac.sum <- f6.lac.sum[order(f6.lac.sum$group, decreasing = T),]

names(f6.lac.sum)[c(1,6)] <- c("lab", "group")

f6.lac.plot <- ggplot(data = f6.lac.sum, 
                      aes(x = factor(group, level = c("group 6", "group 5", 
                                                      "group 4", "group 3")),
                          y = emmean, fill = group))+
geom_col(fill = c("brown1", "deepskyblue3", "chartreuse4","purple"), width = 0.5, 
           color = "black") +
  scale_x_discrete(label = glab) +
  scale_y_continuous(limits = c(0,3), expand = expansion(mult = c(0,.10))) +
  theme_pubr() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title.y = element_text(size = 16))+
  geom_errorbar(data = f6.lac.sum, aes(y = emmean, ymin = emmean, ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  geom_point(data = hff6, aes(x = group, y = lactate), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 4),
             size = 1, stroke = 1) +
  annotate("text", x = .6, y = 3, fontface = 3, label = "AHE (AHE > Ref)",
           hjust = 0, , size = 5) +
  annotate("text", x = 1, y = 1.9, label = "A", size = 5) +
  annotate("text", x = 3, y = 1.9, label = "A", size = 5) +
  annotate("text", x = 2, y = 2.2, label = "B", size = 5) +
  coord_cartesian(clip = "off") +
  labs(title = "Lactate",
       subtitle = "Females 6-Months",
       y ="[Lactate] (mg/mL)",
       x = "Diet")
f6.lac.plot

#### lactate 12 Months ####
hf.lactate.f12.lm <- lm(lactate ~ trt * protein, data = hff12)
anova(hf.lactate.f12.lm)
summary(hf.lactate.f12.lm)

ggqqplot(hff12, x = "lactate", color = "group")

f12.lac.res <- residuals(hf.lactate.f12.lm)
ggqqplot(f12.lac.res)

f12.lactate.em <- emmeans(hf.lactate.f12.lm, specs = c("trt", "protein"))
f12.lactate.em

#f12.lactate.dt <- data.table(summary(f12.lactate.em))
#f12.lactate.dt

# Proper/Simple Post-Hoc

emmip(hf.lactate.f12.lm, trt ~ protein)
emmeans(hf.lactate.f12.lm, pairwise ~ protein)

pairs(f12.lactate.em, simple = "protein")
pairs(f12.lactate.em, simple = c("trt", "protein"))

# Graphing labels ()
#f12.lactate.dt$group <- c("group 4", "group 3", "group 6", "group 5")

### Customize Here ###
f12.lac.sum <- hff12 %>% group_by(group) %>%
  summarise(emmean = mean(lactate),
            SD = sd(lactate),
            SE = sd(lactate)/sqrt(length(lactate)),
            n = length(lactate))

f12.lac.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
#f12.lactate.dt <- f12.lactate.dt[order(f12.lactate.dt$group, decreasing = T),]
f12.lac.sum <- f12.lac.sum[order(f12.lac.sum$group, decreasing = T),]

names(f12.lac.sum)[c(1,6)] <- c("lab", "group")

f12.lac.plot <- ggplot(data = f12.lac.sum, 
                       aes(x = factor(group, level = c("group 6", "group 5", 
                                                       "group 4", "group 3")),
                           y = emmean, fill = group))+
geom_col(fill = c("brown1", "deepskyblue3", "chartreuse4","purple"), width = 0.5, 
           color = "black") +
  scale_y_continuous(limits = c(0,3), expand = expansion(mult = c(0,.10))) +
  scale_x_discrete(label = glab) +
  theme_pubr() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24),
        plot.subtitle = element_text(hjust = 0.5))+
  geom_errorbar(data = f12.lac.sum, aes(y = emmean, ymin = emmean, ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  geom_point(data = hff12, aes(x = group, y = lactate), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 4),
             size = 1, stroke = 1) +
  annotate("segment", x = 1, xend = 2, y = 2, yend = 2) +
  annotate("text", x = 1.5, y = 2.2, label = "A", size = 5) +
  annotate("segment", x = 3, xend = 4, y = 2.1, yend = 2.1) +
  annotate("text", x = 3.5, y = 2.3, label = "B", size = 5) +
  annotate("text", x = .6, y = 3, fontface = 3, label = "DPS (Beef > Casein)",
           hjust = 0, size = 5) +
  coord_cartesian(clip = "off") +
  labs(title = "Lactate",
       subtitle = "Females 12-Months",
       y ="Concentration (mg/mL)",
       x = "Diet")
f12.lac.plot

#### lactate 18 Months ####
hf.lactate.f18.lm <- lm(lactate ~ trt * protein, data = hff18)
anova(hf.lactate.f18.lm)
summary(hf.lactate.f18.lm)

ggqqplot(hff18, x = "lactate", color = "group")

f18.lac.res <- residuals(hf.lactate.f18.lm)
ggqqplot(f18.lac.res)

f18.lactate.em <- emmeans(hf.lactate.f18.lm, specs = c("trt", "protein"))
f18.lactate.em

#f18.lactate.dt <- data.table(summary(f18.lactate.em))
#f18.lactate.dt

# Proper/Simple Post-Hoc

emmip(hf.lactate.f18.lm, trt ~ protein)
emmeans(hf.lactate.f18.lm, pairwise ~ protein)

pairs(f18.lactate.em, simple = "protein")
pairs(f18.lactate.em, simple = c("trt", "protein"))

# Graphing labels ()
#f18.lactate.dt$group <- c("group 4", "group 3", "group 6", "group 5")

### Customize Here ###
f18.lac.sum <- hff18 %>% group_by(group) %>%
  summarise(emmean = mean(lactate),
            SD = sd(lactate),
            SE = sd(lactate)/sqrt(length(lactate)),
            n = length(lactate))

f18.lac.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
#f18.lactate.dt <- f18.lactate.dt[order(f18.lactate.dt$group, decreasing = T),]
f18.lac.sum <- f18.lac.sum[order(f18.lac.sum$group, decreasing = T),]

names(f18.lac.sum)[c(1,6)] <- c("lab", "group")

f18.lac.plot <- ggplot(data = f18.lac.sum, 
                       aes(x = factor(group, level = c("group 6", "group 5", 
                                                       "group 4", "group 3")),
                           y = emmean, fill = group))+
geom_col(fill = c("brown1", "deepskyblue3", "chartreuse4","purple"), width = 0.5, 
           color = "black") +
  scale_y_continuous(limits = c(0, 3), expand = expansion(mult = c(0,.10))) +
  scale_x_discrete(label = glab) +
  theme_pubr() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24),
        plot.subtitle = element_text(hjust = 0.5))+
  geom_errorbar(data = f18.lac.sum, aes(y = emmean, ymin = emmean, ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  geom_point(data = hff18, aes(x = group, y = lactate), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 4),
             size = 1, stroke = 1) +
  coord_cartesian(clip = "off") +
  labs(title = "Lactate",
       subtitle = "Females 18-Months",
       y ="Concentration (mg/mL)",
       x = "Diet")
f18.lac.plot

# Step 1: Create the main plot grid for LAC
f.lac.plot <- plot_grid(
  f6.lac.plot + 
    theme(
      plot.subtitle = element_blank(),
      plot.title = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    ),
  f12.lac.plot +
    theme(
      plot.subtitle = element_blank(),
      plot.title = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    ),
  f18.lac.plot +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    ),
  nrow = 1,
  align = "h",
  axis = "tblr",
  rel_widths = c(1, 1, 1)
)

f.lac.plot
lac.label.plot <- ggplot() + 
  theme_void() +
  annotate("text", 
           x = 0, y = 0.5, 
           label = "[Lactate]\n(mg/mL)",  # mg/mL
           angle = 90, 
           size = 5, 
           fontface = "bold",
           hjust = 0.5)

# Step 3: Combine the label and plot using plot_grid again
f.lac.plot.labeled <- plot_grid(
  lac.label.plot, f.lac.plot,
  ncol = 2,
  rel_widths = c(0.07, 1))  # Adjust label space vs plot


f.lac.plot.labeled

## Alanine ##
ggqqplot(hff.met$alanine)

hf.alanine.f6.lm <- lm(alanine ~ trt * protein, data = hff6)
anova(hf.alanine.f6.lm)
summary(hf.alanine.f6.lm)

ggqqplot(hff6, x = "alanine", color = "group")

f6.ala.res <- residuals(hf.alanine.f6.lm)
ggqqplot(f6.ala.res)

f6.alanine.em <- emmeans(hf.alanine.f6.lm, specs = c("trt", "protein"))
f6.alanine.em

#f6.alanine.dt <- data.table(summary(f6.alanine.em))
#f6.alanine.dt

# Main Effects Only #

hf.ala.ahe.f6.lm <- lm(alanine ~ trt, data = hff6)
anova(hf.ala.ahe.f6.lm)
summary(hf.ala.ahe.f6.lm)
f6.ala.ahe.em <- emmeans(hf.ala.ahe.f6.lm, specs = "trt")
emmeans(f6.ala.ahe.em, pairwise ~ trt)
hff6 %>% group_by(trt) %>%
  summarise(emmean = mean(alanine),
            SD = sd(alanine),
            SE = sd(alanine)/sqrt(length(alanine)),
            n = length(alanine))


hf.ala.p.f6.lm <- lm(alanine ~ protein, data = hff6)
anova(hf.ala.p.f6.lm)
summary(hf.ala.p.f6.lm)
f6.ala.p.em <- emmeans(hf.ala.p.f6.lm, specs = "protein")
emmeans(f6.ala.p.em, pairwise ~ protein)
hff6 %>% group_by(protein) %>%
  summarise(emmean = mean(alanine),
            SD = sd(alanine),
            SE = sd(alanine)/sqrt(length(alanine)),
            n = length(alanine))

f6.ala.me <- rbind(anova(hf.ala.ahe.f6.lm),anova(hf.ala.p.f6.lm))

# No Significance #

# Plot
f6.ala.sum <- hff6 %>% group_by(group) %>%
  summarise(emmean = mean(alanine),
            SD = sd(alanine),
            SE = sd(alanine)/sqrt(length(alanine)),
            n = length(alanine))

g3 <- data.frame(group = "group 3", emmean = 0, SD = 0, SE = 0, n = 0)
f6.ala.sum <- rbind(g3, f6.ala.sum)

f6g.lab$y <- f6.ala.sum$emmean + f6.ala.sum$SE + 0.25

f6.ala.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
#f6.alanine.dt <- f6.alanine.dt[order(f6.alanine.dt$group, decreasing = T),]
f6.ala.sum <- f6.ala.sum[order(f6.ala.sum$group, decreasing = T),]

names(f6.ala.sum)[c(1,6)] <- c("lab", "group")

f6.ala.plot <- ggplot(data = f6.ala.sum, 
                      aes(x = factor(group, level = c("group 6", "group 5", 
                                                      "group 4", "group 3")),
                          y = emmean, fill = group))+
geom_col(fill = c("brown1", "deepskyblue3", "chartreuse4","purple"), width = 0.5, 
           color = "black") +
  scale_x_discrete(label = glab) +
  scale_y_continuous(limits = c(0, 900), expand = expansion(mult = c(0,.10))) +
  theme_pubr() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title.y = element_text(size = 16))+
  geom_errorbar(data = f6.ala.sum, aes(y = emmean, ymin = emmean, ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  geom_point(data = hff6, aes(x = group, y = alanine), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 4),
             size = 1, stroke = 1) +
  coord_cartesian(clip = "off") +
  labs(title = "Alanine",
       subtitle = "Females 6-Months",
       y ="[Alanine] (\U00B5g/mL)",
       x = "Diet")
f6.ala.plot

#### alanine 12 Months ####
hf.alanine.f12.lm <- lm(alanine ~ trt * protein, data = hff12)
anova(hf.alanine.f12.lm)
summary(hf.alanine.f12.lm)

ggqqplot(hff12, x = "alanine", color = "group")

f12.ala.res <- residuals(hf.alanine.f12.lm)
ggqqplot(f12.ala.res)

f12.alanine.em <- emmeans(hf.alanine.f12.lm, specs = c("trt", "protein"))
f12.alanine.em

#f12.alanine.dt <- data.table(summary(f12.alanine.em))
#f12.alanine.dt

# Proper/Simple Post-Hoc

emmip(hf.alanine.f12.lm, trt ~ protein)
emmeans(hf.alanine.f12.lm, pairwise ~ protein)

pairs(f12.alanine.em, simple = "protein")
pairs(f12.alanine.em, simple = c("trt", "protein"))
# Graphing labels ()
#f12.alanine.dt$group <- c("group 4", "group 3", "group 6", "group 5")

# Within Group Indications #

### Customize Here ###
f12.ala.sum <- hff12 %>% group_by(group) %>%
  summarise(emmean = mean(alanine),
            SD = sd(alanine),
            SE = sd(alanine)/sqrt(length(alanine)),
            n = length(alanine))

f12.ala.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
#f12.alanine.dt <- f12.alanine.dt[order(f12.alanine.dt$group, decreasing = T),]
f12.ala.sum <- f12.ala.sum[order(f12.ala.sum$group, decreasing = T),]

names(f12.ala.sum)[c(1,6)] <- c("lab", "group")

f12.ala.plot <- ggplot(data = f12.ala.sum, 
                       aes(x = factor(group, level = c("group 6", "group 5", 
                                                       "group 4", "group 3")),
                           y = emmean, fill = group))+
geom_col(fill = c("brown1", "deepskyblue3", "chartreuse4","purple"), width = 0.5, 
           color = "black") +
  scale_y_continuous(limits = c(0, 900), expand = expansion(mult = c(0,.10))) +
  scale_x_discrete(label = glab) +
  theme_pubr() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24),
        plot.subtitle = element_text(hjust = 0.5))+
  geom_errorbar(data = f12.ala.sum, aes(y = emmean, ymin = emmean, ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  geom_point(data = hff12, aes(x = group, y = alanine), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 4),
             size = 1, stroke = 1) +
  coord_cartesian(clip = "off") +
  labs(title = "Alanine",
       subtitle = "Females 12-Months",
       y ="Concentration (\U00B5g/mL)",
       x = "Diet")
f12.ala.plot

#### Alanine 18 Months ####
hf.alanine.f18.lm <- lm(alanine ~ trt * protein, data = hff18)
anova(hf.alanine.f18.lm)
summary(hf.alanine.f18.lm)

ggqqplot(hff18, x = "alanine", color = "group")

f18.ala.res <- residuals(hf.alanine.f18.lm)
ggqqplot(f18.ala.res)

f18.alanine.em <- emmeans(hf.alanine.f18.lm, specs = c("trt", "protein"))
f18.alanine.em

#f18.alanine.dt <- data.table(summary(f18.alanine.em))
#f18.alanine.dt

# Proper/Simple Post-Hoc

emmip(hf.alanine.f18.lm, trt ~ protein)
emmeans(hf.alanine.f18.lm, pairwise ~ protein)

pairs(f18.alanine.em, simple = "protein")
pairs(f18.alanine.em, simple = c("trt", "protein"))

# Graphing labels ()
#f18.alanine.dt$group <- c("group 4", "group 3", "group 6", "group 5")

### Customize Here ###
f18.ala.sum <- hff18 %>% group_by(group) %>%
  summarise(emmean = mean(alanine),
            SD = sd(alanine),
            SE = sd(alanine)/sqrt(length(alanine)),
            n = length(alanine))

f18.ala.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
#f18.alanine.dt <- f18.alanine.dt[order(f18.alanine.dt$group, decreasing = T),]
f18.ala.sum <- f18.ala.sum[order(f18.ala.sum$group, decreasing = T),]

names(f18.ala.sum)[c(1,6)] <- c("lab", "group")

f18.ala.plot <- ggplot(data = f18.ala.sum, 
                       aes(x = factor(group, level = c("group 6", "group 5", 
                                                       "group 4", "group 3")),
                           y = emmean, fill = group))+
geom_col(fill = c("brown1", "deepskyblue3", "chartreuse4","purple"), width = 0.5, 
           color = "black") +
  scale_y_continuous(limits = c(0, 900), expand = expansion(mult = c(0,.10))) +
  scale_x_discrete(label = glab) +
  theme_pubr() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24),
        plot.subtitle = element_text(hjust = 0.5))+
  geom_errorbar(data = f18.ala.sum, aes(y = emmean, ymin = emmean, ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  geom_point(data = hff18, aes(x = group, y = alanine), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 4),
             size = 1, stroke = 1) +
  coord_cartesian(clip = "off") +
  labs(title = "Alanine",
       subtitle = "Females 18-Months",
       y ="Concentration (\U00B5g/mL)",
       x = "Diet")
f18.ala.plot

f.ala.plot <- plot_grid(
  f6.ala.plot + 
    theme(
      plot.subtitle = element_blank(),
      plot.title = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    ),
  f12.ala.plot +
    theme(
      plot.subtitle = element_blank(),
      plot.title = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    ),
  f18.ala.plot +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    ),
  nrow = 1,
  align = "h",
  axis = "tblr",
  rel_widths = c(1, 1, 1))

f.ala.plot
ala.label.plot <- ggplot() + 
  theme_void() +
  annotate("text", 
           x = 0, y = 0.5, 
           label = "[Alanine]\n(\U00B5g/mL)",  # ug/mL
           angle = 90, 
           size = 5, 
           fontface = "bold",
           hjust = 0.5)

# Step 3: Combine the label and plot using plot_grid again
f.ala.plot.labeled <- plot_grid(
  ala.label.plot, f.ala.plot,
  ncol = 2,
  rel_widths = c(0.07, 1))  # Adjust label space vs plot


f.ala.plot.labeled

## Glutamine ##
hf.glutamine.f6.lm <- lm(glutamine ~ trt * protein, data = hff6)
anova(hf.glutamine.f6.lm)
summary(hf.glutamine.f6.lm)

ggqqplot(hff6, x = "glutamine", color = "group")

f6.gln.res <- residuals(hf.glutamine.f6.lm)
ggqqplot(f6.gln.res)

f6.glutamine.em <- emmeans(hf.glutamine.f6.lm, specs = c("trt", "protein"))
f6.glutamine.em

#f6.glutamine.dt <- data.table(summary(f6.glutamine.em))
#f6.glutamine.dt

# Main Effects Only #

hf.gln.ahe.f6.lm <- lm(glutamine ~ trt, data = hff6)
anova(hf.gln.ahe.f6.lm)
summary(hf.gln.ahe.f6.lm)
f6.gln.ahe.em <- emmeans(hf.gln.ahe.f6.lm, specs = "trt")
emmeans(f6.gln.ahe.em, pairwise ~ trt)
hff6 %>% group_by(trt) %>%
  summarise(emmean = mean(glutamine),
            SD = sd(glutamine),
            SE = sd(glutamine)/sqrt(length(glutamine)),
            n = length(glutamine))


hf.gln.p.f6.lm <- lm(glutamine ~ protein, data = hff6)
anova(hf.gln.p.f6.lm)
summary(hf.gln.p.f6.lm)
f6.gln.p.em <- emmeans(hf.gln.p.f6.lm, specs = "protein")
emmeans(f6.gln.p.em, pairwise ~ protein)
hff6 %>% group_by(protein) %>%
  summarise(emmean = mean(glutamine),
            SD = sd(glutamine),
            SE = sd(glutamine)/sqrt(length(glutamine)),
            n = length(glutamine))

f6.gln.me <- rbind(anova(hf.gln.ahe.f6.lm),anova(hf.gln.p.f6.lm))

#Plot
f6.gln.sum <- hff6 %>% group_by(group) %>%
  summarise(emmean = mean(glutamine),
            SD = sd(glutamine),
            SE = sd(glutamine)/sqrt(length(glutamine)),
            n = length(glutamine))

g3 <- data.frame(group = "group 3", emmean = 0, SD = 0, SE = 0, n = 0)
f6.gln.sum <- rbind(g3, f6.gln.sum)

f6g.lab$y <- f6.gln.sum$emmean + f6.gln.sum$SE + 0.25

f6.gln.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
#f6.glutamine.dt <- f6.glutamine.dt[order(f6.glutamine.dt$group, decreasing = T),]
f6.gln.sum <- f6.gln.sum[order(f6.gln.sum$group, decreasing = T),]

names(f6.gln.sum)[c(1,6)] <- c("lab", "group")

f6.gln.plot <- ggplot(data = f6.gln.sum, 
                      aes(x = factor(group, level = c("group 6", "group 5", 
                                                      "group 4", "group 3")),
                          y = emmean, fill = group))+
geom_col(fill = c("brown1", "deepskyblue3", "chartreuse4","purple"), width = 0.5, 
           color = "black") +
  scale_x_discrete(label = glab) +
  scale_y_continuous(limits = c(0, 350), expand = expansion(mult = c(0,.10))) +
  theme_pubr() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title.y = element_text(size = 16))+
  geom_errorbar(data = f6.gln.sum, aes(y = emmean, ymin = emmean, ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  geom_point(data = hff6, aes(x = group, y = glutamine), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 4),
             size = 1, stroke = 1) +
  coord_cartesian(clip = "off") +
  labs(title = "Glutamine",
       subtitle = "Females 6-Months",
       y ="[Glutamine] (\U00B5g/mL)",
       x = "Diet")
f6.gln.plot

#### glutamine 12 Months ####
hf.glutamine.f12.lm <- lm(glutamine ~ trt * protein, data = hff12)
anova(hf.glutamine.f12.lm)
summary(hf.glutamine.f12.lm)

ggqqplot(hff12, x = "glutamine", color = "group")

f12.gln.res <- residuals(hf.glutamine.f12.lm)
ggqqplot(f12.gln.res)

f12.glutamine.em <- emmeans(hf.glutamine.f12.lm, specs = c("trt", "protein"))
f12.glutamine.em

#f12.glutamine.dt <- data.table(summary(f12.glutamine.em))
#f12.glutamine.dt

# Proper/Simple Post-Hoc

emmip(hf.glutamine.f12.lm, trt ~ protein)
emmeans(hf.glutamine.f12.lm, pairwise ~ protein)

pairs(f12.glutamine.em, simple = "protein")
pairs(f12.glutamine.em, simple = c("trt", "protein"))

# Graphing labels ()
#f12.glutamine.dt$group <- c("group 4", "group 3", "group 6", "group 5")

# Within Group Indications #

### Customize Here ###
f12.gln.sum <- hff12 %>% group_by(group) %>%
  summarise(emmean = mean(glutamine),
            SD = sd(glutamine),
            SE = sd(glutamine)/sqrt(length(glutamine)),
            n = length(glutamine))

f12.gln.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
#f12.glutamine.dt <- f12.glutamine.dt[order(f12.glutamine.dt$group, decreasing = T),]
f12.gln.sum <- f12.gln.sum[order(f12.gln.sum$group, decreasing = T),]

names(f12.gln.sum)[c(1,6)] <- c("lab", "group")

f12.gln.plot <- ggplot(data = f12.gln.sum, 
                       aes(x = factor(group, level = c("group 6", "group 5", 
                                                       "group 4", "group 3")),
                           y = emmean, fill = group))+
geom_col(fill = c("brown1", "deepskyblue3", "chartreuse4","purple"), width = 0.5, 
           color = "black") +
  scale_y_continuous(limits = c(0, 350), expand = expansion(mult = c(0,.10))) +
  scale_x_discrete(label = glab) +
  theme_pubr() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24),
        plot.subtitle = element_text(hjust = 0.5))+
  geom_errorbar(data = f12.gln.sum, aes(y = emmean, ymin = emmean, ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  geom_point(data = hff12, aes(x = group, y = glutamine), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 4),
             size = 1, stroke = 1) +
  coord_cartesian(clip = "off") +
  labs(title = "Glutamine",
       subtitle = "Females 12-Months",
       y ="Concentration (\U00B5g/mL)",
       x = "Diet")
f12.gln.plot

#### glutamine 18 Months ####
hf.glutamine.f18.lm <- lm(glutamine ~ trt * protein, data = hff18)
anova(hf.glutamine.f18.lm)
summary(hf.glutamine.f18.lm)

ggqqplot(hff18, x = "glutamine", color = "group")

f18.gln.res <- residuals(hf.glutamine.f18.lm)
ggqqplot(f18.gln.res)

f18.glutamine.em <- emmeans(hf.glutamine.f18.lm, specs = c("trt", "protein"))
f18.glutamine.em

#f18.glutamine.dt <- data.table(summary(f18.glutamine.em))
#f18.glutamine.dt

# Proper/Simple Post-Hoc

emmip(hf.glutamine.f18.lm, trt ~ protein * trt)
emmeans(hf.glutamine.f18.lm, pairwise ~ protein * trt)
emmeans(hf.glutamine.f18.lm, pairwise ~ trt)

pairs(f18.glutamine.em, simple = "trt")
pairs(f18.glutamine.em, simple = c("trt", "protein"))

# Graphing labels ()
#f18.glutamine.dt$group <- c("group 4", "group 3", "group 6", "group 5")

# Within Group Indications #

f18gln.lab <- paste("group", c(3:6))
f18gln.lab <- data.frame(f18gln.lab)
f18gln.lab$lab <- c("A", "AB", "B", "A")
names(f18gln.lab)[1] <- "group1"
f18gln.lab$group2 <- NA
f18gln.lab$group <- f18gln.lab$group1

### Customize Here ###
f18.gln.sum <- hff18 %>% group_by(group) %>%
  summarise(emmean = mean(glutamine),
            SD = sd(glutamine),
            SE = sd(glutamine)/sqrt(length(glutamine)),
            n = length(glutamine))

f18gln.lab$y <- f18.gln.sum$emmean + f18.gln.sum$SD + 40
f18gln.lab$y[4] <- 280

f18.gln.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
#f18.glutamine.dt <- f18.glutamine.dt[order(f18.glutamine.dt$group, decreasing = T),]
f18.gln.sum <- f18.gln.sum[order(f18.gln.sum$group, decreasing = T),]

names(f18.gln.sum)[c(1,6)] <- c("lab", "group")

f18.gln.plot <- ggplot(data = f18.gln.sum, 
                       aes(x = factor(group, level = c("group 6", "group 5", 
                                                       "group 4", "group 3")),
                           y = emmean, fill = group))+
geom_col(fill = c("brown1", "deepskyblue3", "chartreuse4","purple"), width = 0.5, 
           color = "black") +
  scale_y_continuous(limits = c(0, 350), expand = expansion(mult = c(0,.10))) +
  scale_x_discrete(label = glab) +
  theme_pubr() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24),
        plot.subtitle = element_text(hjust = 0.5))+
  geom_errorbar(data = f18.gln.sum, aes(y = emmean, ymin = emmean, ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  geom_point(data = hff18, aes(x = group, y = glutamine), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 5),
             size = 1, stroke = 1) +
  stat_pvalue_manual(data = f18gln.lab,
                     label = "lab",
                     x = "group",
                     y.position = "y",
                     label.size = 4) +
  annotate("text", x = .6, y = 350, fontface = 3, label = "DPS * AHE", hjust = 0, , size = 5) +
  coord_cartesian(clip = "off") +
  labs(title = "Glutamine",
       subtitle = "Females 18-Months",
       y ="Concentration (\U00B5g/mL)",
       x = "Diet")
f18.gln.plot

f.gln.plot <- plot_grid(
  f6.gln.plot + 
    theme(plot.subtitle = element_blank(),
          plot.title = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank()),
  f12.gln.plot +
    theme(plot.subtitle = element_blank(),
          plot.title = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank()),
  f18.gln.plot +
    theme(plot.title = element_blank(),
          plot.subtitle = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank()),
  nrow = 1,
  align = "v",
  axis = "tblr",
  rel_widths = c(1, 1, 1)
)

f.gln.plot
gln.label.plot <- ggplot() + 
  theme_void() +
  annotate("text", 
           x = 0, y = 0.5, 
           label = "[Glutamine]\n(\U00B5g/mL)",  # ug/mL
           angle = 90, 
           size = 5, 
           fontface = "bold",
           hjust = 0.5)

# Step 3: Combine the label and plot using plot_grid again
f.gln.plot.labeled <- plot_grid(
  gln.label.plot, f.gln.plot,
  ncol = 2,
  rel_widths = c(0.07, 1))  # Adjust label space vs plot


f.gln.plot.labeled

# Carnitine #
hf.carnitine.f6.lm <- lm(carnitine ~ trt * protein, data = hff6)
anova(hf.carnitine.f6.lm)
summary(hf.carnitine.f6.lm)

ggqqplot(hff6, x = "carnitine", color = "group")

f6.crn.res <- residuals(hf.carnitine.f6.lm)
ggqqplot(f6.crn.res)

f6.carnitine.em <- emmeans(hf.carnitine.f6.lm, specs = c("trt", "protein"))
f6.carnitine.em

#f6.carnitine.dt <- data.table(summary(f6.carnitine.em))
#f6.carnitine.dt


pt <- emmeans(hf.carnitine.f6.lm, pairwise ~ protein)
pairs(f6.carnitine.em, simple = "protein")

# Contrast "pairwise comparisions"

hf.crn.ahe.f6.lm <- lm(carnitine ~ trt, data = hff6)
anova(hf.crn.ahe.f6.lm)
summary(hf.crn.ahe.f6.lm)
f6.crn.ahe.em <- emmeans(hf.crn.ahe.f6.lm, specs = "trt")
emmeans(f6.crn.ahe.em, pairwise ~ trt)
hff6 %>% group_by(trt) %>%
  summarise(emmean = mean(carnitine),
            SD = sd(carnitine),
            SE = sd(carnitine)/sqrt(length(carnitine)),
            n = length(carnitine))


hf.crn.p.f6.lm <- lm(carnitine ~ protein, data = hff6)
anova(hf.crn.p.f6.lm)
summary(hf.crn.p.f6.lm)
f6.crn.p.em <- emmeans(hf.crn.p.f6.lm, specs = "protein")
emmeans(f6.crn.p.em, pairwise ~ protein)
hff6 %>% group_by(protein) %>%
  summarise(emmean = mean(carnitine),
            SD = sd(carnitine),
            SE = sd(carnitine)/sqrt(length(carnitine)),
            n = length(carnitine))

f6.crn.me <- rbind(anova(hf.crn.ahe.f6.lm),anova(hf.crn.p.f6.lm))

# Identify Significance #
f6crn.lab <- paste("group", c(3:6))
f6crn.lab <- data.frame(f6crn.lab)
f6crn.lab$lab <- c("", "B", "A", "A")
names(f6crn.lab)[1] <- "group1"
f6crn.lab$group2 <- NA
f6crn.lab$group <- f6crn.lab$group1

#Plot
f6.crn.sum <- hff6 %>% group_by(group) %>%
  summarise(emmean = mean(carnitine),
            SD = sd(carnitine),
            SE = SD/sqrt(length(carnitine)),
            n = length(carnitine))

g3 <- data.frame(group = "group 3", emmean = 0, SD = 0, SE = 0, n = 0)
f6.crn.sum <- rbind(g3, f6.crn.sum)

f6crn.lab$y <- f6.crn.sum$emmean + f6.crn.sum$SD + 15

f6.crn.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
#f6.carnitine.dt <- f6.carnitine.dt[order(f6.carnitine.dt$group, decreasing = T),]
f6.crn.sum <- f6.crn.sum[order(f6.crn.sum$group, decreasing = T),]

names(f6.crn.sum)[c(1,6)] <- c("lab", "group")

f6.crn.plot <- ggplot(data = f6.crn.sum, 
                      aes(x = factor(group, level = c("group 6", "group 5", 
                                                      "group 4", "group 3")),
                          y = emmean, fill = group))+
geom_col(fill = c("brown1", "deepskyblue3", "chartreuse4","purple"), width = 0.5, 
           color = "black") +
  scale_x_discrete(label = glab) +
  scale_y_continuous(limits = c(0, 130), expand = expansion(mult = c(0,.10))) +
  theme_pubr() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title.y = element_text(size = 16))+
  geom_errorbar(data = f6.crn.sum, aes(y = emmean, ymin = emmean, ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  geom_point(data = hff6, aes(x = group, y = carnitine), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 4),
             size = 1, stroke = 1) +
  annotate("segment", x = 1, xend = 2, y = 90, yend = 90) +
  annotate("text", x = 1.5, y = 99, label = "A", size = 5) +
  annotate("text", x = 3, y = 105, label = "B", size = 5) +
  annotate("text", x = .6, y = 130, fontface = 3, label = "DPS (Beef > Casein)",
           hjust = 0, , size = 5) +
  coord_cartesian(clip = "off") +
  labs(title = "Carnitine",
       subtitle = "Females 6-Months",
       y ="[Carnitine] (\U00B5g/mL)",
       x = "Diet")
f6.crn.plot

#### carnitine 12 Months ####
hf.carnitine.f12.lm <- lm(carnitine ~ trt * protein, data = hff12)
anova(hf.carnitine.f12.lm)
summary(hf.carnitine.f12.lm)

ggqqplot(hff12, x = "carnitine", color = "group")

f12.crn.res <- residuals(hf.carnitine.f12.lm)
ggqqplot(f12.crn.res)

f12.carnitine.em <- emmeans(hf.carnitine.f12.lm, specs = c("trt", "protein"))
f12.carnitine.em

#f12.carnitine.dt <- data.table(summary(f12.carnitine.em))
#f12.carnitine.dt

# Proper/Simple Post-Hoc

emmip(hf.carnitine.f12.lm, trt ~ protein)
emmeans(hf.carnitine.f12.lm, pairwise ~ protein)

pairs(f12.carnitine.em, simple = "protein")
pairs(f12.carnitine.em, simple = c("trt", "protein"))

# Graphing labels ()
#f12.carnitine.dt$group <- c("group 4", "group 3", "group 6", "group 5")

# Within Group Indications #

### Customize Here ###
f12.crn.sum <- hff12 %>% group_by(group) %>%
  summarise(emmean = mean(carnitine),
            SD = sd(carnitine),
            SE = sd(carnitine)/sqrt(length(carnitine)),
            n = length(carnitine))

f12.crn.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
#f12.carnitine.dt <- f12.carnitine.dt[order(f12.carnitine.dt$group, decreasing = T),]
f12.crn.sum <- f12.crn.sum[order(f12.crn.sum$group, decreasing = T),]

names(f12.crn.sum)[c(1,6)] <- c("lab", "group")

f12.crn.plot <- ggplot(data = f12.crn.sum, 
                       aes(x = factor(group, level = c("group 6", "group 5", 
                                                       "group 4", "group 3")),
                           y = emmean, fill = group))+
geom_col(fill = c("brown1", "deepskyblue3", "chartreuse4","purple"), width = 0.5, 
           color = "black") +
  scale_y_continuous(limits = c(0, 130), expand = expansion(mult = c(0,.10))) +
  scale_x_discrete(label = glab) +
  theme_pubr() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24),
        plot.subtitle = element_text(hjust = 0.5))+
  geom_errorbar(data = f12.crn.sum, aes(y = emmean, ymin = emmean, ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  geom_point(data = hff12, aes(x = group, y = carnitine), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 4),
             size = 1, stroke = 1) +
  annotate("segment", x = "group 4", xend = "group 3", y = 105, yend = 105) +
  annotate("text", x = 3.5, y = 114, label = "B", size = 5) +
  annotate("segment", x = "group 6", xend = "group 5", y = 81, yend = 81) +
  annotate("text", x = 1.5, y = 90, label = "A", size = 5) +
  annotate("text", x = .6, y = 130, fontface = 3, label = "DPS (Beef > Casein)",
           hjust = 0, size = 5) +
  coord_cartesian(clip = "off") +
  labs(title = "Carnitine",
       subtitle = "Females 12-Months",
       y ="Concentration (\U00B5g/mL)",
       x = "Diet")
f12.crn.plot

#### carnitine 18 Months ####
hf.carnitine.f18.lm <- lm(carnitine ~ trt * protein, data = hff18)
anova(hf.carnitine.f18.lm)
summary(hf.carnitine.f18.lm)

ggqqplot(hff18, x = "carnitine", color = "group")

f18.crn.res <- residuals(hf.carnitine.f18.lm)
ggqqplot(f18.crn.res)

f18.carnitine.em <- emmeans(hf.carnitine.f18.lm, specs = c("trt", "protein"))
f18.carnitine.em

#f18.carnitine.dt <- data.table(summary(f18.carnitine.em))
#f18.carnitine.dt

# Proper/Simple Post-Hoc

emmip(hf.carnitine.f18.lm, trt ~ protein)
emmeans(hf.carnitine.f18.lm, pairwise ~ protein)

pairs(f18.carnitine.em, simple = "protein")
pairs(f18.carnitine.em, simple = c("trt", "protein"))

# Graphing labels ()
#f18.carnitine.dt$group <- c("group 4", "group 3", "group 6", "group 5")

# Within Group Indications #

f18crn.lab <- paste("group", c(3:6))
f18crn.lab <- data.frame(f18crn.lab)
f18crn.lab$lab <- c("A", "AB", "B", "A")
names(f18crn.lab)[1] <- "group1"
f18crn.lab$group2 <- NA
f18crn.lab$group <- f18crn.lab$group1

### Customize Here ###
f18.crn.sum <- hff18 %>% group_by(group) %>%
  summarise(emmean = mean(carnitine),
            SD = sd(carnitine),
            SE = sd(carnitine)/sqrt(length(carnitine)),
            n = length(carnitine))

f18crn.lab$y <- f18.crn.sum$emmean + f18.crn.sum$SD + 20

f18.crn.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
#f18.carnitine.dt <- f18.carnitine.dt[order(f18.carnitine.dt$group, decreasing = T),]
f18.crn.sum <- f18.crn.sum[order(f18.crn.sum$group, decreasing = T),]

names(f18.crn.sum)[c(1,6)] <- c("lab", "group")

f18.crn.plot <- ggplot(data = f18.crn.sum, 
                       aes(x = factor(group, level = c("group 6", "group 5", 
                                                       "group 4", "group 3")),
                           y = emmean, fill = group))+
geom_col(fill = c("brown1", "deepskyblue3", "chartreuse4","purple"), width = 0.5, 
           color = "black") +
  scale_y_continuous(limits = c(0, 130), expand = expansion(mult = c(0,.10))) +
  scale_x_discrete(label = glab) +
  theme_pubr() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24),
        plot.subtitle = element_text(hjust = 0.5))+
  geom_errorbar(data = f18.crn.sum, aes(y = emmean, ymin = emmean, ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  geom_point(data = hff18, aes(x = group, y = carnitine), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 4),
             size = 1, stroke = 1) +
  annotate("segment", x = "group 4", xend = "group 3", y = 85, yend = 85) +
  annotate("text", x = 3.5, y = 94, label = "B", size = 5) +
  annotate("segment", x = "group 6", xend = "group 5", y = 70, yend = 70) +
  annotate("text", x = 1.5, y = 79, label = "A", size = 5) +
  annotate("text", x = .6, y = 130, fontface = 3, label = "DPS (Beef > Casein)",
           hjust = 0, size = 5) +
  coord_cartesian(clip = "off") +
  labs(title = "Carnitine",
       subtitle = "Females 18-Months",
       y ="Concentration (\U00B5g/mL)",
       x = "Diet")
f18.crn.plot

f.crn.plot <- plot_grid(
  f6.crn.plot + 
    theme(plot.subtitle = element_blank(),
          plot.title = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()),
  f12.crn.plot +
    theme(plot.subtitle = element_blank(),
          plot.title = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank()),
  f18.crn.plot +
    theme(plot.title = element_blank(),
          plot.subtitle = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank()),
  nrow = 1,
  align = "v",
  axis = "tblr",
  rel_widths = c(1, 1, 1)
)

f.crn.plot
crn.label.plot <- ggplot() + 
  theme_void() +
  annotate("text", 
           x = 0, y = 0.5, 
           label = "[Carnitine]\n(\U00B5g/mL)",  # ug/mL
           angle = 90, 
           size = 5, 
           fontface = "bold",
           hjust = 0.5)

# Step 3: Combine the label and plot using plot_grid again
f.crn.plot.labeled <- plot_grid(
  crn.label.plot, f.crn.plot,
  ncol = 2,
  rel_widths = c(0.07, 1))  # Adjust label space vs plot


f.crn.plot.labeled

plot_grid(f.glu.plot.labeled, f.lac.plot.labeled, 
          f.ala.plot.labeled, f.gln.plot.labeled, f.crn.plot.labeled, nrow = 5)

#########

clean_plot <- function(p) {
  p + theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 10)
  )
}

first_plot <- function(p) {
  p + theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 10)
  )
}

first_bottom_plot <- function(p) {
  p + theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 10))
}

bottom_plot <- function(p) {
  p + theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 10))
}

#Make y label function
make_y_label <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label, angle = 90,
             size = 7, fontface = "bold") +
    theme_void()
}

#Make Col label function
make_col_label <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label,
             size = 7, hjust = 0.5) +
    theme_void()
}

#### USE THIS ONE ###

ylab_col <- plot_grid(
  make_y_label("[Glucose]\n(mg/mL)"),
  make_y_label("[Lactate]\n(mg/mL)"),
  make_y_label("[Alanine]\n(µg/mL)"),
  make_y_label("[Glutamine]\n(µg/mL)"),
  make_y_label("[Carnitine]\n(µg/mL)"),
  ncol = 1,
  rel_heights = rep(1, 5)
)

# 6 months column: first_plot for first 4, first_bottom_plot for carnitine (bottom)
col_6 <- plot_grid(
  first_plot(f6.glu.plot),
  first_plot(f6.lac.plot),
  first_plot(f6.ala.plot),
  first_plot(f6.gln.plot),
  first_bottom_plot(f6.crn.plot),
  ncol = 1,
  align = "v",
  axis = "lr"
)

# 12 months column: clean_plot for first 4, bottom_plot for carnitine
col_12 <- plot_grid(
  clean_plot(f12.glu.plot),
  clean_plot(f12.lac.plot),
  clean_plot(f12.ala.plot),
  clean_plot(f12.gln.plot),
  bottom_plot(f12.crn.plot),
  ncol = 1,
  align = "v",
  axis = "lr"
)

# 18 months column: clean_plot for first 4, bottom_plot for carnitine
col_18 <- plot_grid(
  clean_plot(f18.glu.plot),
  clean_plot(f18.lac.plot),
  clean_plot(f18.ala.plot),
  clean_plot(f18.gln.plot),
  bottom_plot(f18.crn.plot),
  ncol = 1,
  align = "v",
  axis = "lr"
)

# Combine everything into a 5 x 4 layout: y-axis labels + 3 timepoints
final_plot <- plot_grid(
  ylab_col,
  col_6,
  col_12,
  col_18,
  ncol = 4,
  rel_widths = c(0.3, 1, 1, 1),
  align = "h",
  axis = "tb"
)

final_plot



col_titles <- plot_grid(
  NULL,
  ggdraw() + draw_label("6 Months", fontface = "bold", hjust = 0.5),
  ggdraw() + draw_label("12 Months", fontface = "bold", hjust = 0.5),
  ggdraw() + draw_label("18 Months", fontface = "bold", hjust = 0.5),
  ncol = 4,
  rel_widths = c(0.3, 1, 1, 1)
)

# Final with title row
plot_with_titles <- plot_grid(
  col_titles,
  final_plot,
  ncol = 1,
  rel_heights = c(0.05, 1)
)

main_title <- ggdraw() + draw_label("Metabolite Concentrations from Female Livers", fontface = "bold", size = 20)

f.full_figure <- plot_grid(
  main_title,
  plot_with_titles,
  ncol = 1,
  rel_heights = c(0.07, 1)
)

f.full_figure

ggsave("metabolites from females.pdf", f.full_figure, width = 10, height = 12)

########
aovf.list <- list(f6.glu = f6.glu.me,
                 f6.lac = f6.lac.me,
                 f6.ala = f6.ala.me,
                 f6.gln = f6.gln.me,
                 f6.crn = f6.crn.me,
                 f12.glu = anova(hf.glucose.f12.lm), 
                 f12.lac = anova(hf.lactate.f12.lm),
                 f12.ala = anova(hf.alanine.f12.lm),
                 f12.gln = anova(hf.glutamine.f12.lm),
                 f12.crn = anova(hf.carnitine.f12.lm),
                 f18.glu = anova(hf.glucose.f18.lm),
                 f18.lac = anova(hf.lactate.f18.lm),
                 f18.ala = anova(hf.alanine.f18.lm),
                 f18.gln = anova(hf.glutamine.f18.lm),
                 f18.crn = anova(hf.carnitine.f18.lm))

anovaf.combined <- do.call(rbind, lapply(names(aovf.list), function(name) {
  dff <- as.data.frame(aovf.list[[name]])
  dff$model <- name  # Add identifier
  dff$term <- rownames(aovf.list[[name]])  # Preserve rownames as a column
  rownames(dff) <- NULL  # Optional: clean up row names
  dff
}))

anovaf.combined



write.csv(anovaf.combined, "./Data/f.anova.combined.csv")

write.csv(anova(hf.glucose.f.lm), "./Data/f.glucose.anova.csv")
write.csv(anova(hf.lactate.f.lm), "./Data/f.lactate.anova.csv")
write.csv(anova(hf.ala.f.lm), "./Data/f.ala.anova.csv")
write.csv(anova(hf.gln.f.lm), "./Data/f.gln.anova.csv")
write.csv(anova(hf.crn.f.lm), "./Data/f.crn.anova.csv")
