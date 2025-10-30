library(here) # here makes a project transportable
library(janitor) # clean_names
library(readxl) # read excel, duh!
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
library(cowplot) # combine plots
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

library(cowplot)
library(gridExtra)
library(data.table) # magical data frames

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

m.met <- met[which(met$sex == "M"),] # Males
hfm.met <- m.met[which(m.met$group %in% c("G3", "G4", "G5", "G6")),] #HF only

nrow(hfm.met)

hfm.met$group[which(hfm.met$group == "G3")] <- "group 3"
hfm.met$group[which(hfm.met$group == "G4")] <- "group 4"
hfm.met$group[which(hfm.met$group == "G5")] <- "group 5"
hfm.met$group[which(hfm.met$group == "G6")] <- "group 6"



## Males ##
#Age Separation#
hfm6 <- hfm.met[which(hfm.met$age == "H1"),]
hfm12 <- hfm.met[which(hfm.met$age == "H2"),]
hfm18 <- hfm.met[which(hfm.met$age == "H3"),]

#Glucose First#
hf.glucose.m6.lm <- lm(glucose ~ trt * protein, data = hfm6)
anova(hf.glucose.m6.lm)
summary(hf.glucose.m6.lm)

ggqqplot(hfm6, x = "glucose", color = "group")

m6.glu.res <- residuals(hf.glucose.m6.lm)
ggqqplot(m6.glu.res)


m6.glucose.em <- emmeans(hf.glucose.m6.lm, specs = c("trt", "protein"))
m6.glucose.em

# Proper/Simple Post-Hoc

pairs(m6.glucose.em, simple = c("protein", "trt"))

# Graphing labels ()

### Customize Here ###
m6g.lab <- paste("group", c(3:6))
m6g.lab <- data.frame(m6g.lab)
m6g.lab$lab <- c("B", "A", "A", "A")
names(m6g.lab)[1] <- "group1"
m6g.lab$group2 <- NA
m6g.lab$group <- m6g.lab$group1


# Customize for each set of comp
alab <- c(`H1` = "6 Months", `H2` = "12 Months", 
          `H3` = "18 Months")

glab <- c(`group 3` = "HFBN",
          `group 4` = "HFB",
          `group 5` = "HFCN",
          `group 6` = "HFC")

#Plot
g6.glu.sum <- hfm6 %>% group_by(group) %>%
  summarise(emmean = mean(glucose),
            SD = sd(glucose),
            SE = sd(glucose)/sqrt(length(glucose)),
            n = length(glucose))

m6g.lab$y <- g6.glu.sum$emmean + g6.glu.sum$SD + 5

g6.glu.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
g6.glu.sum <- g6.glu.sum[order(g6.glu.sum$group, decreasing = T),]

names(g6.glu.sum)[c(1,6)] <- c("lab", "group")

m6.glu.plot <- ggplot(data = g6.glu.sum, 
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
  geom_errorbar(data = g6.glu.sum, aes(y = emmean, ymin = emmean, ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  annotate("text", x = .6, y = 40, fontface = 3, label = "DPS * AHE", hjust = 0, size = 5) +
  stat_pvalue_manual(data = m6g.lab,
                     label = "lab",
                     x = "group",
                     y.position = "y",
                     label.size = 5) + 
  geom_point(data = hfm6, aes(x = group, y = glucose), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 4),
             size = 1, stroke = 1) +
  labs(title = "Glucose",
       subtitle = "Males 6-Months",
       y ="mg/mL",
       x = "Diet")
m6.glu.plot

#### Glucose 12 Months ####
hf.glucose.m12.lm <- lm(glucose ~ trt * protein, data = hfm12)
anova(hf.glucose.m12.lm)
summary(hf.glucose.m12.lm)

ggqqplot(hfm12, x = "glucose", color = "group")

m12.glu.res <- residuals(hf.glucose.m12.lm)
ggqqplot(m12.glu.res)

m12.glucose.em <- emmeans(hf.glucose.m12.lm, specs = c("trt", "protein"))
m12.glucose.em

# Proper/Simple Post-Hoc

pairs(m12.glucose.em, simple = c("trt", "protein"))

# Graphing labels ()

### Customize Here ###
m12g.lab <- paste("group", c(3:6))
m12g.lab <- data.frame(m12g.lab)
m12g.lab$lab <- c("C", "B", "AB", "A")
names(m12g.lab)[1] <- "group1"
m12g.lab$group2 <- NA
m12g.lab$group <- m12g.lab$group1

#Plot
g12.glu.sum <- hfm12 %>% group_by(group) %>%
  summarise(emmean = mean(glucose),
            SD = sd(glucose),
            SE = sd(glucose)/sqrt(length(glucose)),
            n = length(glucose))

m12g.lab$y <- g12.glu.sum$emmean + g12.glu.sum$SD + 1

g12.glu.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
g12.glu.sum <- g12.glu.sum[order(g12.glu.sum$group, decreasing = T),]

names(g12.glu.sum)[c(1,6)] <- c("lab", "group")

m12.glu.plot <- ggplot(data = g12.glu.sum, 
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
  geom_errorbar(data = g12.glu.sum, aes(y = emmean, ymin = emmean, ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  stat_pvalue_manual(data = m12g.lab,
                     label = "lab",
                     x = "group",
                     y.position = "y",
                     label.size = 5) + 
  geom_point(data = hfm12, aes(x = group, y = glucose), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 4),
             size = 1, stroke = 1) +
  annotate("text", x = .6, y = 40, fontface = 3, label = "DPS * AHE", hjust = 0, size = 5) +
  labs(title = "Glucose",
       subtitle = "Males 12-Months",
       y ="mg/mL",
       x = "Diet")
m12.glu.plot

#### Glucose 18 Months ####
hf.glucose.m18.lm <- lm(glucose ~ trt * protein, data = hfm18)
anova(hf.glucose.m18.lm)
summary(hf.glucose.m18.lm)

ggqqplot(hfm18, x = "glucose", color = "group")

m18.glu.res <- residuals(hf.glucose.m18.lm)
ggqqplot(m18.glu.res)

m18.glucose.em <- emmeans(hf.glucose.m18.lm, specs = c("trt", "protein"))
m18.glucose.em

# Proper/Simple Post-Hoc

pairs(m18.glucose.em, simple = c("trt", "protein"))

# Graphing labels ()

### Customize Here ###
m18g.lab <- paste("group", c(3:6))
m18g.lab <- data.frame(m18g.lab)
m18g.lab$lab <- c("C", "AB", "B", "AC")
names(m18g.lab)[1] <- "group1"
m18g.lab$group2 <- NA
m18g.lab$group <- m18g.lab$group1

#Plot
g18.glu.sum <- hfm18 %>% group_by(group) %>%
  summarise(emmean = mean(glucose),
            SD = sd(glucose),
            SE = sd(glucose)/sqrt(length(glucose)),
            n = length(glucose))

m18g.lab$y <- g18.glu.sum$emmean + g18.glu.sum$SD + 2

g18.glu.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
g18.glu.sum <- g18.glu.sum[order(g18.glu.sum$group, decreasing = T),]

names(g18.glu.sum)[c(1,6)] <- c("lab", "group")


m18.glu.plot <- ggplot(data = g18.glu.sum, 
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
        plot.subtitle = element_text(hjust = 0.5)) +
  geom_errorbar(data = g18.glu.sum, aes(y = emmean, ymin = emmean, ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  stat_pvalue_manual(data = m18g.lab,
                     label = "lab",
                     x = "group",
                     y.position = "y",
                     label.size = 5) +
  geom_point(data = hfm18, aes(x = group, y = glucose), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 8),
             size = 1, stroke = 1) +
  annotate("text", x = .6, y = 40, fontface = 3, label = "DPS * AHE", hjust = 0, size = 5) +
  labs(title = "Glucose",
       subtitle = "Males 18-Months",
       y ="mg/mL",
       x = "Diet")
m18.glu.plot

m.glu.plot <- ggarrange(m6.glu.plot +
                          theme(plot.title = element_blank(),
                                axis.title.x = element_blank()),
                        m12.glu.plot + 
                          theme(plot.title = element_blank(),
                                axis.text.y = element_blank(),
                                axis.ticks.y = element_blank(),
                                axis.title.y = element_blank(),
                                axis.title.x = element_blank()),
                        m18.glu.plot +
                          theme(plot.title = element_blank(),
                                axis.text.y = element_blank(),
                                axis.ticks.y = element_blank(),
                                axis.title.y = element_blank(),
                                axis.title.x = element_blank()),
                        nrow = 1)

m.glu.plot <- annotate_figure(m.glu.plot,
                              top = text_grob("Male Glucose", size = 18),
                              bottom = text_grob("Diet", size = 16))
m.glu.plot

#Lactate Next#
hf.lactate.m6.lm <- lm(lactate ~ trt * protein, data = hfm6)
anova(hf.lactate.m6.lm)
summary(hf.lactate.m6.lm)

ggqqplot(hfm6, x = "lactate", color = "group")

m6.lac.res <- residuals(hf.lactate.m6.lm)
ggqqplot(m6.lac.res)

m6.lactate.em <- emmeans(hf.lactate.m6.lm, specs = c("trt", "protein"))
m6.lactate.em

# Proper/Simple Post-Hoc

pt <- emmeans(hf.lactate.m6.lm, pairwise ~ protein)

pairs(m6.lactate.em, simple = "protein")

hfm6 %>% 
  emmeans_test(
    lactate ~ protein, model = hf.lactate.m6.lm
  )

plot(pt, comparisons = T)

#Plot
m6.lac.sum <- hfm6 %>% group_by(group) %>%
  summarise(emmean = mean(lactate),
            SD = sd(lactate),
            SE = sd(lactate)/sqrt(length(lactate)),
            n = length(lactate))

m6.lac.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
m6.lac.sum <- m6.lac.sum[order(m6.lac.sum$group, decreasing = T),]

names(m6.lac.sum)[c(1,6)] <- c("lab", "group")

m6.lac.plot <- ggplot(data = m6.lac.sum, 
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
  geom_errorbar(data = m6.lac.sum, aes(y = emmean, ymin = emmean, ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  annotate("segment", x = "group 4", xend = "group 3", y = 2.3, yend = 2.3) +
  annotate("text", x = 3.5, y = 2.5, label = "B", size = 5) +
  annotate("segment", x = "group 6", xend = "group 5", y = 2, yend = 2) +
  annotate("text", x = 1.5, y = 2.2, label = "A", size = 5) +
  annotate("text", x = .6, y = 3, fontface = 3, label = "DPS (Beef > Casein)", hjust = 0, size = 5) +
  geom_point(data = hfm6, aes(x = group, y = lactate), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 4),
             size = 1, stroke = 1) +
  labs(title = "Lactate",
       subtitle = "Males 6-Months",
       y ="mg/mL",
       x = "Diet")
m6.lac.plot

#### lactate 12 Months ####
hf.lactate.m12.lm <- lm(lactate ~ trt * protein, data = hfm12)
anova(hf.lactate.m12.lm)
summary(hf.lactate.m12.lm)

ggqqplot(hfm6, x = "lactate", color = "group")

m12.lac.res <- residuals(hf.lactate.m12.lm)
ggqqplot(m12.lac.res)

m12.lactate.em <- emmeans(hf.lactate.m12.lm, specs = c("trt", "protein"))
m12.lactate.em


# Proper/Simple Post-Hoc
# Graphing labels ()

#Plot
m12.lac.sum <- hfm12 %>% group_by(group) %>%
summarise(emmean = mean(lactate),
          SD = sd(lactate),
          SE = sd(lactate)/sqrt(length(lactate)),
          n = length(lactate))

m12.lac.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
m12.lac.sum <- m12.lac.sum[order(m12.lac.sum$group, decreasing = T),]

names(m12.lac.sum)[c(1,6)] <- c("lab", "group")

m12.lac.plot <- ggplot(data = m12.lac.sum, 
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
  geom_errorbar(data = m12.lac.sum, 
                aes(y = emmean, ymin = emmean, 
                    ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  geom_point(data = hfm12, aes(x = group, y = lactate), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 4),
             size = 1, stroke = 1) +
  labs(title = "Lactate",
       subtitle = "Males 12-Months",
       y ="mg/mL",
       x = "Diet")
m12.lac.plot

#### lactate 18 Months ####
hf.lactate.m18.lm <- lm(lactate ~ trt * protein, data = hfm18)
anova(hf.lactate.m18.lm)
summary(hf.lactate.m18.lm)

ggqqplot(hfm18, x = "lactate", color = "group")

m18.lac.res <- residuals(hf.lactate.m18.lm)
ggqqplot(m18.lac.res)

m18.lactate.em <- emmeans(hf.lactate.m18.lm, specs = c("trt", "protein"))
m18.lactate.em

# Proper/Simple Post-Hoc
pairs(m18.lactate.em, simple = c("trt", "protein"))

# Graphing labels ()

### Customize Here ###
m18lac.lab <- paste("group", c(3:6))
m18lac.lab <- data.frame(m18lac.lab)
m18lac.lab$lab <- c("A", "AB", "B", "A")
names(m18lac.lab)[1] <- "group1"
m18lac.lab$group2 <- NA
m18lac.lab$group <- m18lac.lab$group1

# Customize for each set of comp

#Plot
m18.lac.sum <- hfm18 %>% group_by(group) %>%
  summarise(emmean = mean(lactate),
            SD = sd(lactate),
            SE = sd(lactate)/sqrt(length(lactate)),
            n = length(lactate))

m18lac.lab$y <- m18.lac.sum$emmean + m18.lac.sum$SD + .25

m18.lac.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
m18.lac.sum <- m18.lac.sum[order(m18.lac.sum$group, decreasing = T),]

names(m18.lac.sum)[c(1,6)] <- c("lab", "group")

#m18lac.lab[4,5] <- 1.8

m18.lac.plot <- ggplot(data = m18.lac.sum, 
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
        plot.subtitle = element_text(hjust = 0.5)) +
  geom_errorbar(data = m18.lac.sum, aes(y = emmean, ymin = emmean, ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  stat_pvalue_manual(data = m18lac.lab,
                     label = "lab",
                     x = "group",
                     y.position = "y",
                     label.size = 5) +
  geom_point(data = hfm18, aes(x = group, y = lactate), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 8),
             size = 1, stroke = 1) +
  annotate("text", x = .6, y = 3, fontface = 3, label = "DPS * AHE", hjust = 0, size = 5) +
  labs(title = "Lactate",
       subtitle = "Males 18-Months",
       y ="mg/mL",
       x = "Diet")
m18.lac.plot

m.lac.plot <- ggarrange(m6.lac.plot +
                          theme(plot.title = element_blank(),
                                axis.title.x = element_blank()),
                        m12.lac.plot + 
                          theme(plot.title = element_blank(),
                                axis.text.y = element_blank(),
                                axis.ticks.y = element_blank(),
                                axis.title.y = element_blank(),
                                axis.title.x = element_blank()),
                        m18.lac.plot +
                          theme(plot.title = element_blank(),
                                axis.text.y = element_blank(),
                                axis.ticks.y = element_blank(),
                                axis.title.y = element_blank(),
                                axis.title.x = element_blank()),
                        nrow = 1)

m.lac.plot <- annotate_figure(m.lac.plot,
                              top = text_grob("Male Lactate", size = 18),
                              bottom = text_grob("Diet", size = 16))
m.lac.plot

#Alanine Next#
hf.alanine.m6.lm <- lm(alanine ~ trt * protein, data = hfm6)
anova(hf.alanine.m6.lm)
summary(hf.alanine.m6.lm)

ggqqplot(hfm6, x = "alanine", color = "group")

m6.ala.res <- residuals(hf.alanine.m6.lm)
ggqqplot(m6.ala.res)

m6.alanine.em <- emmeans(hf.alanine.m6.lm, specs = c("trt", "protein"))
m6.alanine.em

# Proper/Simple Post-Hoc

emmeans(hf.alanine.m6.lm, pairwise ~ protein)

pairs(m6.alanine.em, simple = "protein")
hfm6 %>% 
  emmeans_test(
    alanine ~ protein, model = hf.alanine.m6.lm
  )

plot(pt, comparisons = T)

#Plot
m6.ala.sum <- hfm6 %>% group_by(group) %>%
  summarise(emmean = mean(alanine),
            SD = sd(alanine),
            SE = sd(alanine)/sqrt(length(alanine)),
            n = length(alanine))

m6.ala.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
m6.ala.sum <- m6.ala.sum[order(m6.ala.sum$group, decreasing = T),]

names(m6.ala.sum)[c(1,6)] <- c("lab", "group")

m6.ala.plot <- ggplot(data = m6.ala.sum, 
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
  geom_errorbar(data = m6.ala.sum, aes(y = emmean, ymin = emmean, ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  annotate("segment", x = "group 4", xend = "group 3", y = 730, yend = 730) +
  annotate("text", x = 3.5, y = 780, label = "B", size = 5) +
  annotate("segment", x = "group 6", xend = "group 5", y = 740, yend = 740) +
  annotate("text", x = 1.5, y = 790, label = "A", size = 5) +
  annotate("text", x = .6, y = 900, fontface = 3, label = "DPS (Beef > Casein)", hjust = 0, size = 5) +
  geom_point(data = hfm6, aes(x = group, y = alanine), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 4),
             size = 1, stroke = 1) +
  labs(title = "Alanine",
       subtitle = "Males 6-Months",
       y ="\U00B5g/mL",
       x = "Diet")
m6.ala.plot

#### alanine 12 Months ####
hf.alanine.m12.lm <- lm(alanine ~ trt * protein, data = hfm12)
anova(hf.alanine.m12.lm)
summary(hf.alanine.m12.lm)

ggqqplot(hfm12, x = "alanine", color = "group")

m12.ala.res <- residuals(hf.alanine.m12.lm)
ggqqplot(m12.ala.res)

m12.alanine.em <- emmeans(hf.alanine.m12.lm, specs = c("trt", "protein"))
m12.alanine.em

# Proper/Simple Post-Hoc

emmip(hf.alanine.m12.lm, protein ~ trt)
emmeans(hf.alanine.m12.lm, pairwise ~ protein)

pairs(m12.alanine.em, simple = "protein")
pairs(m12.alanine.em, simple = c("trt", "protein"))

# Graphing labels ()
#Plot
m12.ala.sum <- hfm12 %>% group_by(group) %>%
  summarise(emmean = mean(alanine),
            SD = sd(alanine),
            SE = sd(alanine)/sqrt(length(alanine)),
            n = length(alanine))

m12.ala.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
m12.ala.sum <- m12.ala.sum[order(m12.ala.sum$group, decreasing = T),]

names(m12.ala.sum)[c(1,6)] <- c("lab", "group")

m12.ala.plot <- ggplot(data = m12.ala.sum, 
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
  annotate("segment", x = "group 4", xend = "group 3", y = 730, yend = 730) +
  annotate("text", x = 3.5, y = 780, label = "B", size = 5) +
  annotate("segment", x = "group 6", xend = "group 5", y = 740, yend = 740) +
  annotate("text", x = 1.5, y = 790, label = "A", size = 5) +
  annotate("text", x = .6, y = 900, fontface = 3, label = "DPS (Beef > Casein)", hjust = 0, size = 5) +
  geom_errorbar(data = m12.ala.sum, 
                aes(y = emmean, ymin = emmean, 
                    ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  geom_point(data = hfm12, aes(x = group, y = alanine), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 4),
             size = 1, stroke = 1) +
  labs(title = "Alanine",
       subtitle = "Males 12-Months",
       y ="\U00B5g/mL",
       x = "Diet")
m12.ala.plot

#### alanine 18 Months ####
hf.alanine.m18.lm <- lm(alanine ~ trt * protein, data = hfm18)
anova(hf.alanine.m18.lm)
summary(hf.alanine.m18.lm)

ggqqplot(hfm18, x = "alanine", color = "group")

m18.ala.res <- residuals(hf.alanine.m18.lm)
ggqqplot(m18.ala.res)

m18.alanine.em <- emmeans(hf.alanine.m18.lm, specs = c("trt", "protein"))
m18.alanine.em

# Proper/Simple Post-Hoc
pairs(m18.alanine.em, simple = c("trt", "protein"))

emmeans(hf.alanine.m18.lm, pairwise ~ trt * protein)

# Graphing labels ()
### Customize Here ###
m18ala.lab <- paste("group", c(3:6))
m18ala.lab <- data.frame(m18ala.lab)
m18ala.lab$lab <- c("AB", "AB", "B", "A")
names(m18ala.lab)[1] <- "group1"
m18ala.lab$group2 <- NA
m18ala.lab$group <- m18ala.lab$group1

# Customize for each set of comp

#Plot
m18.ala.sum <- hfm18 %>% group_by(group) %>%
  summarise(emmean = mean(alanine),
            SD = sd(alanine),
            SE = sd(alanine)/sqrt(length(alanine)),
            n = length(alanine))

m18ala.lab$y <- m18.ala.sum$emmean + m18.ala.sum$SD + 50

m18.ala.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
m18.ala.sum <- m18.ala.sum[order(m18.ala.sum$group, decreasing = T),]

names(m18.ala.sum)[c(1,6)] <- c("lab", "group")

m18.ala.plot <- ggplot(data = m18.ala.sum, 
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
        plot.subtitle = element_text(hjust = 0.5)) +
  geom_errorbar(data = m18.ala.sum, aes(y = emmean, ymin = emmean, ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  stat_pvalue_manual(data = m18ala.lab,
                     label = "lab",
                     x = "group",
                     y.position = "y",
                     label.size = 5) +
  geom_point(data = hfm18, aes(x = group, y = alanine), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 8),
             size = 1, stroke = 1) +
  annotate("text", x = .6, y = 900, fontface = 3, label = "DPS * AHE", hjust = 0, size = 5) +
  labs(title = "Alanine",
       subtitle = "Males 18-Months",
       y ="\U00B5g/mL",
       x = "Diet")
m18.ala.plot

m.ala.plot <- ggarrange(m6.ala.plot +
                          theme(plot.title = element_blank(),
                                axis.title.x = element_blank(),
                                axis.title.y = element_blank()),
                        m12.ala.plot + 
                          theme(plot.title = element_blank(),
                                axis.text.y = element_blank(),
                                axis.ticks.y = element_blank(),
                                axis.title.y = element_blank(),
                                axis.title.x = element_blank()),
                        m18.ala.plot +
                          theme(plot.title = element_blank(),
                                axis.text.y = element_blank(),
                                axis.ticks.y = element_blank(),
                                axis.title.y = element_blank(),
                                axis.title.x = element_blank()),
                        nrow = 1)

m.ala.plot <- annotate_figure(m.ala.plot,
                              top = text_grob("Male Alanine", size = 18),
                              bottom = text_grob("Diet", size = 16),
                              left = text_grob("Concentration\U00B5g/mL", size = 16,
                                               rot = 90))
m.ala.plot
## Glutamine Next##
ggqqplot(hfm.met$glutamine)

hf.glutamine.m6.lm <- lm(glutamine ~ trt * protein, data = hfm6)
anova(hf.glutamine.m6.lm)
summary(hf.glutamine.m6.lm)

ggqqplot(hfm6, x = "glutamine", color = "group")

m6.gln.res <- residuals(hf.glutamine.m6.lm)
ggqqplot(m6.gln.res)

m6.glutamine.em <- emmeans(hf.glutamine.m6.lm, specs = c("trt", "protein"))
m6.glutamine.em

#Plot
m6.gln.sum <- hfm6 %>% group_by(group) %>%
  summarise(emmean = mean(glutamine),
            SD = sd(glutamine),
            SE = sd(glutamine)/sqrt(length(glutamine)),
            n = length(glutamine))

m6.gln.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
m6.gln.sum <- m6.gln.sum[order(m6.gln.sum$group, decreasing = T),]

names(m6.gln.sum)[c(1,6)] <- c("lab", "group")

m6.gln.plot <- ggplot(data = m6.gln.sum, 
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
  geom_errorbar(data = m6.gln.sum, aes(y = emmean, ymin = emmean, ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  geom_point(data = hfm6, aes(x = group, y = glutamine), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 4),
             size = 1, stroke = 1) +
  labs(title = "Glutamine",
       subtitle = "Males 6-Months",
       y ="\U00B5g/mL",
       x = "Diet")
m6.gln.plot

#### glutamine 12 Months ####
hf.glutamine.m12.lm <- lm(glutamine ~ trt * protein, data = hfm12)
anova(hf.glutamine.m12.lm)
summary(hf.glutamine.m12.lm)

ggqqplot(hfm12, x = "glutamine", color = "group")

m12.gln.res <- residuals(hf.glutamine.m12.lm)
ggqqplot(m12.gln.res)

shapiro.test(hfm12$glutamine)

m12.glutamine.em <- emmeans(hf.glutamine.m12.lm, specs = c("trt", "protein"))
m12.glutamine.em

# Proper/Simple Post-Hoc

emmip(hf.glutamine.m12.lm, protein ~ trt)
emmeans(hf.glutamine.m12.lm, pairwise ~ protein)

pairs(m12.glutamine.em, simple = "protein")

# Graphing labels ()
#Plot
m12.gln.sum <- hfm12 %>% group_by(group) %>%
  summarise(emmean = mean(glutamine),
            SD = sd(glutamine),
            SE = sd(glutamine)/sqrt(length(glutamine)),
            n = length(glutamine))

m12.gln.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
m12.gln.sum <- m12.gln.sum[order(m12.gln.sum$group, decreasing = T),]

names(m12.gln.sum)[c(1,6)] <- c("lab", "group")

m12.gln.plot <- ggplot(data = m12.gln.sum, 
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
  annotate("segment", x = "group 4", xend = "group 3", y = 200, yend = 200) +
  annotate("text", x = 3.5, y = 220, label = "B", size = 5) +
  annotate("segment", x = "group 6", xend = "group 5", y = 270, yend = 270) +
  annotate("text", x = 1.5, y = 290, label = "A", size = 5) +
  annotate("text", x = .6, y = 350, fontface = 3, label = "DPS (Casein > Beef)", hjust = 0, size = 5) +
  geom_errorbar(data = m12.gln.sum, 
                aes(y = emmean, ymin = emmean, 
                    ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  geom_point(data = hfm12, aes(x = group, y = glutamine), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 4),
             size = 1, stroke = 1) +
  labs(title = "Glutamine",
       subtitle = "12-Months",
       y ="\U00B5g/mL",
       x = "Diet")
m12.gln.plot

#### glutamine 18 Months ####
nrow(hfm6)
hf.glutamine.m18.lm <- lm(glutamine ~ trt * protein, data = hfm18)
anova(hf.glutamine.m18.lm)
summary(hf.glutamine.m18.lm)

ggqqplot(hfm18, x = "glutamine", color = "group")

m18.gln.res <- residuals(hf.glutamine.m18.lm)
ggqqplot(m18.gln.res)


m18.glutamine.em <- emmeans(hf.glutamine.m18.lm, specs = c("trt", "protein"))
m18.glutamine.em

# Proper/Simple Post-Hoc
pairs(m18.glutamine.em, simple = c("trt", "protein"))

# Graphing labels ()
### Customize Here ###
m18gln.lab <- paste("group", c(3:6))
m18gln.lab <- data.frame(m18gln.lab)
m18gln.lab$lab <- c("B", "AB", "A", "AB")
names(m18gln.lab)[1] <- "group1"
m18gln.lab$group2 <- NA
m18gln.lab$group <- m18gln.lab$group1

# Customize for each set of comp

#Plot
m18.gln.sum <- hfm18 %>% group_by(group) %>%
  summarise(emmean = mean(glutamine),
            SD = sd(glutamine),
            SE = sd(glutamine)/sqrt(length(glutamine)),
            n = length(glutamine))

m18gln.lab$y <- m18.gln.sum$emmean + m18.gln.sum$SD + 10

m18.gln.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
m18.gln.sum <- m18.gln.sum[order(m18.gln.sum$group, decreasing = T),]

names(m18.gln.sum)[c(1,6)] <- c("lab", "group")

#m18gln.lab[4, 5] <- 260

m18.gln.plot <- ggplot(data = m18.gln.sum, 
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
        plot.subtitle = element_text(hjust = 0.5)) +
  geom_errorbar(data = m18.gln.sum, aes(y = emmean, ymin = emmean, ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  stat_pvalue_manual(data = m18gln.lab,
                     label = "lab",
                     x = "group",
                     y.position = "y",
                     label.size = 5) +
  geom_point(data = hfm18, aes(x = group, y = glutamine), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 8),
             size = 1, stroke = 1) +
  annotate("text", x = .6, y = 350, fontface = 3, label = "DPS * AHE", hjust = 0, size = 5) +
  labs(title = "Glutamine",
       subtitle = "18-Months",
       y ="\U00B5g/mL",
       x = "Diet")
m18.gln.plot

m.gln.plot <- ggarrange(m6.gln.plot +
                          theme(plot.title = element_blank(),
                                axis.title.x = element_blank(),
                                axis.title.y = element_blank()),
                        m12.gln.plot + 
                          theme(plot.title = element_blank(),
                                axis.text.y = element_blank(),
                                axis.ticks.y = element_blank(),
                                axis.title.y = element_blank(),
                                axis.title.x = element_blank()),
                        m18.gln.plot +
                          theme(plot.title = element_blank(),
                                axis.text.y = element_blank(),
                                axis.ticks.y = element_blank(),
                                axis.title.y = element_blank(),
                                axis.title.x = element_blank()),
                        nrow = 1)

m.gln.plot <- annotate_figure(m.gln.plot,
                              top = text_grob("Male Glutamine", size = 18),
                              bottom = text_grob("Diet", size = 16),
                              left = text_grob("Concentration\U00B5g/mL", size = 16,
                                               rot = 90))
m.gln.plot


# Carnitine #
hf.carnitine.m6.lm <- lm(carnitine ~ trt * protein, data = hfm6)
anova(hf.carnitine.m6.lm)
summary(hf.carnitine.m6.lm)

ggqqplot(hfm6, x = "carnitine", color = "group")

m6.crn.res <- residuals(hf.carnitine.m6.lm)
ggqqplot(m6.crn.res)

m6.carnitine.em <- emmeans(hf.carnitine.m6.lm, specs = c("trt", "protein"))
m6.carnitine.em

# Proper/Simple Post-Hoc
emmeans(hf.carnitine.m6.lm, pairwise ~ protein)

pairs(m6.carnitine.em, simple = "protein")

hfm6 %>% 
  emmeans_test(
    carnitine ~ protein, model = hf.carnitine.m6.lm
  )

plot(pt, comparisons = T)

#Plot
m6.crn.sum <- hfm6 %>% group_by(group) %>%
  summarise(emmean = mean(carnitine),
            SD = sd(carnitine),
            SE = sd(carnitine)/sqrt(length(carnitine)),
            n = length(carnitine))

#m6crn.lab$y <- m6.crn.sum$emmean + m6.crn.sum$SE

m6.crn.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
m6.crn.sum <- m6.crn.sum[order(m6.crn.sum$group, decreasing = T),]

names(m6.crn.sum)[c(1,6)] <- c("lab", "group")

m6.crn.plot <- ggplot(data = m6.crn.sum, 
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
  geom_errorbar(data = m6.crn.sum, aes(y = emmean, ymin = emmean, ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  geom_point(data = hfm6, aes(x = group, y = carnitine), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 4),
             size = 1, stroke = 1) +
  annotate("segment", x = "group 4", xend = "group 3", y = 63, yend = 63) +
  annotate("text", x = 3.5, y = 70, label = "B", size = 5) +
  annotate("segment", x = "group 6", xend = "group 5", y = 55, yend = 55) +
  annotate("text", x = 1.5, y = 62, label = "A", size = 5) +
  annotate("text", x = .6, y = 130, fontface = 3, label = "DPS (Beef > Casein)", hjust = 0, size = 5) +
  labs(title = "Carnitine",
       subtitle = "Males 6-Months",
       y ="\U00B5g/mL",
       x = "Diet")
m6.crn.plot

#### carnitine 12 Months ####
hf.carnitine.m12.lm <- lm(carnitine ~ trt * protein, data = hfm12)
anova(hf.carnitine.m12.lm)
summary(hf.carnitine.m12.lm)

ggqqplot(hfm12, x = "carnitine", color = "group")

m12.crn.res <- residuals(hf.carnitine.m12.lm)
ggqqplot(m12.crn.res)

m12.carnitine.em <- emmeans(hf.carnitine.m12.lm, specs = c("trt", "protein"))
m12.carnitine.em

# Proper/Simple Post-Hoc

emmip(hf.carnitine.m12.lm, protein ~ trt)
emmeans(hf.carnitine.m12.lm, pairwise ~ protein * trt)

pairs(m12.carnitine.em, simple = "trt")
pairs(m12.carnitine.em, simple = c("trt", "protein"))

# Graphing labels ()
m12crn.lab <- paste("group", c(3:6))
m12crn.lab <- data.frame(m12crn.lab)
m12crn.lab$lab <- c("A", "B", "A", "A")
names(m12crn.lab)[1] <- "group1"
m12crn.lab$group2 <- NA
m12crn.lab$group <- m12crn.lab$group1

#Plot
m12.crn.sum <- hfm12 %>% group_by(group) %>%
  summarise(emmean = mean(carnitine),
            SD = sd(carnitine),
            SE = sd(carnitine)/sqrt(length(carnitine)),
            n = length(carnitine))

m12crn.lab$y <- m12.crn.sum$emmean + m12.crn.sum$SD + 5

m12.crn.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
m12.crn.sum <- m12.crn.sum[order(m12.crn.sum$group, decreasing = T),]

names(m12.crn.sum)[c(1,6)] <- c("lab", "group")

m12.crn.plot <- ggplot(data = m12.crn.sum, 
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
  stat_pvalue_manual(data = m12crn.lab,
                     label = "lab",
                     x = "group",
                     y.position = "y",
                     label.size = 5) +
  geom_errorbar(data = m12.crn.sum, 
                aes(y = emmean, ymin = emmean, 
                    ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  geom_point(data = hfm12, aes(x = group, y = carnitine), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 4),
             size = 1, stroke = 1) +
  annotate("text", x = .6, y = 130, fontface = 3, label = "DPS * AHE",  hjust = 0, size = 5) +
  labs(title = "Carnitine",
       subtitle = "Males 12-Months",
       y ="\U00B5g/mL",
       x = "Diet")
m12.crn.plot

#### carnitine 18 Months ####
hf.carnitine.m18.lm <- lm(carnitine ~ trt * protein, data = hfm18)
anova(hf.carnitine.m18.lm)
summary(hf.carnitine.m18.lm)

ggqqplot(hfm18, x = "carnitine", color = "group")

m18.crn.res <- residuals(hf.carnitine.m18.lm)
ggqqplot(m18.crn.res)

m18.carnitine.em <- emmeans(hf.carnitine.m18.lm, specs = c("trt", "protein"))
m18.carnitine.em

### Customize Here ###
#Plot
m18.crn.sum <- hfm18 %>% group_by(group) %>%
  summarise(emmean = mean(carnitine),
            SD = sd(carnitine),
            SE = sd(carnitine)/sqrt(length(carnitine)),
            n = length(carnitine))

m18.crn.sum$lab <- c("group 3", "group 4", "group 5", "group 6")
m18.crn.sum <- m18.crn.sum[order(m18.crn.sum$group, decreasing = T),]

names(m18.crn.sum)[c(1,6)] <- c("lab", "group")

m18.crn.plot <- ggplot(data = m18.crn.sum, 
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
        plot.subtitle = element_text(hjust = 0.5)) +
  geom_errorbar(data = m18.crn.sum, aes(y = emmean, ymin = emmean, ymax = emmean + SE, x = group), width = 0.1, color = "black") +
  geom_point(data = hfm18, aes(x = group, y = carnitine), shape = 16, 
             position = position_jitterdodge(jitter.width = .3, dodge.width = .1,
                                             seed = 8),
             size = 1, stroke = 1) +
  labs(title = "Carnitine",
       subtitle = "Males 18-Months",
       y ="\U00B5g/mL",
       x = "Diet")
m18.crn.plot

m.crn.plot <- ggarrange(m6.crn.plot +
                          theme(plot.title = element_blank(),
                                axis.title.x = element_blank(),
                                axis.title.y = element_blank()),
                        m12.crn.plot + 
                          theme(plot.title = element_blank(),
                                axis.text.y = element_blank(),
                                axis.ticks.y = element_blank(),
                                axis.title.y = element_blank(),
                                axis.title.x = element_blank()),
                        m18.crn.plot +
                          theme(plot.title = element_blank(),
                                axis.text.y = element_blank(),
                                axis.ticks.y = element_blank(),
                                axis.title.y = element_blank(),
                                axis.title.x = element_blank()),
                        nrow = 1)

m.crn.plot <- annotate_figure(m.crn.plot,
                              top = text_grob("Male Carnitine", size = 18),
                              bottom = text_grob("Diet", size = 16),
                              left = text_grob("Concentration\U00B5g/mL", size = 16,
                                               rot = 90))
m.crn.plot



#### Use Correct Graphs ####
# 6 months column: first_plot for first 4, first_bottom_plot for carnitine (bottom)
col_m_6 <- plot_grid(
  first_plot(m6.glu.plot),
  first_plot(m6.lac.plot),
  first_plot(m6.ala.plot),
  first_plot(m6.gln.plot),
  first_bottom_plot(m6.crn.plot),
  ncol = 1,
  align = "v",
  axis = "lr"
)

# 12 months column: clean_plot for first 4, bottom_plot for carnitine
col_m_12 <- plot_grid(
  clean_plot(m12.glu.plot),
  clean_plot(m12.lac.plot),
  clean_plot(m12.ala.plot),
  clean_plot(m12.gln.plot),
  bottom_plot(m12.crn.plot),
  ncol = 1,
  align = "v",
  axis = "lr"
)

# 18 months column: clean_plot for first 4, bottom_plot for carnitine
col_m_18 <- plot_grid(
  clean_plot(m18.glu.plot),
  clean_plot(m18.lac.plot),
  clean_plot(m18.ala.plot),
  clean_plot(m18.gln.plot),
  bottom_plot(m18.crn.plot),
  ncol = 1,
  align = "v",
  axis = "lr"
)

# Combine everything into a 5 x 4 layout: y-axis labels + 3 timepoints
final_m_plot <- plot_grid(
  ylab_col,
  col_m_6,
  col_m_12,
  col_m_18,
  ncol = 4,
  rel_widths = c(0.3, 1, 1, 1),
  align = "h",
  axis = "tb"
)

final_m_plot



col_titles <- plot_grid(
  NULL,
  ggdraw() + draw_label("6 Months", fontface = "bold", hjust = 0.5),
  ggdraw() + draw_label("12 Months", fontface = "bold", hjust = 0.5),
  ggdraw() + draw_label("18 Months", fontface = "bold", hjust = 0.5),
  ncol = 4,
  rel_widths = c(0.3, 1, 1, 1)
)

# Final with title row
plot_m_with_titles <- plot_grid(
  col_titles,
  final_m_plot,
  ncol = 1,
  rel_heights = c(0.05, 1)
)

main_m_title <- ggdraw() + draw_label("Metabolite Concentrations from Male Livers", fontface = "bold", size = 20)

m.full_figure <- plot_grid(
  main_m_title,
  plot_m_with_titles,
  ncol = 1,
  rel_heights = c(0.07, 1)
)

m.full_figure

ggsave("metabolites from males.pdf", m.full_figure, width = 10, height = 12)

metabolites.full_figure <- f.full_figure + m.full_figure

ggsave("Metabolites From Both Sexes.pdf", metabolites.full_figure, width = 20, height = 12)


aov.list <- list(m6.glu = anova(hf.glucose.m6.lm),
                 m6.lac = anova(hf.lactate.m6.lm),
                 m6.ala = anova(hf.alanine.m6.lm),
                 m6.gln = anova(hf.glutamine.m6.lm),
                 m6.crn = anova(hf.carnitine.m6.lm),
                 m12.glu = anova(hf.glucose.m12.lm),
                 m12.lac = anova(hf.lactate.m12.lm),
                 m12.ala = anova(hf.alanine.m12.lm),
                 m12.gln = anova(hf.glutamine.m12.lm),
                 m12.crn = anova(hf.carnitine.m12.lm),
                 m18.glu = anova(hf.glucose.m18.lm),
                 m18.lac = anova(hf.lactate.m18.lm),
                 m18.ala = anova(hf.alanine.m18.lm),
                 m18.gln = anova(hf.glutamine.m18.lm),
                 m18.crn = anova(hf.carnitine.m18.lm))

anova.combined <- do.call(rbind, lapply(names(aov.list), function(name) {
  df <- as.data.frame(aov.list[[name]])
  df$model <- name  # Add identifier
  df$term <- rownames(aov.list[[name]])  # Preserve rownames as a column
  rownames(df) <- NULL  # Optional: clean up row names
  df
}))

anova.combined

write.csv(anova.combined, "./Data/m.anova.combined.csv")


# Correct this for each new comparison 6, 12, 18 #
write.csv(anova(hf.glucose.m.lm), "./Data/m.glucose.anova.csv")
write.csv(anova(hf.lactate.m.lm), "./Data/m.lactate.anova.csv")
write.csv(anova(hf.ala.m.lm), "./Data/m.ala.anova.csv")
write.csv(anova(hf.gln.m.lm), "./Data/m.gln.anova.csv")
write.csv(anova(hf.crn.m.lm), "./Data/m.crn.anova.csv")
#####