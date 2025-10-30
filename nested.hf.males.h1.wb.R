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
library(cowplot) # combine plots
library(knitr) # kable tables
library(kableExtra) # kable_styling tables

library(insight)
library(lazyWeave)
library(dplyr)

## Data Read-In
getwd()
h1m <- read_xlsx("./raw-R work/Dissertation Work/Male Age and Group WB.xlsx",
                    sheet = "H1 Males")
bcat1m <- h1m[which(names(h1m) %in% c("B-Cat R1","B-Cat R2",
                                 "B-Cat R3", "B-Cat R4"))]
bcat1m
bcat1m$group <- "HFC"
bcat1m$group[1:9] <- rep(c("HFBN", "HFB", "HFCN"), each = 3)

bcat1m$trt <- "no"
bcat1m$trt[c(1:3,7:9)] <- "yes"

bcat1m$protein <- "c"
bcat1m$protein[1:6] <- "beef"

bcat1m <- melt(bcat1m, by = group)
bcat1m

bcat1m$id <- rep(c(1:3), times = 16)

bcat1m.sum <-  bcat1m %>% group_by(group, id, trt, protein) %>%
  summarise(emmean = mean(value),
            SE = sd(value)/sqrt(length(value)),
            n = length(value))

bcat1m.sum

bcat1m.lm <- lm(emmean ~ trt * protein, data = bcat1m.sum)
anova(bcat1m.lm)
summary(bcat1m.lm)
#write.csv(anova(bcat1m.lm), "./Data/bcat1m.lm.csv")

bcat1m.nest <- bcat1m.sum %>% group_by(group) %>%
  summarise(emmean1 = mean(emmean),
            n = length(emmean),
            SD = sd(emmean),
            SE = SD/sqrt(n))

bcat1m.plot <- ggplot(bcat1m.nest, aes(x = factor(group, level = c("HFC", "HFCN", "HFB", "HFBN")),
                      emmean1, fill = group)) +
    geom_col(fill = c("chartreuse4","purple", "brown1", "deepskyblue3"), width = 0.5,
           col = "black") +
  theme_pubr() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 13),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_errorbar(data = bcat1m.nest, 
                aes(y = emmean1, ymin = emmean1, ymax = emmean1 + SE),
                width = 0.1, color = "black") +
  geom_point(data = bcat1m.sum, aes(x = group, y = emmean), fill = "black", 
              size = 2) +
    labs(title = "\u03b2-catenin") +
    scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0))
bcat1m.plot

## CYP

cyp1m <- h1m[which(names(h1m) %in% c("CYP R1","CYP R2",
                                      "CYP R3", "CYP R4"))]
cyp1m
cyp1m$group <- "HFC"
cyp1m$group[1:9] <- rep(c("HFBN", "HFB", "HFCN"), each = 3)

cyp1m$trt <- "no"
cyp1m$trt[c(1:3,7:9)] <- "yes"

cyp1m$protein <- "c"
cyp1m$protein[1:6] <- "beef"

cyp1m <- melt(cyp1m, by = group)
cyp1m$id <- rep(c(1:3), times = 16)

cyp1m <- cyp1m[-c(5, 17, 29, 41),]

cyp1m.sum <-  cyp1m %>% group_by(group, id, trt, protein) %>%
  summarise(emmean = mean(value),
            SE = sd(value)/sqrt(length(value)),
            n = length(value))
cyp1m.sum

cyp1m.lm <- lm(emmean ~ trt * protein, data = cyp1m.sum)
anova(cyp1m.lm)
summary(cyp1m.lm)
#write.csv(anova(cyp1m.lm), "./Data/cyp1m.lm.csv")

#Graph it
cyp1m.nest<-  cyp1m.sum %>% group_by(group) %>%
  summarise(emmean1 = mean(emmean),
            SD = sd(emmean),
            SE = SD/sqrt(length(emmean)),
            n = length(emmean))

cyp1m.plot <- ggplot(cyp1m.nest, aes(x = factor(group, level = c("HFC", "HFCN", "HFB", "HFBN")),
                     emmean1, fill = group)) +
  geom_col(fill = c("chartreuse4","purple",
                    "brown1", "deepskyblue3"), col = "black", width = 0.5) +
  theme_pubr() +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        legend.position = "none",
        axis.title = element_blank()) +
  geom_errorbar(data = cyp1m.nest, 
                aes(y = emmean1, ymin = emmean1, ymax = emmean1 + SE),
                width = 0.1, color = "black") +
  geom_point(data = cyp1m.sum, aes(x = group, y = emmean), fill = "black", 
              size = 2) +  
  labs(title = "CYP3A4") +
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0))
cyp1m.plot

## GS

gs1m <- h1m[which(names(h1m) %in% c("GS R1","GS R2",
                                     "GS R3", "GS R4"))]
gs1m
gs1m$group <- "HFC"
gs1m$group[1:9] <- rep(c("HFBN", "HFB", "HFCN"), each = 3)


gs1m$trt <- "no"
gs1m$trt[c(1:3,7:9)] <- "yes"

gs1m$protein <- "c"
gs1m$protein[1:6] <- "beef"

gs1m <- melt(gs1m, by = group)

gs1m$id <- rep(c(1:3), times = 16)

gs1m.sum <-  gs1m %>% group_by(group, id, trt, protein) %>%
  summarise(emmean = mean(value),
            SE = sd(value)/sqrt(length(value)),
            n = length(value))
gs1m.sum


gs1m.lm <- lm(emmean ~ trt * protein, data = gs1m.sum)
anova(gs1m.lm)
summary(gs1m.lm)
#write.csv(anova(gs1m.lm), "./Data/gs1m.lm.csv")

gs1m.em <- emmeans(gs1m.lm, specs = c("trt", "protein"))
gs1m.em

gs1m.dt <- data.table(summary(gs1m.em))
gs1m.dt

gs1m.dt$group <- c("HFB", "HFBN", "HFC", "HFCN")

#Graph it
gs1m.nest <-  gs1m.sum %>% group_by(group) %>%
  summarise(emmean1 = mean(emmean),
            SD = sd(emmean),
            SE = SD/sqrt(length(emmean)),
            n = length(emmean))

gs1m.plot <- ggplot(gs1m.nest, aes(x = factor(group, level = c("HFC", "HFCN", "HFB", "HFBN")),
                    emmean1, fill = group)) +
  geom_col(fill = c("chartreuse4","purple",
                    "brown1", "deepskyblue3"), col = "black", width = 0.5) +
  theme_pubr() +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        legend.position = "none",
        axis.title = element_blank()) +
  geom_errorbar(data = gs1m.nest, 
                aes(y = emmean1, ymin = emmean1, ymax = emmean1 + SE),
                width = 0.1, color = "black") +
  geom_point(data = gs1m.sum, aes(x = group, y = emmean), fill = "black", 
              size = 2) +
  labs(title = "Glutamine Synthetase") +
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0))
gs1m.plot

m6 <- ggarrange(bcat1m.plot + border(), gs1m.plot + border(), cyp1m.plot + border(), nrow = 1)
m6 <- annotate_figure(m6, 
                top = text_grob("Relative Protein Levels in Males at 6 Months", size = 14),
                bottom = text_grob("Diets", size = 14))
m6
