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


## Data Read-In
getwd()
h2m <- read_xlsx("./raw-R work/Dissertation Work/Male Age and Group WB.xlsx",
                    sheet = "H2 Males")
bcat2m <- h2m[which(names(h2m) %in% c("B-Cat R1","B-Cat R2",
                                 "B-Cat R3", "B-Cat R4"))]
bcat2m
bcat2m$group <- "HFC"
bcat2m$group[1:9] <- rep(c("HFBN", "HFB", "HFCN"), each = 3)

bcat2m$trt <- "no"
bcat2m$trt[c(1:3,7:9)] <- "yes"

bcat2m$protein <- "c"
bcat2m$protein[1:6] <- "beef"

bcat2m <- melt(bcat2m, by = group)
bcat2m
bcat2m$id <- rep(c(1:3), times = 16)

bcat2m.sum <-  bcat2m %>% group_by(group, id, trt, protein) %>%
  summarise(emmean = mean(value),
            SE = sd(value)/sqrt(length(value)),
            n = length(value))

bcat2m.sum

bcat2m.lm <- lm(emmean ~ trt * protein, data = bcat2m.sum)
anova(bcat2m.lm)
summary(bcat2m.lm)
#write.csv(anova(bcat2m.lm), "./Data/bcat2m.lm.csv")

bcat2m.em <- emmeans(bcat2m.lm, specs = c("trt", "protein"))
bcat2m.em

bcat2m.dt <- data.table(summary(bcat2m.em))
bcat2m.dt

#Graph it
bcat2m.nest <- bcat2m.sum %>% group_by(group) %>%
  summarise(emmean1 = mean(emmean),
            n = length(emmean),
            SD = sd(emmean),
            SE = SD/sqrt(n))
bcat2m.nest

bcat2m.plot <- ggplot(bcat2m.nest, aes(x = factor(group, level = c("HFC", "HFCN", "HFB", "HFBN")),
                      emmean1, fill = group)) +
  geom_col(fill = c("chartreuse4","purple",
                    "brown1", "deepskyblue3"), col = "black", width = 0.5) +
  theme_pubr() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 13),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_errorbar(data = bcat2m.nest, 
                aes(y = emmean1, ymin = emmean1, ymax = emmean1 + SE),
                width = 0.1, color = "black") +
  geom_point(data = bcat2m.sum, aes(x = group, y = emmean), fill = "black", 
              size = 2) +
  labs(title = "\u03b2-catenin") +
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0))

bcat2m.plot

## CYP

cyp2m <- h2m[which(names(h2m) %in% c("CYP R1","CYP R2",
                                      "CYP R3", "CYP R4"))]
cyp2m
cyp2m$group <- "HFC"
cyp2m$group[1:9] <- rep(c("HFBN", "HFB", "HFCN"), each = 3)

cyp2m$trt <- "no"
cyp2m$trt[c(1:3,7:9)] <- "yes"

cyp2m$protein <- "c"
cyp2m$protein[1:6] <- "beef"

cyp2m <- melt(cyp2m, by = group)

cyp2m$id <- rep(c(1:3), times = 16)

cyp2m.sum <-  cyp2m %>% group_by(group, id, trt, protein) %>%
  summarise(emmean = mean(value),
            SE = sd(value)/sqrt(length(value)),
            n = length(value))

cyp2m.sum

cyp2m.lm <- lm(emmean ~ trt * protein, data = cyp2m.sum)
anova(cyp2m.lm)
summary(cyp2m.lm)
#write.csv(anova(cyp2m.lm), "./Data/cyp2m.lm.csv")

cyp2m.em <- emmeans(cyp2m.lm, specs = c("trt", "protein"))
cyp2m.em

cyp2m.dt <- data.table(summary(cyp2m.em))
cyp2m.dt

cyp2m.dt$group <- c("HFB", "HFBN", "HFC", "HFCN")

cyp2m.nest <- cyp2m.sum %>% group_by(group) %>%
  summarise(emmean1 = mean(emmean),
            n = length(emmean),
            SD = sd(emmean),
            SE = SD/sqrt(n))
cyp2m.nest

cyp2m.plot <- ggplot(cyp2m.nest, aes(x = factor(group, level = c("HFC", "HFCN", "HFB", "HFBN")),
                     emmean1, fill = group)) +
  geom_col(fill = c("chartreuse4","purple",
                    "brown1", "deepskyblue3"), col = "black", width = 0.5) +
  theme_pubr() +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        legend.position = "none",
        axis.title = element_blank()) +
  geom_errorbar(data = cyp2m.nest, 
                aes(y = emmean1, ymin = emmean1, ymax = emmean1 + SE),
                width = 0.1, color = "black") +
  geom_point(data = cyp2m.sum, aes(x = group, y = emmean), fill = "black", 
              size = 2) +
  labs(title = "CYP3A4") +
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0))
 
cyp2m.plot  
## GS

gs2m <- h2m[which(names(h2m) %in% c("GS R1","GS R2",
                                     "GS R3", "GS R4"))]

gs2m
gs2m$group <- "HFC"
gs2m$group[1:9] <- rep(c("HFBN", "HFB", "HFCN"), each = 3)


gs2m$trt <- "no"
gs2m$trt[c(1:3,7:9)] <- "yes"

gs2m$protein <- "c"
gs2m$protein[1:6] <- "beef"


gs2m <- melt(gs2m, by = group)

gs2m$id <- rep(c(1:3), times = 16)

gs2m <- na.omit(gs2m)

gs2m.sum <-  gs2m %>% group_by(group, id, trt, protein) %>%
  summarise(emmean = mean(value),
            SE = sd(value)/sqrt(length(value)),
            n = length(value))

gs2m.sum

gs2m.lm <- lm(emmean ~ trt * protein, data = gs2m.sum)
anova(gs2m.lm)
summary(gs2m.lm)
#write.csv(anova(gs2m.lm), "./Data/gs2m.lm.csv")

gs2m.em <- emmeans(gs2m.lm, specs = c("trt", "protein"))
gs2m.em

gs2m.dt <- data.table(summary(gs2m.em))
gs2m.dt

gs2m.dt$group <- c("HFB", "HFBN", "HFC", "HFCN")

#Graph 
gs2m.nest <- gs2m.sum %>% group_by(group) %>%
  summarise(emmean1 = mean(emmean),
            n = length(emmean),
            SD = sd(emmean),
            SE = SD/sqrt(n))
gs2m.nest

gs2m.plot <- ggplot(gs2m.nest, aes(x = factor(group, level = c("HFC", "HFCN", "HFB", "HFBN")),
                    emmean1, fill = group)) +
  geom_col(fill = c("chartreuse4","purple",
                    "brown1", "deepskyblue3"), col = "black", , width = 0.5) +
  theme_pubr() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 13),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_errorbar(data = gs2m.nest, 
                aes(y = emmean1, ymin = emmean1, ymax = emmean1 + SE),
                width = 0.1, color = "black") +
  geom_point(data = gs2m.sum, aes(x = group, y = emmean), fill = "black", size = 2) +
  labs(title = "Glutamine Synthetase") +
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0))

gs2m.plot

m12 <- ggarrange(bcat2m.plot + border(), gs2m.plot + border(), cyp2m.plot + border(), nrow = 1)
m12 <- annotate_figure(m12, 
                top = text_grob("Relative Protein Levels in Males at 12 Months", size = 14),
                bottom = text_grob("Diets", size = 14))
m12
