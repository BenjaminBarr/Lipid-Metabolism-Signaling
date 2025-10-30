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
h1f<- read_xlsx("./raw-R work/Dissertation Work/Female Age and Group WB.xlsx",
                sheet = "H1 Females")
h1f
bcat1f<- h1f[which(names(h1f) %in% c("B-Cat R1","B-Cat R2",
                                     "B-Cat R3", "B-Cat R4"))]
bcat1f
bcat1f$group <- "HFC"
bcat1f$group[1:6] <- rep(c("HFB", "HFCN"), each = 3)

bcat1f$trt <- "no"
bcat1f$trt[c(4:6)] <- "yes"

bcat1f$protein <- "c"
bcat1f$protein[1:3] <- "beef"

bcat1f <- bcat1f[-1]
bcat1f$trt[10:11] <- "yes"
bcat1f$group[10:11] <- "HFBN"
bcat1f$protein[10:11] <- "beef"

bcat1f

bcat1f<- reshape2::melt(bcat1f, by = "group")
bcat1f

bcat1f$id <- rep(c(rep(c(1:3), times = 3), 1:2), times = 3)
bcat1f.sum <-  bcat1f %>% group_by(group, id, trt, protein) %>%
  summarise(emmean = mean(value),
            SE = sd(value)/sqrt(length(value)),
            n = length(value))

bcat1f.sum
bcat1f.lm <- lm(emmean ~ trt * protein, data = bcat1f.sum)
lm(emmean ~ trt, data = bcat1f.sum)
bcat1f.ahe.aov <- anova(lm(emmean ~ trt, data = bcat1f.sum))
bcat1f.dps.aov <- anova(lm(emmean ~ protein, data = bcat1f.sum))
#write.csv(bcat1f.ahe.aov, "./Data/bcat1f.ahe.csv")
#write.csv(bcat1f.dps.aov, "./Data/bcat1f.dps.csv")


#anova(bcat1f.lm)
#summary(bcat1f.lm)

#bcat1f.em <- emmeans(bcat1f.lm, specs = c("trt", "protein"))
#bcat1f.em

#is.factor(bcat1f$trt)

#bcat1f.dt <- data.table(summary(bcat1f.em))
#bcat1f.dt
#bcat1f.dt$group <- c("HFB", "HFBN", "HFC", "HFCN")


#Graph it
bcat1f.nest <- bcat1f.sum %>% group_by(group) %>%
  summarise(emmean1 = mean(emmean),
            n = length(emmean),
            SD = sd(emmean),
            SE = SD/sqrt(length(emmean)))
bcat1f.nest


bcat1f.nest[2,c(2, 4,5)] <- 0


bcat1f.nest$group <- factor(bcat1f.nest$group, levels = c("HFC", "HFCN", "HFB", "HFBN"))
bcat1f.sum$group <- factor(bcat1f.sum$group, levels = c("HFC", "HFCN", "HFB", "HFBN"))

bcat1f.plot <- ggplot(bcat1f.nest, aes(x = group, y = emmean1, fill = group)) +
  geom_col(col = "black", width = 0.5) +
  scale_fill_manual(name = "Diets", values = c(
    "HFC" = "brown1",
    "HFCN" = "deepskyblue3",
    "HFB" = "chartreuse4",
    "HFBN" = "purple"
  )) +
  geom_errorbar(aes(ymin = emmean1, ymax = emmean1 + SE),
                width = 0.1, color = "black") +
  geom_point(data = bcat1f.sum, aes(x = group, y = emmean), fill = "black", size = 2) +
  labs(title = expression("beta-Catenin")) +
  theme_pubr() +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0))

bcat1f.plot

## CYP

cyp1f<- h1f[which(names(h1f) %in% c("CYP R1","CYP R2",
                                    "CYP R3", "CYP R4"))]
cyp1f
cyp1f$group <- "HFC"
cyp1f$group[1:6] <- rep(c("HFB", "HFCN"), each = 3)

cyp1f$trt <- "no"
cyp1f$trt[c(4:6)] <- "yes"

cyp1f$protein <- "c"
cyp1f$protein[1:3] <- "beef"

cyp1f$trt[10:11] <- "yes"
cyp1f$group[10:11] <- "HFBN"
cyp1f$protein[10:11] <- "beef"
cyp1f[11, 1] <- NA
cyp1f[1] <- as.numeric(unlist(cyp1f[1]))


cyp1f<- reshape2::melt(cyp1f, by = "group")

cyp1f

cyp1f$id <- rep(c(rep(c(1:3), times = 3), 1:2), times = 4)
cyp1f.sum <-  cyp1f %>% group_by(group, id, trt, protein) %>%
  summarise(emmean = mean(value),
            SE = sd(value)/sqrt(length(value)),
            n = length(value))

cyp1f.sum
cyp1f.lm <- lm(emmean ~ trt * protein, data = cyp1f.sum)
cyp1f.ahe.aov <- anova(lm(emmean ~ trt, data = cyp1f.sum))
cyp1f.dps.aov <- anova(lm(emmean ~ protein, data = cyp1f.sum))
#write.csv(cyp1f.ahe.aov, "./Data/cyp1f.ahe.csv")
#write.csv(cyp1f.dps.aov, "./Data/cyp1f.dps.csv")

anova(cyp1f.lm)
summary(cyp1f.lm)

#cyp1f.em <- emmeans(cyp1f.lm, specs = c("trt", "protein"))
#cyp1f.em

#cyp1f.dt <- data.table(summary(cyp1f.em))
#cyp1f.dt

#cyp1f.dt$group <- c("HFB", "HFBN", "HFC", "HFCN")

#Graph 
cyp1f.nest <- cyp1f.sum %>% group_by(group) %>%
  summarise(emmean1 = mean(emmean),
            n = length(emmean),
            SD = sd(emmean),
            SE = SD/sqrt(length(emmean)))
cyp1f.nest
cyp1f.nest[2, c(2,4,5)] <- 0

cyp1f.nest$group <- factor(cyp1f.nest$group, levels = c("HFC", "HFCN", "HFB", "HFBN"))
cyp1f.sum$group <- factor(cyp1f.sum$group, levels = c("HFC", "HFCN", "HFB", "HFBN"))

cyp1f.plot <- ggplot(cyp1f.nest, aes(x = group, y = emmean1, fill = group)) +
  geom_col(col = "black", width = 0.5) +
  scale_fill_manual(name = "Diets", values = c(
    "HFC" = "brown1",
    "HFCN" = "deepskyblue3",
    "HFB" = "chartreuse4",
    "HFBN" = "purple"
  )) +
  geom_errorbar(aes(ymin = emmean1, ymax = emmean1 + SE),
                width = 0.1, color = "black") +
  geom_point(data = cyp1f.sum, aes(x = group, y = emmean), fill = "black", size = 2) +
  labs(title = "CYP3A4") +
  theme_pubr() +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0))
cyp1f.plot


## GS

gs1f<- h1f[which(names(h1f) %in% c("GS R1","GS R2",
                                   "GS R3", "GS R4"))]
gs1f
gs1f$group <- "HFC"
gs1f$group[1:6] <- rep(c("HFB", "HFCN"), each = 3)

gs1f$trt <- "no"
gs1f$trt[c(4:6)] <- "yes"

gs1f$protein <- "c"
gs1f$protein[1:3] <- "beef"

gs1f <- gs1f[-1]
gs1f$trt[10:11] <- "yes"
gs1f$group[10:11] <- "HFBN"
gs1f$protein[10:11] <- "beef"

gs1f<- reshape2::melt(gs1f, by = "group")

gs1f$id <- rep(c(rep(c(1:3), times = 3), 1:2), times = 3)
gs1f.sum <-  gs1f %>% group_by(group, id, trt, protein) %>%
  summarise(emmean = mean(value),
            SE = sd(value)/sqrt(length(value)),
            n = length(value))

gs1f.sum
gs1f.lm <- lm(emmean ~ trt * protein, data = gs1f.sum)
#anova(gs1f.lm)
#summary(gs1f.lm)

gs1f.ahe.aov <- anova(lm(emmean ~ trt, data = gs1f.sum))
gs1f.dps.aov <- anova(lm(emmean ~ protein, data = gs1f.sum))
#summary(lm(emmean ~ trt, data = gs1f.sum))
#write.csv(gs1f.ahe.aov, "./Data/gs1f.ahe.csv")
#write.csv(gs1f.dps.aov, "./Data/gs1f.dps.csv")

#gs1f.em <- emmeans(gs1f.lm, specs = c("trt", "protein"))
#gs1f.em

#gs1f.dt <- data.table(summary(gs1f.em))
#gs1f.dt


#gs1f.dt$group <- c("HFB", "HFBN", "HFC", "HFCN")

#Graph 
gs1f.nest <- gs1f.sum %>% group_by(group) %>%
  summarise(emmean1 = mean(emmean),
            n = length(emmean),
            SD = sd(emmean),
            SE = SD/sqrt(length(emmean)))
gs1f.nest
gs1f.nest[2, c(2,4,5)] <- 0

gs1f.nest$group <- factor(gs1f.nest$group, levels = c("HFC", "HFCN", "HFB", "HFBN"))
gs1f.sum$group <- factor(gs1f.sum$group, levels = c("HFC", "HFCN", "HFB", "HFBN"))

gs1f.plot <- ggplot(gs1f.nest, aes(x = group,
                    emmean1, fill = group)) +
  geom_col(col = "black", width = 0.5) +
  scale_fill_manual(name = "Diets", values = c("HFC" = "brown1",
                                               "HFCN" = "deepskyblue3",
                                               "HFB" = "chartreuse4",
                                               "HFBN" = "purple")) +
  geom_errorbar(data = gs1f.nest, 
                aes(y = emmean1, ymin = emmean1, ymax = emmean1 + SE),
                width = 0.1, color = "black") +
  geom_point(data = gs1f.sum, aes(x = group, y = emmean), fill = "black", 
              size = 2) +
  labs(title = "Glutamine Synthetase") +
  theme_pubr() +
  theme(legend.position = "none",
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0))

gs1f.plot

f6wb.title <- expression(atop("Relative Protein Levels in Female Livers at 6 Months*"))#, scriptstyle("(Main Effects Only)"

f6 <- ggarrange(first_plot(bcat1f.plot), clean_plot(gs1f.plot), clean_plot(cyp1f.plot), nrow = 1)
                #common.legend = TRUE, legend = "bottom")
f6 <- annotate_figure(f6, 
                top = text_grob("Relative Protein Levels in Female Livers at 6 Months*", size = 22))
f6
