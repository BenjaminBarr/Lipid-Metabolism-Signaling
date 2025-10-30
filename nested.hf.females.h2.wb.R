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
h2f<- read_xlsx("./raw-R work/Dissertation Work/Female Age and Group WB.xlsx",
                    sheet = "H2 Females")
h2f
bcat2f<- h2f[which(names(h2f) %in% c("B-Cat R1","B-Cat R2",
                                 "B-Cat R3", "B-Cat R4"))]
bcat2f
bcat2f$group <- "HFC"
bcat2f$group[1:9] <- rep(c("HFBN", "HFB", "HFCN"), each = 3)

bcat2f$trt <- "no"
bcat2f$trt[c(1:3,7:9)] <- "yes"

bcat2f$protein <- "c"
bcat2f$protein[1:6] <- "beef"

bcat2f <- bcat2f[-c(13,14),]

bcat2f<- reshape2::melt(bcat2f, by = "group")
bcat2f

bcat2f$id <- rep(c(1:3), times = 16)

bcat2f.sum <-  bcat2f %>% group_by(group, id, trt, protein) %>%
  summarise(emmean = mean(value),
            SE = sd(value)/sqrt(length(value)),
            n = length(value))

bcat2f.sum
bcat2f.sum <- bcat2f.sum[-2,]
bcat2f.lm <- lm(emmean ~ trt * protein, data = bcat2f.sum)
anova(bcat2f.lm)
summary(bcat2f.lm)
#write.csv(anova(bcat2f.lm), "./Data/bcat2f.lm.csv")

#Graph it
bcat2f.nest <- bcat2f.sum %>% group_by(group) %>%
  summarise(emmean1 = mean(emmean),
            SD = sd(emmean),
            SE = SD/sqrt(length(emmean)),
            n = length(emmean))

bcat2f.nest$group <- factor(bcat2f.nest$group, levels = c("HFC", "HFCN", "HFB", "HFBN"))
bcat2f.sum$group <- factor(bcat2f.sum$group, levels = c("HFC", "HFCN", "HFB", "HFBN"))

bcat2f.plot <- ggplot(bcat2f.nest, aes(x = group, emmean1, fill = group)) +
  geom_col(col = "black", width = 0.5) +
  scale_fill_manual(name = "Diets", values = c(
    "HFC" = "brown1",
    "HFCN" = "deepskyblue3",
    "HFB" = "chartreuse4",
    "HFBN" = "purple"
  )) +
  geom_errorbar(data = bcat2f.nest, 
                aes(y = emmean1, ymin = emmean1, ymax = emmean1 + SE),
                width = 0.1, color = "black") +
  geom_point(data = bcat2f.sum, aes(x = group, y = emmean), fill = "black", size = 2) +
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

bcat2f.plot

## CYP

cyp2f<- h2f[which(names(h2f) %in% c("CYP R1","CYP R2",
                                      "CYP R3", "CYP R4"))]
cyp2f
cyp2f$group <- "HFC"
cyp2f$group[1:9] <- rep(c("HFBN", "HFB", "HFCN"), each = 3)

cyp2f$trt <- "no"
cyp2f$trt[c(1:3,7:9)] <- "yes"

cyp2f$protein <- "c"
cyp2f$protein[1:6] <- "beef"

cyp2f <- cyp2f[-c(13,14),]

cyp2f<- reshape2::melt(cyp2f, by = "group")

cyp2f$id <- rep(c(1:3), times = 16)

cyp2f.sum <-  cyp2f %>% group_by(group, id, trt, protein) %>%
  summarise(emmean = mean(value),
            SE = sd(value)/sqrt(length(value)),
            n = length(value))

cyp2f.sum
cyp2f.sum <- cyp2f.sum[-2,]
cyp2f.lm <- lm(emmean ~ trt * protein, data = cyp2f.sum)
anova(cyp2f.lm)
summary(cyp2f.lm)
#write.csv(anova(cyp2f.lm), "./Data/cyp2f.lm.csv")

#cyp2f.em <- emmeans(cyp2f.lm, specs = c("trt", "protein"))
#cyp2f.em

#cyp2f.dt <- data.table(summary(cyp2f.em))
#cyp2f.dt

#cyp2f.dt$group <- c("HFB", "HFBN", "HFC", "HFCN")

#Graph 
cyp2f.nest <- cyp2f.sum %>% group_by(group) %>%
  summarise(emmean1 = mean(emmean),
            SD = sd(emmean),
            SE = SD/sqrt(length(emmean)),
            n = length(emmean))
cyp2f.sum

cyp2f.nest$group <- factor(cyp2f.nest$group, levels = c("HFC", "HFCN", "HFB", "HFBN"))
cyp2f.sum$group <- factor(cyp2f.sum$group, levels = c("HFC", "HFCN", "HFB", "HFBN"))

cyp2f.plot <- ggplot(cyp2f.nest, aes(x = group, emmean1, fill = group)) +
  geom_col(col = "black", width = 0.5) +
  scale_fill_manual(name = "Diets", values = c(
    "HFC" = "brown1",
    "HFCN" = "deepskyblue3",
    "HFB" = "chartreuse4",
    "HFBN" = "purple"
  )) +
  geom_errorbar(data = cyp2f.nest, 
                aes(y = emmean1, ymin = emmean1, ymax = emmean1 + SE),
                width = 0.1, color = "black") +
  geom_point(data = cyp2f.sum, aes(x = group, y = emmean), fill = "black", size = 2) +
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

cyp2f.plot

## GS

gs2f<- h2f[which(names(h2f) %in% c("GS R1","GS R2",
                                     "GS R3", "GS R4"))]

gs2f
gs2f$group <- "HFC"
gs2f$group[1:9] <- rep(c("HFBN", "HFB", "HFCN"), each = 3)


gs2f$trt <- "no"
gs2f$trt[c(1:3,7:9)] <- "yes"

gs2f$protein <- "c"
gs2f$protein[1:6] <- "beef"

gs2f <- gs2f[-c(13,14),]

gs2f<- reshape2::melt(gs2f, by = "group")

gs2f$id <- rep(c(1:3), times = 16)

gs2f.sum <-  gs2f %>% group_by(group, id, trt, protein) %>%
  summarise(emmean = mean(value),
            SE = sd(value)/sqrt(length(value)),
            n = length(value))

gs2f.sum
gs2f.sum <- gs2f.sum[-2,]

gs2f.lm <- lm(emmean ~ trt * protein, data = gs2f.sum)
anova(gs2f.lm)
summary(gs2f.lm)

#write.csv(anova(gs2f.lm), "./Data/gs2f.lm.csv")


#gs2f.em <- emmeans(gs2f.lm, specs = c("trt", "protein"))
#gs2f.em

#gs2f.dt <- data.table(summary(gs2f.em))
#gs2f.dt

#gs2f.dt$group <- c("HFB", "HFBN", "HFC", "HFCN")

#Graph 
gs2f.nest <- gs2f.sum %>% group_by(group) %>%
  summarise(emmean1 = mean(emmean),
            SD = sd(emmean),
            SE = SD/sqrt(length(emmean)),
            n = length(emmean))
gs2f.nest

gs2f.nest$group <- factor(gs2f.nest$group, levels = c("HFC", "HFCN", "HFB", "HFBN"))
gs2f.sum$group <- factor(gs2f.sum$group, levels = c("HFC", "HFCN", "HFB", "HFBN"))

gs2f.plot <- ggplot(gs2f.nest, aes(x = group, emmean1, fill = group)) +
  geom_col(col = "black", width = 0.5) +
  scale_fill_manual(name = "Diets", values = c(
    "HFC" = "brown1",
    "HFCN" = "deepskyblue3",
    "HFB" = "chartreuse4",
    "HFBN" = "purple"
  )) +
  geom_errorbar(data = gs2f.nest, 
                aes(y = emmean1, ymin = emmean1, ymax = emmean1 + SE),
                width = 0.1, color = "black") +
  geom_point(data = gs2f.sum, aes(x = group, y = emmean), fill = "black", size = 2) +
  labs(title = "Glutamine Synthetase") +
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

gs2f.plot

f12 <- ggarrange(first_plot(bcat2f.plot), clean_plot(gs2f.plot), clean_plot(cyp2f.plot), nrow = 1)
                 #common.legend = TRUE, legend = "bottom")
f12 <- annotate_figure(f12, 
                top = text_grob("Relative Protein Levels in Females at 12 Months", size = 22))
f12
