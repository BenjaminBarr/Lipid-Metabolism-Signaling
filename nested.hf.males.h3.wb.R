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
h3m <- read_xlsx("./raw-R work/Dissertation Work/Male Age and Group WB.xlsx",
                    sheet = "H3 Males")
h3m
bcat3m <- h3m[which(names(h3m) %in% c("B-Cat R1","B-Cat R2",
                                 "B-Cat R3", "B-Cat R4"))]
bcat3m
bcat3m$group <- "HFC"
bcat3m$group[1:9] <- rep(c("HFBN", "HFB", "HFCN"), each = 3)

bcat3m$trt <- "no"
bcat3m$trt[c(1:3,7:9)] <- "yes"

bcat3m$protein <- "c"
bcat3m$protein[1:6] <- "beef"

# Remove Tumors

bcat3m <- bcat3m[1:11,]

bcat3m <- melt(bcat3m, by = group)
bcat3m$id <- rep(c(rep(c(1:3), times = 3), 1:2), times = 4)


bcat3m.sum <-  bcat3m %>% group_by(group, id, trt, protein) %>%
  summarise(emmean = mean(value),
            SE = sd(value)/sqrt(length(value)),
            n = length(value))

bcat3m.sum
bcat3m.lm <- lm(emmean ~ trt * protein, data = bcat3m.sum)
anova(bcat3m.lm)
summary(bcat3m.lm)
#write.csv(anova(bcat3m.lm), "./Data/bcat3m.lm.csv")

bcat3m.em <- emmeans(bcat3m.lm, specs = c("trt", "protein"))
bcat3m.em

bcat3m.dt <- data.table(summary(bcat3m.em))
bcat3m.dt

pairs(bcat3m.em, simple = "trt")
pairs(bcat3m.em, simple = c("protein", "trt"))
emmeans(bcat3m.lm, pairwise ~ trt)

bcat3m.dt$group <- c("HFB", "HFBN", "HFC", "HFCN")

#Graph it
bcat3m.nest <- bcat3m.sum %>% group_by(group) %>%
  summarise(emmean1 = mean(emmean),
            n = length(emmean),
            SD = sd(emmean),
            SE = SD/sqrt(n))
bcat3m.nest
bcat3m.nest$sig <- c("A", "AB", "A", "B")
bcat3m.nest$y <- bcat3m.nest$emmean1 + bcat3m.nest$SD + .05
bcat3m.nest$group1 <- bcat3m.nest$group
bcat3m.nest$group2 <- NA

bcat3m.plot <- ggplot(bcat3m.nest, aes(x = factor(group, level = c("HFC", "HFCN", "HFB", "HFBN")),
                      emmean1, fill = group)) +
  geom_col(fill = c("chartreuse4","purple",
                    "brown1", "deepskyblue3"), col = "black", width = 0.5) +
  theme_pubr() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 13),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_errorbar(data = bcat3m.nest, 
                aes(y = emmean1, ymin = emmean1, ymax = emmean1 + SE),
                width = 0.1, color = "black") +
  geom_point(data = bcat3m.sum, aes(x = group, y = emmean), fill = "black", 
              size = 2) +
  annotate("text", x = 1, y = 1.3, label = "A") +
  annotate("text", x = 3, y = 1.45, label = "A") +
  annotate("text", x = 2, y = .8, label = "B") +
  annotate("text", x = 4, y = 1, label = "B") +
  annotate("text", x = .6, y = 2.3, fontface = 3, 
           label = "AHE (Ref > AHE)", hjust = 0, size = 5) +
  labs(title = "\u03b2-catenin") +
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0))
bcat3m.plot

## CYP

cyp3m <- h3m[which(names(h3m) %in% c("CYP R1","CYP R2",
                                      "CYP R3", "CYP R4"))]
cyp3m
cyp3m$group <- "HFC"
cyp3m$group[1:9] <- rep(c("HFBN", "HFB", "HFCN"), each = 3)

cyp3m$trt <- "no"
cyp3m$trt[c(1:3,7:9)] <- "yes"

cyp3m$protein <- "c"
cyp3m$protein[1:6] <- "beef"

# Remove Tumors
cyp3m <- cyp3m[1:11,]

cyp3m <- melt(cyp3m, by = group)
cyp3m$id <- rep(c(rep(c(1:3), times = 3), 1:2), times = 4)

cyp3m.sum <-  cyp3m %>% group_by(group, id, trt, protein) %>%
  summarise(emmean = mean(value),
            SE = sd(value)/sqrt(length(value)),
            n = length(value))

cyp3m.sum
cyp3m.lm <- lm(emmean ~ trt * protein, data = cyp3m.sum)
anova(cyp3m.lm)
summary(cyp3m.lm)
#write.csv(anova(cyp3m.lm), "./Data/cyp3m.lm.csv")


cyp3m.em <- emmeans(cyp3m.lm, specs = c("trt", "protein"))
cyp3m.em

cyp3m.dt <- data.table(summary(cyp3m.em))
cyp3m.dt

pairs(cyp3m.em, simple = "trt")
emmeans(cyp3m.lm, pairwise ~ trt)

cyp3m.dt$group <- c("HFB", "HFBN", "HFC", "HFCN")

#Graph it
cyp3m.nest <- cyp3m.sum %>% group_by(group) %>%
  summarise(emmean1 = mean(emmean),
            n = length(emmean),
            SD = sd(emmean),
            SE = SD/sqrt(n))
cyp3m.nest

cyp3m.plot <- ggplot(cyp3m.nest, aes(x = factor(group, level = c("HFC", "HFCN", "HFB", "HFBN")),
                     emmean1, fill = group)) +
  geom_col(fill = c("chartreuse4","purple",
                    "brown1", "deepskyblue3"), col = "black", width = 0.5) +
  theme_pubr() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 13),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_errorbar(data = cyp3m.nest, 
                aes(y = emmean1, ymin = emmean1, ymax = emmean1 + SE),
                width = 0.1, color = "black") +
  geom_point(data = cyp3m.sum, aes(x = group, y = emmean), fill = "black", 
              size = 2) +
  annotate("text", x = .6, y = 2.3, fontface = 3, 
           label = "AHE (Ref > AHE)", hjust = 0, size = 5) +
  annotate("text", x = 1, y = 1.5, label = "A") +
  annotate("text", x = 3, y = 1.3, label = "A") +
  annotate("text", x = 2, y = 1.1, label = "B") +
  annotate("text", x = 4, y = 1.1, label = "B") +
  labs(title = "CYP3A4") +
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0))
cyp3m.plot

## GS

gs3m <- h3m[which(names(h3m) %in% c("GS R1","GS R2",
                                     "GS R3", "GS R4"))]

gs3m
gs3m$group <- "HFC"
gs3m$group[1:9] <- rep(c("HFBN", "HFB", "HFCN"), each = 3)


gs3m$trt <- "no"
gs3m$trt[c(1:3,7:9)] <- "yes"

gs3m$protein <- "c"
gs3m$protein[1:6] <- "beef"

# Remove Tumors
gs3m <- gs3m[1:11,]

gs3m <- melt(gs3m, by = group)
gs3m$id <- rep(c(rep(c(1:3), times = 3), 1:2), times = 4)

gs3m.sum <-  gs3m %>% group_by(group, id, trt, protein) %>%
  summarise(emmean = mean(value),
            SE = sd(value)/sqrt(length(value)),
            n = length(value))

gs3m.sum
gs3m.lm <- lm(emmean ~ trt * protein, data = gs3m.sum)
anova(gs3m.lm)
summary(gs3m.lm)

#write.csv(anova(gs3m.lm), "./Data/gs3m.lm.csv")

gs3m.em <- emmeans(gs3m.lm, specs = c("trt", "protein"))
gs3m.em

gs3m.dt <- data.table(summary(gs3m.em))
gs3m.dt

gs3m.dt$group <- c("HFB", "HFBN", "HFC", "HFCN")

#Graph 
gs3m.nest <- gs3m.sum %>% group_by(group) %>%
  summarise(emmean1 = mean(emmean),
            n = length(emmean),
            SD = sd(emmean),
            SE = SD/sqrt(n))
gs3m.nest

gs3m.plot <- ggplot(gs3m.nest, aes(x = factor(group, level = c("HFC", "HFCN", "HFB", "HFBN")),
                    emmean1, fill = group)) +
  geom_col(fill = c("chartreuse4","purple",
                    "brown1", "deepskyblue3"), col = "black",, width = 0.5) +
  theme_pubr() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 13),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_errorbar(data = gs3m.nest, 
                aes(y = emmean1, ymin = emmean1, ymax = emmean1 + SE),
                width = 0.1, color = "black") +
  geom_point(data = gs3m.sum, aes(x = group, y = emmean), fill = "black", 
              size = 2) +
  labs(title = "Glutamine Synthetase") +
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0))
gs3m.plot

col_m_bcat <- plot_grid(
  first_plot(bcat1m.plot),
  first_plot(bcat2m.plot),
  first_bottom_plot(bcat3m.plot),
  ncol = 1,
  align = "v",
  axis = "lr"
)

col_m_gs <- plot_grid(
  clean_plot(gs1m.plot),
  clean_plot(gs2m.plot),
  bottom_plot(gs3m.plot),
  ncol = 1,
  align = "v",
  axis = "lr"
)

col_m_cyp <- plot_grid(
  clean_plot(cyp1m.plot),
  clean_plot(cyp2m.plot),
  bottom_plot(cyp3m.plot),
  ncol = 1,
  align = "v",
  axis = "lr"
)


final_m_wb_plot <- plot_grid(
  age.lab,
  col_m_bcat,
  col_m_gs,
  col_m_cyp,
  ncol = 4,
  rel_widths = c(0.2, 1, 1, 1),
  align = "h",
  axis = "tb"
)

wb_col_titles <- plot_grid(
  NULL,
  ggdraw() + draw_label("\u03B2-Catenin", fontface = "bold", hjust = 0.5),
  ggdraw() + draw_label("GS", fontface = "bold", hjust = 0.5),
  ggdraw() + draw_label("CYP3A4", fontface = "bold", hjust = 0.5),
  ncol = 4,
  rel_widths = c(0.3, 1, 1, 1))

# Final with title row
mwb_plot_with_titles <- plot_grid(
  wb_col_titles,
  final_m_wb_plot,
  ncol = 1,
  rel_heights = c(0.05, 1)
)

wb_m_title <- ggdraw() + draw_label("Relative Hepatic Protein Levels in Males", fontface = "bold", size = 20)

wb.m.full_figure <- plot_grid(
  wb_m_title,
  mwb_plot_with_titles,
  ncol = 1,
  rel_heights = c(0.07, 1)
)

wb.m.full_figure

ggsave("Hepatic Protein Levels in Males.pdf", wb.m.full_figure, width = 8, height = 9)

wb.f.full_figure + wb.m.full_figure
ggsave("Hepatic Protein Levels in Females and Males.pdf", 
       wb.f.full_figure + wb.m.full_figure, width = 16, height = 9)
