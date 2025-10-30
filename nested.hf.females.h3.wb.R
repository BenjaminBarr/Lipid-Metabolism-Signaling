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
library(grid)

## Data Read-In
getwd()
h3f<- read_xlsx("./raw-R work/Dissertation Work/Female Age and Group WB.xlsx",
                    sheet = "H3 Females")
h3f
bcat3f<- h3f[which(names(h3f) %in% c("B-Cat R1","B-Cat R2",
                                 "B-Cat R3", "B-Cat R4"))]
bcat3f
bcat3f$group <- "HFC"
bcat3f$group[1:9] <- rep(c("HFBN", "HFB", "HFCN"), each = 3)

bcat3f$trt <- "no"
bcat3f$trt[c(1:3,7:9)] <- "yes"

bcat3f$protein <- "c"
bcat3f$protein[1:6] <- "beef"


bcat3f<- reshape2::melt(bcat3f, by = "group")
bcat3f

bcat3f$id <- rep(c(1:3), times = 16)

bcat3f.sum <-  bcat3f %>% group_by(group, id, trt, protein) %>%
  summarise(emmean = mean(value),
            SE = sd(value)/sqrt(length(value)),
            n = length(value))

bcat3f.sum
bcat3f.lm <- lm(emmean ~ trt * protein, data = bcat3f.sum)
anova(bcat3f.lm)
summary(bcat3f.lm)

#write.csv(anova(bcat3f.lm), "./Data/bcat3f.lm.csv")


#bcat3f.em <- emmeans(bcat3f.lm, specs = c("trt", "protein"))
#bcat3f.em

#bcat3f.dt <- data.table(summary(bcat3f.em))
#bcat3f.dt

### Graph 
bcat3f.nest <- bcat3f.sum %>% group_by(group) %>%
  summarise(emmean1 = mean(emmean),
            SD = sd(emmean),
            SE = SD/sqrt(length(emmean)),
            n = length(emmean))
bcat3f.nest

bcat3f.nest$group <- factor(bcat3f.nest$group, levels = c("HFC", "HFCN", "HFB", "HFBN"))
bcat3f.sum$group <- factor(bcat3f.sum$group, levels = c("HFC", "HFCN", "HFB", "HFBN"))

bcat3f.plot <- ggplot(bcat3f.nest, aes(x = group, emmean1, fill = group)) +
  geom_col(col = "black", width = 0.5) +
  scale_fill_manual(name = "Diets", values = c(
    "HFC" = "brown1",
    "HFCN" = "deepskyblue3",
    "HFB" = "chartreuse4",
    "HFBN" = "purple"
  )) +
  geom_errorbar(data = bcat3f.nest, 
                aes(y = emmean1, ymin = emmean1, ymax = emmean1 + SE),
                width = 0.1, color = "black") +
  geom_point(data = bcat3f.sum, aes(x = group, y = emmean), fill = "black", size = 2) +
  labs(title = expression("beta-Catenin")) +
  theme_pubr() +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0))

bcat3f.plot


## CYP

cyp3f<- h3f[which(names(h3f) %in% c("CYP R1","CYP R2",
                                      "CYP R3", "CYP R4"))]
cyp3f
cyp3f$group <- "HFC"
cyp3f$group[1:9] <- rep(c("HFBN", "HFB", "HFCN"), each = 3)

cyp3f$trt <- "no"
cyp3f$trt[c(1:3,7:9)] <- "yes"

cyp3f$protein <- "c"
cyp3f$protein[1:6] <- "beef"

cyp3f<- reshape2::melt(cyp3f, by = "group")

cyp3f$id <- rep(c(1:3), times = 16)

cyp3f.sum <-  cyp3f %>% group_by(group, id, trt, protein) %>%
  summarise(emmean = mean(value),
            SE = sd(value)/sqrt(length(value)),
            n = length(value))

cyp3f.sum
cyp3f.lm <- lm(emmean ~ trt * protein, data = cyp3f.sum)
anova(cyp3f.lm)
summary(cyp3f.lm)

#write.csv(anova(cyp3f.lm), "./Data/cyp3f.lm.csv")

cyp3f.em <- emmeans(cyp3f.lm, specs = c("trt", "protein"))
cyp3f.em

emmip(cyp3f.lm, protein ~ trt)
emmeans(cyp3f.lm, pairwise ~ protein)

pairs(cyp3f.em, simple = "protein")

#cyp3f.dt <- data.table(summary(cyp3f.em))
#cyp3f.dt

#cyp3f.dt$group <- c("HFB", "HFBN", "HFC", "HFCN")

#Graph 
cyp3f.nest <- cyp3f.sum %>% group_by(group) %>%
  summarise(emmean1 = mean(emmean),
            SD = sd(emmean),
            SE = SD/sqrt(length(emmean)),
            n = length(emmean))
cyp3f.nest
cyp3f.nest$group <- factor(cyp3f.nest$group, levels = c("HFC", "HFCN", "HFB", "HFBN"))
cyp3f.sum$group <- factor(cyp3f.sum$group, levels = c("HFC", "HFCN", "HFB", "HFBN"))

cyp3f.plot <- ggplot(cyp3f.nest, aes(x = group, emmean1, fill = group)) +
  geom_col(col = "black", width = 0.5) +
  scale_fill_manual(name = "Diets", values = c(
    "HFC" = "brown1",
    "HFCN" = "deepskyblue3",
    "HFB" = "chartreuse4",
    "HFBN" = "purple"
  )) +
  geom_errorbar(data = cyp3f.nest, 
                aes(y = emmean1, ymin = emmean1, ymax = emmean1 + SE),
                width = 0.1, color = "black") +
  geom_point(data = cyp3f.sum, aes(x = group, y = emmean), fill = "black", size = 2) +
  annotate("segment", x = 1, xend = 2, y = 1.2, yend = 1.2) +
  annotate("text", x = 1.5, y = 1.4, label = "A", size = 5) +
  annotate("segment", x = 3, xend = 4, y = 1.4, yend = 1.4) +
  annotate("text", x = 3.5, y = 1.6, label = "B", size = 5) +
  annotate("text", x = .6, y = 2.3, fontface = 3, label = "DPS (Beef > Casein)", hjust = 0, size = 5) +
  labs(title = "CYP3A4") +
  theme_pubr() +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0))

cyp3f.plot

## GS

gs3f<- h3f[which(names(h3f) %in% c("GS R1","GS R2",
                                     "GS R3", "GS R4"))]
#gs3f<- gs3m[-c(13,26,39,52),]
gs3f
gs3f$group <- "HFC"
gs3f$group[1:9] <- rep(c("HFBN", "HFB", "HFCN"), each = 3)

gs3f$trt <- "no"
gs3f$trt[c(1:3,7:9)] <- "yes"

gs3f$protein <- "c"
gs3f$protein[1:6] <- "beef"

gs3f<- reshape2::melt(gs3f, by = "group")

gs3f$id <- rep(c(1:3), times = 16)

gs3f.sum <-  gs3f %>% group_by(group, id, trt, protein) %>%
  summarise(emmean = mean(value),
            SE = sd(value)/sqrt(length(value)),
            n = length(value))

gs3f.sum
gs3f.lm <- lm(emmean ~ trt * protein, data = gs3f.sum)
anova(gs3f.lm)
summary(gs3f.lm)

#write.csv(anova(gs3f.lm), "./Data/gs3f.lm.csv")

#gs3f.em <- emmeans(gs3f.lm, specs = c("trt", "protein"))
#gs3f.em

#gs3f.dt <- data.table(summary(gs3f.em))
#gs3f.dt

#gs3f.dt$group <- c("HFB", "HFBN", "HFC", "HFCN")

#Graph 
gs3f.nest <- gs3f.sum %>% group_by(group) %>%
  summarise(emmean1 = mean(emmean),
            SD = sd(emmean),
            SE = sd(emmean)/sqrt(length(emmean)),
            n = length(emmean))
gs3f.nest

gs3f.nest$group <- factor(gs3f.nest$group, levels = c("HFC", "HFCN", "HFB", "HFBN"))
gs3f.sum$group <- factor(gs3f.sum$group, levels = c("HFC", "HFCN", "HFB", "HFBN"))

gs3f.plot <- ggplot(gs3f.nest, aes(x = group, emmean1, fill = group)) +
  geom_col(col = "black", width = 0.5) +
  scale_fill_manual(name = "Diets", values = c(
    "HFC" = "brown1",
    "HFCN" = "deepskyblue3",
    "HFB" = "chartreuse4",
    "HFBN" = "purple"
  )) +
  geom_errorbar(data = gs3f.nest, 
                aes(y = emmean1, ymin = emmean1, ymax = emmean1 + SE),
                width = 0.1, color = "black") +
  geom_point(data = gs3f.sum, aes(x = group, y = emmean), fill = "black", size = 2) +
  labs(title = "Glutamine Synthetase") +
  theme_pubr() +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0))

gs3f.plot



######### Adjust functions appropriately

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
             size = 6, fontface = "bold") +
    theme_void()
}

#Make Col label function
make_col_label <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label,
             size = 7, hjust = 0.5) +
    theme_void()
}
age.lab <- plot_grid(make_y_label("6 Months"),
  make_y_label("12 Months"),
  make_y_label("18 Months"),
  ncol = 1,
  rel_heights = rep(1, 5)
)

col_bcat <- plot_grid(
  first_plot(bcat1f.plot),
  first_plot(bcat2f.plot),
  first_bottom_plot(bcat3f.plot),
  ncol = 1,
  align = "v",
  axis = "lr"
)

col_gs <- plot_grid(
  clean_plot(gs1f.plot),
  clean_plot(gs2f.plot),
  bottom_plot(gs3f.plot),
  ncol = 1,
  align = "v",
  axis = "lr"
)

col_cyp <- plot_grid(
  clean_plot(cyp1f.plot),
  clean_plot(cyp2f.plot),
  bottom_plot(cyp3f.plot),
  ncol = 1,
  align = "v",
  axis = "lr"
)


final_f_wb_plot <- plot_grid(
  age.lab,
  col_bcat,
  col_gs,
  col_cyp,
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
  rel_widths = c(0.2, 1, 1, 1))

# Final with title row
fwb_plot_with_titles <- plot_grid(
  wb_col_titles,
  final_f_wb_plot,
  ncol = 1,
  rel_heights = c(0.05, 1)
)

wb_f_title <- ggdraw() + draw_label("Relative Hepatic Protein Levels in Females", fontface = "bold", size = 20)

wb.f.full_figure <- plot_grid(
  wb_f_title,
  fwb_plot_with_titles,
  ncol = 1,
  rel_heights = c(0.07, 1)
)

wb.f.full_figure

ggsave("Hepatic Protein Levels in Females.pdf", wb.f.full_figure, width = 8, height = 9)




