#====================================================================================
#========== Final analysis
#========== Lupine x Subsidy impacts on Arthropods
#====================================================================================

# load packages
library(tidyverse)
library(lubridate)
library(vegan)
library(lme4)
library(car)
library(AICcmodavg)
library(cowplot)
library(ggplot2)
library(tidyverse)

# read data
pitfall <- read_csv("lupine_pitfall.csv")
soil <- read_csv("lupine_soil_arthropods.csv")
plants <- read_csv("lupine_plant.csv")

#Data notes on treatments
#Location = Border (Lupine-invaded), Heath (Control)
#Treatment = Control, KNO3, and Midge addition
#Taxa = Individual taxa representing the arthropods OR plants

# set theme
theme_set(theme_bw() %+replace%
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  plot.margin = margin(t = 1,
                                       b = 1,
                                       l = 1,
                                       r = 1),
                  axis.text = element_text(size = 10, 
                                           color="black", 
                                           family = "sans"),
                  axis.title = element_text(size =10),
                  axis.title.y = element_text(angle = 90,
                                              margin=margin(t = 0,
                                                            b = 0,
                                                            l = 0,
                                                            r = 5)),
                  axis.title.x = element_text(margin=margin(t = 5,
                                                            b = 0,
                                                            l = 0,
                                                            r = 0)),
                  panel.spacing = unit(0.1, "lines"),
                  axis.ticks = element_line(size = 0.25)))

#====================================================================================









#====================================================================================
#========== Pitfall traps (ground-active arthropods)
#====================================================================================

# clean data
# remove Lupine plot and springtails; group spiders and beetles
pitfall_filtered <- pitfall %>%
  filter(treatment != "Lupine",
         !(taxa %in% c("clmb","colj")))
pitfall_clean <- pitfall_filtered %>%
  mutate(taxa2 = ifelse(taxa %in% c("cara", "cole", "cryp", "curc", "stap"),
                         "coleoptera",
                         ifelse(taxa %in% c("aran", "gnap", "lyco","thom"),
                                "araneae", taxa)),
         year = year(coldate)) %>%
  group_by(block, plot, treatment, location, year, coldate, daysout, taxa2) %>%
  summarize(count = sum(count)) %>%
  mutate(catch_rate = count / daysout) %>%
  ungroup()

# examine data
summary(pitfall_clean)
pitfall_clean %>% group_by(taxa2) %>%
  summarize(mean = mean(count),
            zeros = mean(ifelse(count == 0, 1, 0)))

# preliminary plot (counts)
ggplot(data = pitfall_clean,
       aes(x = year,
           y = count,
           group = plot,
           color = treatment))+
  facet_grid(location~taxa2)+
  geom_line(size = 0.2)+
  scale_x_continuous(breaks = c(2014,2016,2018))+
  scale_y_continuous(trans = "log1p")

# preliminary plot (catch rate)
ggplot(data = pitfall_clean,
       aes(x = year,
           y = catch_rate,
           group = plot,
           color = treatment))+
  facet_grid(location~taxa2)+
  geom_line(size = 0.2)+
  scale_x_continuous(breaks = c(2014,2016,2018))+
  scale_y_continuous(trans = "log1p")

# prepare data for analysis
# transform response variable
# log1p transform and z-score (NOT by taxa)
# define factors
pitfall_prep <- pitfall_clean %>%
  mutate(zlog = log1p(catch_rate),
         zlog = (zlog - mean(zlog)) / sd(zlog),
         year_f = factor(year),
         treat_f = factor(treatment),
         plot_f = factor(plot),
         block_f = factor(block),
         loc_f = factor(location))

# fit LMM
# treat taxa as fixed effect (no partial pooling, but easy)
# plot as random effect; similar to repeated measures ANOVA
# interaction between year and taxon to account for variation among years
# NO interactions between year and treatment / location, since this is not really of primary interest
pitfall_lmm <- lmer(zlog ~ year_f * taxa2 + treat_f * taxa2 * loc_f + (1 | plot_f),
                    data = pitfall_prep)

# p-values for effects on overall catch rate
# use F-test (reasonable for this design)
Anova(pitfall_lmm, type = 3, test.statistic = "F") # full
Anova(update(pitfall_lmm, .~. - treat_f:taxa2:loc_f), 
      test.statistic = "F") # 3-way interaction dropped
Anova(update(pitfall_lmm, .~. - year_f:taxa2  - 
                                taxa2:treat_f - 
                                treat_f:loc_f - 
                                taxa2:loc_f - 
                                treat_f:taxa2:loc_f), 
      test.statistic = "F") # all interactions dropped

# These results indicate that the taxa differed in their responses to the treatments.
# Furthermore, the taxa differed among locations
# The taxa x loc x treat interaction is weak; this is probably just noise
# Taxa differed in their trends through time
# No obvious effect of traetmnt consistent across taxa (but there are common year effects)

# standardize observations to 2014, using coefficient estimates from model
summary(pitfall_lmm)$coefficients
coef_pos <- c(2,3,11:18)
coef_mat <- summary(pitfall_lmm)$coefficients[coef_pos,1]
mod_mat <- model.matrix(pitfall_lmm)[,coef_pos]
pitfall_prep$year_cor <- as.numeric(mod_mat %*% coef_mat)
pitfall_prep$zlog_cor = pitfall_prep$zlog - pitfall_prep$year_cor 

# generate predicted values
pitfall_nd <- tidyr::expand(data = pitfall_prep,
                            treat_f, year_f, taxa2, loc_f, plot_f = "dummy") %>%
  filter(year_f == "2014")
pitfall_preds <- predictSE.merMod(pitfall_lmm, 
                                  newdata = pitfall_nd, 
                                  REForm = NA)
pitfall_nd$fit <- pitfall_preds$fit
pitfall_nd$se <- pitfall_preds$se.fit

# offset treatment values by location for plotting
pitfall_prep <- pitfall_prep %>%
  mutate(treat_p = as.numeric(treat_f) + 0.4 * (as.numeric(loc_f) - 1.5))
pitfall_nd <- pitfall_nd %>%
  mutate(treat_p = as.numeric(treat_f) + 0.4 * (as.numeric(loc_f) - 1.5))

# plot
ggplot(data = pitfall_prep,
       aes(x = treat_p,
           y = zlog_cor,
           color = loc_f))+
  geom_hline(yintercept = 0, size = 0.2, color = "gray50")+
  facet_wrap(~taxa2, ncol = 1,strip.position = "right")+
  geom_rect(data = pitfall_nd,
            aes(xmin = treat_p - 0.2,
                xmax = treat_p + 0.2,
                ymin = fit - 1.96 * se,
                ymax = fit + 1.96* se,
                y = fit,
                fill = loc_f),
            alpha = 0.2,
            size = 0,
            linetype = 0)+
  geom_segment(data = pitfall_nd,
               aes(x = treat_p - 0.2,
                   xend = treat_p + 0.2,
                   y = fit,
                   yend = fit),
               size = 0.5)+
  geom_jitter(width = 0.2,
              alpha = 0.5,
              size = 1)+
  scale_color_manual("",values = c("firebrick",
                                   "dodgerblue"))+
  scale_fill_manual("",values = c("firebrick",
                                  "dodgerblue"))+
  scale_y_continuous("Log catch rate (z-scored)",
                     limits = c(-3.5,3.5), 
                     breaks = c(-2,0,2))+
  scale_x_continuous("Treatment",breaks = c(1,2,3),
                     labels = c("Control","KNO3","Midges"))+
  theme(legend.position = "top")
# ggsave(file = "figures/pit.pdf",width = 3.5, height = 7)



# PERMANOVA
pitfall_wide <- pitfall_prep %>%
  select(plot,year_f,taxa2, treat_f, loc_f, zlog) %>%
  pivot_wider(names_from = taxa2,
              values_from = zlog) 
pitfall_matrix <- pitfall_wide %>%
  select(unique(pitfall_prep$taxa2))
set.seed(1e2)
pitfall_euc <- adonis(pitfall_matrix ~ year_f + treat_f * loc_f,
                      data = pitfall_wide,
                      method = "euclidean",
                      permutations = 10000)
pitfall_euc$aov.tab

pitfall_wide <- pitfall_prep %>%
  select(plot,year_f,taxa2, treat_f, loc_f, count) %>%
  pivot_wider(names_from = taxa2,
              values_from = count) 
pitfall_matrix <- pitfall_wide %>%
  select(unique(pitfall_prep$taxa2))
set.seed(1e2)
pitfall_euc <- adonis(pitfall_matrix ~ year_f + treat_f * loc_f,
                      data = pitfall_wide,
                      method = "bray",
                      permutations = 10000)
pitfall_euc$aov.tab

# significant effects of year, treatment, and location; non-significant treat x location 
# confirms LMM for individual taxa responses

#====================================================================================











#====================================================================================
#========== Soil core sampling (soil-dwelling arthropods)
#====================================================================================

# prepare data
soil_prep <- soil %>% 
  filter(!is.na(location),
         location != "Lupine") %>%
  mutate(treat_f = factor(treatment),
         loc_f = factor(location),
         year_f = factor(year(coldate)),
         zlog = log1p(count),
         zlog = (zlog - mean(zlog, na.rm = T))/sd(zlog, na.rm = T))
soil_acar <- soil_prep %>% filter(taxa == "acar") 
soil_clmb <- soil_prep %>% filter(taxa == "clmb")

# plot acar
ggplot(data = soil_acar,
       aes(x = year_f,
           y = zlog,
           color = treat_f))+
  facet_wrap(~loc_f)+
  geom_jitter(width = 0.2)

# plot acar
ggplot(data = soil_clmb,
       aes(x = year_f,
           y = zlog,
           color = treat_f))+
  facet_wrap(~loc_f)+
  geom_jitter(width = 0.2)

# heteroscedasticity is a bit concerning

# fit model : acar
soil_acar_lmm <- lmer(zlog ~ year_f + treat_f * loc_f + (1 | plot),
                       data = soil_acar)

# p-values for effects on overall catch rate
# use F-test (reasonable for this design)
Anova(soil_acar_lmm, type = 3, test.statistic = "F") # full
Anova(update(soil_acar_lmm, .~. - treat_f:loc_f), type = 3, 
      test.statistic = "F") # drop interaction
# no clear patterns here

# fit model : clmb
soil_clmb_lmm <- lmer(zlog ~ year_f + treat_f * loc_f + (1 | plot),
                       data = soil_clmb) # plot variance is zero; probably not important for other results

# p-values for effects on overall catch rate
# use F-test (reasonable for this design)
Anova(soil_clmb_lmm, type = 3, test.statistic = "F") # full
Anova(update(soil_clmb_lmm, .~. - treat_f:loc_f), type = 3, 
      test.statistic = "F") # drop interaction
# year effects only; heteroscedasticity may contribute

# try without plot effect
# fit model : clmb
soil_clmb_lmm2 <- lm(zlog ~ year_f + treat_f * loc_f,
                      data = soil_clmb)

# p-values for effects on overall catch rate
# use F-test (reasonable for this design)
Anova(soil_clmb_lmm2, type = 3, test.statistic = "F") # full
Anova(update(soil_clmb_lmm2, .~. - treat_f:loc_f), type = 3, 
      test.statistic = "F") # drop interaction
# inference is not affected, so I would stick with th LMM for consistency

# standardize observations to 2014, using coefficient estimates from model
# acar
summary(soil_acar_lmm)$coefficients
acar_2018 <-  summary(soil_acar_lmm)$coefficients[2,1]
soil_acar <- soil_acar %>%
  mutate(zlog_cor = zlog - acar_2018 * (as.numeric(factor(year, levels = c(2014,2018))) - 1))
# clmb
summary(soil_clmb_lmm)$coefficients
clmb_2018 <-  summary(soil_clmb_lmm)$coefficients[2,1]
soil_clmb <- soil_clmb %>%
  mutate(zlog_cor = zlog - clmb_2018 * (as.numeric(factor(year, levels = c(2014,2018))) - 1))

# predictd values : acar
soil_acar_nd <- tidyr::expand(data = soil_acar,
                            treat_f, year_f, taxa, loc_f, plot_f = "dummy") %>%
  filter(year_f == "2014")
soil_acar_preds <- predictSE.merMod(soil_acar_lmm, 
                                  newdata = soil_acar_nd, 
                                  REForm = NA)
soil_acar_nd$fit <- soil_acar_preds$fit
soil_acar_nd$se <- soil_acar_preds$se.fit

# predictd values : clmb
soil_clmb_nd <- tidyr::expand(data = soil_clmb,
                              treat_f, year_f, taxa, loc_f, plot_f = "dummy") %>%
  filter(year_f == "2014")
soil_clmb_preds <- predictSE.merMod(soil_clmb_lmm, 
                                    newdata = soil_clmb_nd, 
                                    REForm = NA)
soil_clmb_nd$fit <- soil_clmb_preds$fit
soil_clmb_nd$se <- soil_clmb_preds$se.fit

# combine fitted values and plot
soil_nd <- bind_rows(soil_acar_nd,
                     soil_clmb_nd)

# offset treatment values by location for plotting
soil_plot <- bind_rows(soil_acar,
                       soil_clmb) %>%
  mutate(treat_p = as.numeric(treat_f) + 0.4 * (as.numeric(loc_f) - 1.5))
soil_nd <- soil_nd %>%
  mutate(treat_p = as.numeric(treat_f) + 0.4 * (as.numeric(loc_f) - 1.5))

# plot
ggplot(data = soil_plot,
       aes(x = treat_p,
           y = zlog_cor,
           color = loc_f))+
  geom_hline(yintercept = 0, size = 0.2, color = "gray50")+
  facet_wrap(~taxa, ncol = 1,strip.position = "right")+
  geom_rect(data = soil_nd,
            aes(xmin = treat_p - 0.2,
                xmax = treat_p + 0.2,
                ymin = fit - 1.96 * se,
                ymax = fit + 1.96* se,
                y = fit,
                fill = loc_f),
            alpha = 0.2,
            size = 0,
            linetype = 0)+
  geom_segment(data = soil_nd,
               aes(x = treat_p - 0.2,
                   xend = treat_p + 0.2,
                   y = fit,
                   yend = fit),
               size = 0.5)+
  geom_jitter(width = 0.2,
              alpha = 0.5,
              size = 1)+
  scale_color_manual("",values = c("firebrick",
                                   "dodgerblue"))+
  scale_fill_manual("",values = c("firebrick",
                                  "dodgerblue"))+
  scale_y_continuous("Log catch (z-scored)",
                     limits = c(-3.5,3.5), 
                     breaks = c(-2,0,2))+
  scale_x_continuous("Treatment",breaks = c(1,2,3),
                     labels = c("Control","KNO3","Midges"))+
  theme(legend.position = "top")
# ggsave(file = "figures/soil.pdf",width = 3.5, height = 5)


#====================================================================================










#====================================================================================
#========== Plants
#====================================================================================

# examine data
plants %>%
  group_by(taxa) %>%
  summarize(val = median(perc.cover)) %>%
  pivot_wider(names_from = taxa, values_from = val)

# extract taxa with median cover greater than zero (i.e., fewer than half of observations are zero)
plant_nonzero <- {plants %>%
  group_by(taxa) %>%
  summarize(val = median(perc.cover)) %>%
  filter(val > 0)}$taxa

# prepare data
plants_prep <- plants %>%
  filter(taxa %in% plant_nonzero) %>%
  rename(cover = perc.cover,
         n_int = no.intersections)%>%
  mutate(treat_f = factor(treatment),
         loc_f = factor(location),
         year_f = factor(year(coldate)))

# initial cover (2013)
plants_13 <- plants_prep %>%
  filter(year == 2013)

# plot
ggplot(data = plants_13,
       aes(x = treat_f,
           y = cover,
           color = loc_f))+
  facet_wrap(~taxa, ncol = 1)+
  geom_jitter(width = 0.2)+
  scale_color_manual("",values = c("firebrick",
                                   "dodgerblue"))

# fit model for initial cover
plants_13 <- plants_13 %>% 
  mutate(ztrans = asin(sqrt(cover)),
         ztrans = (ztrans - mean(ztrans, na.rm = T)) / sd(ztrans, na.rm = T))
plants_13_lm <- lm(ztrans ~ treat_f * loc_f * taxa,
                  data = plants_13)

# p-values for effects on overall catch rate
# use F-test (reasonable for this design)
Anova(plants_13_lm, type = 3, test.statistic = "F") # full
Anova(update(plants_13_lm, .~. - treat_f:loc_f:taxa), type = 3, 
      test.statistic = "F") # drop 3-way interaction
Anova(update(plants_13_lm, .~. - treat_f:loc_f:taxa - 
               treat_f:loc_f - treat_f:taxa - loc_f:taxa), 
      type = 3, test.statistic = "F") # drop all interactions
# taxa differ overall and among locations

# try glm
plants_13_glm <- glm(cbind(count, n_int - count) ~ treat_f * loc_f * taxa,
                    data = plants_13, family = "quasibinomial")
Anova(plants_13_glm, type = 3, test.statistic = "F") # full
Anova(update(plants_13_glm, .~. - treat_f:loc_f:taxa), type = 3, 
      test.statistic = "F") # drop 3-way interaction
Anova(update(plants_13_glm, .~. - treat_f:loc_f:taxa - 
               treat_f:loc_f - treat_f:taxa - loc_f:taxa), 
      test.statistic = "F",type = 3) # drop all interactions
# same result as before



# change in cover relative to 2013 and change between years
plants_change <- plants_prep %>%
  group_by(taxa, plot) %>%
  mutate(cover_0 = cover[1],
         delt_cover = cover - cover_0,
         diff_cover = c(0, diff(cover))) %>%
  ungroup()

# plot change relative to 2013 
ggplot(data = plants_change,
       aes(x = year,
           y = delt_cover,
           color = treat_f,
           group = plot))+
  ylab("Change in Plant Cover by Year")+
  theme(legend.title=element_blank())+
  facet_grid(taxa~loc_f)+
  geom_line()

# plot change between years
ggplot(data = plants_change,
       aes(x = year,
           y = diff_cover,
           color = treat_f,
           group = plot))+
  theme(legend.title=element_blank())+
  facet_grid(taxa~loc_f)+
  geom_line()

# uneven sampling between years complicates yearly analysis
# plots look quite consistent through time
# to me, cleanest route would be to analyze 2018 cover, with cover 2013 as a covariate
plants_18 <- plants_change %>%
  filter(year == 2018) %>%
  mutate(zcover = (cover - mean(cover, na.rm = T)) / sd(cover, na.rm = T),
         zcover0 = (cover_0 - mean(cover_0, na.rm = T)) / sd(cover_0, na.rm = T),
         id = row_number())

# fit model for initial cover
plants_18_glmm <- glmer(cbind(count, n_int - count) 
                          ~  zcover0 +  treat_f * loc_f * taxa + (1 | id),
                          data = plants_18,
                          family = "binomial",
                        control=glmerControl(
                          optimizer="bobyqa", 
                          optCtrl=list(maxfun=2e6)))
summary(plants_18_glmm)
# convergence warning; try quasibinomial

plants_18_glm <- glm(cbind(count, n_int - count) 
                        ~ zcover0 +  treat_f * loc_f * taxa,
                        data = plants_18,
                        family = "quasibinomial")
summary(plants_18_glm)

# p-values with F-test
Anova(plants_18_glm, type = 3, test.statistic = "F") # full
Anova(update(plants_18_glm, .~. - treat_f:taxa:loc_f), 
      test.statistic = "F") # 3-way interaction dropped
Anova(update(plants_18_glm, .~. - year_f:taxa  - 
               taxa:treat_f - 
               treat_f:loc_f - 
               taxa:loc_f - 
               treat_f:taxa:loc_f), 
      test.statistic = "F") # all interactions dropped
# taxa differ in responses to treatment and among locations, but not by their interaction



# remove effect of initial cover
summary(plants_18_glm)$coefficients
zcover_e <- summary(plants_18_glm)$coefficients["zcover0",1]
plants_18 <- plants_18 %>%
  mutate(cover_cor = log(cover / (1 - cover)) - zcover_e * zcover0,
         cover_cor = exp(cover_cor) / (1 + exp(cover_cor)))

# generate predicted values
plants_18_nd <- tidyr::expand(data = plants_18,
                              treat_f, taxa, loc_f, zcover0 = 0)
plants_18_preds <- predict(plants_18_glm,
                           newdata = plants_18_nd,
                           se.fit = T,
                           type = "link")
plants_18_nd <- plants_18_nd %>%
  mutate(fit = plants_18_preds$fit,
         se = plants_18_preds$se.fit,
         lo = fit - 1.96 * se,
         hi = fit + 1.96 * se,
         fit = exp(fit) / (1 + exp(fit)),
         lo = exp(lo) / (1 + exp(lo)),
         hi = exp(hi) / (1 + exp(hi)))


plants_18_nd$se_link <- plants_18_preds$se.fit
plants_18_nd$fit


# offset treatment values by location for plotting
plants_18 <- plants_18 %>%
  mutate(treat_p = as.numeric(treat_f) + 0.4 * (as.numeric(loc_f) - 1.5))
plants_18_nd <- plants_18_nd %>%
  mutate(treat_p = as.numeric(treat_f) + 0.4 * (as.numeric(loc_f) - 1.5))

# plot
ggplot(data = plants_18,
       aes(x = treat_p,
           y = cover,
           color = loc_f))+
  geom_hline(yintercept = 0, size = 0.2, color = "gray50")+
  facet_wrap(~taxa, ncol = 1,strip.position = "right")+
  geom_rect(data = plants_18_nd,
            aes(xmin = treat_p - 0.2,
                xmax = treat_p + 0.2,
                ymin = lo,
                ymax = hi,
                y = fit,
                fill = loc_f),
            alpha = 0.2,
            size = 0,
            linetype = 0)+
  geom_segment(data = plants_18_nd,
               aes(x = treat_p - 0.2,
                   xend = treat_p + 0.2,
                   y = fit,
                   yend = fit),
               size = 0.5)+
  geom_jitter(width = 0.2,
              alpha = 0.5,
              size = 1)+
  scale_color_manual("",values = c("firebrick",
                                   "dodgerblue"))+
  scale_fill_manual("",values = c("firebrick",
                                  "dodgerblue"))+
  scale_y_continuous("Proportion cover (z-scored)",
                     limits = c(0,1),
                     breaks = c(0.2,0.5,0.8))+
  scale_x_continuous("Treatment",breaks = c(1,2,3),
                     labels = c("Control","KNO3","Midges"))+
  theme(legend.position = "top")
# ggsave(file = "figures/plant.pdf",width = 3.5, height = 7)

# PERMANOVA
plants_18_wide <- plants_18 %>%
  mutate(zcover = asin(sqrt(cover)),
         zcover = (zcover - mean(zcover, na.rm = T)) / sd(zcover, na.rm = T)) %>%
  select(plot,taxa, treat_f, loc_f, zcover) %>%
  pivot_wider(names_from = taxa,
              values_from = zcover) 
plants_18_matrix <- plants_18_wide %>%
  select(unique(plants_18$taxa))
set.seed(1e2)
plants_18_euc <- adonis(plants_18_matrix ~  treat_f * loc_f,
                      data = plants_18_wide,
                      method = "euclidean",
                      permutations = 10000)
plants_18_euc$aov.tab

#====================================================================================











#====================================================================================
#========== Pitfall & Plants
#====================================================================================

plants_agg <- plants_prep %>%
  mutate(zcover = asin(sqrt(cover)),
         zcover = (cover - mean(cover)) / sd(cover)) %>%
  group_by(block, plot, treat_f, loc_f, taxa) %>%
  summarize(zcover = mean(zcover)) %>%
  pivot_wider(names_from = taxa,
              values_from = zcover) %>%
  ungroup()

pca <- prcomp(plants_agg[,c("BENA","GRASS","LUNO","SAPH","VAUL")])
pca
summary(pca)

plants_agg$pc1 <- pca$x[,1]
plants_agg$pc2 <- pca$x[,2]
plants_agg$pc3 <- pca$x[,3]

pitfall_plants <- pitfall_prep %>%
  group_by(block, plot, treat_f, loc_f, taxa2) %>%
  summarize(zlog = mean(zlog)) %>%
  ungroup() %>%
  left_join(plants_agg)

pitfall_plant_mod <- lm(zlog ~ taxa2 + taxa2*pc1 + taxa2*pc2 + taxa2*pc3,
                        data = pitfall_plants)
Anova(pitfall_plant_mod, test.statistic = "F", type = 3)
Anova(update(pitfall_plant_mod, .~. - taxa2:pc1 - taxa2:pc2  - taxa2:pc3), 
      test.statistic = "F", type = 3)
summary(pitfall_plant_mod)

pitfall_plants_nd_a <- tidyr::expand(data = pitfall_plants,
                                   taxa2 = taxa2,
                                   pc1 = seq(min(pc1), max(pc1), length.out = 10),
                                   pc2 = mean(pc2),
                                   pc3 = mean(pc3),
                                   type = "PC1")
pred_a <- predict(pitfall_plant_mod, newdata = pitfall_plants_nd_a, se.fit = T)
pitfall_plants_nd_a$fit <- pred_a$fit
pitfall_plants_nd_a$se <- pred_a$se.fit
pitfall_plants_nd_b <- tidyr::expand(data = pitfall_plants,
                                     taxa2 = taxa2,
                                     pc2 = seq(min(pc2), max(pc2), length.out = 10),
                                     pc1 = mean(pc1),
                                     pc3 = mean(pc3),
                                     type = "PC2")
pred_b <- predict(pitfall_plant_mod, newdata = pitfall_plants_nd_b, se.fit = T)
pitfall_plants_nd_b$fit <- pred_b$fit
pitfall_plants_nd_b$se <- pred_b$se.fit
pitfall_plants_nd_c <- tidyr::expand(data = pitfall_plants,
                                     taxa2 = taxa2,
                                     pc3 = seq(min(pc3), max(pc3), length.out = 10),
                                     pc2 = mean(pc2),
                                     pc1 = mean(pc1),
                                     type = "PC3")
pred_c <- predict(pitfall_plant_mod, newdata = pitfall_plants_nd_c, se.fit = T)
pitfall_plants_nd_c$fit <- pred_c$fit
pitfall_plants_nd_c$se <- pred_c$se.fit

pitfall_plants_nd <- bind_rows(pitfall_plants_nd_a, 
                               pitfall_plants_nd_b,
                               pitfall_plants_nd_c) 


p1<- ggplot(data = pitfall_plants,
            aes(x = pc1,
                y = zlog,
                color = taxa2))+
  geom_ribbon(data = pitfall_plants_nd %>% filter(type == "PC1"),
              aes(y = fit,
                  ymin = fit - 1.96 * se,
                  ymax = fit + 1.96 * se,
                  fill = taxa2),
              alpha = 0.1,
              linetype = 0)+
  geom_point(size = 1)+
  geom_line(data = pitfall_plants_nd %>% filter(type == "PC1"),
            aes(y = fit))+
  scale_color_viridis_d("", guide = F)+
  scale_fill_viridis_d("", guide = F)+
  scale_y_continuous("")+
  scale_x_continuous("PC1")+
  theme(plot.margin = margin(t = 1,
                             b = 10,
                             l = 1,
                             r = 1))
p1


p2 <- ggplot(data = pitfall_plants,
            aes(x = pc2,
                y = zlog,
                color = taxa2))+
  geom_ribbon(data = pitfall_plants_nd %>% filter(type == "PC2"),
              aes(y = fit,
                  ymin = fit - 1.96 * se,
                  ymax = fit + 1.96 * se,
                  fill = taxa2),
              alpha = 0.1,
              linetype = 0)+
  geom_point(size = 1)+
  geom_line(data = pitfall_plants_nd %>% filter(type == "PC2"),
            aes(y = fit))+
  scale_color_viridis_d("")+
  scale_fill_viridis_d("")+
  scale_y_continuous("Catch rate (z-scored)")+
  scale_x_continuous("PC2")+
  theme(plot.margin = margin(t = 1,
                             b = 10,
                             l = 1,
                             r = 1))
p2

p3 <-ggplot(data = pitfall_plants,
            aes(x = pc3,
                y = zlog,
                color = taxa2))+
  geom_ribbon(data = pitfall_plants_nd %>% filter(type == "PC3"),
              aes(y = fit,
                  ymin = fit - 1.96 * se,
                  ymax = fit + 1.96 * se,
                  fill = taxa2),
              alpha = 0.1,
              linetype = 0)+
  geom_point(size = 1)+
  geom_line(data = pitfall_plants_nd %>% filter(type == "PC3"),
            aes(y = fit),
            size = 0.4)+
  scale_color_viridis_d("", guide = F)+
  scale_fill_viridis_d("", guide = F)+
  scale_y_continuous("")+
  scale_x_continuous("PC3")+
  theme(plot.margin = margin(t = 1,
                             b = 10,
                             l = 1,
                             r = 1))
p3


# combine
p4 <- plot_grid(p1, p2, p3,
                nrow = 3,
                align = "v",
                axis = "tblr")

p4
# ggsave(file = "figures/pit_plants.pdf",width = 3.5, height = 5)
