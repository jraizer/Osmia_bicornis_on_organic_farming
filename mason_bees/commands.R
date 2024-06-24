# Header----
# Manuscript draft
# Provisory tittle: Conventional agriculture interferes with sex communication in a wild bee with negative impacts on local population size

# Authors: Samuel Boff1*, Sara Olberz1, Irem Gülsoy2, Marvin Preuß1,  Josué Raizer3, Manfred Ayasse1

# Affiliations:
#1 Institute of Evolutionary Ecology and Conservation Genomics, Ulm University, Ulm, Germany
#2 Department of Molecular Biology and Genetics, İhsan Doğramacı Bilkent University, Ankara, Turkey
#3 Graduate Program in Entomology and Biodiversity Conservation, Federal University of Grande Dourados, Dourados, Brazil

# Correspondence
#*Email: samuel.boff@uni-ulm.de
#Phone number: +49 (0)731 50 22665

# updated 2024, Jun, 24
# R version 4.4.1

# Required package----
library(readxl)
library(lme4)
library(glmmTMB)
library(car)
library(emmeans)
library(effects)
library(DHARMa)
library(pwr)
library(vegan)
library(kableExtra)
library(ggplot2)

## Data----
### Import
bees <- readxl::read_excel("data/dataset.xlsx") #data per bee
bees <- bees[-which(is.na(bees$Organicity_rev) == T), ] #sites excluded

reprod <- readxl::read_excel("data/reprod.xlsx") #additional reproduction data
meadow <- readxl::read_excel("data/meadow.xlsx") #flower abundance in meadows
exp_1 <- readxl::read_excel("data/experiment_1.xlsx") #data from experiment 1

### Inspect
str(bees)
str(reprod)
str(meadow)
str(exp_1)

unique(bees$Farm) # 11 farms
unique(bees$Site) # 15 sites
unique(bees$Station) # 33 stations

table(bees$Farm)
table(bees$Site)
table(bees$Station)

rowSums(table(bees$Farm, bees$Site) > 0)
rowSums(table(bees$Farm, bees$Station) > 0)
rowSums(table(bees$Site, bees$Station) > 0)

### Subset: traits of bees and sampling information 
envi <- bees[, c(2:9)] 
str(envi)
#Station          - station identifier
#Site             - site identifier;
#Size             - bee size (width between eyes in mm);
#sex              - bee sex: "Male", "Female";
#Farm             - farm identifier;
#Farm_system      - farming management: "Conventional", "Organic";
#Organicity_rev   - extent of organic farming (buffer 500 m radius);
#Forest           - forest cover (buffer 500 m radius).

### Subset: CHC compounds
compounds <- bees[, 10:30] # 21 cuticular hydrocarbons
compounds_names <- colnames(compounds)
colnames(compounds) <- #simplify names to insert in the figures
  c("I", "Henei", "Do", "Tri", "9_Tetra", "Tetra", "11_Penta", "9_Penta", 
    "7_Penta", "Penta", "Hexa", "11_Hepta", "9_Hepta", "7_Hepta", "Hepta",
    "Octa", "9_Nona", "7_Nona", "Nona", "Hentri", "Dotri")

comp_rel <- vegan::decostand(compounds, "total") #relative amount of compounds

### Subsets females
envi_f <- envi[which(envi$sex == "Female"), ]
female <- comp_rel[which(envi$sex == "Female"), ]
### Subsets males
envi_m <- envi[which(envi$sex == "Male"), ]
male <- comp_rel[which(envi$sex == "Male"), ]

# Flower abundance in the near of bee hotels----
## Inspect
str(meadow)
stripchart(meadow$flower_abund, method = "stack")
stripchart(meadow$flower_abund ~ meadow$Farm_system, vertical = T,
           method = "stack")
boxplot(meadow$flower_abund ~meadow$Farm_system)

### Basic stats
aggregate(meadow$flower_abund, list(meadow$Farm_system), mean, na.rm = T)
aggregate(meadow$flower_abund, list(meadow$Farm_system), sd, na.rm = T)
aggregate(meadow$flower_abund, list(meadow$Farm_system), quantile, na.rm = T)

## Meadows in bubble map----

# Sort the dataset by Site_1
meadow <- meadow[order(meadow$Site_1), ]

# Ensure the Site_1 column is treated as a factor with levels in the reverse order
meadow$Site_1 <- factor(meadow$Site_1, levels = rev(unique(meadow$Site_1)))

# Plot the bubble map
ggplot(meadow, aes(x = Farm_system, y = Site_1, size = flower_abund, 
                 fill = Farm_system)) +
  geom_point(alpha = 0.6, shape = 21, color = "black") +
  scale_size(range = c(3, 15)) +  # Adjust the range of bubble sizes
  scale_fill_manual(values = c("Conventional" = "orange", 
                               "Organic" = "skyblue")) +
  labs(title = "Flower abundance in conventional and organic meadows",
       x = "Farming management",
       y = "Farming sites",
       size = "Flower abundance") +
  theme_minimal() +
  theme(legend.position = "right") +
  guides(fill = "none")  + # Remove a legenda de fill (Farm_system)
  coord_cartesian(xlim = c(0.5, 
                           length(unique(meadow$Farm_system)) + 0.5),
                  ylim = c(0.5, 
                           length(unique(meadow$Site_1)) + 0.5))

## Save as jpeg:
ggsave("figures/figure_1.jpg", width = 17, height = 13, units = "cm", dpi = 300)
## Save as pdf:
ggsave("figures/figure_1.pdf", width = 17, height = 13, units = "cm", dpi = 300)
## Save as tiff:
ggsave("figures/figure_1.tiff", width = 17, height = 13, units = "cm", dpi = 300)
## Save as png:
ggsave("figures/figure_1.png", width = 17, height = 13, units = "cm", dpi = 300)

## Linear model----
flower_model <- lm(flower_abund ~ Farm_system, data = meadow)
anova(flower_model)

# Nest occupation and brood cells----
reprod <- reprod[-which(is.na(reprod$Organicity_rev) == T), ] #three sites excluded
str(reprod)

## Nest frequency----
reprod$nest_corrected == reprod$Nest_per_station / reprod$Total_Reeds
stripchart(round(reprod$nest_corrected, 3), method = "stack")

# Test of zero proportion
zero_norm <- dnorm(0, mean(reprod$nest_corrected), 
                   sd(reprod$nest_corrected))/100
prop_test <- prop.test(sum(reprod$nest_corrected == 0), nrow(reprod), 
                       p = zero_norm, alternative = "greater")
prop_test

## Power test (Type II Error)
obs_prop <- sum(reprod$nest_corrected == 0) / nrow(reprod)

# Effect size
h <- ES.h(obs_prop, p2 = zero_norm)

# Power
pwr.p.test(h = h, n = nrow(reprod), sig.level = 0.05,
                           alternative = "greater")

#Although the test did not indicate a significant difference at the 0.05 significance level (prop_test: p = 0.06), the low power of the test (0.49) suggests that we may not have enough data to detect zero inflation, if it exists. Given the potential consequences of failing to account for zero inflation in subsequent modeling, we consider it prudent to proceed under the assumption that the data may be zero inflated.

### Generalized linear mixed-effect model----
nest_freq_model <- glmmTMB::glmmTMB(nest_corrected ~ Farm_system + 
                                      Organicity_rev +  Forest  +  
                                      (1 | Site), 
                                    data = reprod, family = tweedie())

#### Assumptions test (DHARMa) of the model----
# Residuals simulation
resid_sim_freq <- DHARMa::simulateResiduals(fittedModel = nest_freq_model, 
                                            n = 1000) 
plot(resid_sim_freq) #no violations

#### Model selection----

# The competing models were compared using AIC, AICc, BIC, and LogLik.
source("functions/backward_selection.R")
freq_results <- backward_selection(nest_freq_model)
freq_results
#AIC differences < 4: stay initial model

### Global significance of the model----
freq_null_mod <- glmmTMB::glmmTMB(nest_corrected ~ 1 +
                                    (1 | Site), 
                                  data = reprod, family = tweedie())
anova(freq_null_mod, nest_freq_model) 

### Hypothesis test----
summary(nest_freq_model)
car::Anova(nest_freq_model)

### Plot----
scatter <- ggplot(reprod_no_na, 
                  aes(Organicity_rev, nest_corrected)) + 
  theme_classic(30) + 
  geom_point(size = 3)
scatter + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  labs(x = "Extent of organic farming (%)", y = "Nesting frequency")


## Brood cells----

# Inspect
plot(reprod$nest_corrected, reprod$broodcells)
cor.test(reprod$broodcells, reprod$nest_corrected) #r = 0.96

stripchart(reprod$broodcells, method = "stack")

# Test of zero proportion
reprod_no_na <- reprod[complete.cases(reprod$broodcells), ]
zero_pois <- exp(-mean(reprod_no_na$broodcells))
prop_test <- prop.test(sum(reprod_no_na$broodcells == 0), nrow(reprod_no_na), 
                       p = zero_pois, alternative = "greater")
prop_test

#zero inflated

### Generalized linear mixed-effect model----
brood_model <- glmmTMB::glmmTMB(broodcells ~ Farm_system + 
                                  Organicity_rev +  
                                  Forest  +  
                                  (1 | Site ), ziformula = ~1,
                                data = reprod_no_na, 
                                family = tweedie)

#### Assumptions test (DHARMa) of the model----
# Residuals simulation
resid_sim_brood <- DHARMa::simulateResiduals(fittedModel = brood_model, 
                                             n = 1000) 
plot(resid_sim_brood) #no violations

#### Model selection----

# The competing models were compared using AIC, AICc, BIC, and LogLik.
brood_results <- backward_selection(brood_model)
brood_results
#AIC differences < 4: stay initial model

#### Global significance of the model----
brood_null_mod <- glmmTMB::glmmTMB(broodcells ~ 1  +  
                                     (1 | Site ), ziformula = ~1,
                                   data = reprod_no_na, 
                                   family = tweedie)
anova(brood_null_mod, brood_model) 

### Hypothesis test----
summary(brood_model)
car::Anova(brood_model)

### Plot----
scatter <- ggplot(reprod_no_na, 
                  aes(Organicity_rev, broodcells)) + 
  theme_classic(30) + 
  geom_point(size = 3)
scatter + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  labs(x = "Extent of organic farming (%)", y = "Number of broodcells")


# Size of bees----
## Inspect
stripchart(round(envi$Size, 1), at = 0,
           method = "stack", offset = .25, xpd = T)
stripchart(round(envi$Size, 1) ~
             envi$sex, at = c(1, 2),
           method = "stack", offset = .25, xpd = T)
stripchart(round(envi$Size, 1) ~
             envi$Farm_system, at = c(1, 2),
           method = "stack", offset = .25, xpd = T)
stripchart(round(envi$Size, 1) ~ 
             paste(envi$sex, envi$Farm_system), 
           vertical = T,
           method = "stack", xpd = T, cex = 1.25,
           pch = c("\u2640", "\u2640", "\u2642", "\u2642"),
           frame.plot = F, 
           ylab = "Bee size (width between eyes in mm)",
           xaxt = "n", col = c("orange", "skyblue"),
           family = "Arial Unicode MS")
axis(1, at = c(1.5, 3.5), labels = c("Female", "Male"), )
mtext("Sex", 1, at = 2.5, line = 2)
legend(3.5, 2.6, c("Conventional", "Organic"), 
       pch = 21, bty = "n", col = c("orange", "skyblue"), 
       title = "Faming system", pt.cex = 1.5)

### Basic stats
aggregate(envi$Size, list(envi$Farm_system, envi$sex), mean, na.rm = T)
aggregate(envi$Size, list(envi$Farm_system, envi$sex), sd, na.rm = T)
aggregate(envi$Size, list(envi$Farm_system, envi$sex), quantile, na.rm = T)

### Correlations organicity and farm system
stripchart(envi$Organicity_rev ~ envi$Farm_system, vertical = T,
           method = "stack", offset = 0.25, xpd = T)
boxplot(envi$Organicity_rev ~ envi$Farm_system)

### Correlations forest and farm system
stripchart(envi$Forest ~ envi$Farm_system, vertical = T,
           method = "stack", offset = 0.25, xpd = T)
boxplot(envi$Forest ~ envi$Farm_system)

### Correlations organicity and forest
org_site <- aggregate(envi$Organicity_rev, list(envi$Site), mean)$x
for_site <- aggregate(envi$Forest, list(envi$Site), mean)$x
plot(org_site ~ for_site, xlab = "Forest by site", 
     ylab = "Organicity by site")
cor.test(org_site, for_site) #uncorrelated

## Linear mixed-effect model----
envi_non_na <- envi[complete.cases(envi$Size), ] #exclude NA cases on bee size
size_model <- lme4::lmer(Size ~ Farm_system * sex + Organicity_rev + Forest +
                     (1 | Site/Station), 
                   data = envi_non_na)

# Model assumptions teste (DHARMa)
# Residuals simulation
resid_sim <- DHARMa::simulateResiduals(fittedModel = size_model, n = 1000) 

plot(resid_sim) #no violations of model assumptions
DHARMa::outliers(resid_sim)

### Model selection----

# The competing models were compared using AIC, AICc, BIC, and LogLik.
size_results <- backward_selection(size_model)
size_sel_mod <- size_results$final_model #selected model
size_info <- size_results$model_info
size_info

# Export size_info as HTML
size_info |>
  kableExtra::kable("html", caption = "Model selection") |>
  kableExtra::kable_styling(bootstrap_options = c("striped", "hover",
                                                  "condensed")) |>
  kableExtra::save_kable("size_models.html")
# See competing models in "size_models.html" 

#### Assumptions test (DHARMa) of the selected model----
# Residuals simulation
resid_sim_sel <- DHARMa::simulateResiduals(fittedModel = size_sel_mod, 
                                           n = 1000) 

DHARMa::plotQQunif(resid_sim_sel) #assumptions okay
DHARMa::plotResiduals(resid_sim_sel) #heterocedasticity detected but:
car::leveneTest(residuals(size_sel_mod) ~ Farm_system * sex, 
                data = envi_non_na) 
DHARMa::testDispersion(resid_sim_sel)
DHARMa::outliers(resid_sim_sel)


### Global significance of the selected model----
size_null_mod <- lme4::lmer(Size ~ 1 +
                        (1 | Site/Station), 
                      data = envi_non_na)
anova(size_null_mod, size_sel_mod)

### Hypothesis test----
summary(size_sel_mod)
car::Anova(size_sel_mod)

# Estimated marginal means (Least-squares means)
(emm <- emmeans::emmeans(size_sel_mod, ~ Farm_system + sex))

# Contrast
pairs(emm)
plot(emm)

### Plot----
farm_eff <- effects::allEffects(size_sel_mod)[["Farm_system"]]
farm_eff$x
farm_eff$fit
farm_eff$lower
farm_eff$upper
plot(farm_eff)

sex_eff <- effects::allEffects(size_sel_mod)[["sex"]]
sex_eff$x
sex_eff$fit
sex_eff$lower
sex_eff$upper
plot(sex_eff)

### See "figure_2_alt.R"
### or "figure_2.R" (Figure 2 in the paper)


# Cuticular hydrocarbon composition----
## CHC amount----
envi_non_na$amountCHC <- rowSums(compounds[complete.cases(envi$Size), ])
plot(envi_non_na$Size, envi_non_na$amountCHC)
abline(lm(amountCHC ~ Size, data = envi_non_na))
cor.test(envi_non_na$Size, envi_non_na$amountCHC)
### ratio amountCHC by size
chc_by_size <- envi_non_na$amountCHC / envi_non_na$Size

### Inspect
stripchart(round(chc_by_size), 
           method = "stack", xpd = T)
stripchart(round(chc_by_size) ~
             envi_non_na$sex, at = c(1, 2),
           method = "stack", xpd = T)
stripchart(round(chc_by_size) ~
             envi_non_na$Farm_system, at = c(1, 2),
           method = "stack", xpd = T)
op <- par()
par(mar = c(5, 6, 2, 1))
stripchart(round(chc_by_size) ~ 
             paste(envi_non_na$sex, envi_non_na$Farm_system), 
           vertical = T,
           method = "stack", xpd = T, cex = 1.25,
           pch = c("\u2640", "\u2640", "\u2642", "\u2642"),
           frame.plot = F, 
           ylab = expression(
             frac("Amount of cuticular hidrocarbon production (ng)", 
                  "Bee size (width between eyes in mm)")),
           xaxt = "n", col = c("orange", "skyblue"))
axis(1, at = c(1.5, 3.5), labels = c("Female", "Male"), )
mtext("Sex", 1, at = 2.5, line = 2)
legend(3.5, 35, c("Conventional", "Organic"), 
       pch = 21, bty = "n", col = c("orange", "skyblue"), 
       title = "Faming system", pt.cex = 1.5)
par(op)

### Basic stats
aggregate(envi_non_na$amountCHC, 
          list(envi_non_na$Farm_system, envi_non_na$sex), 
          mean, na.rm = T)
aggregate(envi_non_na$amountCHC, 
          list(envi_non_na$Farm_system, envi_non_na$sex), 
          sd, na.rm = T)
aggregate(envi_non_na$amountCHC, 
          list(envi_non_na$Farm_system, envi_non_na$sex), 
          quantile, na.rm = T)

#### Linear mixed-effect model----
# The response variable, which is the ratio of amountCHC to Size, inherently reflects the significant correlation between amountCHC and Size (r = 0.56, p < 0.001). Given this relationship, it is methodologically sound to utilize the predictor variables identified in the previously selected model for Size when modeling this ratio. This approach ensures that the established predictor structure is maintained and eliminates the need for redundant model selection procedures. By including Farm_system and sex as predictors, the model can adequately capture the relationship between the amountCHC/Size ratio and the influencing factors, leveraging the robustness of the previously validated model. This strategy is both efficient and scientifically justified, as it builds upon the existing understanding of the relationships among the variables.

amountCHC_model <- lme4::lmer(chc_by_size ~ Farm_system + sex +
                           (1 | Site/Station), 
                         data = envi_non_na)
# Model assumptions teste(DHARMa)
# Residuals simulation
resid_chc_sim <- DHARMa::simulateResiduals(fittedModel = amountCHC_model, 
                                           n = 1000) 

plot(resid_chc_sim) #no violations of model assumptions

#### Global significance of the model----
chc_null_mod <- lme4::lmer(chc_by_size ~ 1 +
                        (1 | Site/Station), 
                      data = envi_non_na)
anova(chc_null_mod, amountCHC_model) #non significant

## NMDS to CHC----
### All bees----
(ord <- vegan::metaMDS(comp_rel, trymax = 500)) #NMDS from vegan package
vegan::stressplot(ord) #Shepard diagram 
loads_tot <- ord$species #loadings

#### Permutational test----
vegan::envfit(ord ~ sex + Farm_system,
              data = envi, na.rm = T, perm = 1000)

#### Plot----
plot(ord, type = "n", bty = "n", ylim = c(min(ord$points[, 2]) * 1.1, 
                                          max(ord$points[, 2]) * 1.1),
     xlim = c(min(ord$points[, 1]) * 1.1, 
              max(ord$points[, 1]) * 1.1))
with(envi, 
     points(ord, 
            display = "sites", 
            pch = ifelse(sex == "Female", -0x2640L, -0x2642L), 
            cex = 1.5,
            col = ifelse(Farm_system == "Organic", 
                         "skyblue", 
                         "orange")))
with(envi, ordispider(ord$points, 
                      Farm_system, col = c("orange", "skyblue")))
with(envi, 
     ordiellipse(ord, 
                 Farm_system, 
                 kind = "ehull", 
                 lwd = 2,
                 draw = "polygon", 
                 col = c("orange", "skyblue"), 
                 border = c("orange", "skyblue"),
                 alpha = 23,
                 xpd = T))
with(comp_rel,
     text(ord,
          display = "species",
          cex = .75,
          col = "gray20",
          xpd = T))
abline(h = 0, v = 0, lty = 2, col = "gray")
legend("bottomright", xpd = T, 
       legend = c("Organic", "Conventional"),
       lwd = 1.5,
       col = c("skyblue", "orange"),
       bty = "n",
       title = "Farming System")

##### jpeg file----
jpeg("figures/fig3_nmds_total.jpg", width = 15, height = 15, units = "cm", 
     res = 300)
plot(ord$points, 
     type = "n", 
     bty = "n", 
     xlim = c(min(ord$points[, 1]) * 1.1, 
              max(ord$points[, 1]) * 1.1),
     ylim = c(min(ord$points[, 2]) * 1.1, 
              max(ord$points[, 2]) * 1.1),
     xlab = "NMDS 1",
     ylab = "NMDS 2")
with(envi, 
     points(ord, 
            display = "sites", 
            pch = ifelse(sex == "Female", -0x2640L, -0x2642L), 
            cex = 1.5,
            col = ifelse(Farm_system == "Organic", 
                         "skyblue", 
                         "orange"),
            xpd = T))
with(envi, ordispider(ord$points, 
                      Farm_system, col = c("orange", "skyblue")))
with(envi, 
     ordiellipse(ord, 
                 Farm_system, 
                 kind = "ehull", 
                 lwd = 2,
                 draw = "polygon", 
                 col = c("orange", "skyblue"), 
                 border = c("orange", "skyblue"),
                 alpha = 23,
                 xpd = T))
with(comp_rel,
     text(ord,
          display = "species",
          cex = .75,
          col = "gray20",
          xpd = T))
legend(min(ord$points[, 1]) * 1.5, max(ord$points[, 2]) * 1.5, xpd = T, 
       legend = c("Organic", "Conventional"),
       lwd = 1.5,
       col = c("skyblue", "orange"),
       bty = "n",
       title = "Farming System",
       cex = .9)
dev.off()


### By sex----
#### Females NMDS----
ord_f <- vegan::metaMDS(female, trymax = 500)
vegan::stressplot(ord_f)
loads_f <- ord_f$species

##### Permutational test----
vegan::envfit(ord_f ~ Farm_system, 
              data = envi_f, na.rm = T, perm = 5000)

##### Plot----
plot(ord_f, type = "n", bty = "n")
with(envi_f, 
     points(ord_f, 
            display = "sites", 
            pch = -0x2640L, 
            cex = 1.5,
            col = ifelse(Farm_system == "Organic", 
                         "skyblue", 
                         "orange")))
with(envi_f, ordispider(ord_f$points, 
                        Farm_system, col = c("orange", "skyblue")))
with(envi_f, ordiellipse(ord_f$points, 
                         Farm_system, 
                         kind = "ehull", 
                         lwd = 2,
                         draw = "polygon", 
                         col = c("orange", "skyblue"), 
                         border = c("orange", "skyblue"),
                         alpha = 23,
                         xpd = T))
with(female,
     text(ord_f,
          display = "species",
          cex = .75,
          col = "gray20",
          xpd = T))
legend(0.2, 0.6, xpd = T, 
       legend = c("Organic", "Conventional"),
       lwd = 1.5,
       col = c("skyblue", "orange"),
       bty = "n",
       title = "Farming System")

###### jpeg figure----
jpeg("figures/fig4_nmds_female.jpg", width = 15, height = 15, units = "cm", 
     res = 300)
plot(ord_f$points, 
     type = "n", 
     bty = "n", 
     xlim = c(min(ord_f$points[, 1]) * 1.1, 
              max(ord_f$points[, 1]) * 1.1), 
     ylim = c(min(ord_f$points[, 2]) * 1.25, 
              max(ord_f$points[, 2])), 
     xlab = "NMDS 1",
     ylab = "NMDS 2")
with(envi_f, 
     points(ord_f, 
            display = "sites", 
            pch = -0x2640L, 
            cex = 1.5,
            col = ifelse(Farm_system == "Organic", 
                         "skyblue", 
                         "orange")))
with(envi_f, ordispider(ord_f$points, 
                        Farm_system, col = c("orange", "skyblue")))
with(envi_f, ordiellipse(ord_f, 
                         Farm_system, 
                         kind = "ehull", 
                         lwd = 2,
                         draw = "polygon", 
                         col = c("orange", "skyblue"), 
                         border = c("orange", "skyblue"),
                         alpha = 23,
                         xpd = T))
with(female,
     text(ord_f,
          display = "species",
          cex = .75,
          col = "gray20",
          xpd = T))
legend("topleft", 
       legend = c("Organic", "Conventional"),
       lwd = 1.5,
       col = c("skyblue", "orange"),
       bty = "n",
       title = "Farming System",
       cex = 0.75)
dev.off()


#### Males NMDS----
(ord_m <- vegan::metaMDS(male, trymax = 500))
vegan::stressplot(ord_m)
plot(ord_m) #Here, we see a bee that differs from others due to a specific cuticular hydrocarbon. Consequently, we will redo the NMDS analysis excluding this particular bee.

which(ord_m$points[, 1] == max(ord_m$points[, 1])) #bee identified
male <- male[-which(ord_m$points[, 1] == max(ord_m$points[, 1])), ]
envi_m <- envi_m[-which(ord_m$points[, 1] == max(ord_m$points[, 1])), ]
(ord_m <- vegan::metaMDS(male, trymax = 500))
vegan::stressplot(ord_m)
vegan::envfit(ord_m ~ Farm_system, 
       data = envi_m,
       na.rm = T, perm = 5000)

#### Plot----
plot(ord_m, type = "n")
with(envi_m, 
     points(ord_m, 
            display = "sites", 
            col = ifelse(Farm_system == "Organic", "skyblue", "orange"), 
            cex = 1.5, 
            pch = -0x2642L))
with(envi_m, ordispider(ord_m, 
                        Farm_system, col = c("orange", "skyblue")))
with(envi_m, ordiellipse(ord_m, 
                         Farm_system, 
                         kind = "ehull", 
                         lwd = 2,
                         draw = "polygon", 
                         col = c("orange", "skyblue"), 
                         border = c("orange", "skyblue"),
                         alpha = 23,
                         xpd = T))
with(envi_m, 
     text(ord_m, 
          display = "species", 
          col = "gray20", 
          cex = 1,
          xpd = T))
legend("bottomleft", 
       legend = c("Conventional", "Organic"),
       lwd = 1.5,
       col = c("orange", "skyblue"),
       bty = "n",
       title = "Farming System")

##### jpeg figure----
jpeg("figures/figS5_nmds_male.jpg", width = 15, height = 15, units = "cm", 
     res = 300)
plot(ord_m$points, type = "n", bty = "n", 
     xlim = c(-.4, .5), 
     ylim = c(-.5, 1.25),
     xlab = "NMDS 1",
     ylab = "NMDS 2")
with(envi_m, 
     points(ord_m, 
            display = "sites", 
            col = ifelse(Farm_system == "Organic", "skyblue", "orange"), 
            cex = 1.5, 
            pch = -0x2642L))
with(envi_m, ordispider(ord_m, 
                        Farm_system, col = c("orange", "skyblue")))
with(envi_m, ordiellipse(ord_m, 
                         Farm_system, 
                         kind = "ehull", 
                         lwd = 2,
                         draw = "polygon", 
                         col = c("orange", "skyblue"), 
                         border = c("orange", "skyblue"),
                         alpha = 23,
                         xpd = T))
with(envi_m, 
     text(ord_m, 
          display = "species", 
          col = "gray20", 
          cex = .75,
          xpd = T))
legend("topright", 
       legend = c("Conventional", "Organic"),
       lwd = 1.5,
       col = c("orange", "skyblue"),
       bty = "n",
       title = "Farming System",
       cex = .75)
dev.off()


# Male choice----
#Exposing males of Osmia bicornis to females originally
#from organic and conventional farms. In the flight cage there is 5 males and 
#two frezzed-killied females from these two different farming management.

str(exp_1)

## Male latency----
#Inspect
stripchart(exp_1$Male_latency, method = "stack")
stripchart(exp_1$Male_latency ~ exp_1$Farm_system, method = "stack")
boxplot(exp_1$Male_latency ~ exp_1$Farm_system)


male_lat_mod <- lmer(Male_latency ~ Farm_system + 
                        (1 | Female_ID), 
                       data = exp_1)

#### Assumptions test (DHARMa) of the selected model----
# Residuals simulation
resid_male_lat <- DHARMa::simulateResiduals(fittedModel = male_lat_mod, 
                                           n = 1000) 
plot(resid_male_lat) 
DHARMa::plotQQunif(resid_male_lat) #
DHARMa::plotResiduals(resid_male_lat) #
DHARMa::testDispersion(resid_male_lat)
DHARMa::outliers(resid_male_lat)


### Global significance of the selected model----
size_null_mod <- lme4::lmer(Size ~ 1 +
                              (1 | Site/Station), 
                            data = envi_non_na)
anova(size_null_mod, size_sel_mod)

### Hypothesis test----
summary(male_lat_mod)
car::Anova(male_lat_mod)

## Time latency, time interval until males interact with females
A1st_touching <- glmer(st_touch ~ Farm_system + (1|Female_ID), 
                       data = exp_1, family = "binomial")
summary(A1st_touching)                     
Anova(A1st_touching)

simulationOutput <- simulateResiduals(fittedModel = A1st_touching, 
                                      n = 250)
simulationOutput
plot(simulationOutput)

summary(Time_1st_touching)                     
Anova(Time_1st_touching)
boxplot(Dataset$Touching_Time_Second~Dataset$Farm)
#Response: Touching_Time_Second
#Chisq Df Pr(>Chisq)  
#Farm 4.9837  1    0.02559 *

# Define the order of variables
variable_order <- c("Organic", "Conventional")

# Define colors for fill
variable_colors <- c("Organic" = "skyblue", "Conventional" = "orange")
# Define color for data points
data_point_color <- "#FF69B4"  # Pink color for data points

# Create a factor variable with the desired order
Dataset$Farm <- factor(Dataset$Farm, levels = variable_order)

# Plot the figure
TIME<-ggplot(Dataset, aes(x = Farm, y =Touching_Time_Second, fill = Farm)) +
  geom_boxplot(color = "black") +
  # geom_jitter(color = "black", fill = data_point_color, width = 0.2, height = 0, shape = 21, size = 3) +  # Add data points with black border
  # Customize other plot aesthetics as needed
  scale_fill_manual(values = variable_colors) +
  xlab("Female origin") +
  ylab("Time latency (s)") +
  geom_signif(y_position=750, xmin=0.7, xmax= 2.35,
              annotation="*", tip_length=0)+
  theme_classic()+theme(axis.text=element_text(size=17), axis.title=element_text(size=18,face="bold"))+
  labs(title = "") +  # Add the title of the plot if needed
  theme_classic() +  # Apply the classic theme
  theme(axis.text = element_text(size = 12),  # Increase the size of axis labels
        axis.title = element_text(size = 14))  # Increase the size of axis titles
# Add significance symbols for multiple comparisons

