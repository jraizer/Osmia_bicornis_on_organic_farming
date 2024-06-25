# Header----
# Manuscript draft
# Conventional agriculture interferes with sex communication and impacts on local population size in a wild bee

# Authors: Samuel Boff1*, Sara Olberz1, Irem Gülsoy2, Marvin Preuß1,  Josué Raizer3, Manfred Ayasse1

# Affiliations:
#1 Institute of Evolutionary Ecology and Conservation Genomics, Ulm University, Ulm, Germany
#2 Department of Molecular Biology and Genetics, İhsan Doğramacı Bilkent University, Ankara, Turkey
#3 Graduate Program in Entomology and Biodiversity Conservation, Federal University of Grande Dourados, Dourados, Brazil

# Correspondence
#*Email: samuel.boff@uni-ulm.de
#Phone number: +49 (0)731 50 22665

# Figure 2
# updated 2024, Jun, 10


# Required package----
library(effects)

# Figure 2----
# Figure 2. Estimates of bee size (95% IC) from farm system (A) and sex (B) effects based on the selected mixed model.

jpeg("figures/fig_2_alt.jpg", width = 20, height = 15, units = "cm", res = 300)
par(mfrow = c(1, 2))
plot(1:2, farm_eff$fit, 
     xaxt = "n", 
     xlab = "Farming system", 
     xlim = c(0.5, 2.5),
     ylab = "Bee size (width between eyes in mm)",
     ylim = c(min(farm_eff$lower), max(farm_eff$upper)), 
     col = c("orange", "skyblue"), pch = 19, cex = 1.5, 
     bty = "l")
mtext("A", at = .5, line = 1.5, cex = 2)
segments(1, farm_eff$lower[1, 1], 1, farm_eff$upper[1, 1], 
         col = "orange", lwd = 3)
segments(2, farm_eff$lower[2, 1], 2, farm_eff$upper[2, 1], 
         col = "skyblue", lwd = 3)
axis(1, at = 1:2, labels = c("Conventional", "Organic"))

plot(1:2, sex_eff$fit, 
     xaxt = "n", 
     xlab = "Sex", 
     xlim = c(0.5, 2.5),
     ylab = "Bee size (width between eyes in mm)",
     ylim = c(min(sex_eff$lower), max(sex_eff$upper)), 
     pch = 19, cex = 1.5, col = "gray",
     bty = "l")
mtext("B", at = .5, line = 1.5, cex = 2)
segments(1, sex_eff$lower[1, 1], 1, sex_eff$upper[1, 1], 
         col = "gray", lwd = 3)
segments(2, sex_eff$lower[2, 1], 2, sex_eff$upper[2, 1], 
         col = "gray", lwd = 3)
axis(1, at = 1:2, labels = c("Female", "Male"))

dev.off()
