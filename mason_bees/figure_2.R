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
# updated 2024, Jun, 11


# Required package----
library(ggplot2)
library(ggsignif)

# Figure 2----
## Figure 2. Boxplot shows the size of Osmia bicornis sampled in conventional and organic farming sites. The box in the plot represents the interquartile range (IQR), which spans from the first quartile (25th percentile) to the third quartile (75th percentile). The line inside the box represents the median (50th percentile), and vertical interval is the range of variation. Points outside these intervals are outliers. 

ggplot(envi, 
       aes(x = sex, 
           y = Size, 
           fill = Farm_system)) +
  geom_boxplot() + 
  scale_fill_manual(values = c("Organic" = "skyblue", 
                               "Conventional" = "orange"), 
                    name = "Farming system") +
  xlab("Sex") +
  ylab("Bee size (width between eyes in mm)") +
  theme_classic() +
  theme(axis.text = element_text(size = 12),  # Size of axis labels
        axis.title = element_text(size = 14)) +  
  scale_x_discrete(breaks = c("Female", "Male"),
                   labels = c("Female", "Male")) + # Change tick labels
  geom_signif(y_position = max(envi$Size, na.rm = T) + 
                max(envi$Size, na.rm = T) * 0.02, 
              xmin = 0.8, xmax = 1.2,
              annotation = "p = 0.042", 
              tip_length = 0.009) + 
  geom_signif(y_position = max(envi$Size, na.rm = T) + 
                max(envi$Size, na.rm = T) * 0.02, 
              xmin = 1.8, xmax = 2.2, 
              annotation = "p = 0.042", 
              tip_length = 0.009) + 
  geom_signif(y_position = min(envi$Size, na.rm = T), 
              xmin = 1.2, xmax = 1.8,
              annotation = "p = 0.001", 
              tip_length = -0.009) + 
  geom_signif(y_position = min(envi$Size, na.rm = T) - 
                min(envi$Size, na.rm = T) * 0.05, 
              xmin = 0.8, xmax = 2.2, 
              annotation = "p < 0.001", 
              tip_length = -0.009) + 
  geom_signif(y_position = max(envi$Size, na.rm = T) + 
                max(envi$Size, na.rm = T) * 0.05, 
              xmin = 1, xmax = 2,
              annotation = "p < 0.001", 
              tip_length = 0.009)

## Save as jpeg:
ggsave("figures/figure_2.jpg", width = 17, height = 13, units = "cm", dpi = 300)
## Save as pdf:
ggsave("figures/figure_2.pdf", width = 17, height = 13, units = "cm", dpi = 300)
## Save as tiff:
ggsave("figures/figure_2.tiff", width = 17, height = 13, units = "cm", dpi = 300)
## Save as png
ggsave("figures/figure_2.png", width = 17, height = 13, units = "cm", dpi = 300)

