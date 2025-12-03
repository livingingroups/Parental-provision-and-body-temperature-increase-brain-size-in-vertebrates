#############################
###20250826_by_Zitan#########
#############################
#
# Load required libraries
library(tidyr)         # Data tidying
library(ggplot2)       # Data visualization
library(dplyr)         # Data manipulation
library(intrval)       # Interval-based calculations
library(magick)        # Image processing
library(grid)          # Grid graphics
library(ggimage)       # Image annotations for ggplot2
library(ggridges)      # Ridgeline plots
library(rphylopic)     # Phylogenetic silhouette images
library(RCurl)         # HTTP requests
library(png)           # PNG image handling
library(ggforce)  


# Data input
#ParProv_Song_2025 <- readxl::read_excel( "~/ParProv_Song_etal2025_Data_submit.xlsx", sheet = 1)
rownames(ParProv_Song_2025)<-ParProv_Song_2025$Species # Set row names to the 'Species' column

# Reorder factor levels for plotting
ParProv_Song_2025$EggCarePlot <- factor(
  ParProv_Song_2025$EggCarePlot, 
  levels = c("Abandon", "Guard", "Bear", "Pre-hat. prov.", "Pre-&post-hat. prov.")
)
ParProv_Song_2025$Class <- factor(
  ParProv_Song_2025$Class, 
  levels = c("Hagfishes", "Lampreys", "Cartilaginous fishes", "Ray-finned fishes", "Coelacanths", "Lungfishes", "Amphibians", "Reptiles", "Birds", "Mammals")
)


##Figure 1, Allometry: brain ~ body#####
# Fit phylogenetic linear models for brain mass as a function of body mass across several vertebrate classes.
library(phylolm)
# Model for Lampreys
Lamprey_brain_body <- phylolm(log10(BrainMass.BrM..g.)~log10(Bodymass.BdM..g.), data=ParProv_Song_2025[which(ParProv_Song_2025$Class=="Lampreys"),], phy=jawlessTree, model = "lambda")
# Model for Cartilaginous fishes
shark_brain_body <- phylolm(log10(BrainMass.BrM..g.)~log10(Bodymass.BdM..g.), data=ParProv_Song_2025[which(ParProv_Song_2025$Class=="Cartilaginous fishes"),], phy=SharkTree, model = "lambda")
# Model for Ray-finned fishes
fish_brain_body <- phylolm(log10(BrainMass.BrM..g.)~log10(Bodymass.BdM..g.), data=ParProv_Song_2025[which(ParProv_Song_2025$Class=="Ray-finned fishes"),], phy=fishtree, model = "lambda")
# Model for Amphibians
Amphi_brain_body <- phylolm(log10(BrainMass.BrM..g.)~log10(Bodymass.BdM..g.), data=ParProv_Song_2025[which(ParProv_Song_2025$Class=="Amphibians"),], phy=AmphibiaTree, model = "lambda")
# Model for Reptiles
Reptil_brain_body <- phylolm(log10(BrainMass.BrM..g.)~log10(Bodymass.BdM..g.), data=ParProv_Song_2025[which(ParProv_Song_2025$Class=="Reptiles"),], phy=CrocoTurtleMCCTree, model = "lambda")
# Model for Birds
Bird_brain_body <- phylolm(log10(BrainMass.BrM..g.)~log10(Bodymass.BdM..g.), data=ParProv_Song_2025[which(ParProv_Song_2025$Class=="Birds"),], phy=Birdtree, model = "lambda")
# Model for Mammals
Mammal_brain_body <- phylolm(log10(BrainMass.BrM..g.)~log10(Bodymass.BdM..g.), data=ParProv_Song_2025[which(ParProv_Song_2025$Class=="Mammals"),], phy=MammalTree, model = "lambda")

colorblind_friendly_colors <- c(
  "Hagfishes" = "#00008B",   
  "Lampreys" = "#404040",
  "Cartilaginous fishes" = "#394A92",
  "Ray-finned fishes" = "#4F8095",
  "Coelacanths" = "#0d0d0d",
  "Lungfishes" = "#E69F00",   
  "Amphibians" = "#9D87C6",    
  "Reptiles" = "#83B24D", 
  "Birds" = "#D2691E",
  "Mammals" = "#9B3A4D")

# Filter the main dataset for relevant classes and non-missing categories, preparing data for further analysis and visualization
ParProv_Song_2025Ect <- ParProv_Song_2025 %>%
  filter(Class %in% c("Hagfishes", "Lampreys", "Cartilaginous fishes", "Lungfishes", "Coelacanths", "Ray-finned fishes", "Amphibians", "Reptiles", "Birds", "Mammals") & !is.na(EggCare5Cat))

# Calculate the range of log-transformed body mass and brain mass for Lampreys, which will be used to set up axes ranges or plot regression lines in a visual output
x_range_Lamprey <- log10(range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class == "Lampreys"),]$Bodymass.BdM..g., na.rm = TRUE))
y_range_Lamprey <- log10(range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class == "Lampreys"),]$BrainMass.BrM..g., na.rm = TRUE))

# Calculate the start and end points of the regression line
slope_Lamprey <- Lamprey_brain_body$coefficients[2]
intercept_Lamprey <- Lamprey_brain_body$coefficients[1]
start_y_Lamprey <- slope_Lamprey * x_range_Lamprey[1] + intercept_Lamprey
end_y_Lamprey <- slope_Lamprey * x_range_Lamprey[2] + intercept_Lamprey

x_range_shark <- log10(range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class=="Cartilaginous fishes"),]$Bodymass.BdM..g., na.rm = TRUE))
y_range_shark <- log10(range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class=="Cartilaginous fishes"),]$BrainMass.BrM..g., na.rm = TRUE))

# Calculate the start and end points of the regression line
slope_shark <- shark_brain_body$coefficients[2]
intercept_shark <- shark_brain_body$coefficients[1]
start_y_shark <- slope_shark * x_range_shark[1] + intercept_shark
end_y_shark <- slope_shark * x_range_shark[2] + intercept_shark

x_range_fish <- log10(range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class=="Ray-finned fishes"),]$Bodymass.BdM..g., na.rm = TRUE))
y_range_fish <- log10(range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class=="Ray-finned fishes"),]$BrainMass.BrM..g., na.rm = TRUE))

# Calculate the start and end points of the regression line
slope_fish <- fish_brain_body$coefficients[2]
intercept_fish <- fish_brain_body$coefficients[1]
start_y_fish <- slope_fish * x_range_fish[1] + intercept_fish
end_y_fish <- slope_fish * x_range_fish[2] + intercept_fish

x_range_Amphi<- log10(range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class=="Amphibians"),]$Bodymass.BdM..g., na.rm = TRUE))
y_range_Amphi <- log10(range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class=="Amphibians"),]$BrainMass.BrM..g., na.rm = TRUE))

# Calculate the start and end points of the regression line
slope_Amphi <- Amphi_brain_body$coefficients[2]
intercept_Amphi <-Amphi_brain_body$coefficients[1]
start_y_Amphi <- slope_Amphi * x_range_Amphi[1] + intercept_Amphi
end_y_Amphi <- slope_Amphi * x_range_Amphi[2] + intercept_Amphi

x_range_Reptil<- log10(range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class=="Reptiles"),]$Bodymass.BdM..g., na.rm = TRUE))
y_range_Reptil <- log10(range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class=="Reptiles"),]$BrainMass.BrM..g., na.rm = TRUE))

# Calculate the start and end points of the regression line
slope_Reptil <- Reptil_brain_body$coefficients[2]
intercept_Reptil <- Reptil_brain_body$coefficients[1]
start_y_Reptil <- slope_Reptil * x_range_Reptil[1] + intercept_Reptil
end_y_Reptil <- slope_Reptil * x_range_Reptil[2] + intercept_Reptil

x_range_Bird<- log10(range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class=="Birds"),]$Bodymass.BdM..g., na.rm = TRUE))
y_range_Bird <- log10(range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class=="Birds"),]$BrainMass.BrM..g., na.rm = TRUE))

# Calculate the start and end points of the regression line
slope_Bird <- Bird_brain_body$coefficients[2]
intercept_Bird <- Bird_brain_body$coefficients[1]
start_y_Bird <- slope_Bird * x_range_Bird[1] + intercept_Bird
end_y_Bird <- slope_Bird * x_range_Bird[2] + intercept_Bird

x_range_Mammal<- log10(range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class=="Mammals"),]$Bodymass.BdM..g., na.rm = TRUE))
y_range_Mammal <- log10(range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class=="Mammals"),]$BrainMass.BrM..g., na.rm = TRUE))

# Calculate the start and end points of the regression line
slope_Mammal <- Mammal_brain_body$coefficients[2]
intercept_Mammal <- Mammal_brain_body$coefficients[1]
start_y_Mammal <- slope_Mammal * x_range_Mammal[1] + intercept_Mammal
end_y_Mammal <- slope_Mammal * x_range_Mammal[2] + intercept_Mammal

ParProv_Song_2025Ect$Class<-factor(ParProv_Song_2025Ect$Class, levels = c("Hagfishes", "Lampreys", "Cartilaginous fishes", "Ray-finned fishes", "Coelacanths", "Lungfishes", "Amphibians", "Reptiles", "Birds", "Mammals"))



p<-NA
p<-ggplot(ParProv_Song_2025Ect, 
       aes(x =  log10(Bodymass.BdM..g.), y =log10(BrainMass.BrM..g.), group = Class, colour = Class , fill = Class)) +
  geom_point( data = subset(ParProv_Song_2025Ect, Class %in% c("Lampreys","Cartilaginous fishes", "Ray-finned fishes","Amphibians","Reptiles", "Birds", "Mammals")),  alpha = 0.4, size = 1, stroke= 0) +
  
  geom_mark_hull( data = subset(ParProv_Song_2025Ect, Class %in% c("Lampreys","Cartilaginous fishes", "Ray-finned fishes","Amphibians","Reptiles", "Birds", "Mammals")), 
                  alpha = 0.2, size = 0.5,  concavity = 10,expand=0,radius=0) +
  
  geom_segment(aes(x = x_range_Lamprey[1], xend = x_range_Lamprey[2], y = start_y_Lamprey, yend = end_y_Lamprey),  colour = "#404040", size = 1) + 
  geom_segment(aes(x = x_range_Amphi[1], xend = x_range_Amphi[2], y = start_y_Amphi, yend = end_y_Amphi ), colour = "#9D87C6" , size = 1) +
  geom_segment(aes(x = x_range_shark[1], xend = x_range_shark[2], y = start_y_shark, yend = end_y_shark ), colour = "#394A92", size = 1) +
  geom_segment(aes(x = x_range_Reptil[1], xend = x_range_Reptil[2], y = start_y_Reptil, yend = end_y_Reptil) , colour = "#83B24D"  , size = 1) +
  geom_segment(aes(x = x_range_fish[1], xend = x_range_fish[2], y = start_y_fish, yend = end_y_fish ), colour =  "#4F8095" , size = 1) +
  geom_segment(aes(x = x_range_Bird[1], xend = x_range_Bird[2], y = start_y_Bird, yend = end_y_Bird ), colour =  "#D2691E"   , size = 1) +
  geom_segment(aes(x = x_range_Mammal[1], xend = x_range_Mammal[2], y = start_y_Mammal, yend = end_y_Mammal ), colour = "#9B3A4D", size = 1) +
  
  geom_point( data = subset(ParProv_Song_2025Ect, Class %in% c("Hagfishes", "Lungfishes","Coelacanths")),  alpha = 0.8, size = 0.6) +
  geom_image(data = subset(ParProv_Song_2025Ect, Class == "Hagfishes"), aes(image = "~/Hagfish.svg"), size = 0.08, alpha = 0.9, color = "#00008B") + 
  geom_image(data = subset(ParProv_Song_2025Ect, Class == "Lungfishes"), aes(image = "~/Lungfish.svg"), size = 0.08, alpha = 0.9, color = "#E69F00") +
  geom_image(data = subset(ParProv_Song_2025Ect, Class == "Coelacanths"), aes(image = "~/Coelacanth.svg"), size = 0.08, alpha = 0.9, color = "#0d0d0d")+ 
  add_phylopic(uuid = "071ee517-a0f1-4d19-aa29-812b9f86cb53", alpha = 1, x = 7, y = 3.8 , ysize = 1, fill = "#9B3A4D"  )+ # mammals
  add_phylopic(uuid = "63953094-ac64-42c3-920e-53ee87ab188f", alpha = 1, x = 3, y = (-1.6) , ysize = 0.2, fill = "#404040" )+ # Lampreys
  add_phylopic(uuid = "b65312ae-91b5-45b1-9553-c192f1000aba", alpha = 1, x = 7, y = (1.8) , ysize = 0.5, fill = "#394A92" )+ # Cartilaginous fishes
  add_phylopic(uuid = "ea7ecd77-8c84-4da6-bca4-d7cfcc466558", alpha = 1, x = (-1), y = (-2) , ysize = 0.5, fill =  "#4F8095")+ # Ray-finned fishes
  add_phylopic(uuid = "43497e8a-45e7-4fa2-a8a0-ffadac8401fc", alpha = 1, x = (-1.2), y = (-2.7) , ysize = 0.35, fill =  "#9D87C6")+ # Amphibians
  add_phylopic(uuid = "bf7d9c5f-83c0-435a-b09f-dc6111ece257", alpha = 1, x = 5, y = (0.8) , ysize = 0.5, fill =  "#83B24D" )+ # Reptiles
  add_phylopic(uuid = "92589388-08e3-422f-b452-aa7454411a9c", alpha = 1, x = 5.4, y = (1.8) , ysize = 1, fill =  "#D2691E"  )+ # Bird
  
  scale_colour_manual(values = colorblind_friendly_colors, 
                      breaks =  rev(c("Hagfishes", "Lampreys", "Cartilaginous fishes", "Ray-finned fishes", "Coelacanths", "Lungfishes", "Amphibians", "Reptiles", "Birds", "Mammals"))) +
  scale_fill_manual(values = colorblind_friendly_colors, 
                    breaks =  rev(c("Hagfishes", "Lampreys", "Cartilaginous fishes", "Ray-finned fishes", "Coelacanths", "Lungfishes", "Amphibians", "Reptiles", "Birds", "Mammals"))) +
  
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size = 15),
        legend.title = element_text(colour = "steelblue", face = "bold.italic", size = 14   ),
        legend.text = element_text(face = "italic", colour = "steelblue4", size =10),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  labs(color= "Class", 
       x = "Body mass (log-10)",
       y = "Brain mass (log-10)") +
  guides( colour = guide_legend(override.aes = list(   linetype = "solid" , size = 3 )) ) + 
  coord_fixed(ratio = 1) 



pdf(file= paste(savepath, "Fig.1_2" , ".pdf", sep = "" ), width = 22/2.54, height = 18/2.54)
print(p)
dev.off()



#####

# Figure 2, Offspring mass####

# Function to add comparative annotations to ggplot figures
# It uses geom_segment and annotate functions to draw horizontal lines and text markers between groups to signify the statistical importance or differences.
# The function dynamically adjusts the y-positions of the lines and text based on the maximum values from the sample_sizes data frame, allowing for a flexible application across different data sets.
annotate_comparisons <- function(sample_sizes, a, b, c) {
  list(
    geom_segment(aes(x = 1, xend = 3, y = max(sample_sizes$max, na.rm = TRUE) * 1.45, yend = max(sample_sizes$max, na.rm = TRUE) * 1.45), color = "black"),
    geom_segment(aes(x = 1, xend = 1, y = max(sample_sizes$max, na.rm = TRUE) * 1.42, yend = max(sample_sizes$max, na.rm = TRUE) * 1.45), color = "black"),
    geom_segment(aes(x = 3, xend = 3, y = max(sample_sizes$max, na.rm = TRUE) * 1.42, yend = max(sample_sizes$max, na.rm = TRUE) * 1.45), color = "black"),
    annotate("text", x = 2, y = max(sample_sizes$max, na.rm = TRUE) * 1.46, label = a, size = 6, vjust = 0),
    geom_segment(aes(x = 1, xend = 2, y = max(sample_sizes$max, na.rm = TRUE) * 1.35, yend = max(sample_sizes$max, na.rm = TRUE) * 1.35), color = "black"),
    geom_segment(aes(x = 1, xend = 1, y = max(sample_sizes$max, na.rm = TRUE) * 1.32, yend = max(sample_sizes$max, na.rm = TRUE) * 1.35), color = "black"),
    geom_segment(aes(x = 2, xend = 2, y = max(sample_sizes$max, na.rm = TRUE) * 1.32, yend = max(sample_sizes$max, na.rm = TRUE) * 1.35), color = "black"),
    annotate("text", x = 1.5, y = max(sample_sizes$max, na.rm = TRUE) * 1.32, label = b, size = 6, vjust = 0),
    geom_segment(aes(x = 2, xend = 3, y = max(sample_sizes$max, na.rm = TRUE) * 1.25, yend = max(sample_sizes$max, na.rm = TRUE) * 1.25), color = "black"),
    geom_segment(aes(x = 2, xend = 2, y = max(sample_sizes$max, na.rm = TRUE) * 1.22, yend = max(sample_sizes$max, na.rm = TRUE) * 1.25), color = "black"),
    geom_segment(aes(x = 3, xend = 3, y = max(sample_sizes$max, na.rm = TRUE) * 1.22, yend = max(sample_sizes$max, na.rm = TRUE) * 1.25), color = "black"),
    annotate("text", x = 2.5, y = max(sample_sizes$max, na.rm = TRUE) * 1.26, label = c, size = 6, vjust = 0)
  )
}

# Filter and summarize data for Cartilaginous fishes
# Calculate sample sizes, maximum, and minimum residual offspring size
sample_sizes <- ParProv_Song_2025 %>% 
  filter(Class %in% c("Cartilaginous fishes"), !is.na(Residual_OS)) %>%
  group_by(EggCarePlot) %>%
  summarize(
    count = n(), # Number of samples
    max = max(Residual_OS), # Maximum residual offspring size
    min = min(Residual_OS)  # Minimum residual offspring size
  ) 

# Reorder factor levels for plotting in sample_sizes
sample_sizes$EggCarePlot <- factor(
  sample_sizes$EggCarePlot, 
  levels = c("Abandon", "Guard", "Bear", "Pre-hat. prov.", "Pre-&post-hat. prov.")
)

# Display the sample_sizes data frame
sample_sizes


# Plot data for Cartilaginous fishes
p<-NA
p<-ggplot(
  ParProv_Song_2025[which(
    ParProv_Song_2025$Class %in% c("Cartilaginous fishes") & 
      !is.na(ParProv_Song_2025$Residual_OS) & 
      ParProv_Song_2025$EggCarePlot %in% c("Abandon", "Bear", "Pre-hat. prov.")
  ),], 
  aes(x = EggCarePlot, y = Residual_OS, fill = EggCarePlot)
) +
  # Add violin plot
  geom_violin(width = 1.1) +
  # Add boxplot inside violin plot
  geom_boxplot(width = 0.25, position = position_dodge(1.2), color = "black", alpha = 0.2) +
  # Set manual fill colors for the plots
  scale_fill_manual(values = c("#1E466E", "#FFE6B7", "#E76254")) +
  # Set x-axis limits and ensure factors are not dropped
  scale_x_discrete(drop = FALSE, limits = c("Abandon", "Bear", "Pre-hat. prov.")) +
  # Customize the theme
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.position = "none"
  ) +
  # Set y-axis limits based on sample sizes
  ylim(
    min(sample_sizes$min, na.rm = TRUE), 
    max(sample_sizes$max, na.rm = TRUE) * 1.5
  ) +
  # Label y-axis
  ylab("Newborn mass (residual)") +
  # Remove x-axis label
  xlab("") +
  # Add a phylogenetic silhouette image for Cartilaginous fishes
  add_phylopic(
    uuid = "b65312ae-91b5-45b1-9553-c192f1000aba", 
    alpha = 1, 
    x = 3, 
    y = -1.2, 
    ysize = 0.5, 
    color = "black"
  ) +
  # Add sample sizes as text on the plot
  geom_text(
    data = sample_sizes[which(!is.na(sample_sizes$max)),], 
    aes(x = EggCarePlot, y = max, label = count), 
    size = 5, 
    position = position_dodge(width = 0.3), 
    vjust = -0.5
  ) +
  annotate_comparisons(sample_sizes, "**", " ", "." )  # Apply the annotation function here - Add significance indicators# Plot data for Cartilaginous fishes

pdf(file= paste(savepath, "Fig.2a" , ".pdf", sep = "" ), width = 11/2.54, height = 11/2.54)
print(p)
dev.off()



# Fish Data
sample_sizes <- ParProv_Song_2025 %>% 
  filter(Class == "Ray-finned fishes", !is.na(EggCarePlot)) %>%
  group_by(EggCarePlot) %>%
  summarize(
    count = n(),
    max = max(Residual_OS, na.rm = TRUE),
    min = min(Residual_OS, na.rm = TRUE)
  )
sample_sizes$EggCarePlot <- factor(sample_sizes$EggCarePlot, levels = c("Abandon", "Guard", "Bear", "Pre-hat. prov.", "Pre-&post-hat. prov."))

# Display the sample_sizes data frame
sample_sizes

p<-NA
p<-ggplot(ParProv_Song_2025[which(
  ParProv_Song_2025$Class %in% c("Ray-finned fishes") &
    !is.na(ParProv_Song_2025$Residual_OS) &
    ParProv_Song_2025$EggCarePlot %in% c("Abandon", "Guard", "Bear")
),], 
aes(x = EggCarePlot, y = Residual_OS, fill = EggCarePlot)
) +
  geom_violin(width = 1.1) +
  geom_boxplot(width = 0.25, position = position_dodge(1.2), color = "black", alpha = 0.2) +
  scale_fill_manual(values = c("#1E466E","#AADCE0","#FFE6B7")) +
  scale_x_discrete(drop = FALSE, limits = c("Abandon", "Guard", "Bear")) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.position = "none"
  ) +
  ylim(min(sample_sizes$min, na.rm = TRUE), max(sample_sizes$max, na.rm = TRUE) * 1.5) +
  ylab("Newborn mass (residual)") +
  xlab("") +
  add_phylopic(uuid = "ea7ecd77-8c84-4da6-bca4-d7cfcc466558", alpha = 1, x = 3, y = -1.5, ysize = 1, color = "black") +  # Ray-finned fishes
  geom_text(data = sample_sizes, aes(x = EggCarePlot, y = max, label = count), size = 5, position = position_dodge(width = 0.3), vjust = -0.5) +
  annotate_comparisons(sample_sizes, "**", " ", "*")  # Apply the annotation function here


pdf(file= paste(savepath, "Fig.2b" , ".pdf", sep = "" ), width = 11/2.54, height = 11/2.54)
print(p)
dev.off()

# Amphibian Data
sample_sizes <- ParProv_Song_2025 %>% 
  filter(Class == "Amphibians", !is.na(EggCarePlot)) %>%
  group_by(EggCarePlot) %>%
  summarize(
    count = n(),
    max = max(Residual_OS, na.rm = TRUE),
    min = min(Residual_OS, na.rm = TRUE)
  ) 

sample_sizes$EggCarePlot <- factor(sample_sizes$EggCarePlot, levels = c("Abandon", "Guard", "Bear", "Pre-hat.\n provision", "Post-hat.\n provision"))
sample_sizes

p<-NA
p<-ggplot(ParProv_Song_2025[which(
  ParProv_Song_2025$Class %in% c("Amphibians") &
    ParProv_Song_2025$EggCarePlot %in% c("Abandon", "Guard", "Bear")
),], 
aes(x = EggCarePlot, y = Residual_OS, fill = EggCarePlot)) +
  geom_violin(width = 1.1) +
  geom_boxplot(width = 0.25, position = position_dodge(1.2), color = "black", alpha = 0.2) +
  scale_fill_manual(values = c("#1E466E","#AADCE0","#FFE6B7")) +
  scale_x_discrete(drop = FALSE, limits = c("Abandon", "Guard", "Bear")) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.position = "none"
  ) +
  ylim(min(sample_sizes$min, na.rm = TRUE), max(sample_sizes$max, na.rm = TRUE) * 1.5) +
  ylab("Newborn mass (residual)") +
  xlab("") +
  add_phylopic(uuid = "43497e8a-45e7-4fa2-a8a0-ffadac8401fc", alpha = 1, x = 3, y = -0.9, ysize = 0.7, color = "black") +  # Amphibians
  geom_text(data = sample_sizes, aes(x = EggCarePlot, y = max, label = count), size = 5, position = position_dodge(width = 0.3), vjust = -0.5) +
  annotate_comparisons(sample_sizes, ".", "*", " ")  # Apply the annotation function here

pdf(file= paste(savepath, "Fig.2c" , ".pdf", sep = "" ), width = 11/2.54, height = 11/2.54)
print(p)
dev.off()

# Reptile Data
sample_sizes <- ParProv_Song_2025 %>% 
  filter(Class == "Reptiles", !is.na(EggCarePlot)) %>%
  group_by(EggCarePlot) %>%
  summarize(
    count = n(),
    max = max(Residual_OS, na.rm = TRUE),
    min = min(Residual_OS, na.rm = TRUE)
  )

sample_sizes$EggCarePlot <- factor(sample_sizes$EggCarePlot, levels = c("Abandon", "Guard","Bear", "Pre-hat. prov.", "Pre-&post-hat. prov."))
sample_sizes

p<-NA
p<-ggplot(ParProv_Song_2025[which(
  ParProv_Song_2025$Class %in% c("Reptiles")  & 
    ParProv_Song_2025$EggCarePlot %in% c("Abandon", "Guard", "Bear")),], 
  aes(x = EggCarePlot, y = Residual_OS, fill = EggCarePlot)
) +
  geom_violin(width = 1.0) +
  geom_boxplot(width = 0.25, position = position_dodge(1.2), color = "black", alpha = 0.2) +
  scale_fill_manual(values = c("#1E466E","#AADCE0","#FFE6B7")) +
  scale_x_discrete(drop = FALSE, limits = c("Abandon", "Guard", "Bear")) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.position = "none"
  ) +
  ylim(min(sample_sizes$min, na.rm = TRUE), max(sample_sizes$max, na.rm = TRUE) * 1.5) + 
  ylab("Newborn mass (residual)") +
  xlab("") +
  add_phylopic(uuid = "bf7d9c5f-83c0-435a-b09f-dc6111ece257", alpha = 1, x = 3, y = (-0.6) , ysize = 0.3, color = "black") + # Reptiles
  geom_text(
    data = sample_sizes[which(!is.na(sample_sizes$max)),], 
    aes(x = EggCarePlot, y = max, label = count), 
    size = 5, 
    position = position_dodge(width = 0.3), 
    vjust = -0.5
  ) +
  annotate_comparisons(sample_sizes, " ", " ", " ")  # Apply the annotation function here

pdf(file= paste(savepath, "Fig.2d" , ".pdf", sep = "" ), width = 11/2.54, height = 11/2.54)
print(p)
dev.off()

#####


## Figure 3a, Forest Plot ####

# Load and preprocess data

EstmateplotScale <- read.csv("/ForestPlot_log.csv")  # Load data from a specified path
colnames(EstmateplotScale) <- c("Var", "mean", "lb", "ub", "EffSample", "p", "Taxon", "Group")

# Convert columns to numeric where necessary
EstmateplotScale$mean <- as.numeric(EstmateplotScale$mean)
EstmateplotScale$lb <- as.numeric(EstmateplotScale$lb)
EstmateplotScale$ub <- as.numeric(EstmateplotScale$ub)

# Factorize the 'Taxon' column with a specific order for plotting
EstmateplotScale$Taxon <- factor(EstmateplotScale$Taxon, levels = rev(c("Cartilaginous fishes", "Ray-finned fishes", "Amphibians", "Reptiles", "Birds", "Mammals")))


# Standardize group names
EstmateplotScale[which(EstmateplotScale$Group == "Egg\\n abandoning"),]$Group <- "Abandon"
EstmateplotScale[which(EstmateplotScale$Group == "Egg\\n bearing"),]$Group <- "Bear"
EstmateplotScale[which(EstmateplotScale$Group == "Egg\\n guarding"),]$Group <- "Guard"
EstmateplotScale[which(EstmateplotScale$Group == "Pre-hatching\\n provision"),]$Group <- "Pre-hat.\n provision"

# Factorize the 'Group' column with a specific order
EstmateplotScale$Group <- factor(EstmateplotScale$Group, levels = rev(c("Pre-hat.\n provision", "Bear", "Guard", "Abandon")))

# Assign significance labels based on p-values
EstmateplotScale$significance <- NA
EstmateplotScale[which(EstmateplotScale$p == 0.0001),]$significance <- "***"
EstmateplotScale[which(EstmateplotScale$p > 0.001 & EstmateplotScale$p <= 0.01),]$significance <- "**"
EstmateplotScale[which(EstmateplotScale$p > 0.01 & EstmateplotScale$p <= 0.05),]$significance <- "*"
EstmateplotScale[which(EstmateplotScale$p > 0.05 & EstmateplotScale$p <= 0.1),]$significance <- "+"


# Create the forest plot using ggplot2
p<-ggplot(EstmateplotScale[which(EstmateplotScale$Var == "Offspring (log-10)"),], 
       aes(x = Taxon, y = mean, ymin = lb, ymax = ub)) +
  geom_point(aes(group = Group, color = Group), position = position_dodge(width = 0.25), 
             size = 1.5, shape = 21, colour = "black", stroke = 0.6) +
  geom_linerange(aes(group = Group, color = Group), position = position_dodge(width = 0.25), size = 0.6) +
  geom_hline(yintercept = 0.0, linetype = "dashed", color = "black", alpha = 0.2) +
  geom_text(aes(label = significance, group = Group), position = position_dodge(width = 0.25), vjust = 1.8, size = 2) +
  facet_wrap(~ Group, ncol = 4, scales = "free_x") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 0, hjust = 1, size = 6),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 8),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 8, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey", size = 0.2, linetype = "solid"),
    panel.grid.minor = element_blank()
    #axis.line = element_line(colour = "black")
  ) +
  coord_flip() +  # Flip coordinates to make it a horizontal plot
  scale_colour_manual(values = c("#1E466E", "#AADCE0", "#FFE6B7", "#E76254")) +  # Set colors for different groups
  scale_y_continuous(name = "Effect estimate", limits = c(-0.35, 0.72), breaks = c(-0.2, 0.2, 0.6))  # Set y-axis properties


pdf(file= paste(savepath, "Fig.3a" , ".pdf", sep = "" ), width = 8.9/2.54, height = 8.9/2.54)
print(p)
dev.off()

 
#####

##Figure 3b, Pie plot####
size1 <-2.5
sample_sizes <- ParProv_Song_2025 %>% 
  filter(Class %in% c("Cartilaginous fishes"), !is.na(Residual_OS)) %>%
  group_by(EggCarePlot) %>%
  summarize(
    count = n(),
    max = max(Residual_OS),
    min = min(Residual_OS)
  ) 
sample_sizes$EggCarePlot<-factor(sample_sizes$EggCarePlot, levels = c("Abandon", "Guard","Bear", "Pre-hat. prov.", "Pre-&post-hat. prov."))
sample_sizes

pieShark<-ggplot(sample_sizes[which(!is.na(sample_sizes$max)),], aes(x = "", y = count, fill = EggCarePlot)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") + 
  theme_void() +  
  theme(legend.position = "none") +   
  scale_fill_manual(values=c("Abandon" = "#1E466E",
                             "Bear"="#FFE6B7",
                             "Pre-hat. prov." = "#E76254"))  +
  geom_text(aes(label = count, color = EggCarePlot), position = position_stack(vjust = 0.5), size = size1) +
  scale_color_manual(values = c("white", "black", "black", "black")) 

sample_sizes <- ParProv_Song_2025 %>% 
  filter(Class %in% c("Ray-finned fishes"), !is.na(Residual_OS), !is.na(EggCarePlot)) %>%
  group_by(EggCarePlot) %>%
  summarize(
    count = n(),
    max = max(Residual_OS),
    min = min(Residual_OS)
  ) 
sample_sizes$EggCarePlot<-factor(sample_sizes$EggCarePlot, levels = c("Abandon", "Guard","Bear", "Pre-hat. prov.", "Pre-&post-hat. prov."))
sample_sizes


piefish<-ggplot(sample_sizes[which(!is.na(sample_sizes$max)),], aes(x = "", y = count, fill = EggCarePlot)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") + 
  theme_void() +   
  theme(legend.position = "none") +   
  scale_fill_manual(values=c("Abandon" = "#1E466E",
                             "Guard" = "#AADCE0",
                             "Bear"="#FFE6B7",
                             "Pre-hat. prov." = "#E76254"))   +
  geom_text(aes(label = count, color = EggCarePlot), position = position_stack(vjust = 0.5), size = size1) +
  scale_color_manual(values = c("white", "black", "black", "black")) 

sample_sizes <- ParProv_Song_2025 %>% 
  filter(Class %in% c("Amphibians"), !is.na(Residual_OS), !is.na(EggCarePlot)) %>%
  group_by(EggCarePlot) %>%
  summarize(
    count = n(),
    max = max(Residual_OS),
    min = min(Residual_OS)
  ) 

sample_sizes$EggCarePlot<-factor(sample_sizes$EggCarePlot, levels = c("Abandon", "Guard","Bear","Pre-hat.\n provision", "Post-hat.\n provision"))
sample_sizes


pieAmphi<-ggplot(sample_sizes[which(!is.na(sample_sizes$max)),], aes(x = "", y = count, fill = EggCarePlot)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") + 
  theme_void() +   
  theme(legend.position = "none") +   
  scale_fill_manual(values=c("Abandon" = "#1E466E",
                             "Guard" = "#AADCE0",
                             "Bear"="#FFE6B7",
                             "Pre-hat. prov." = "#E76254"))   +
  geom_text(aes(label = count, color = EggCarePlot), position = position_stack(vjust = 0.5), size = size1) +
  scale_color_manual(values = c("white", "black", "black", "black")) 


sample_sizes <- ParProv_Song_2025 %>% 
  filter(Class %in% c("Reptiles"), !is.na(Residual_OS), !is.na(EggCarePlot)) %>%
  group_by(EggCarePlot) %>%
  summarize(
    count = n(),
    max = max(Residual_OS),
    min = min(Residual_OS)
  ) 
sample_sizes$EggCarePlot<-factor(sample_sizes$EggCarePlot, levels = c("Abandon", "Guard","Bear", "Pre-hat. prov.", "Pre-&post-hat. prov."))
sample_sizes


pieReptil<-ggplot(sample_sizes[which(!is.na(sample_sizes$max)),], aes(x = "", y = count, fill = EggCarePlot)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") + 
  theme_void() +  
  theme(legend.position = "none") +   
  scale_fill_manual(values=c("Abandon" = "#1E466E",
                             "Guard" = "#AADCE0",
                             "Bear"="#FFE6B7",
                             "Pre-hat. prov." = "#E76254"))   +
  geom_text(aes(label = count, color = EggCarePlot), position = position_stack(vjust = 0.5), size = size1) +
  scale_color_manual(values = c("white", "black", "black", "black")) 


sample_sizes <- ParProv_Song_2025 %>% 
  filter(Class %in% c("Birds"), !is.na(Residual_OS)) %>%
  group_by(EggCarePlot) %>%
  summarize(
    count = n(),
    max = max(Residual_OS),
    min = min(Residual_OS)
  ) 
sample_sizes$EggCarePlot<-factor(sample_sizes$EggCarePlot, levels = c("Abandon", "Guard","Bear", "Pre-hat. prov.", "Pre-&post-hat. prov."))
sample_sizes

sample_sizes[which(sample_sizes$EggCarePlot=="Pre-hat. prov."),]$count<- I(299+833)
sample_sizes <- sample_sizes[-which(sample_sizes$EggCarePlot=="Pre-&post-hat. prov."),]
pieBird<-ggplot(sample_sizes[which(!is.na(sample_sizes$max)),], aes(x = "", y = count, fill = EggCarePlot)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") + 
  theme_void() +   
  theme(legend.position = "none") +   
  scale_fill_manual(values=c("Abandon" = "#1E466E",
                             "Guard" = "#AADCE0",
                             "Bear"="#FFE6B7",
                             "Pre-hat. prov." = "#E76254"))   +
  geom_text(aes(label = count, color = EggCarePlot), position = position_stack(vjust = 0.5), size = size1) +
  scale_color_manual(values = c("white", "black", "black", "black")) 


sample_sizes <- ParProv_Song_2025 %>% 
  filter(Class %in% c("Mammals"), !is.na(Residual_OS)) %>%
  group_by(EggCarePlot) %>%
  summarize(
    count = n(),
    max = max(Residual_OS),
    min = min(Residual_OS)
  ) 
sample_sizes
sample_sizes$EggCarePlot<-factor(sample_sizes$EggCarePlot, levels = c("Abandon", "Guard","Bear", "Pre-hat. prov.", "Pre-&post-hat. prov."))
sample_sizes

sample_sizes[which(sample_sizes$EggCarePlot=="Pre-&post-hat. prov."),]$EggCarePlot<-"Pre-hat. prov."

pieMamm<-ggplot(sample_sizes[which(sample_sizes$count!=0),], aes(x = "", y = count, fill = EggCarePlot)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") + 
  theme_void() +  
  theme(legend.position = "none")+  
  scale_fill_manual(values=c("Abandon" = "#1E466E",
                             "Guard" = "#AADCE0",
                             "Bear"="#FFE6B7",
                             "Pre-hat. prov." = "#E76254"))+
  geom_text(aes(label = count, color = EggCarePlot), position = position_stack(vjust = 0.5), size = size1) +
  scale_color_manual(values = c("white")) 


library(patchwork)

combined_plot <- (pieShark / piefish / pieAmphi / pieReptil / pieBird / pieMamm)  
combined_plot

pdf(file= paste(savepath, "Fig.3b" , ".pdf", sep = "" ), width = 8.9/2.54, height = 8.9/2.54)
print(combined_plot)
dev.off()
 

combined_plot2<-p+combined_plot +plot_layout(ncol = 2, widths = c(8, 2))

pdf(file= paste(savepath, "Fig.3a&b" , ".pdf", sep = "" ), width = 22/2.54, height = 18/2.54)
print(combined_plot2)
dev.off()


#####


##Figure 3c, Allometry: brain ~ Offspring#####
# Each phylolm call fits a phylogenetic linear model to log-transformed brain mass as a function of log-transformed offspring mass for each specified class. The models account for evolutionary relationships using phylogenetic trees specific to each class.

# Fitting models for different vertebrate classes
Lamprey_Brain_Off <- phylolm(log10(AdultBrain)~log10(Offspring), data=ParProv_Song_2025[which(ParProv_Song_2025$Class=="Lampreys"),], phy=jawlessTree, model = "lambda")
shark_Brain_Off <- phylolm(log10(AdultBrain)~log10(Offspring), data=ParProv_Song_2025[which(ParProv_Song_2025$Class=="Cartilaginous fishes"),], phy=SharkTree, model = "lambda")
fish_Brain_Off <- phylolm(log10(AdultBrain)~log10(Offspring), data=ParProv_Song_2025[which(ParProv_Song_2025$Class=="Ray-finned fishes"  ),], phy=fishtree, model = "lambda")
Amphi_Brain_Off <- phylolm(log10(AdultBrain)~log10(Offspring),data=ParProv_Song_2025[which(ParProv_Song_2025$Class=="Amphibians"  ),], phy=AmphibiaTree, model = "lambda")
Reptil_Brain_Off <- phylolm(log10(AdultBrain)~log10(Offspring), data=ParProv_Song_2025[which(ParProv_Song_2025$Class=="Reptiles"  ),], phy=CrocoTurtleMCCTree, model = "lambda")
Bird_Brain_Off <- phylolm(log10(AdultBrain)~log10(Offspring), data=ParProv_Song_2025[which(ParProv_Song_2025$Class=="Birds"  ),], phy=Birdtree, model = "lambda")
Mammal_Brain_Off <- phylolm(log10(AdultBrain)~log10(Offspring), data=ParProv_Song_2025[which(ParProv_Song_2025$Class=="Mammals"  ),], phy=MammalTree, model = "lambda")

# Setting up color palette for better visibility for colorblind users
colorblind_friendly_colors <- c(
  "Hagfishes" = "#00008B",   
  "Lampreys" = "#404040",
  "Cartilaginous fishes" = "#394A92",
  "Ray-finned fishes" = "#4F8095",
  "Coelacanths" = "#0d0d0d",
  "Lungfishes" = "#E69F00",  
  "Amphibians" = "#9D87C6",    
  "Reptiles" = "#83B24D", 
  "Birds" = "#D2691E",
  "Mammals" = "#9B3A4D")

# Filter data for relevant classes and non-missing categories
ParProv_Song_2025Ect <- ParProv_Song_2025 %>%
  filter(Class %in% c("Hagfishes", "Lampreys", "Cartilaginous fishes", "Lungfishes", "Coelacanths", "Ray-finned fishes", "Amphibians", "Reptiles", "Birds", "Mammals") & !is.na(EggCare5Cat))

# Calculating x and y ranges for regression line plotting for each class based on the range of log-transformed offspring and brain masses
x_range_Lamprey <- log10(range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class=="Lampreys"),]$Offspring, na.rm = TRUE))
y_range_Lamprey <- log10(range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class=="Lampreys"),]$AdultBrain, na.rm = TRUE))

# Calculate the start and end points of the regression line
slope_Lamprey <- Lamprey_Brain_Off$coefficients[2]
intercept_Lamprey <- Lamprey_Brain_Off$coefficients[1]
start_y_Lamprey <- slope_Lamprey * x_range_Lamprey[1] + intercept_Lamprey
end_y_Lamprey <- slope_Lamprey * x_range_Lamprey[2] + intercept_Lamprey

x_range_shark <- log10(range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class=="Cartilaginous fishes"),]$Offspring, na.rm = TRUE))
y_range_shark <- log10(range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class=="Cartilaginous fishes"),]$AdultBrain, na.rm = TRUE))

# Calculate the start and end points of the regression line
slope_shark <- shark_Brain_Off$coefficients[2]
intercept_shark <- shark_Brain_Off$coefficients[1]
start_y_shark <- slope_shark * x_range_shark[1] + intercept_shark
end_y_shark <- slope_shark * x_range_shark[2] + intercept_shark

x_range_fish <- log10(range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class=="Ray-finned fishes"),]$Offspring, na.rm = TRUE))
y_range_fish <- log10(range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class=="Ray-finned fishes"),]$AdultBrain, na.rm = TRUE))

# Calculate the start and end points of the regression line
slope_fish <- fish_Brain_Off$coefficients[2]
intercept_fish <- fish_Brain_Off$coefficients[1]
start_y_fish <- slope_fish * x_range_fish[1] + intercept_fish
end_y_fish <- slope_fish * x_range_fish[2] + intercept_fish

x_range_Amphi<- log10(range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class=="Amphibians"),]$Offspring, na.rm = TRUE))
y_range_Amphi <- log10(range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class=="Amphibians"),]$AdultBrain, na.rm = TRUE))

# Calculate the start and end points of the regression line
slope_Amphi <- Amphi_Brain_Off$coefficients[2]
intercept_Amphi <-Amphi_Brain_Off$coefficients[1]
start_y_Amphi <- slope_Amphi * x_range_Amphi[1] + intercept_Amphi
end_y_Amphi <- slope_Amphi * x_range_Amphi[2] + intercept_Amphi

x_range_Reptil<- log10(range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class=="Reptiles"),]$Offspring, na.rm = TRUE))
y_range_Reptil <- log10(range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class=="Reptiles"),]$AdultBrain, na.rm = TRUE))

# Calculate the start and end points of the regression line
slope_Reptil <- Reptil_Brain_Off$coefficients[2]
intercept_Reptil <- Reptil_Brain_Off$coefficients[1]
start_y_Reptil <- slope_Reptil * x_range_Reptil[1] + intercept_Reptil
end_y_Reptil <- slope_Reptil * x_range_Reptil[2] + intercept_Reptil

x_range_Bird<- log10(range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class=="Birds"),]$Offspring, na.rm = TRUE))
y_range_Bird <- log10(range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class=="Birds"),]$AdultBrain, na.rm = TRUE))

# Calculate the start and end points of the regression line
slope_Bird <- Bird_Brain_Off$coefficients[2]
intercept_Bird <- Bird_Brain_Off$coefficients[1]
start_y_Bird <- slope_Bird * x_range_Bird[1] + intercept_Bird
end_y_Bird <- slope_Bird * x_range_Bird[2] + intercept_Bird

x_range_Mammal<- log10(range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class=="Mammals"),]$Offspring, na.rm = TRUE))
y_range_Mammal <- log10(range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class=="Mammals"),]$AdultBrain, na.rm = TRUE))

# Calculate the start and end points of the regression line
slope_Mammal <- Mammal_Brain_Off$coefficients[2]
intercept_Mammal <- Mammal_Brain_Off$coefficients[1]
start_y_Mammal <- slope_Mammal * x_range_Mammal[1] + intercept_Mammal
end_y_Mammal <- slope_Mammal * x_range_Mammal[2] + intercept_Mammal


ParProv_Song_2025Ect$Class <- factor(ParProv_Song_2025Ect$Class, levels = c("Hagfishes", "Lampreys", "Cartilaginous fishes", "Ray-finned fishes", "Coelacanths", "Lungfishes", "Amphibians", "Reptiles", "Birds", "Mammals"))

p<-NA
p<-ggplot(ParProv_Song_2025Ect, aes(x = log10(Offspring), y = log10(AdultBrain), group = Class, colour = Class, fill = Class)) +
  geom_point(data = subset(ParProv_Song_2025Ect, Class %in% c("Lampreys", "Cartilaginous fishes", "Ray-finned fishes", "Amphibians", "Reptiles", "Birds", "Mammals")), alpha = 0.4, size = 1, stroke= 0) +
  #geom_smooth(method = "lm", se = F, alpha = 0 , size = 0,show.legend = T ) +
  
  geom_mark_hull( data = subset(ParProv_Song_2025Ect, Class %in% c("Lampreys","Cartilaginous fishes", "Ray-finned fishes","Amphibians","Reptiles", "Birds", "Mammals")), 
                  alpha = 0.2,  concavity = 10,expand=0,radius=0,aes(fill=Class)) +
  
  geom_segment(aes(x = x_range_Lamprey[1], xend = x_range_Lamprey[2], y = start_y_Lamprey, yend = end_y_Lamprey), colour = "#404040", size = 1) +
  geom_segment(aes(x = x_range_Amphi[1], xend = x_range_Amphi[2], y = start_y_Amphi, yend = end_y_Amphi), colour = "#BAAFD1", size = 1) +
  geom_segment(aes(x = x_range_shark[1], xend = x_range_shark[2], y = start_y_shark, yend = end_y_shark), colour = "#394A92", size = 1) +
  geom_segment(aes(x = x_range_Reptil[1], xend = x_range_Reptil[2], y = start_y_Reptil, yend = end_y_Reptil), colour = "#68AC57", size = 1) +
  geom_segment(aes(x = x_range_fish[1], xend = x_range_fish[2], y = start_y_fish, yend = end_y_fish), colour = "#70A0AC", size = 1) +
  geom_segment(aes(x = x_range_Bird[1], xend = x_range_Bird[2], y = start_y_Bird, yend = end_y_Bird), colour = "#F4C28F", size = 1) +
  geom_segment(aes(x = x_range_Mammal[1], xend = x_range_Mammal[2], y = start_y_Mammal, yend = end_y_Mammal), colour = "#9B3A4D", size = 1) +
  
  geom_point( data = subset(ParProv_Song_2025Ect, Class %in% c("Hagfishes", "Lungfishes","Coelacanths")),  alpha = 0.8, size = 0.6) +
  geom_image(data = subset(ParProv_Song_2025Ect, Class == "Hagfishes"), aes(image = "~/Hagfish.svg"), size = 0.06, alpha = 0.8, color = "#00008B") + 
  geom_image(data = subset(ParProv_Song_2025Ect, Class == "Lungfishes"), aes(image = "~/Lungfish.svg"), size = 0.04, alpha = 0.8, color = "#E69F00") +
  geom_image(data = subset(ParProv_Song_2025Ect, Class == "Coelacanths"), aes(image = "~/Coelacanth.svg"), size = 0.06, alpha = 0.8, color = "#0d0d0d")+ 
  
  add_phylopic(uuid = "071ee517-a0f1-4d19-aa29-812b9f86cb53", alpha = 1, x = 6, y = 5.2 , ysize = 1, fill = "#9B3A4D")+ # mammals
  add_phylopic(uuid = "63953094-ac64-42c3-920e-53ee87ab188f", alpha = 1, x = (-4.8), y = (-1.8) , ysize = 0.2, fill = "#404040")+ # Lampreys
  add_phylopic(uuid = "b65312ae-91b5-45b1-9553-c192f1000aba", alpha = 1, x = 5.2, y = (2.4) , ysize = 0.5, fill = "#394A92")+ # Cartilaginous fishes
  add_phylopic(uuid = "ea7ecd77-8c84-4da6-bca4-d7cfcc466558", alpha = 1, x = (-5), y = (-0.8) , ysize = 0.5, fill = "#4F8095")+ # Ray-finned fishes
  add_phylopic(uuid = "43497e8a-45e7-4fa2-a8a0-ffadac8401fc", alpha = 1, x = (-1.8), y = (-2.3) , ysize = 0.35, fill = "#9D87C6")+ # Amphibians
  add_phylopic(uuid = "bf7d9c5f-83c0-435a-b09f-dc6111ece257", alpha = 1, x = 2, y = (0.6) , ysize = 0.5, fill = "#83B24D")+ # Reptiles
  add_phylopic(uuid = "92589388-08e3-422f-b452-aa7454411a9c", alpha = 1, x = (-0.5), y = (-0.5) , ysize = 1, fill = "#D2691E")+ # Bird
  scale_colour_manual(values = colorblind_friendly_colors, 
                      breaks =  rev(c("Hagfishes", "Lampreys", "Cartilaginous fishes", "Ray-finned fishes", "Coelacanths", "Lungfishes", "Amphibians", "Reptiles", "Birds", "Mammals"))) +
  scale_fill_manual(values = colorblind_friendly_colors, 
                    breaks =  rev(c("Hagfishes", "Lampreys", "Cartilaginous fishes", "Ray-finned fishes", "Coelacanths", "Lungfishes", "Amphibians", "Reptiles", "Birds", "Mammals"))) +
  
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size = 15),
        legend.title = element_text(colour = "steelblue", face = "bold.italic", size = 14   ),
        legend.text = element_text(face = "italic", colour = "steelblue4", size =10),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  labs(color= "Class", 
       x = "Newborn mass (log-10)",
       y = "Brain mass (log-10)") +
  guides(colour = guide_legend(override.aes = list(size = 3))) + 
  coord_fixed(ratio = 1)

pdf(file= paste(savepath, "Fig.3c_2" , ".pdf", sep = "" ), width = 22/2.54, height = 18/2.54)
print(p)
dev.off()

#####


##### Figure 4, Residual_BS_all ~ Tb #####

# Perform phylogenetic linear modeling on log-transformed offspring mass as a function of adult mass for various classes
# Each line models the relationship within a specific taxonomic group, adjusting for phylogenetic relatedness.

# Filter dataset for the classes of interest, ensuring relevant data points are not missing
ParProv_Song_2025Ect <- ParProv_Song_2025 %>%
  filter(Class %in% c( "Lampreys", "Cartilaginous fishes", "Ray-finned fishes", "Amphibians", "Reptiles", "Birds", "Mammals") & !is.na(Tb1))

# Model for Lampreys
Lamprey_BS_Tb <- phylolm(Residual_BS_all ~ Tb, data = ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class == "Lampreys"),], phy = jawlessTree, model = "lambda")
# Model for Cartilaginous fishes
shark_BS_Tb <- phylolm(Residual_BS_all ~ Tb, data = ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class == "Cartilaginous fishes"),], phy = SharkTree, model = "lambda")
# Model for Ray-finned fishes
fish_BS_Tb <- phylolm(Residual_BS_all~ Tb, data = ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class == "Ray-finned fishes"),], phy = fishtree, model = "lambda")
# Model for Amphibians
Amphi_BS_Tb <- phylolm(Residual_BS_all~ Tb, data = ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class == "Amphibians"),], phy = AmphibiaTree, model = "lambda")
# Model for Reptiles
Reptil_BS_Tb <- phylolm(Residual_BS_all ~ Tb, data = ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class == "Reptiles"),], phy = CrocoTurtleMCCTree, model = "lambda")
# Model for Birds
Bird_BS_Tb <- phylolm(Residual_BS_all~ Tb, data = ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class == "Birds"),], phy = Birdtree, model = "lambda")
# Model for Mammals
Mammal_BS_Tb <- phylolm(Residual_BS_all ~ Tb, data = ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class == "Mammals"),], phy = MammalTree, model = "lambda")

# Define color palette that is friendly for colorblind individuals
colorblind_friendly_colors <- c(
  
  "Lampreys" = "#404040",
  "Cartilaginous fishes" = "#394A92",
  "Ray-finned fishes" = "#4F8095",
  
  "Amphibians" = "#9D87C6",    
  "Reptiles" = "#83B24D", 
  "Birds" = "#D2691E",
  "Mammals" = "#9B3A4D")


# For Lampreys
# Calculate log10 transformed ranges for adult mass and offspring mass.
x_range_Lamprey <- range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class == "Lampreys"),]$Tb1, na.rm = TRUE)
y_range_Lamprey <- range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class == "Lampreys"),]$Residual_BS_all, na.rm = TRUE)

# Determine the start and end points of the regression line based on computed slopes and intercepts from the phylolm models.
slope_Lamprey <- Lamprey_BS_Tb1$coefficients[2]
intercept_Lamprey <- Lamprey_BS_Tb1$coefficients[1]
start_y_Lamprey <- slope_Lamprey * x_range_Lamprey[1] + intercept_Lamprey
end_y_Lamprey <- slope_Lamprey * x_range_Lamprey[2] + intercept_Lamprey

# Similar calculations for Cartilaginous fishes
x_range_shark <- range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class == "Cartilaginous fishes"),]$Tb1, na.rm = TRUE)
y_range_shark <-  range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class == "Cartilaginous fishes"),]$Residual_BS_all, na.rm = TRUE)

slope_shark <- shark_BS_Tb1$coefficients[2]
intercept_shark <- shark_BS_Tb1$coefficients[1]
start_y_shark <- slope_shark * x_range_shark[1] + intercept_shark
end_y_shark <- slope_shark * x_range_shark[2] + intercept_shark

# Similar calculations for Ray-finned fishes
x_range_fish <- range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class == "Ray-finned fishes"),]$Tb1, na.rm = TRUE)
y_range_fish <- range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class == "Ray-finned fishes"),]$Residual_BS_all, na.rm = TRUE)

slope_fish <- fish_BS_Tb1$coefficients[2]
intercept_fish <- fish_BS_Tb1$coefficients[1]
start_y_fish <- slope_fish * x_range_fish[1] + intercept_fish
end_y_fish <- slope_fish * x_range_fish[2] + intercept_fish

# Similar calculations for Amphibians
x_range_Amphi <- range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class == "Amphibians"),]$Tb1, na.rm = TRUE)
y_range_Amphi <- range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class == "Amphibians"),]$Residual_BS_all, na.rm = TRUE)

slope_Amphi <- Amphi_BS_Tb1$coefficients[2]
intercept_Amphi <- Amphi_BS_Tb1$coefficients[1]
start_y_Amphi <- slope_Amphi * x_range_Amphi[1] + intercept_Amphi
end_y_Amphi <- slope_Amphi * x_range_Amphi[2] + intercept_Amphi

# Similar calculations for Reptiles
x_range_Reptil <- range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class == "Reptiles"),]$Tb1, na.rm = TRUE)
y_range_Reptil <- range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class == "Reptiles"),]$Residual_BS_all, na.rm = TRUE)

slope_Reptil <- Reptil_BS_Tb1$coefficients[2]
intercept_Reptil <- Reptil_BS_Tb1$coefficients[1]
start_y_Reptil <- slope_Reptil * x_range_Reptil[1] + intercept_Reptil
end_y_Reptil <- slope_Reptil * x_range_Reptil[2] + intercept_Reptil

# Similar calculations for Birds
x_range_Bird <- range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class == "Birds"),]$Tb1, na.rm = TRUE)
y_range_Bird <- range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class == "Birds"),]$Residual_BS_all, na.rm = TRUE)

slope_Bird <- Bird_BS_Tb1$coefficients[2]
intercept_Bird <- Bird_BS_Tb1$coefficients[1]
start_y_Bird <- slope_Bird * x_range_Bird[1] + intercept_Bird
end_y_Bird <- slope_Bird * x_range_Bird[2] + intercept_Bird

# Similar calculations for Mammals
x_range_Mammal <-  range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class == "Mammals"),]$Tb1, na.rm = TRUE)
y_range_Mammal <- range(ParProv_Song_2025Ect[which(ParProv_Song_2025Ect$Class == "Mammals"),]$Residual_BS_all, na.rm = TRUE)

slope_Mammal <- Mammal_BS_Tb1$coefficients[2]
intercept_Mammal <- Mammal_BS_Tb1$coefficients[1]
start_y_Mammal <- slope_Mammal * x_range_Mammal[1] + intercept_Mammal
end_y_Mammal <- slope_Mammal * x_range_Mammal[2] + intercept_Mammal


# Set the 'Class' column as a factor with specific order for plotting purposes
ParProv_Song_2025Ect$Class <- factor(ParProv_Song_2025Ect$Class, levels = c( "Lampreys", "Cartilaginous fishes", "Ray-finned fishes", "Amphibians", "Reptiles", "Birds", "Mammals"))


# Create a ggplot
p<-NA
p<-ggplot(ParProv_Song_2025Ect, aes(x = Tb1 , y =Residual_BS_all, group = Class, colour = Class, fill = Class)) +
  geom_point(data = ParProv_Song_2025Ect, alpha = 0.4, size = 1, stroke= 0) +
  #geom_smooth(method = "lm", se = FALSE, alpha = 0 ) +
  
  geom_mark_hull( data = ParProv_Song_2025Ect , 
                  alpha = 0.2,  concavity = 10,expand=0,radius=0,aes(fill=Class)) +
  
  # Add regression line segments for each class, colored and sized appropriately
  geom_segment(aes(x = x_range_Lamprey[1], xend = x_range_Lamprey[2], y = start_y_Lamprey, yend = end_y_Lamprey, colour = "Lampreys"), size = 1, show_guide=T) +
  geom_segment(aes(x = x_range_Amphi[1], xend = x_range_Amphi[2], y = start_y_Amphi, yend = end_y_Amphi, colour = "Amphibians"), size = 1,show_guide=T) +
  geom_segment(aes(x = x_range_shark[1], xend = x_range_shark[2], y = start_y_shark, yend = end_y_shark, colour = "Cartilaginous fishes"),  size = 1,show_guide=T) +
  geom_segment(aes(x = x_range_Reptil[1], xend = x_range_Reptil[2], y = start_y_Reptil, yend = end_y_Reptil, colour = "Reptiles" ), size = 1,show_guide=T) +
  geom_segment(aes(x = x_range_fish[1], xend = x_range_fish[2], y = start_y_fish, yend = end_y_fish, colour = "Ray-finned fishes" ),  size = 1,show_guide=T) +
  geom_segment(aes(x = x_range_Bird[1], xend = x_range_Bird[2], y = start_y_Bird, yend = end_y_Bird, colour = "Birds" ), size = 1,show_guide=T) +
  geom_segment(aes(x = x_range_Mammal[1], xend = x_range_Mammal[2], y = start_y_Mammal, yend = end_y_Mammal, colour = "Mammals", ),  size = 1,show_guide=T) +
  # Adjust the color scale according to a predefined palette for accessibility
  scale_colour_manual(values = colorblind_friendly_colors, 
                      breaks =  rev(c("Hagfishes", "Lampreys", "Cartilaginous fishes", "Ray-finned fishes", "Coelacanths", "Lungfishes", "Amphibians", "Reptiles", "Birds", "Mammals"))) +
  scale_fill_manual(values = colorblind_friendly_colors, 
                    breaks =  rev(c("Hagfishes", "Lampreys", "Cartilaginous fishes", "Ray-finned fishes", "Coelacanths", "Lungfishes", "Amphibians", "Reptiles", "Birds", "Mammals"))) +
  
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size = 12),
        legend.title = element_text(colour = "steelblue", face = "bold.italic", size = 8   ),
        legend.text = element_text(face = "italic", colour = "steelblue4", size =6),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  labs(color= "Class", linetype= "Class", 
       x = "Body Temperature (°C) ",
       y = "Brain mass (residual)") +
  guides(colour = guide_legend(override.aes = list( size = 3)))  +
  #coord_fixed(ratio = 1) + # Fixed aspect ratio for scale consistency
  # Custom scales for x and y axes
  scale_x_continuous(breaks = seq(from = floor(0), to = ceiling(44), by = 5)) +
  scale_y_continuous(breaks = seq(from = floor(-0.85), to = ceiling(1.85), by = 0.5))

pdf(file= paste(savepath, "Fig.4a_2" , ".pdf", sep = "" ), width = 11/2.54, height = 11/2.54)
print(p)
dev.off()


#####


  

##Figure 4B####
ParProv_Song_2025Plot<-ParProv_Song_2025[which(ParProv_Song_2025$Class %in% c(  "Lampreys", "Cartilaginous fishes",  "Ray-finned fishes", "Amphibians", "Reptiles", "Birds", "Mammals") & !is.na(ParProv_Song_2025$Tb1)),]

sample_sizes1 <- ParProv_Song_2025Plot %>%
  filter() %>%
  group_by(Class) %>%
  summarize(
    Residual_BS_all = mean(Residual_BS_all ), 
    Residual_OS_all = median(Residual_OS_all ), 
    Tb1 = median(Tb1 ), 
    .groups = "drop" 
  )
sample_sizes1



quartile_data <- ParProv_Song_2025Plot %>%
  group_by(Class) %>%
  summarize(
    Q1_OS = quantile(Residual_OS_all, 0.25), 
    Q3_OS = quantile(Residual_OS_all, 0.75), 
    
    OS_Upper_Whisker = min(max(Residual_OS_all), Q3_OS + 1.5 * (Q3_OS - Q1_OS)), 
    OS_Lower_Whisker = max(min(Residual_OS_all), Q1_OS - 1.5 * (Q3_OS - Q1_OS)), 
    
    Q1_Tb1 = quantile(Tb1, 0.25),           
    Q3_Tb1 = quantile(Tb1, 0.75),          
    Tb1_Upper_Whisker = min(max(Tb1), Q3_Tb1 + 1.5 * (Q3_Tb1 - Q1_Tb1)), 
    Tb1_Lower_Whisker = max(min(Tb1), Q1_Tb1 - 1.5 * (Q3_Tb1 - Q1_Tb1)), 
    
    Median_OS = median(Residual_OS_all),   
    Median_Tb1 = median(Tb1),               
    .groups = "drop"
  )

quartile_data

p<-NA
p<-ggplot() +
  #geom_smooth(data = quartile_data, aes(x = Q1_OS, y= Q3_OS,   fill = Class,colour = Class), size = 0, alpha=0)+

  geom_rect(
    data = quartile_data,
    aes(
      xmin = Q1_Tb1,
      xmax = Q3_Tb1,
      ymin = Q1_OS,
      ymax = Q3_OS,
      fill = Class,
      colour = Class
    ),
    size=0.2,
    alpha = 0.3
  ) +
  
 
  geom_segment(
    data = quartile_data,
    aes(
      x = Median_Tb1,
      xend = Median_Tb1,
      y = Q1_OS,
      yend = Q3_OS,
      colour = Class
  
    ),
    size = 0.5
  ) +
  geom_segment(
    data = quartile_data,
    aes(
      x = Q1_Tb1,
      xend = Q3_Tb1,
      y = Median_OS,
      yend = Median_OS,
      colour = Class
    ),
    size = 0.5
  ) +
  
  geom_segment(
    data = quartile_data,
    aes(
      x = Q1_Tb1,  
      xend = Tb1_Lower_Whisker,
      y = (Q1_OS+ Q3_OS)/2,
      yend = (Q1_OS+ Q3_OS)/2,
      colour = Class
    ),
    size = 0.2
  ) +
  
  geom_segment(
    data = quartile_data,
    aes(
      x = Q3_Tb1,  
      xend = Tb1_Upper_Whisker,
      y = (Q1_OS+ Q3_OS)/2,
      yend = (Q1_OS+ Q3_OS)/2,
      colour = Class
    ),
    size = 0.2
  ) +
  
  #
  geom_segment(
    data = quartile_data,
    aes(
      y = Q1_OS,  
      yend = OS_Lower_Whisker,
      x = (Q1_Tb1+ Q3_Tb1)/2,
      xend = (Q1_Tb1+ Q3_Tb1)/2,
      colour = Class
    ),
    size = 0.2
  ) +
  
  geom_segment(
    data = quartile_data,
    aes(
      y = Q3_OS,  
      yend = OS_Upper_Whisker,
      x = (Q1_Tb1+ Q3_Tb1)/2,
      xend = (Q1_Tb1+ Q3_Tb1)/2,
      colour = Class
    ),
    size = 0.2
  ) +
  

  geom_point(
    data = sample_sizes1,
    aes(x = Tb1, y = Residual_OS_all, size = I(Residual_BS_all +1)*4, colour = Class ,fill = Class),
    alpha = 0.8
  ) +
  scale_colour_manual(values = colorblind_friendly_colors, 
                      breaks =  rev(c("Hagfishes", "Lampreys", "Cartilaginous fishes", "Ray-finned fishes", "Coelacanths", "Lungfishes", "Amphibians", "Reptiles", "Birds", "Mammals"))) +
  scale_fill_manual(values = colorblind_friendly_colors, 
                    breaks =  rev(c("Hagfishes", "Lampreys", "Cartilaginous fishes", "Ray-finned fishes", "Coelacanths", "Lungfishes", "Amphibians", "Reptiles", "Birds", "Mammals"))) +
  
 
  labs(
    x = "Body Temperature °C",
    y = "Newborn mass (residual)",
    colour = "Class"
  ) +
  
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size = 12),
        legend.title = element_text(colour = "steelblue", face = "bold.italic", size = 8),
        legend.text = element_text(face = "italic", colour = "steelblue4", size = 6),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  guides(colour = guide_legend(override.aes = list( size = 3)))+
  scale_x_continuous(breaks = seq(from = floor(0), to = ceiling(44), by = 5)) 
  

pdf(file= paste(savepath, "Fig.4b_2" , ".pdf", sep = "" ), width = 11/2.54, height = 11/2.54)
print(p)
dev.off()

#####


## Figure S1B ####

p<-NA
p<-ggplot(ParProv_Song_2025[which(ParProv_Song_2025$Class %ni% c("Hagfishes", "Lampreys", "Lungfishes", "Coelacanths") & !is.na(ParProv_Song_2025$AdultMass)),], 
       aes(x = log(Offspring/AdultMass), fill = Class)) +  # Mapping log-transformed ratio with Taxon plot fill
  geom_density(size = 0.2, alpha = 0.7, adjust = 1, kernel = "gaussian") +  # Density plot with Gaussian smoothing
  # Text and segment annotations for different classes, pointing out mean log values with labels and arrows
  geom_text(label = "H", x = mean(log(ParProv_Song_2025[which(ParProv_Song_2025$Class %in% c("Hagfishes") & !is.na(ParProv_Song_2025$AdultMass)), ]$Offspring/ParProv_Song_2025[which(ParProv_Song_2025$Class %in% c("Hagfishes") & !is.na(ParProv_Song_2025$AdultMass)), ]$AdultMass)),
            y = 0.5, size = 5) +
  geom_segment(aes(x = mean(log(ParProv_Song_2025[which(ParProv_Song_2025$Class %in% c("Hagfishes") & !is.na(ParProv_Song_2025$AdultMass)), ]$Offspring/ParProv_Song_2025[which(ParProv_Song_2025$Class %in% c("Hagfishes") & !is.na(ParProv_Song_2025$AdultMass)), ]$AdultMass)), xend = mean(log(ParProv_Song_2025[which(ParProv_Song_2025$Class %in% c("Hagfishes") & !is.na(ParProv_Song_2025$AdultMass)), ]$Offspring/ParProv_Song_2025[which(ParProv_Song_2025$Class %in% c("Hagfishes") & !is.na(ParProv_Song_2025$AdultMass)), ]$AdultMass)), 
                   y = 0.43, yend = 0.48), color = "black", 
               arrow = arrow(length = unit(0.2, "cm"), ends = "first")) +
  geom_text(label = "La", x = mean(log(ParProv_Song_2025[which(ParProv_Song_2025$Class %in% c("Lampreys") & !is.na(ParProv_Song_2025$AdultMass)), ]$Offspring/ParProv_Song_2025[which(ParProv_Song_2025$Class %in% c("Lampreys") & !is.na(ParProv_Song_2025$AdultMass)),]$AdultMass)),
            y = 0.5, size = 5) +
  geom_segment(aes(x = mean(log(ParProv_Song_2025[which(   ParProv_Song_2025$Class %in% c( "Lampreys") & !is.na(ParProv_Song_2025$AdultMass)), ]$Offspring/ParProv_Song_2025[which(ParProv_Song_2025$Class %in% c("Lampreys") & !is.na(ParProv_Song_2025$AdultMass)), ]$AdultMass)), xend = mean(log(ParProv_Song_2025[which(ParProv_Song_2025$Class %in% c("Lampreys") & !is.na(ParProv_Song_2025$AdultMass)), ]$Offspring/ParProv_Song_2025[which(ParProv_Song_2025$Class %in% c("Lampreys") & !is.na(ParProv_Song_2025$AdultMass)), ]$AdultMass)), 
                   y = 0.43, yend = 0.48), color = "black", 
               arrow = arrow(length = unit(0.2, "cm"), ends = "first")) +
  
  geom_text(label = "Lu" ,x = mean(log(ParProv_Song_2025[which(   ParProv_Song_2025$Class %in% c("Lungfishes") & !is.na(ParProv_Song_2025$AdultMass)), ]$Offspring/ParProv_Song_2025[which(ParProv_Song_2025$Class %in% c("Lungfishes") & !is.na(ParProv_Song_2025$AdultMass)),]$AdultMass)),
            y = 0.5, size = 5) +
  geom_segment(aes(x = mean(log(ParProv_Song_2025[which(   ParProv_Song_2025$Class %in% c("Lungfishes") & !is.na(ParProv_Song_2025$AdultMass)), ]$Offspring/ParProv_Song_2025[which(ParProv_Song_2025$Class %in% c("Lungfishes") & !is.na(ParProv_Song_2025$AdultMass)), ]$AdultMass)), xend = mean(log(ParProv_Song_2025[which(ParProv_Song_2025$Class %in% c("Lungfishes") & !is.na(ParProv_Song_2025$AdultMass)), ]$Offspring/ParProv_Song_2025[which(ParProv_Song_2025$Class %in% c("Lungfishes") & !is.na(ParProv_Song_2025$AdultMass)), ]$AdultMass)), 
                   y = 0.43, yend = 0.48), color = "black", 
               arrow = arrow(length = unit(0.2, "cm"), ends = "first")) +
  
  geom_text(label = "C" ,x = mean(log(ParProv_Song_2025[which(   ParProv_Song_2025$Class %in% c("Coelacanths") & !is.na(ParProv_Song_2025$AdultMass)), ]$Offspring/ParProv_Song_2025[which(ParProv_Song_2025$Class %in% c("Coelacanths") & !is.na(ParProv_Song_2025$AdultMass)),]$AdultMass)),
            y = 0.5, size = 5) +
  geom_segment(aes(x = mean(log(ParProv_Song_2025[which(   ParProv_Song_2025$Class %in% c("Coelacanths") & !is.na(ParProv_Song_2025$AdultMass)), ]$Offspring/ParProv_Song_2025[which(ParProv_Song_2025$Class %in% c("Coelacanths") & !is.na(ParProv_Song_2025$AdultMass)), ]$AdultMass)), xend = mean(log(ParProv_Song_2025[which(ParProv_Song_2025$Class %in% c("Coelacanths") & !is.na(ParProv_Song_2025$AdultMass)), ]$Offspring/ParProv_Song_2025[which(ParProv_Song_2025$Class %in% c("Coelacanths") & !is.na(ParProv_Song_2025$AdultMass)), ]$AdultMass)), 
                   y = 0.43, yend = 0.48), color = "black", 
               arrow = arrow(length = unit(0.2, "cm"), ends = "first")) +
  
  # Custom color scale for visual distinction of different taxonomic groups
  scale_fill_manual(values = c(
    "Cartilaginous fishes" = "#394A92",
    "Ray-finned fishes" = "#70A0AC",
    "Amphibians" = "#BAAFD1",
    "Reptiles" = "#68AC57",
    "Birds" = "#F4C28F",
    "Mammals" = "#9B3A4D")) +
  # Custom theme settings for axis text, titles, and panel background
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        legend.position = "right",
        axis.text.x = element_text(angle = 0, hjust = 1)) +
  labs(x = "log (newborn/adult body mass)", y = "Distribution", title = "") +  # Labels for axes and title
  theme(legend.title = element_blank())  # Removes the legend title

pdf(file= paste(savepath, "Fig.S1a" , ".pdf", sep = "" ))
print(p)
dev.off()

#####

##Figure S1A####

p<-NA
p<-ggplot(ParProv_Song_2025[which( ParProv_Song_2025$Class %ni% c( "Hagfishes", "Lampreys", "Lungfishes", "Coelacanths")),], 
       aes(x = log10(Offspring), fill = Class)) +
  #geom_histogram(bins = 55, alpha = 0.5, position = "identity") + 
  geom_density( size = 0.2, alpha = 0.7,adjust = 1, kernel = "gaussian") +
  geom_text(label = "H" ,x = mean(log10(ParProv_Song_2025[which(   ParProv_Song_2025$Class %in% c( "Hagfishes") ),]$Offspring)),
            y = 0.6, size = 5) +
  geom_segment(aes(x = mean(log10(ParProv_Song_2025[which(ParProv_Song_2025$Class %in% c( "Hagfishes") ),]$Offspring)), xend =  mean(log10(ParProv_Song_2025[which(   ParProv_Song_2025$Class %in% c( "Hagfishes") ),]$Offspring)), 
                   y = 0.5, yend = 0.58), color = "black", 
               arrow = arrow(length = unit(0.2, "cm"), ends = "first")) +
  
  geom_text(label = "La" ,x = mean(log10(ParProv_Song_2025[which(   ParProv_Song_2025$Class %in% c( "Lampreys") ),]$Offspring)),
            y = 0.6, size = 5) +
  geom_segment(aes(x = mean(log10(ParProv_Song_2025[which(ParProv_Song_2025$Class %in% c( "Lampreys") ),]$Offspring)), xend =  mean(log10(ParProv_Song_2025[which(   ParProv_Song_2025$Class %in% c( "Lampreys") ),]$Offspring)), 
                   y = 0.5, yend = 0.58), color = "black", 
               arrow = arrow(length = unit(0.2, "cm"), ends = "first")) +
  geom_text(label = "Lu" ,x = mean(log10(ParProv_Song_2025[which(   ParProv_Song_2025$Class %in% c( "Lungfishes") ),]$Offspring)),
            y = 0.6, size = 5) +
  geom_segment(aes(x = mean(log10(ParProv_Song_2025[which(ParProv_Song_2025$Class %in% c( "Lungfishes") ),]$Offspring)), xend =  mean(log10(ParProv_Song_2025[which(   ParProv_Song_2025$Class %in% c( "Lungfishes") ),]$Offspring)), 
                   y = 0.5, yend = 0.58), color = "black", 
               arrow = arrow(length = unit(0.2, "cm"), ends = "first")) +
  geom_text(label = "C" ,x = mean(log10(ParProv_Song_2025[which(   ParProv_Song_2025$Class %in% c( "Coelacanths")),]$Offspring)),
            y = 0.6, size = 5) +
  geom_segment(aes(x = mean(log10(ParProv_Song_2025[which(ParProv_Song_2025$Class %in% c( "Coelacanths") ),]$Offspring)), xend =  mean(log10(ParProv_Song_2025[which(   ParProv_Song_2025$Class %in% c( "Coelacanths") ),]$Offspring)), 
                   y = 0.5, yend = 0.58), color = "black", 
               arrow = arrow(length = unit(0.2, "cm"), ends = "first")) +
  scale_fill_manual(values = c(
    "Cartilaginous fishes" = "#394A92",
    "Ray-finned fishes" = "#70A0AC",
    "Amphibians" = "#BAAFD1",
    "Reptiles" = "#68AC57",
    "Birds" = "#F4C28F",
    "Mammals" = "#9B3A4D"),
    breaks = c( "Cartilaginous fishes", 
                "Ray-finned fishes", "Amphibians", "Reptiles",
                "Birds", "Mammals")
  ) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        legend.position = "right",
        axis.text.x = element_text(angle = 0, hjust = 1)) +
  labs(x = "Newborn mass (log-10)", y = "Distribution", title = "") +
  theme(legend.title = element_blank()) 
pdf(file= paste(savepath, "Fig.S1b" , ".pdf", sep = "" ))
print(p)
dev.off()

#####


##Figure S2####

sample_sizes <- ParProv_Song_2025 %>% 
  filter(Class %in% c("Ray-finned fishes")& !is.na(Fertilization)) %>%
  group_by(Fertilization ) %>%
  summarize(
    count = n(),
    max = max(Residual_OS),
    min = min(Residual_OS)
  ) 

sample_sizes$Fertilization<-factor(sample_sizes$Fertilization, levels = c("external", "internal"))
sample_sizes

ggplot(ParProv_Song_2025[which(ParProv_Song_2025$Class %in% c( "Ray-finned fishes") & !is.na(ParProv_Song_2025$Fertilization)),], 
       aes(x = Fertilization, y = Residual_OS, fill = Fertilization)) +
  geom_violin(width = 1.4) +
  geom_boxplot(width = 0.3, color = "black", alpha = 0.2) +
  scale_fill_manual(values = c("#1E466E", "#FFE6B7")) +
  scale_x_discrete(drop = FALSE, limits = c("external", "internal")) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        legend.position = "none") +
  ylim(min(sample_sizes$min, na.rm = TRUE), max(sample_sizes$max, na.rm = TRUE)+0.2*( max(sample_sizes$max, na.rm = TRUE) - min(sample_sizes$min, na.rm = TRUE))  ) + 
  ylab("New-born mass (residual)") +
  xlab("") +
  geom_text(data = sample_sizes[which(!is.na(sample_sizes$max)),], aes(x = Fertilization, y = max, label = count), size = 8, position = position_dodge(width = 0.3), vjust = -0.5) +
  geom_segment(aes(x = 1, xend = 2, y = max(sample_sizes$max, na.rm = T) + 0.6, yend = max(sample_sizes$max, na.rm = T) + 0.6), color = "black") +
  geom_segment(aes(x = 1, xend = 1, y = max(sample_sizes$max, na.rm = T) + 0.5, yend = max(sample_sizes$max, na.rm = T) + 0.6), color = "black") +
  geom_segment(aes(x = 2, xend = 2, y = max(sample_sizes$max, na.rm = T) + 0.5, yend = max(sample_sizes$max, na.rm = T) + 0.6), color = "black") +
  annotate("text", x = 1.5, y = max(sample_sizes$max, na.rm = T) + 0.6, label = " ", size = 10, vjust = 0) 

#
sample_sizes <- ParProv_Song_2025 %>% 
  filter(Class %in% c("Ray-finned fishes")& !is.na(Fertilization)) %>%
  group_by(Fertilization ) %>%
  summarize(
    count = n(),
    max = max(Residual_BS),
    min = min(Residual_BS)
  ) 
sample_sizes$Fertilization<-factor(sample_sizes$Fertilization, levels = c("external", "internal"))
sample_sizes

ggplot(ParProv_Song_2025[which(ParProv_Song_2025$Class %in% c( "Ray-finned fishes") & !is.na(ParProv_Song_2025$Fertilization)),], 
       aes(x = Fertilization, y = Residual_BS, fill = Fertilization)) +
  geom_violin(width = 1.4) +
  geom_boxplot(width = 0.3, color = "black", alpha = 0.2) +
  scale_fill_manual(values = c("#1E466E", "#FFE6B7")) +
  scale_x_discrete(drop = FALSE, limits = c("external", "internal")) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        legend.position = "none") +
  ylim(min(sample_sizes$min, na.rm = TRUE), max(sample_sizes$max, na.rm = TRUE)+0.2*( max(sample_sizes$max, na.rm = TRUE) - min(sample_sizes$min, na.rm = TRUE))  ) + 
  ylab("Brain mass (residual)") +
  xlab("") +
  geom_text(data = sample_sizes[which(!is.na(sample_sizes$max)),], aes(x = Fertilization, y = max, label = count), size = 8, position = position_dodge(width = 0.3), vjust = -0.5) +
  geom_segment(aes(x = 1, xend = 2, y = max(sample_sizes$max, na.rm = T) +0.2, yend = max(sample_sizes$max, na.rm = T) +0.2), color = "black") +
  geom_segment(aes(x = 1, xend = 1, y = max(sample_sizes$max, na.rm = T)  +0.16, yend = max(sample_sizes$max, na.rm = T) +0.2), color = "black") +
  geom_segment(aes(x = 2, xend = 2, y = max(sample_sizes$max, na.rm = T) +0.16, yend = max(sample_sizes$max, na.rm = T)+0.2), color = "black") +
  annotate("text", x = 1.5, y = max(sample_sizes$max, na.rm = T) +0.2, label = "*", size = 10, vjust = 0) 



sample_sizes <- ParProv_Song_2025 %>% 
  filter(Class %in% c("Ray-finned fishes")& !is.na(Fertilization)) %>%
  group_by(Fertilization ) %>%
  summarize(
    count = n(),
    max = max(log10(AdultMass)),
    min = min(log10(AdultMass))
  ) 
sample_sizes$Fertilization<-factor(sample_sizes$Fertilization, levels = c("external", "internal"))
sample_sizes

ggplot(ParProv_Song_2025[which(ParProv_Song_2025$Class %in% c( "Ray-finned fishes") & !is.na(ParProv_Song_2025$Fertilization)),], 
       aes(x = Fertilization, y = log10(AdultMass), fill = Fertilization)) +
  geom_violin(width = 1.4) +
  geom_boxplot(width = 0.3, color = "black", alpha = 0.2) +
  scale_fill_manual(values = c("#1E466E", "#FFE6B7")) +
  scale_x_discrete(drop = FALSE, limits = c("external", "internal")) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        legend.position = "none") +
  ylim(min(sample_sizes$min, na.rm = TRUE), max(sample_sizes$max, na.rm = TRUE)+0.2*( max(sample_sizes$max, na.rm = TRUE) - min(sample_sizes$min, na.rm = TRUE))  ) + 
  ylab("Body mass (log-10)") +
  xlab("") +
  geom_text(data = sample_sizes[which(!is.na(sample_sizes$max)),], aes(x = Fertilization, y = max, label = count), size = 8, position = position_dodge(width = 0.3), vjust = -0.5) +
  geom_segment(aes(x = 1, xend = 2, y = max(sample_sizes$max, na.rm = T) + 0.8, yend = max(sample_sizes$max, na.rm = T) +0.8), color = "black") +
  geom_segment(aes(x = 1, xend = 1, y = max(sample_sizes$max, na.rm = T)  + 0.6, yend = max(sample_sizes$max, na.rm = T) +0.8), color = "black") +
  geom_segment(aes(x = 2, xend = 2, y = max(sample_sizes$max, na.rm = T) + 0.6, yend = max(sample_sizes$max, na.rm = T)+0.8), color = "black") +
  annotate("text", x = 1.5, y = max(sample_sizes$max, na.rm = T) + 0.88, label = ".", size = 12, vjust = 0) 



#####

##Figure S3####

ggplot(ParProv_Song_2025[which(   ParProv_Song_2025$Class %in% c(  "Cartilaginous fishes" ,  "Ray-finned fishes" ) & !is.na(ParProv_Song_2025$AdultMass)),], aes(x = log10(AdultMass), fill = Class)) +
  #geom_histogram(bins = 55, alpha = 0.5, position = "identity") + 
  geom_density( size = 0.2, alpha = 0.7,adjust = 1, kernel = "gaussian") +
  
  scale_fill_manual(values = c(
    "Cartilaginous fishes" = "#394A92",
    "Ray-finned fishes" = "#70A0AC"
  )  
  ) +
  #scale_color_manual(values = c("Cartilaginous fishes" = "#1E466E", "Ray-finned fishes" = "#F7AA58")) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        legend.position = "right",
        legend.text = element_text(size = 16), 
        axis.text.x = element_text(angle = 0, hjust = 1)) +
  labs(x = "Adult body mass (log-10)", y = "Distribution", title = "") +
  theme(legend.title = element_blank())




ggplot(ParProv_Song_2025[which(   ParProv_Song_2025$Class %in% c(  "Cartilaginous fishes" ,  "Ray-finned fishes" ) & !is.na(ParProv_Song_2025$Tb)),],
       aes(x = Class,y = Tb, fill = Class)) +
  geom_violin(width = 0.8) +
  geom_boxplot(width = 0.3, color = "black", alpha = 0.2) +
  scale_fill_manual(values = c("Cartilaginous fishes" = "#394A92" , "Ray-finned fishes" =   "#70A0AC"   )) +
  scale_x_discrete(drop = FALSE, limits = c( "Cartilaginous fishes", "Ray-finned fishes")) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        legend.position = "none",
        axis.text.x = element_text()) + 
  geom_text(aes(x = 1, y = 16.45, label = "16.5"), size = 6,  vjust = -0.5, col = "white") +
  geom_text(aes(x = 2, y = 13.30, label = "13.3"), size = 6,  vjust = -0.5, col = "white") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") + 
  ylab("Ambient temperature (°C)") +
  xlab("")

####




