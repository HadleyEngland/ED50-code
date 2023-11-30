##### Essential libraries ####
library(drc)
library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(stringr)
library(ggh4x)

### read data into global environment ### 
setwd("~/Desktop")
input <- read.table("Phenoplate Data.txt", sep = "\t", header = T)

#change temperature and PAM values to numeric values
input$Temperature=as.numeric(input$Temperature) 
input$PAM=as.numeric(input$PAM) 

#remove rows with missing data if you have any.
input<-input[complete.cases(input), ]

#view input
View(input)

#make a new column with a unique sample ID. Drc will use it to group all PAM values from one colony (biological replicate) per species together and calculate the ED50. 
input$Sample=as.factor(paste(input$Species,input$Colony, sep = "_"))
levels(input$Sample)
View(input)

#Now we use the drc package to calculate the ED50 for each colony. 
#Demo to one colony with limits (upper and lower temperature limits can need to be adjusted to data)
drm(PAM ~ Temperature, data=input[input$Sample=="Acro_1",],
    fct = LL.3(names = c('Slope', 'Max', 'ED50')), 
    upperl = c(120, 0.72, 42), lowerl = c(10, 0.55, 20))


#For more colonies, it makes sense to run this function as a loop and calculate the ED50s for each replicate in one go. This will return a list of models saved under the variable mod1.
mod1<-lapply(unique(input$Sample), 
             function(x) drm(PAM ~ Temperature, data=input[input$Sample==x,],
                             fct = LL.3(names = c('Slope', 'Max', 'ED50')),
                             upperl = c(120, 0.72, 42), lowerl = c(10, 0.55, 20)))
View(mod1)

#For comparison: model without limits - usually helps to fit the datapoints to the curve. However, here it results in very high ED50s. Only run this command with the rest of the script if no limits should be applied. 
mod1 <- lapply(unique(input$Sample), 
               function(x) drm(PAM ~ Temp, data = input[input$Sample == x,],
                               fct = LL.3(names = c('Slope', 'Max', 'ED50'))))



#Now we need to extract the ED50s out of the list of models we just created.
ed50_list <- lapply(c(1:length(mod1)), function(x) mod1[[x]][["coefficients"]][["ED50:(Intercept)"]])
View(ed50_list)

ed50_df <- as.data.frame(do.call(rbind, ed50_list))              #make a new dataframe with the extracted ED50 values.
ed50_df$Sample = unique(input$Sample)                            #add sample information to new dataframe.
ed50_df = tidyr::separate(ed50_df,Sample,into = c("Species", "Colony" ),sep = "_",remove = FALSE,extra = "merge") #seperate sample information into species and colony and add them as new columns to the dataframe.
ed50_df$Colony = factor(ed50_df$Colony) #make variable colony into a factor. Levels can be reordered. Default is alphabetic. 
colnames(ed50_df)[1] = "ED50" 
levels(ed50_df$Colony)
View(ed50_df)



#### Model ED50 regression curves and make a plot ####
temp_x <- seq(22, 40, length = 100) #print 100 values from 27 to 39.

#prediction of the fitted values corresponding to the range of temperatures 
pred1 <- lapply(mod1, function(x) predict(x, data.frame(Temperature = temp_x))) #use list of models from before, predict() predicts all values from x/y to form our curve.
pred_df <- as.data.frame(do.call(rbind, pred1)) #create a new dataframe.
colnames(pred_df) = round(temp_x, digits = 2)   #rename column names.
pred_df$Sample = unique(input$Sample)  #add column with sample names.          
pred_df_long = reshape2::melt(pred_df, id.vars = c("Sample")) #rearrange from short (length)/long(width) to long (length)/short (width) table with melt() for plotting.
colnames(pred_df_long)[2:3] <- c("Temperature", "Fv/Fm")      #rename column names
pred_df_long = tidyr::separate(pred_df_long,Sample,into =c("Species", "Colony"),sep = "_",remove = FALSE,extra = "merge")  #split column into variables. 
pred_df_long$Temperature = as.numeric(as.character(pred_df_long$Temperature)) #make temperature as numeric 
pred_df_long$group = paste(pred_df_long$Species, pred_df_long$Colony)     

View(pred_df_long)

#calculate population ED50s per System
ED50_means<-ed50_df %>% 
  group_by(Species, Colony) %>%
  summarise(mean=mean(ED50), sd=sd(ED50)) %>%
  unite(Group, c(Colony), sep = "-", remove = FALSE)
ED50_means$group = paste(ED50_means$Species, ED50_means$Colony)
pred_df_long$meanED50=round(ED50_means$mean[match(pred_df_long$group,ED50_means$group)], 2)
curve_input <- pred_df_long #save pred_df_long as curve_input for plotting.
input$Colony <- as.factor(input$Colony)
curve_input_meanED50s <- curve_input[c(1:6),]

View(curve_input_meanED50s)

#change species names for plotting 
input$Species <- gsub("Acro", "Acropora Tenuis", input$Species)
input$Species <- gsub("Poci", "Pocillopora Damicornis", input$Species)
View(input$Species)

curve_input$Species <- gsub("Acro", "Acropora Tenuis", curve_input$Species)
curve_input$Species <- gsub("Poci", "Pocillopora Damicornis", curve_input$Species)
View(curve_input$Species)

curve_input_meanED50s$Species <- gsub("Acro", "Acropora Tenuis", curve_input_meanED50s$Species)
curve_input_meanED50s$Species <- gsub("Poci", "Pocillopora Damicornis", curve_input_meanED50s$Species)
View(curve_input_meanED50s$Species)


####plot ED50 regression curves.####
ED50curve <-  ggplot() +
              geom_line(data = curve_input, aes(x=Temperature, y=`Fv/Fm`, group=Colony, color = factor(Colony)), size = 1, show.legend = F) +
              geom_segment(data = curve_input, aes(x = meanED50, y = 0, xend = meanED50, yend = 0.65, color = Colony), linetype=3, size=1, show.legend = F) +
              ggrepel::geom_text_repel(aes(x = meanED50, y = 0.78, label = meanED50, color = Colony, fontface = 2), data = curve_input_meanED50s, size= 6, angle = 90, hjust = 0, max.overlaps = 10, direction = "x", point.size = NA, segment.color = NA, show.legend = F) +
              scale_x_continuous(breaks=c(24, 28, 32, 36, 40), limits = c(22, 42), expand = c(0, 0)) +
              scale_y_continuous(limits = c(-0.02,0.78), expand = c(0, 0)) + labs(color='') +
              geom_jitter(data = input, aes(x = Temperature, y = PAM, group=Colony, color=Colony, shape=Colony), size = 3, width = 0.25) +
              theme_classic() +
              scale_shape_manual(name = "Colony",
                     labels = c("1","2","3"),
                     values = c(15, 15, 15)) +
              scale_color_manual(name = "Colony",
                     labels = c("1","2","3"),
                     values = c("#FDA550", "#A799B7", "#9ADEFE")) +
              labs(y = "Photosynthetic Efficiency (Fv/Fm)", x = "Temperature (°C)") +
              facet_nested(~ Species, scales ="free_x" )


ED50curve_theme <- ED50curve + theme (legend.position= "bottom", 
                                      legend.title = element_text(colour="black", size=13,face="bold"),
                                      legend.text=element_text(size = 13),
                                      line = element_line(size = 0.8),
                                      axis.line = element_line(colour = "black"),
                                      axis.ticks = element_line(colour = "black"),
                                      axis.ticks.length = unit(0.2 , "cm"),
                                      axis.text = element_text(size = 13, colour = "black"),
                                      text = element_text(size = 13, colour = "black"),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      panel.background = element_blank(),
                                      strip.background = element_blank(), 
                                      strip.text.x = element_text(color = "black", size = 12, angle = 0, hjust = 0.5, vjust = 0.5, face = "italic"))
ED50curve_theme #view your finished plot.

#safe plots as a pdf
pdf(file="~/Desktop/Phenoplate ED50",
    10, 8)

ED50curve_theme
dev.off()



#### boxplot ####
species_colors<-c("#CE6479", "#708599") 

ed50_df$Species <- gsub("Acro", "Acropora Tenuis", ed50_df$Species)
ed50_df$Species <- gsub("Poci", "Pocillopora Damicornis", ed50_df$Species)
View(ed50_df)

b1 <- ggplot(ed50_df, aes(x=Species, y=ED50, fill = Species)) +
  stat_boxplot(geom = "errorbar", width = 0.1) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, lwd=0.8) + theme_bw() +
  scale_fill_manual(values=species_colors) 

b1

b1_plot <- b1+ theme (legend.position= "none", 
                      legend.title = element_text(colour="black", size=13,face="bold"),
                      legend.text=element_text(size = 13),
                      line = element_line(size = 0.8),
                      axis.line = element_line(colour = "black"),
                      axis.ticks = element_line(colour = "black"),
                      axis.ticks.length = unit(0.2 , "cm"),
                      axis.text.y = element_text(size = 13, colour = "black"),
                      axis.text.x = element_text(size = 13, colour = "black", face = "italic"),
                      text = element_text(size = 13, colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank())
b1_plot

pdf(file="~/Desktop/Phenoplate ED50",
    6, 8)
b1_plot
dev.off()
