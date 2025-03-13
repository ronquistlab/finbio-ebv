library(janitor)
library(tidyverse)
library(scales)
# data ----------------------------------------------------------------------------------------


# Populations and community EBVs
ebv_data <- read_delim("~/Downloads/example_cs_EBVs.csv" , delim = ";" , ) |> clean_names()
ebv_data$spatial_scale <- factor(ebv_data$spatial_scale , levels = c("Low" , "Medium", "High"))
ebv_data$taxon <- factor(ebv_data$taxon , levels = c("Birds" , "Tetrapods" , "Invertebrates",
                                                     "Plants" , "Fungi"))
ebv_data$taxonomic_scale <- factor(ebv_data$taxonomic_scale , levels = c("Very Low", "Low", "Medium", "High"))
ebv_data$survey_method <- factor(ebv_data$survey_method , 
                                 levels = c("Human observation" , "eDNA" , "Terrestrial image sensors" , "PAM" , "Earth observation"))

# Ecosystem EBVs
eco_data <- read_delim("~/Downloads/ecosystem_ebvs.csv"  , delim=";") |> clean_names()
eco_data$spatial_scale <- factor(eco_data$spatial_scale , levels = c("Low" ,"Medium", "High"))
eco_data$survey_method <- factor(eco_data$survey_method , 
                                 levels = c("Human observation" , "eDNA" , "Terrestrial image sensors" , "PAM" , "Earth observation"))


# plot ----------------------------------------------------------------------------------------


# Pop & com
ramp <- colour_ramp(c("lightblue" , "red3")) 
pal  <- ramp(seq(0,1,l=3)) 
p1 <- ggplot(ebv_data , aes(taxon, ebv))+
  geom_point(aes(fill= taxonomic_scale,size=spatial_scale) ,pch=21,stroke = 1)+
  scale_size_discrete(range = c(2,10) , na.translate = FALSE)+
  scale_fill_manual(values = pal , na.translate = FALSE)+
  theme_linedraw()+
  geom_hline(yintercept = c(0.5, 1.5 , 2.5 , 3.5, 4.5 , 5.5 , 6.5) , lty = 3)+
  geom_vline(xintercept = c(0.5, 1.5 , 2.5 , 3.5, 4.5, 5.5 , 6.5) , lty = 3)+
  theme(axis.text.x = element_text(angle = 65,hjust=0.95),
        axis.title.x = element_text(vjust = -2.5), 
        axis.title =element_text(size = 15),
        axis.text = element_text(size = 20),
        legend.position = "bottom" , 
        legend.title.position = "top",
        legend.title.align=0.5,
        legend.text = element_text(size = 20),
        strip.text = element_text(size = 20),
        panel.grid.major = element_blank())+
  facet_grid(class~survey_method ,scales = "free_y" , space = "free",switch="y")+
  guides(fill = guide_legend(override.aes = list(size = 7) ) )+
  labs(x= "Taxon" , y = "Essential biodiversity variable" , size = "Spatial scale" , fill = "Taxonomic coverage")



# Ecosystem
ggplot(eco_data , aes(" ",ebv_name))+
  geom_point(aes(size=spatial_scale , fill = ecosystem_coverage),pch=21, colour = "black" , stroke = 1)+
  scale_fill_manual(values = pal , na.translate = FALSE)+
  scale_fill_discrete(na.value = 'white' , na.translate  = FALSE)+
  scale_size_discrete(range = c(4,8) , na.translate = FALSE)+
  theme_linedraw()+
  geom_hline(yintercept = c(0.5, 1.5 , 2.5 , 3.5, 4.5 , 5.5 , 6.5) , lty = 3)+
  geom_vline(xintercept = c(0.5, 1.5 , 2.5 , 3.5, 4.5, 5.5 , 6.5) , lty = 3)+
  theme(axis.text.x = element_text(angle = 65,hjust=0.95),
        axis.title.x = element_text(vjust = -2.5), 
        axis.title =element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.position = "bottom" , 
        legend.title.position = "top",
        legend.title.align=0.5,
        legend.text = element_text(size = 20),
        strip.text = element_text(size = 10),
        panel.grid.major = element_blank())+
  facet_grid(ebv_class~survey_method ,scales = "free" , space = "free",switch="y")+
  guides(fill = guide_legend(override.aes = list(size = 7) ) )+
  labs(x= " " , y = "Essential biodiversity variable" , size = "Spatial scale" , fill = "Ecosystem coverage")


library(patchwork)
summary <- ebv_data |> 
  group_by(ebv , taxon,class) |> 
  summarise(spatial_scale = sum(as.numeric(spatial_scale) , na.rm = TRUE),
            taxonomic_scale = sum(as.numeric(taxonomic_scale) , na.rm=TRUE))


summary |> 
  ggplot(aes(taxon , ebv))+
  geom_point(aes(fill = taxonomic_scale,size=spatial_scale),pch=21, colour = "black")+
  scale_fill_viridis_c(option = "turbo" , direction = 1)+
  scale_size_continuous(range = c(1,10))+
  theme_linedraw()+
  theme(axis.text.x = element_text(size =10, angle = 90))+
  facet_grid(class~.,scales = "free_y" , space = "free")



# plot ----------------------------------------------------------------------------------------


tiff("ebv_cs.tiff" , width = 1500 , height = 1000 , compression = 'lzw')
p1
dev.off()

browseURL("ebv_cs.tiff")
