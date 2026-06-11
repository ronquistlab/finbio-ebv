library(tidyverse)
library(xgboost)
library(patchwork)
library(ggh4x)

source("R/functions.R")

# Load data -----------------------------------------------------------------------------------
full_M <- readRDS("data/tidydata/full_M.rds")
resFiles <- list.files("results/",full.names = TRUE)
resList  <- resFiles |> purrr::map(~readRDS(.x), .id = "source")
names(resList) <- resFiles


# Get feature importance 
get_fimp <- function(model_list){
  xgboost::xgb.importance(model=model_list$model_fits$full)
}



# Run over results list , tidy, and summarise for plot
fimp_list <- lapply(X = resList , FUN = get_fimp) %>% 
             bind_rows(.id = "model") %>% 
             mutate(model = str_remove_all(model , "results//results_|.rds"),
                    country = str_extract(model , "sweden|madagascar"),
                    ebv = str_extract(model,"SR|LCBD|FD|FE|GSH|PD")) %>% 
            filter(ebv != "PD",Feature!= "sample_time")

# Rename Features
nm_key <- c("Low LAI" , "High LAI" ,"Temp max" , "Temp min", "Precipitaion","Photoperiod",
            "Grass cover","Other forest cover","Shrub cover",
            "Sin(week)", "DBF cover" , "Cos(week)",
            "EBF cover" ,"Crop cover", "ENF cover", 
            "Moss cover","Urban cover" , "Mixed forest cover", "Water cover")

names(nm_key) <- unique(fimp_list$Feature)
fimp_list <- fimp_list %>% mutate(Feature   = recode_factor(factor(Feature) , !!!nm_key))
fimp_list$Feature <- factor(fimp_list$Feature , 
                            levels = c("Temp max" , "Temp min","Precipitaion", "Photoperiod",
                                       "Cos(week)" , "Sin(week)",
                                       "Low LAI" , "High LAI",
                                       "DBF cover",  "EBF cover", "ENF cover"    , "Crop cover","Mixed forest cover","Other forest cover",
                                       "Grass cover" , "Moss cover", "Shrub cover" , "Urban cover" , "Water cover"))

# Hacky solution to get bar charts as segments (rather than summed via stacking)
fimp_sum <- fimp_list %>% 
  group_by(ebv,Feature) %>% 
  arrange(-Gain) %>% 
  summarise(ymax = max(Gain), 
            ymin = min(Gain),
            cmax = first(country),
            cmin = last(country)) %>% 
            ungroup() 
  

# Plot
p1 <- ggplot(fimp_sum , aes(x = Feature)) +
  geom_segment(aes(yend = ymin , y = ymax,colour = cmax), lwd = 8, size = 3)+
  geom_segment(aes(yend = 0 , y = ymin,colour=cmin), lwd = 8, size = 3)+
  coord_flip()+
  theme_linedraw(base_size = 35)+
  facet_wrap(~ebv,scales = "free_x" , nrow = 1)+
  theme(legend.position = "bottom")+
  scale_colour_manual(values = c("darkred" , "steelblue3"))+
  labs(y = "Feature importance" , colour = "Country")
  scale_y_continuous(limits = c(0,0.175))

p2 <- ggplot(fimp_list %>% filter(country == "sweden"))+
  geom_col(stat =  "identity", aes(x = Feature , y = Gain , fill = country) , width = 0.5)+
  coord_flip()+
  theme_linedraw()+
  facet_wrap(~ebv, nrow = 1)+
  theme(axis.text.y = element_text(size = 10),
        legend.position = "bottom")+
  scale_fill_manual(values = c("steelblue3"))+
  labs(y = "Feature importance")+
  scale_y_continuous(limits = c(0,0.3))

p3 <- ggplot(fimp_list %>% filter(country == "madagascar"))+
  geom_col(stat =  "identity", aes(x = Feature , y = Gain , fill = country) , width = 0.5)+
  coord_flip()+
  theme_linedraw()+
  facet_wrap(~ebv, nrow = 1)+
  theme(axis.text.y = element_text(size = 10),
        legend.position = "bottom")+
  scale_fill_manual(values = c("darkred"))+
  labs(y = "Feature importance")+
  scale_y_continuous(limits = c(0,0.3))


tiff("plots/feat_imp.tif" , width = 2000,height=1000)
p1 
dev.off()

browseURL("plots/feat_imp.tif")
