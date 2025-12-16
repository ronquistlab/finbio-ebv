library(tidyverse)
library(xgboost)
library(patchwork)
library(ggh4x)

source("R/functions.R")

# Load data -----------------------------------------------------------------------------------
resFiles <- list.files("results/",full.names = TRUE)
resList  <- resFiles |> purrr::map(~readRDS(.x), .id = "source")

# plot ----------------------------------------------------------------------------------------

# loop over results
outList <- list()
for(i in seq_along(resList)){
  
  results     <- resList[[i]]
  y           <- results$data$y
  tmod_loso   <- results$model_fits$loso
  tmod_kfold   <- results$model_fits$kfold
  tmod_strat  <- results$model_fits$strat
  loss_loso   <- get_fold_loss(tmod_loso,y)
  loss_kfold  <- get_fold_loss(tmod_kfold,y) |> mutate(fold = as.character(fold))
  loss_strat  <- get_fold_loss(tmod_strat ,fit_data = results$data$full , stratified = TRUE ,y) |>
                 mutate(fold = as.character(fold))
  
  
  outList[[i]] <- bind_rows(list("K-fold"=loss_kfold,"New site"=loss_loso , "Within site" = loss_strat),.id="cv_type") 
  
}

names(outList) <- list.files("results")

loss_df_all <- bind_rows(outList,.id="model_fit") |> 
              mutate(country = str_extract(model_fit , "sweden|madagascar"),
                     ebv = str_extract(model_fit , "(?<=_)[^_.]+(?=\\.)")) |> 
              dplyr::select(-model_fit) |> 
              filter(cv_type != "K-fold") |> 
              mutate(cv_type = fct_relevel(cv_type, "Within site" , "New site"))

  
loss_df_all$country <- factor(loss_df_all$country ,
                            levels = c("sweden","madagascar"),
                            labels = c("Sweden","Madagascar"))

loss_df_all$ebv <- factor(loss_df_all$ebv , levels = c("SR","LCBD","FD","FE","GSH"))

loss_sum <- loss_df_all |> 
  pivot_longer(!matches("cv|fold|country|ebv")  , names_to = "metric" , values_to = "value") |> 
  summarise(mean_value = median(value) , 
            qi_80_upper = quantile(value , probs=c(0.1,0.9))[2] ,
            qi_80_lower = quantile(value , probs=c(0.1,0.9))[1] ,
            qi_50_upper = quantile(value , probs=c(0.25,0.75))[2] ,
            qi_50_lower = quantile(value , probs=c(0.25,0.75))[1] ,
            
            .by=c("country","ebv","metric","cv_type")) |> 
      filter(metric == "MAE")


# plot ----------------------------------------------------------------------------------------


loss_tmp <- loss_df_all |> 
            filter(ebv != "PD") %>% 
            pivot_longer(!matches("cv|fold|country|ebv")  , names_to = "metric" , values_to = "value") |> 
            filter( metric=="MAE") 
            

loss_tmp_sum <- loss_sum |> 
  filter(ebv != "PD") %>% 
  
                  mutate(pc_error_change = ((first(mean_value) - last(mean_value)) / last(mean_value)) * 100,.by=c("ebv" ,"country" ))


ann_dat <- loss_tmp_sum |>
  filter(ebv != "PD") %>% 
  
            dplyr::select(ebv,cv_type ,pc_error_change,mean_value,country) |> 
            filter(cv_type == "New site")



p1 <- loss_tmp |> 
  ggplot()+
  #geom_jitter(aes(cv_type , value),width = 0.05 , alpha = .3,size=5,colour="grey")+
  geom_errorbar(data=loss_tmp_sum  , aes(x = cv_type, ymin = qi_80_lower , ymax = qi_80_upper ),width=0 ,size=4, colour="black")+
  geom_errorbar(data=loss_tmp_sum  , aes(x = cv_type, ymin = qi_50_lower , ymax = qi_50_upper ),width=0 ,size=7, colour="darkred")+
  geom_line(data=loss_tmp_sum,aes(cv_type , mean_value,group=1),lty=3 , size=2)+
  geom_point(data=loss_tmp_sum , aes(cv_type, mean_value),size=12)+
  geom_point(data=loss_tmp_sum , aes(cv_type, mean_value),size=8, colour="white")+
  theme_linedraw(base_size = 40)+
  theme(axis.text.x = element_text(angle = 90))+
  geom_text(data=ann_dat , aes(x=2 , y = mean_value , label = round(pc_error_change)),size=10, fontface="bold" , hjust = -1)+
  facet_grid(ebv~country,scales="free_y")+
  labs(x = "CV type" , y = "Error")

tiff("plots/prelim_results" , width = 1200,height=1800,compression = "lzw")
p1 
dev.off()

browseURL("plots/prelim_results")
