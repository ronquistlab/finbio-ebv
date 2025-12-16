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

# function ----------------------------------------------------------------

get_rank_change <- function(resDat,list_name){
resdf <- data.frame(
          obs         = resDat$data$full[,2],
          preds_strat = resDat$model_fits$strat$pred,
          preds_loso  = resDat$model_fits$loso$pred,
          ebv         = colnames(resDat$data$full)[2],
          country     = str_extract(list_name,"sweden|madagascar"),
          trap_id     = resDat$data$full[,1]
  ) %>% 
  summarise(
    across(
    obs:preds_loso,
    mean,
    .names = "mean_{.col}"), .by = c("trap_id","country","ebv")) %>% 
  mutate(
    across(
      mean_obs:mean_preds_loso, 
      rank, 
      .names = "rank_{.col}")
  ) %>% 
  mutate(
    across(
      rank_mean_obs:rank_mean_preds_loso, 
      ~ ntile(.x, 4),
      .names = "group_{.col}") , .by = ebv)
 }


# get rank changes --------------------------------------------------------

dd <- mapply(get_rank_change , resList,resFiles,SIMPLIFY = FALSE) %>%
    bind_rows() %>% 
    filter(ebv != "PD") 

ebv_d_rank <- dd %>% 
  mutate(c_rank_obs = rank_mean_obs - rank_mean_obs , 
         c_rank_strat = rank_mean_preds_strat - rank_mean_obs,
         c_rank_loso = rank_mean_preds_loso - rank_mean_obs,
        ) %>% 
  pivot_longer( cols = starts_with("c_"), names_to = "observation" , values_to = "rank_value") %>% 
mutate(observation = fct_relevel(observation , "c_rank_obs" , "c_rank_strat" , "c_rank_loso")) %>% 
mutate(observation = fct_recode(observation , "Observed"="c_rank_obs" , "Within-site"="c_rank_strat" , "New site" = "c_rank_loso")) %>% 
  select(country,ebv,trap_id , observation,rank_value)


ebv_d_group <- dd %>% 
  mutate(c_group_obs = group_rank_mean_obs - group_rank_mean_obs , 
         c_group_strat = group_rank_mean_preds_strat - group_rank_mean_obs,
         c_group_loso = group_rank_mean_preds_loso - group_rank_mean_obs,
  ) %>% 
  pivot_longer( cols = starts_with("c_"), names_to = "observation" , values_to = "group_value") %>% 
  mutate(observation = fct_relevel(observation , "c_group_obs" , "c_group_strat" , "c_group_loso")) %>% 
  mutate(observation = fct_recode(observation , "Observed"="c_group_obs" , "Within-site"="c_group_strat" , "New site" = "c_group_loso")) %>% 
  select(country,ebv,trap_id , observation,group_value) %>% 
  group_by(ebv,country,observation,group_value) %>% 
  mutate(count=n()) %>% 
  ungroup() 
  

# plot --------------------------------------------------------------------


ebv_names <- list(
  'FDis'="Fun Dispersal",
  'FEve'="Fun Eveness",
  'lcbd'="LCBD",
  'nOTU'="Species richness",
  'mean_shn' = "Genetic diversity"
)

ebv_labeller <- function(variable,value){
  return(ebv_names[value])
}


p1 <- ebv_d_rank %>% filter(country == "sweden") %>% 
  ggplot(aes(observation , rank_value,group=trap_id))+
  geom_line(lty=1,alpha=.3)+
  geom_point(size=3, aes(fill = rank_value))+
  theme_linedraw(base_size = 25)+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~ebv,scales="free_y",labeller = ebv_labeller,nrow = 1)+
  labs(x = "Observation" , y = "Change in rank")


p2 <- ebv_d_rank %>% filter(country == "sweden") %>% 
  ggplot(aes(observation , rank_value,group=trap_id))+
  geom_line(lty=1,alpha=.3)+
  geom_point(size=3, aes(fill = rank_value))+
  theme_linedraw(base_size = 25)+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~ebv,scales="free_y",labeller = ebv_labeller,nrow = 1)+
  labs(x = "Observation" , y = "Change in rank")


# group change ------------------------------------------------------------

ebv_d_group$country <- factor(ebv_d_group$country , 
                              levels = c("sweden","madagascar"),
                              labels = c("Sweden" , "Madagascar"))

ebv_d_group$ebv <- factor(ebv_d_group$ebv , 
                              levels = c("nOTU","lcbd" , "FDis" , "FEve" , "mean_shn"),
                              labels = c("SR" , "LCBD" , "FD","FE","GSH"))



ebv_pr_rank_g  <- ebv_d_group %>% 
                  select(-trap_id) %>% distinct() %>% 
                  group_by(country , ebv, observation) %>% 
                  mutate(count_s =  count / sum(count))

ebv_g_sum <-  ebv_d_group %>% 
  select(-trap_id) %>% 
  distinct() %>% 
  mutate(sum_count=sum(count),
         pc = (count / sum_count)*100,
         .by = c("ebv","country","observation"))

 ebv_seg <- ebv_d_group %>% select(-count)  %>% 
  mutate(pattern = paste(group_value,collapse = " "), 
         .by = c("country" , "ebv" , "trap_id")) %>% 
  select(-trap_id) %>% 
  group_by(ebv , country) %>% 
  count(pattern) %>% 
   mutate(id = row_number()) %>% 
  separate(pattern, into = c("Observed" , "Within-site" , "New site") , sep = " ") %>%
  pivot_longer(c("Observed" , "Within-site" , "New site") , 
               names_to = "observation",
               values_to ="group_value") %>% 
   mutate(group_value = as.numeric(group_value)) %>% 
   ungroup() %>% 
   mutate(n = n/sum(n) , .by = "country")
          

  
p1 <- ebv_pr_rank_g %>% 
ggplot(aes(observation , group_value))+
  geom_line(data =ebv_seg,aes(group = id,lwd = n))+
  geom_point(pch=21, aes(fill = factor(group_value),size=count_s))+
  geom_text(data = ebv_g_sum ,
            aes(observation , group_value , label = round(pc), group=NULL),
            size=8 , colour = "white")+
  theme_linedraw(base_size = 40)+
  theme(axis.text.x = element_text(angle = 90),
        legend.position="none",aspect.ratio = 1.1)+
  scale_fill_viridis_d(option="F",end=0.9)+
  scale_size(range = c(10, 30)) +
  facet_grid(country~ebv)+
  scale_y_continuous(limits = c(-3.2, 3.2))+
  labs(x = "Observation" , y = "Change in priority tier")

# plot --------------------------------------------------------------------

library(patchwork)
tiff("plots/rank_changes.tif" , width = 2000,height=1000)
p1 + theme(panel.spacing = unit(2, "lines"))
dev.off()

browseURL("plots/rank_changes.tif")


# sankey diagrams ---------------------------------------------------------

library(ggsankey)

df_w <-dd %>% 
  mutate(c_group_obs = group_rank_mean_obs - group_rank_mean_obs , 
         c_group_strat = group_rank_mean_preds_strat - group_rank_mean_obs,
         c_group_loso = group_rank_mean_preds_loso - group_rank_mean_obs
  ) 

df_split <- dd %>%
  mutate(c_group_obs = group_rank_mean_obs - group_rank_mean_obs , 
                          c_group_strat = group_rank_mean_preds_strat - group_rank_mean_obs,
                          c_group_loso = group_rank_mean_preds_loso - group_rank_mean_obs) %>%
  split(.$ebv)

df_long_all <- purrr::imap_dfr(
  .x = df_split,
  .f = ~ make_long(.x, c_group_obs:c_group_loso) %>%
    mutate(ebv = .y)  # .y is the ebv group name
)

 ggplot(df_long_all, aes(x = x, 
               next_x = next_x, 
               node = node, 
               next_node = next_node,
               fill = factor(node),
               label = node)) +
  geom_sankey(flow.alpha = 0.6, node.color = 1) +
  geom_sankey_label(size = 3.5, color = 1, fill = "white") +
  scale_fill_viridis_d(option='A')+
  theme_sankey(base_size = 16) +
  guides(fill = guide_legend(title = "Title"))+
  theme(legend.position = "none")+
   facet_wrap(~ebv,nrow=1)

  
df_loso <- dd %>%
  filter(ebv == "nOTU",country=="sweden") %>% 
  make_long(rank_group_rank_mean_obs,rank_group_rank_mean_preds_loso)

p2 <- ggplot(df_loso, aes(x = x, 
               next_x = next_x, 
               node = node, 
               next_node = next_node,
               fill = factor(node),
               label = node)) +
  geom_sankey(flow.alpha = 0.6, node.color = 1) +
  geom_sankey_label(size = 3.5, color = 1, fill = "white") +
  scale_fill_viridis_d(option='A') +
  theme_sankey(base_size = 16) +
  guides(fill = guide_legend(title = "Title"))+
  theme(legend.position = "none")


p1 + p2
