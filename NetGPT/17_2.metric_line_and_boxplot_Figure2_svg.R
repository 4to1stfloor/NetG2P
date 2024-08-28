library(readxl)

metric = read_xlsx("~/nas/04.Results/ml_metric/Total_for_figure_only_mcc.xlsx")

metric_tes = metric %>%
  filter(cancer_type %in% c("Average", "SD")) %>% 
  t()
colnames(metric_tes) = c("MCC", "SD")
metric_tes = metric_tes[-1,]
metric_tes = metric_tes %>% as.data.frame() %>%
  mutate(across(everything(), as.numeric))

metric_tes_filt = metric_tes %>% 
  mutate(data_type = rownames(.))

rownames(metric_tes_filt) = NULL
metric_tes_filt$data_type = factor(metric_tes_filt$data_type , levels = metric_tes_filt$data_type)
original = metric_tes_filt %>%
  ggplot(aes(x = data_type, y = MCC, group = 1)) +
  geom_point(color = "grey", size = 5)+
  geom_line(color = "grey") +
  geom_errorbar(aes(ymin = MCC - SD , ymax = MCC + SD),
                width=.2,
                position=position_dodge(.9) , 
                color = "grey") +
  theme_classic() + 
  ylim(c(0,0.8))
setwd("~/nas/04.Results/ml_metric/")

ggsave(file= "original_metric.svg", plot=original,width = 8,height = 4)
####
metric_box <- metric %>%
  filter(!cancer_type %in% c("Average", "SD")) %>% t()
colnames(metric_box) = metric_box[1,]
metric_box = metric_box[-1,]
metric_box = metric_box %>% as.data.frame() %>%
  mutate(across(everything(), as.numeric))
metric_box = metric_box %>%
  mutate(data_type = rownames(.))
rownames(metric_box) = NULL

library(reshape2)
library(ggpubr)
library(rstatix)
metric_box_re = melt(metric_box , id.vals = c("data_type"), measure.vars = colnames(metric_box)[-ncol(metric_box)])
colnames(metric_box_re) = c("data_type", "cancer_type", "mcc")
metric_box_re$data_type = factor(metric_box_re$data_type , levels = unique(metric_box_re$data_type ))

mean <- metric_box_re %>%  
  group_by(data_type) %>%  
  summarize(average = mean(mcc)) %>% 
  ungroup() 
mean$data_type = factor(mean$data_type , levels = unique(metric_box_re$data_type ))

stat_test <- metric_box_re %>%
  select(-cancer_type) %>%
  t_test(mcc ~ data_type) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

stat_test <- stat_test %>%
  filter(group2 == "multiomics_combine") %>%
  mutate(data_type = group1) %>% 
  select(data_type, everything()) %>%
  add_xy_position(fun = "mean_sd", x = "data_type", dodge = 0.8) 

stat_test = stat_test %>% arrange(data_type)

stat_test = stat_test %>% 
  mutate(x = row_number()) %>% 
  ungroup() %>% 
  mutate(xmin = x - 0.2,
         xmax = x + 0.2) 

stat_test_add_signif = stat_test %>% mutate(p_signif = case_when(p > 0.05 ~ "ns",
                                                                 p <= 0.05 & p > 0.01 ~ "*",
                                                                 p <= 0.01 & p > 0.001 ~ "**" ,
                                                                 p <= 0.001 ~ "***", 
                                                                 .default = NA))
library(readxl)

saveRDS(stat_test_add_signif, "~/nas/04.Results/ml_metric/stat_for_figure2.rds")

stat_test_add_signif = stat_test_add_signif %>% arrange(data_type)
stat_test_add_signif$y.position = stat_test_add_signif$y.position / 2.22

stat_test_add_signif$data_type = factor(stat_test_add_signif$data_type , levels = unique(metric_box_re$data_type ))
stat_test_add_signif$x = c(1.5:8.5)
stat_test_add_signif$xmin = stat_test_add_signif$x -0.5
stat_test_add_signif$xmax = stat_test_add_signif$x +0.5

box_with_line = metric_box_re %>%
  ggplot(aes(x = data_type, y = mcc)) +
  geom_boxplot() +
  geom_line(mean,mapping = aes(x = data_type, y = average , group = 1))+
  stat_pvalue_manual(
    stat_test_add_signif,  
    label = "p_signif", 
    tip.length = 0.02,
    bracket.nudge.y = 0.2
  )  + 
  ylim(c(0,0.8))+
  theme_classic()

setwd("~/nas/04.Results/ml_metric/")

ggsave(file= "box_with_line_metric.svg", plot=box_with_line,width = 8,height = 4)
