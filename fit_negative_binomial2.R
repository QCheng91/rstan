library(rstan)
library(rlang)
library(tidyverse)
library(doParallel)
library(doMC)

doMC::registerDoMC(cores = 14)
rstan_options(auto_write = TRUE)

# options(mc.cores = parallel::detectCores())
load("Data_K2.rda") # gated_data

tukey_lim = function(x, k = 2){
  q1_q3 = quantile(x, c(0.75, 0.9))
  upper_lim = q1_q3[2] + k*(q1_q3[2]-q1_q3[1])
  upper_lim
}

# prot_data %>% 
#   ggplot(aes(x = prot_count))  + geom_bar(aes(y = ..count../sum(..count..))) + facet_wrap(~Day) #+ xlim(-1, 400) 
# # + geom_histogram(stat = "density", bins = 100) 

fit_neg_binom = function(x, y,  day_index, d, cmp_model){
  num_chains = 4
  data_for_stan = list(N = length(x), x = as.integer(x), num_proteins = ncol(y), y = y, day_index = data_day, d=day_number)
  # run_fit = stan("model_poisson_gamma.stan", iter = 1000, data = data_for_stan, chains = num_chains, thin = 2, cores = 1)
  run_fit = sampling(cmp_model, iter = 1000, data = data_for_stan, chains = num_chains, thin = 2, cores = 4)
  tbl_summary = summary(run_fit)$summary %>% as_tibble() 
  if(any(tbl_summary$Rhat < 0.9 | tbl_summary$Rhat > 1.1)) print("convergence_problem")
  #tbl_summary %>% select(`2.5%`, `50%`, `97.5%`) %>% mutate(type = c("mu", "phi", "a", "lp__"))
  run_fit
}

cmp_model = stan_model(file = "Research/RBC/Data Analysis/Data_October/model_poisson_gamma2.stan")


all_prot = names(gated_data)[-1*(1:2)]
pdf("Neg_Bin_fit.pdf", 14, 9)
## Loop through proteins

for(protein in all_prot){
  
  protein = "TAL1"
  
  prot_data = gated_data %>% rename(prot_count = !!!sym(protein)) %>% mutate(prot_count = ceiling(prot_count)) %>%  
    group_by(Day) %>% sample_n(5000) %>% filter(prot_count<tukey_lim(prot_count, 3)) 
  
  day_index = prot_data %>% group_indices(.) 
  
  prot_data_0 = prot_data %>% ungroup() %>% select(-Event_Length, -sum_prot, -mean_Ir, -mean_Pt, -Ir191, -Ir193, -Pt194, -Pt195) %>% bind_cols(Day_index = data.frame(day_index), .) 
  
  prot_data1 = prot_data_0 %>% select(Day, day_index, prot_count)
  
  prot_data2 = prot_data_0 %>% select(Day, day_index, RUNX1, IKZF1, ATRX, CD90, GATA2, CD49f, KLF1, KAT3B) %>% select(-Day) %>% ceiling(.) %>% bind_cols(Day = prot_data$Day, .) %>% 
    group_by(Day) %>% ungroup()
  
  #prot_data2 = prot_data %>% select(-starts_with(protein), -Cell_id, -prot_count) %>% select(-Day) %>% ceiling(.) %>% bind_cols(Day = prot_data$Day, .) %>% 
  #  group_by(Day) %>% ungroup()
  
  # prot_data1 = gated_data %>% select(Day, Cell_id, !!!sym(protein)) %>% rename(prot_count = !!!sym(protein)) %>% mutate(prot_count = ceiling(prot_count)) %>%  
  #group_by(Day) %>% sample_n(5000) %>% filter(prot_count<tukey_lim(prot_count, 3)) %>% ungroup()
  
  #prot_data2 = gated_data %>% select(-starts_with(protein)) %>% select(Cell_id = pr) select(-Day) %>% ceiling(.) %>% bind_cols(Day = gated_data$Day, .) %>% 
  #  group_by(Day) %>% sample_n(5000) %>% ungroup()
  
  ## Parallel loop through all days
  all_fits = foreach(i = unique(prot_data$Day))%dopar%{
    
    #i = "2"
    data_x = prot_data1 %>% select(prot_count) %>% as_vector() #%>% filter(Day == i)
    data_y = prot_data2 %>% select(-Day, -day_index) %>% as.matrix() #%>% filter(Day == i)
    data_day = prot_data2 %>% select(day_index) %>% as_vector()
    day_number = length(table(data_day))
    
    length(data_x)
    ncol(data_y)
    colnames(data_y)
    
    fit_out = fit_neg_binom(data_x, data_y, data_day, day_number, cmp_model) #%>% mutate(Day = i)
    
    print(fit_out, digits=8)
    plot(fit_out, pars = "mu")
    
    tbl_pred = tibble(x = 0:max(data_x), y = dnbinom(x, mu = fit_out$`50%`[1], size = fit_out$`50%`[2]))
    pp_out = tibble(x = data_x) %>% 
      ggplot(aes(x)) + geom_bar(aes(y = ..count../sum(..count..))) + geom_point(data = tbl_pred, aes(x = x, y =y), col = "red", alpha= 0.9) +
      ggtitle(paste("Day", i)) + ylab("Proportion of counts") #+  xlim(-1, 50)
    list(pp = pp_out, fit_out = fit_out)
  }
  
  pp_mu = bind_rows(map(all_fits, 2)) %>% filter(type =="mu") %>% ggplot(aes(x = Day, y = `50%`)) + geom_point(size = 2) + 
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%` ), width = 0.5) + facet_wrap(~type, scale = "free_y")
  
  pp_phi = bind_rows(map(all_fits, 2)) %>% filter(type =="phi") %>% ggplot(aes(x = Day, y = `50%`)) + geom_point(size = 2) + 
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%` ), width = 0.5) + facet_wrap(~type, scale = "free_y")
  
  plots_list = map(all_fits, 1)
  
  gridExtra::grid.arrange(grobs = c(plots_list, list(pp_mu), list(pp_phi)),  top=grid::textGrob(protein, gp=grid::gpar(fontsize=15,font=8)))
  
  print(paste0("Done", protein))
  
}

dev.off()



