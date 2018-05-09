library(rstan)
library(rlang)
library(tidyverse)
library(doParallel)
library(stringr)
doMC::registerDoMC(cores = 14)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

tukey_lim = function(x, k = 2){
  q1_q3 = quantile(x, c(0.75, 0.9))
  upper_lim = q1_q3[2] + k*(q1_q3[2]-q1_q3[1])
  upper_lim
}

qplot(seq(0, 1, length.out = 100), dbeta(seq(0, 1, length.out = 100), 1, 3))

load("Desktop/Apr2018/May_K1.rda") # gated_data

cmp_model = stan_model(file = "Research/RBC/Data Analysis/mixture_models_stan/mix_nb_model_V3.stan")

fit_neg_binom = function(xx, day_number, cmp_model){
  
  num_chains = 4
  
  data_for_stan = list(N = nrow(xx), x = xx$prot_count, day_index=xx$day_index, d=day_number) 
  
  run_fit = sampling(cmp_model, iter = 500, warmup = 250, data = data_for_stan, chains = num_chains, thin = 2, cores = 1)
  
  run_fit
}

all_prot = names(gated_data1)[-1*(1:2)] 

all_prot = all_prot[-c(21, 22, 27, 28, 32, 33)]

gated_data1 = gated_data %>% gather(prot_name, Intensity, ATRX:mean_Ir) %>% filter(Day == "MNCs1"| Day == "MNCs2"| Day == "MNCs3") %>%  
  spread(prot_name, Intensity)

gated_data1 = gated_data1 %>% mutate(Day="MNCs")

pdf("Desktop/Apr2018/K1_neg_bino/May_K1_MNCs.pdf", 14, 9)
## Loop through proteins
foreach(protein = all_prot)%dopar%{
#for(protein in all_prot){ 
  protein = "CD235ab"
   pp_out = list()
prot_data = gated_data1 %>% select(Day, !!!sym(protein)) %>% rename(prot_count = !!!sym(protein)) %>% mutate(prot_count = ceiling(prot_count)) %>%
  group_by(Day) %>% sample_n(2000) #%>% filter(prot_count<tukey_lim(prot_count, 6)) #%>% filter(Day== "2" | Day== "8" | Day == "12" | Day == "14" | Day == "16" | Day == "20")

day_index = prot_data %>% group_indices(.)

data_x1 = prot_data %>% ungroup() %>% select(Day, prot_count) %>% bind_cols(day_index = data.frame(day_index), .) %>% as.tibble()

day_number = length(table(data_x1$day_index))

median_n = prot_data %>% group_by(Day) %>% summarise(median = median(prot_count, na.rm = TRUE))

median_n_min = min(median_n$median)

direct1 = "Desktop/Apr2018/K1_neg_bino/"
direct2 = "_K1.csv"
direct3 = "_K1.rda"
direct4 = "_K1.pdf"
direct0 = paste(direct1, protein, direct2, sep="")
direct = paste(direct1, protein, direct3, sep="")
directp = paste(direct1, protein, direct4, sep="") 

pdf(file = directp, 14, 9)

run_fit = fit_neg_binom(data_x1, day_number, cmp_model) #, median_n_min)

traceplot(run_fit, inc_warmup = TRUE, pars = "mu_n")

save(run_fit, file = direct)

#load(file ="Desktop/Apr2018/K1_neg_bino/CD235ab_K1.rda")

tbl_summary = summary(run_fit)$summary %>% as.data.frame() %>% rownames_to_column() %>%  as_tibble()

if(any(tbl_summary$Rhat < 0.9 | tbl_summary$Rhat > 1.1)) print("convergence_problem")

#write.csv(tbl_summary, file = direct0)

fit_out = tbl_summary %>% select(rowname,`2.5%`, `50%`, `97.5%`) 

day_vals =c("MNCs")# "2", "4", "6", "8", "10", "12", "14", "16", "18", "20") 

# c("2", "8", "12", "14", "16", "20") c("0", "2", "4", "6", "8", "10", "11", "12", "14", "16", "18", "20", "MNCs", "Jurkats")

for(d in day_vals){
  
  i = d
  
  kk = which(day_vals == d)
  print(kk)
  data_x = data_x1 %>% filter(Day == i) %>% select(prot_count)
  
  # tbl_summary1 = tbl_summary[1:20,] %>% mutate(mu = str_replace_all(rowname, "mu_h\\[|\\]", "")) %>% 
  # separate(mu, c("obstime", "state"), ",") %>% mutate(obstime = as.numeric(obstime), state = as.numeric(state))
  # 
  # tbl_summary2 = tbl_summary[22:31,] %>% mutate(phi = str_replace_all(rowname, "phi\\[|\\]", "")) %>% 
  #   separate(phi, c("obstime"), ",") %>% mutate(obstime = as.numeric(obstime))
  
  tbl_pred = tibble(x = c(0, 1:max(data_x)), 
                    y1 = c((1-fit_out$`50%`[23+kk-1])*(fit_out$`50%`[33+kk-1] + (1-fit_out$`50%`[33+kk-1])*dnbinom(0, mu = fit_out$`50%`[1], size = fit_out$`50%`[12])), 
                           (1-fit_out$`50%`[23+kk-1])*((1-fit_out$`50%`[33+kk-1])*dnbinom(x[-1], mu = fit_out$`50%`[1], size = fit_out$`50%`[12]))), 
                    y2 = (fit_out$`50%`[23+kk-1])*dnbinom(x, mu = fit_out$`50%`[2+kk-1], size = fit_out$`50%`[13+kk-1]))
  
 
  #fit_out$`50%`[11+kk-1]*
    
  #mm = length(which(data_x$prot_count == 0)) 
 # length_data_x = length(data_x$prot_count)
  #tbl_pred_0 = tibble(y3 = rbernoulli(mm, p= fit_out$`50%`[15]), length_data_x = length_data_x)
  
 # tbl_pred_0$y3[FALSE] <- 0
  # tbl_pred = tibble(x = 0:max(data_x), 
  #                   y1 = fit_out$`50%`[7+3*kk-3]*dnbinom(x, mu = fit_out$`50%`[1], size = fit_out$`50%`[4]), 
  #                   y2 = (fit_out$`50%`[8+3*kk-3])*dnbinom(x, mu = fit_out$`50%`[2], size = fit_out$`50%`[5]),
  #                   y3 = (fit_out$`50%`[9+3*kk-3])*dnbinom(x, mu = fit_out$`50%`[3], size = fit_out$`50%`[6]))
 
  #mm = length(which(data_x$prot_count == 0))/length(data_x$prot_count)
  pp_out[[kk]] =  data_x %>% 
    ggplot(aes(x=prot_count)) + geom_bar(aes(y = ..count../sum(..count..))) + geom_point(data = tbl_pred, aes(x = x, y =y1), col = "red", alpha= 0.9) + 
    geom_point(data = tbl_pred, aes(x = x, y =y2), col = "blue", alpha= 0.9)  +
    ggtitle(paste("Day", i)) + ylab("Proportion of counts") +  xlab("Counts")
  #+ geom_bar(data = tbl_pred_0, aes(x = y3, y=..count..), col = "green", alpha= 0.9)
  print(paste0("Done", d))
}

pp_out
# tbl_summary$`50%`
# 
# theta1 = c(9.436588e-01, 9.567756e-01, 9.796791e-01, 9.819549e-01, 9.949175e-01, 9.829168e-01, 9.430644e-01, 6.862853e-01,
#            1.011087e-01, 2.075781e-02, 3.842676e-02, 1.578665e-02, 9.914345e-01, 9.928245e-01)
# 
# theta2 = c(5.485633e-02, 4.161768e-02, 1.856077e-02, 1.626686e-02, 4.225237e-03, 1.611090e-02, 5.354240e-02, 1.792807e-01,
#            5.596186e-02, 8.544897e-02, 5.343530e-02, 4.235841e-03, 8.285950e-03, 6.099456e-03)
# 
# theta3 = c(1.155183e-03, 1.079145e-03, 1.127803e-03, 4.819867e-04, 4.752151e-04, 4.759796e-04, 3.189323e-03, 1.796458e-01,
#            8.638063e-01, 9.005260e-01, 9.177685e-01, 9.796872e-01, 4.619559e-04, 4.814225e-04)
# 
# theta = data.frame(theta1 = theta1, theta2 = theta2, theta3 = theta3) %>% mutate(Day = c(1:14)) %>% gather(theta, value, -Day)
# 
# theta %>%  ggplot(aes(x=Day, y=value)) + geom_point() + facet_wrap(~theta)

# #num_chains = 4
# 
# #data_for_stan = list(N = nrow(data_x1), x = data_x1$prot_count, day_index=data_x1$day_index, d=day_number)
# 
# run_fit = stan("Research/RBC/Data Analysis/mixture_models_stan/mix_gamma_model.stan", data = data_for_stan, chains = num_chains, iter = 1500, warmup = 500)
# 
# #traceplot(run_fit, pars = "theta")
# 
# save(x, y, file = "xy.RData")
# 
# tbl_summary = summary(run_fit)$summary %>% as.data.frame() %>% rownames_to_column() %>%  as_tibble()
# 
# if(any(tbl_summary$Rhat < 0.9 | tbl_summary$Rhat > 1.1)) print("convergence_problem")
# 
# direct1 = "Desktop/Nov9/Mix_model/"
# direct2 = "_K1.csv"
# direct = paste(direct1, protein, direct2, sep="")
# 
# write.csv(tbl_summary, file = direct)
# 
# fit_out = tbl_summary %>% select(rowname,`2.5%`, `50%`, `97.5%`) 
# 
# day_vals = c("0", "2", "4", "6", "8", "10", "11", "12", "14", "16", "18", "20", "MNCs", "Jurkats")
# 
# 
# 
# for(d in day_vals) {
#   
#   i = d
#   
#   kk = which(day_vals == d)
# 
#   data_x = data_x1 %>% filter(Day == i) %>% select(prot_count)
# 
#  tbl_pred = tibble(x = 0:max(data_x), y1 =fit_out$`50%`[30+kk]*dnbinom(x, mu = fit_out$`50%`[1], size = fit_out$`50%`[16]), y2 = (1-fit_out$`50%`[30+kk])*dnbinom(x, mu = fit_out$`50%`[1+kk], size = fit_out$`50%`[16+kk]))
# 
#  pp_out[[kk]] =  data_x %>% 
#w   ggplot(aes(x=prot_count)) + geom_bar(aes(y = ..count../sum(..count..))) + geom_point(data = tbl_pred, aes(x = x, y =y1), col = "red", alpha= 0.9) + 
#   geom_point(data = tbl_pred, aes(x = x, y =y2), col = "blue", alpha= 0.9) +
#   ggtitle(paste("Day", i)) + ylab("Proportion of counts") #+  xlim(-1, 50)
#  
#  print(paste0("Done", d))
# }

newpic = gridExtra::grid.arrange(grobs = pp_out[1:10], ncol = 4, top=grid::textGrob(protein, gp=grid::gpar(fontsize=15,font=8)))
print(paste0("Done", protein))
dev.off()
}
dev.off()
#########################################################################
load("Desktop/Apr2018/May_K1.rda") # gated_data

all_prot = names(gated_data)[-1*(1:2)] 

all_prot = all_prot[-c(15, 21, 22)]

pdf("Desktop/Apr2018/May_K1.pdf", 14, 9)

all_prob <- data.frame()

for(protein in all_prot){
  
  #protein = "CD34"
  
  prot_data = gated_data %>% select(Day, !!!sym(protein)) %>% rename(prot_count = !!!sym(protein)) %>% mutate(prot_count = ceiling(prot_count)) %>%
    group_by(Day) %>% sample_n(1499) %>% filter(prot_count<tukey_lim(prot_count, 6)) #%>% filter(Day== "2" | Day== "8" | Day == "12" | Day == "14" | Day == "16" | Day == "20")
  
  day_index = prot_data %>% group_indices(.)
  
  data_x1 = prot_data %>% ungroup() %>% select(Day, prot_count) %>% bind_cols(day_index = data.frame(day_index), .) %>% as.tibble()
  
  direct1 = "Desktop/Apr2018/K1_neg_bino/"
  direct2 = "_K1.rda"
  direct0 = paste(direct1, protein, direct2, sep="")
  
  load(file =direct0)
  
  tbl_summary = summary(run_fit)$summary %>% as.data.frame() %>% rownames_to_column() %>%  as_tibble()
  
  fit_out = tbl_summary %>% select(rowname,`2.5%`, `50%`, `97.5%`) 
  
  day_vals =c("2", "4", "6", "8", "10", "12", "14", "16", "18", "20") 
  
  prob_vals <- data.frame()
  
  for(d in day_vals){
    
    i = d
    kk = which(day_vals == d)
    print(kk)
    
    data_x = data_x1 %>% filter(Day == i) %>% select(prot_count)
    data_max = max(data_x)#data_x1 %>% select(prot_count) %>% max()
    
    theta_vals = fit_out$`50%`[23+kk-1]

    tbl_pred = tibble(Day = d,
                      
                      prot_name = protein,
                      
                      x = c(0, 1:data_max), 
                      
                      prob_n = c((fit_out$`50%`[33+kk-1] + (1-fit_out$`50%`[33+kk-1])*dnbinom(0, mu = fit_out$`50%`[1], size = fit_out$`50%`[12])), 
                              ((1-fit_out$`50%`[33+kk-1])*dnbinom(x[-1], mu = fit_out$`50%`[1], size = fit_out$`50%`[12]))), 
                      
                      prob_s = dnbinom(x, mu = fit_out$`50%`[2+kk-1], size = fit_out$`50%`[13+kk-1]),
    
                      prob_th = theta_vals*prob_s/(theta_vals*prob_s + (1-theta_vals)*prob_n))
    
    tbl_pred_prob = tbl_pred %>% select(prot_name, Day, x, prob_th)
    
    prob_vals <- rbind(prob_vals, tbl_pred_prob)
}
  
  all_prob = rbind(all_prob, prob_vals)
  
  
# pp = prob_vals %>% mutate(Day = factor(Day, levels = day_vals)) %>% ggplot(aes(x=x, y = prob_th)) +
#     geom_point()+
#     facet_wrap(~Day, scales = "free_x") +
#     ggtitle(protein) + ylab("Probability of ON") + xlab("Counts")
#     #filter(prot_name == "CD235ab") %>% 
# plot(pp)  
}

dev.off()
    
    th_ = all_prob %>% filter(prob_th>=0.45 & prob_th <= 0.6)
    
    th_comp = data.frame(prot_name = c(rep("ATRX", )),
                         Day = c(),
                         x = c())



##########################################################################
pdf("Desktop/Apr2018/K1_neg_bino/fit_mixed_K1.pdf", 14, 9)

## Loop through proteins
for(protein in all_prot){
  
  prot_data = gated_data %>% select(Day, !!!sym(protein)) %>% rename(prot_count = !!!sym(protein)) %>% mutate(prot_count = ceiling(prot_count)) %>%  
    group_by(Day) %>% sample_n(1000) %>% filter(prot_count<tukey_lim(prot_count, 3)) 
  
  day_index = prot_data %>% group_indices(.) 
  
  data_x1 = prot_data %>% ungroup() %>% select(Day, prot_count) %>% bind_cols(day_index = data.frame(day_index), .) %>% as.tibble()
  
  day_number = length(table(data_x1$day_index))
  
  direct1 = "Desktop/Apr2018/K1_neg_bino/"
  direct2 = "_K1.csv"
  direct3 = "_K1.RData"
  direct4 = "_K1.pdf"
  direct0 = paste(direct1, protein, direct2, sep="")
  direct = paste(direct1, protein, direct3, sep="")
  directp = paste(direct1, protein, direct4, sep="") 
  
  #pdf(directp, 14, 9)
  
  #load(direct) # gated_data
 load("Desktop/Dec06/2+1/K1/CD38_K1.rda")
 
  tbl_summary = summary(run_fit)$summary %>% as.data.frame() %>% rownames_to_column() %>%  as_tibble()
  
  if(any(tbl_summary$Rhat < 0.9 | tbl_summary$Rhat > 1.1)) print("convergence_problem")
  
  fit_out = tbl_summary %>% select(rowname,`2.5%`, `50%`, `97.5%`) 
  
  day_vals = c("2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
  
  pp_out = list()
  
  for(d in day_vals){
    
    i = d
    
    kk = which(day_vals == d)
    
    data_x = data_x1 %>% filter(Day == i) %>% select(prot_count)
    
    tbl_pred = tibble(x = 0:max(data_x), y1 =fit_out$`50%`[30+kk]*dnbinom(x, mu = fit_out$`50%`[1], size = fit_out$`50%`[16]), y2 = (1-fit_out$`50%`[30+kk])*dnbinom(x, mu = fit_out$`50%`[1+kk], size = fit_out$`50%`[16+kk]))
    
    pp_out[[kk]] =  data_x %>% 
      ggplot(aes(x=prot_count)) + geom_bar(aes(y = ..count../sum(..count..))) + geom_point(data = tbl_pred, aes(x = x, y =y1), col = "red", alpha= 0.9) + 
      geom_point(data = tbl_pred, aes(x = x, y =y2), col = "blue", alpha= 0.9) +
      ggtitle(paste("Day", i)) + ylab("Proportion of counts") #+  xlim(-1, 50)
    
    print(paste0("Done", d))
  }

  newpic = gridExtra::grid.arrange(grobs = pp_out[1:10], ncol = 4, top=grid::textGrob(protein, gp=grid::gpar(fontsize=15,font=8)))

  print(paste0("Done", protein))
}

dev.off()
# 
# #num_chains = 4
# 
# #data_for_stan = list(N = nrow(data_x1), x = data_x1$prot_count, day_index=data_x1$day_index, d=day_number)
# 
# run_fit = stan("Research/RBC/Data Analysis/mixture_models_stan/mix_gamma_model.stan", data = data_for_stan, chains = num_chains, iter = 1500, warmup = 500)
# 
# #traceplot(run_fit, pars = "theta")
# 
# save(x, y, file = "xy.RData")
# 
# tbl_summary = summary(run_fit)$summary %>% as.data.frame() %>% rownames_to_column() %>%  as_tibble()
# 
# if(any(tbl_summary$Rhat < 0.9 | tbl_summary$Rhat > 1.1)) print("convergence_problem")
# 
# direct1 = "Desktop/Nov9/Mix_model/"
# direct2 = "_K1.csv"
# direct = paste(direct1, protein, direct2, sep="")
# 
# write.csv(tbl_summary, file = direct)
# 
# fit_out = tbl_summary %>% select(rowname,`2.5%`, `50%`, `97.5%`) 
# 
# day_vals = c("0", "2", "4", "6", "8", "10", "11", "12", "14", "16", "18", "20", "MNCs", "Jurkats")
# 
# 
# 
# for(d in day_vals) {
#   
#   i = d
#   
#   kk = which(day_vals == d)
# 
#   data_x = data_x1 %>% filter(Day == i) %>% select(prot_count)
# 
#  tbl_pred = tibble(x = 0:max(data_x), y1 =fit_out$`50%`[30+kk]*dnbinom(x, mu = fit_out$`50%`[1], size = fit_out$`50%`[16]), y2 = (1-fit_out$`50%`[30+kk])*dnbinom(x, mu = fit_out$`50%`[1+kk], size = fit_out$`50%`[16+kk]))
# 
#  pp_out[[kk]] =  data_x %>% 
#   ggplot(aes(x=prot_count)) + geom_bar(aes(y = ..count../sum(..count..))) + geom_point(data = tbl_pred, aes(x = x, y =y1), col = "red", alpha= 0.9) + 
#   geom_point(data = tbl_pred, aes(x = x, y =y2), col = "blue", alpha= 0.9) +
#   ggtitle(paste("Day", i)) + ylab("Proportion of counts") #+  xlim(-1, 50)
#  
#  print(paste0("Done", d))
# }


# protein = "GATA1"
# 
# prot_data = gated_data_k1 %>% select(Day, !!!sym(protein)) %>% rename(prot_count = !!!sym(protein)) %>% mutate(prot_count = ceiling(prot_count)) %>%  
#   group_by(Day) %>% sample_n(1000) %>% filter(prot_count<tukey_lim(prot_count, 3)) 
# 
# day_index = prot_data %>% group_indices(.) 
# 
# data_x1 = prot_data %>% ungroup() %>% select(Day, prot_count) %>% bind_cols(day_index = data.frame(day_index), .) %>% as.tibble()
# 
# data_x1 %>% 
#   ggplot(aes(x = prot_count)) + geom_bar(aes(y = ..count../sum(..count..))) + facet_wrap(~Day, scales = "free")
# 
# day_number = length(table(data_x1$day_index))
# 
# num_chains = 4
# 
# data_for_stan    = list(N = nrow(data_x1), x = data_x1$prot_count, day_index=data_x1$day_index, d=day_number)
# 
# is.atomic(data_x1$prot_count)
# # initialise_means = list(m1 = c(1, 2, 4))
# # init_par         = replicate(num_chains, initialise_means, simplify=FALSE)
# run_fit = stan("Research/RBC/Data Analysis/mixture_models_stan/mix_gamma_model.stan", data = data_for_stan, chains = num_chains, iter = 2000, warmup = 1000)
# 
# traceplot(run_fit, pars = "theta")
# 
# tbl_summary = summary(run_fit)$summary %>% as.data.frame() %>% rownames_to_column() %>%  as_tibble()
# 
# fit_out = tbl_summary %>% select(rowname,`2.5%`, `50%`, `97.5%`) #%>% mutate(type = c("mu1", "mu2", "phi1", "phi2", "theta", "lp__"))
# 
# aa = rstan:::extract(run_fit)
# 
# mu_pre = data.frame(mu = aa$mu)
# phi_pre = data.frame(phi = aa$phi)
# theta_pre = data.frame(theta = aa$theta)
# 
# mcmc = bind_cols(mu_pre, phi_pre, theta_pre)
# 
# mcmc %>% gather(var, val) %>% ggplot(aes(x = val))+geom_histogram()+facet_wrap(~var, scales = "free")
# 
# day_vals = c("0", "2", "4", "6", "8", "10", "11", "12", "14", "16", "18", "20", "MNCs", "Jurkats")
# 
# pp_out = list()
# 
# for(d in day_vals) {
#   
#   i = d
#   kk = which(day_vals == d)
#   
#   data_x = data_x1 %>% filter(Day == i) %>% select(prot_count)
#   
#   # mu_n = fit_out$`50%`[1]
#   # phi_n = fit_out$`50%`[16]
#   # mu = fit_out$`50%`[2:15]
#   # phi = fit_out$`50%`[17:30]
#   
#   tbl_pred = tibble(x = 0:max(data_x), y1 =fit_out$`50%`[30+kk]*dnbinom(x, mu = fit_out$`50%`[1], size = fit_out$`50%`[16]), y2 = (1-fit_out$`50%`[30+kk])*dnbinom(x, mu = fit_out$`50%`[1+kk], size = fit_out$`50%`[16+kk]))
#   
#   pp_out[[kk]] =  data_x %>% 
#     ggplot(aes(x=prot_count)) + geom_bar(aes(y = ..count../sum(..count..))) + geom_point(data = tbl_pred, aes(x = x, y =y1), col = "red", alpha= 0.9) + 
#     geom_point(data = tbl_pred, aes(x = x, y =y2), col = "blue", alpha= 0.9) +
#     ggtitle(paste("Day", i)) + ylab("Proportion of counts") #+  xlim(-1, 50)
#   
#   print(paste0("Done", d))
# }
# 
# length(pp_out)
# 
# newpic = gridExtra::grid.arrange(grobs = pp_out[1:14], ncol =4)
# 
# dev.off()

#################################################
# protein = "GATA1"
# 
# prot_data = gated_data_k1 %>% select(Day, !!!sym(protein)) %>% rename(prot_count = !!!sym(protein)) %>% mutate(prot_count = ceiling(prot_count)) %>%  
#   group_by(Day) %>% sample_n(5000) %>% filter(prot_count<tukey_lim(prot_count, 3)) %>% ungroup() 
# 
# i = "Jurkats"
# 
# data_x1 = prot_data %>% filter(Day == i) %>% select(prot_count) %>% as.tibble()
# 
# # data_sim = bind_rows(tibble(x1 = rgamma(400, 1, 0.5))   %>% mutate(cluster  = "Cluster1"),
# #                      tibble(x1 = rgamma(100, 10, 1)) %>% mutate(cluster  = "Cluster2"))
# 
# data_x1 %>% 
#   ggplot(aes(x = prot_count)) + geom_histogram(binwidth = 1.0)
# 
# num_chains = 4
# data_for_stan    = list(N = nrow(data_x1), x = data_x1$prot_count)
# 
# is.atomic(data_x1$prot_count)
# # initialise_means = list(m1 = c(1, 2, 4))
# # init_par         = replicate(num_chains, initialise_means, simplify=FALSE)
# run_fit = stan("Research/RBC/Data Analysis/mixture_models_stan/mix_gamma_model.stan", data = data_for_stan, chains = num_chains, iter = 4000, warmup = 2000)
# 
# traceplot(run_fit)
# 
# tbl_summary = summary(run_fit)$summary %>%  as_tibble()
# 
# fit_out = tbl_summary %>% select(`2.5%`, `50%`, `97.5%`) %>% mutate(type = c("mu1", "mu2", "phi1", "phi2", "theta", "lp__"))
# 
# aa = rstan:::extract(run_fit)
# 
# mu_pre = data.frame(mu = aa$mu)
# phi_pre = data.frame(phi = aa$phi)
# theta_pre = data.frame(theta = aa$theta)
# 
# mcmc = bind_cols(mu_pre, phi_pre, theta_pre)
# 
# mcmc %>% gather(var, val) %>% ggplot(aes(x = val))+geom_histogram()+facet_wrap(~var, scales = "free")
# 
# tbl_pred = tibble(x = 0:max(data_x1), y1 = dnbinom(x, mu = fit_out$`50%`[1], size = fit_out$`50%`[3]), y2 = dnbinom(x, mu = fit_out$`50%`[2], size = fit_out$`50%`[4]))
# 
# 
# pp_out =  data_x1 %>% 
#   ggplot(aes(x=prot_count)) + geom_bar(aes(y = ..count../sum(..count..))) + geom_point(data = tbl_pred, aes(x = x, y =y1), col = "red", alpha= 0.9) + 
#   geom_point(data = tbl_pred, aes(x = x, y =y2), col = "blue", alpha= 0.9) +
#   ggtitle(paste("Day", i)) + ylab("Proportion of counts") #+  xlim(-1, 50)


######################20171212
## Loop through proteins
# foreach(protein = all_prot)%dopar%{
#   #for(protein in all_prot){ 
#   protein = "CD235ab"
#   pp_out = list()
#   prot_data = gated_data_k1 %>% select(Day, !!!sym(protein)) %>% rename(prot_count = !!!sym(protein)) %>% mutate(prot_count = ceiling(prot_count)) %>%
#     group_by(Day) %>% sample_n(1000) %>% filter(prot_count<tukey_lim(prot_count, 3)) %>% filter(Day== "0" | Day== "10" | Day == "16" | Day == "20")
#   
#   day_index = prot_data %>% group_indices(.)
#   
#   data_x1 = prot_data %>% ungroup() %>% select(Day, prot_count) %>% bind_cols(day_index = data.frame(day_index), .) %>% as.tibble()
#   
#   day_number = length(table(data_x1$day_index))
#   
#   direct1 = "Desktop/Dec06/"
#   direct2 = "_K1.csv"
#   direct3 = "_K1.rda"
#   direct4 = "_K1.pdf"
#   direct0 = paste(direct1, protein, direct2, sep="")
#   direct = paste(direct1, protein, direct3, sep="")
#   directp = paste(direct1, protein, direct4, sep="") 
#   
#   run_fit = fit_neg_binom(data_x1, day_number, cmp_model)
#   
#   traceplot(run_fit, par="theta")
#   
#   save(run_fit, file = direct)
#   
#   load(direct)
#   
#   tbl_summary = summary(run_fit)$summary %>% as.data.frame() %>% rownames_to_column() %>%  as_tibble()
#   
#   if(any(tbl_summary$Rhat < 0.9 | tbl_summary$Rhat > 1.1)) print("convergence_problem")
#   
#   write.csv(tbl_summary, file = direct0)
#   
#   fit_out = tbl_summary %>% select(rowname,`2.5%`, `50%`, `97.5%`) 
#   
#   day_vals = c("0",  "10", "16", "20") #c("0", "2", "4", "6", "8", "10", "11", "12", "14", "16", "18", "20", "MNCs", "Jurkats")
#   
#   for(d in day_vals){
#     
#     i = d
#     
#     kk = which(day_vals == d)
#     print(kk)
#     data_x = data_x1 %>% filter(Day == i) %>% select(prot_count)
#     
#     tbl_pred = tibble(x = 0:max(data_x), 
#                       y1 = fit_out$`50%`[7+3*kk-3]*dnbinom(x, mu = fit_out$`50%`[1], size = fit_out$`50%`[4]), 
#                       y2 = (fit_out$`50%`[8+3*kk-3])*dnbinom(x, mu = fit_out$`50%`[2], size = fit_out$`50%`[5]),
#                       y3 = (fit_out$`50%`[9+3*kk-3])*dnbinom(x, mu = fit_out$`50%`[3], size = fit_out$`50%`[6]))
#     
#     pp_out[[kk]] =  data_x %>% 
#       ggplot(aes(x=prot_count)) + geom_bar(aes(y = ..count../sum(..count..))) + geom_point(data = tbl_pred, aes(x = x, y =y1), col = "red", alpha= 0.9) + 
#       geom_point(data = tbl_pred, aes(x = x, y =y2), col = "blue", alpha= 0.9) + geom_point(data = tbl_pred, aes(x = x, y =y3), col = "green", alpha= 0.9) +
#       ggtitle(paste("Day", i)) + ylab("Proportion of counts") +  xlab("Counts")
#     
#     print(paste0("Done", d))
#   }
#   
#   tbl_summary$`50%`
#   
#   theta1 = c(9.436588e-01, 9.567756e-01, 9.796791e-01, 9.819549e-01, 9.949175e-01, 9.829168e-01, 9.430644e-01, 6.862853e-01,
#              1.011087e-01, 2.075781e-02, 3.842676e-02, 1.578665e-02, 9.914345e-01, 9.928245e-01)
#   
#   theta2 = c(5.485633e-02, 4.161768e-02, 1.856077e-02, 1.626686e-02, 4.225237e-03, 1.611090e-02, 5.354240e-02, 1.792807e-01,
#              5.596186e-02, 8.544897e-02, 5.343530e-02, 4.235841e-03, 8.285950e-03, 6.099456e-03)
#   
#   theta3 = c(1.155183e-03, 1.079145e-03, 1.127803e-03, 4.819867e-04, 4.752151e-04, 4.759796e-04, 3.189323e-03, 1.796458e-01,
#              8.638063e-01, 9.005260e-01, 9.177685e-01, 9.796872e-01, 4.619559e-04, 4.814225e-04)
#   
#   theta = data.frame(theta1 = theta1, theta2 = theta2, theta3 = theta3) %>% mutate(Day = c(1:14)) %>% gather(theta, value, -Day)
#   
#   theta %>%  ggplot(aes(x=Day, y=value)) + geom_point() + facet_wrap(~theta)
#   
#   # #num_chains = 4
#   # 
#   # #data_for_stan = list(N = nrow(data_x1), x = data_x1$prot_count, day_index=data_x1$day_index, d=day_number)
#   # 
#   # run_fit = stan("Research/RBC/Data Analysis/mixture_models_stan/mix_gamma_model.stan", data = data_for_stan, chains = num_chains, iter = 1500, warmup = 500)
#   # 
#   # #traceplot(run_fit, pars = "theta")
#   # 
#   # save(x, y, file = "xy.RData")
#   # 
#   # tbl_summary = summary(run_fit)$summary %>% as.data.frame() %>% rownames_to_column() %>%  as_tibble()
#   # 
#   # if(any(tbl_summary$Rhat < 0.9 | tbl_summary$Rhat > 1.1)) print("convergence_problem")
#   # 
#   # direct1 = "Desktop/Nov9/Mix_model/"
#   # direct2 = "_K1.csv"
#   # direct = paste(direct1, protein, direct2, sep="")
#   # 
#   # write.csv(tbl_summary, file = direct)
#   # 
#   # fit_out = tbl_summary %>% select(rowname,`2.5%`, `50%`, `97.5%`) 
#   # 
#   # day_vals = c("0", "2", "4", "6", "8", "10", "11", "12", "14", "16", "18", "20", "MNCs", "Jurkats")
#   # 
#   # 
#   # 
#   # for(d in day_vals) {
#   #   
#   #   i = d
#   #   
#   #   kk = which(day_vals == d)
#   # 
#   #   data_x = data_x1 %>% filter(Day == i) %>% select(prot_count)
#   # 
#   #  tbl_pred = tibble(x = 0:max(data_x), y1 =fit_out$`50%`[30+kk]*dnbinom(x, mu = fit_out$`50%`[1], size = fit_out$`50%`[16]), y2 = (1-fit_out$`50%`[30+kk])*dnbinom(x, mu = fit_out$`50%`[1+kk], size = fit_out$`50%`[16+kk]))
#   # 
#   #  pp_out[[kk]] =  data_x %>% 
#   #   ggplot(aes(x=prot_count)) + geom_bar(aes(y = ..count../sum(..count..))) + geom_point(data = tbl_pred, aes(x = x, y =y1), col = "red", alpha= 0.9) + 
#   #   geom_point(data = tbl_pred, aes(x = x, y =y2), col = "blue", alpha= 0.9) +
#   #   ggtitle(paste("Day", i)) + ylab("Proportion of counts") #+  xlim(-1, 50)
#   #  
#   #  print(paste0("Done", d))
#   # }
#   
#   newpic = gridExtra::grid.arrange(grobs = pp_out[1:4], ncol = 2, top=grid::textGrob(protein, gp=grid::gpar(fontsize=15,font=8)))
#   print(paste0("Done", protein))
# }
# dev.off()

