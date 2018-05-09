library(rstan)
library(rlang)
library(tidyverse)
library(doParallel)

doMC::registerDoMC(cores = 14)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

tukey_lim = function(x, k = 2){
  q1_q3 = quantile(x, c(0.75, 0.9))
  upper_lim = q1_q3[2] + k*(q1_q3[2]-q1_q3[1])
  upper_lim
}

Yule_func = function(t_i, lambda, n0, n_max){
  n_vec = 1:n_max
  vec_out = choose(n_vec-1, n_vec-n0)*exp(-lambda*n0*t_i)*(1-exp(-lambda*t_i))^(n_vec-n0)
  vec_out[n_max] = 1- sum(vec_out[-n_max])
  vec_out
}

state_data = count3 %>% filter(comb_code != "0 0 1" & comb_code != "1 1 0") %>% mutate(state = c(2,3,4,1)) %>%
  select(-comb_code) %>% arrange(state) %>% select(-`11`) 

ham_dist = day_levels %>% as.tibble() %>% mutate(state=c(2,3,4,1)) %>% arrange(state)

save(state_data, file = "Desktop/Mar2018/state_data.rda")
save(ham_dist, file = "Desktop/Mar2018/ham_dist.rda")

load(file = "Desktop/Mar2018/state_data.rda")
load(file = "Desktop/Mar2018/ham_dist.rda")

data_t = state_data %>% select(-state) %>% as.matrix() %>% t() #`10`, `12`, `14`, `16`, `18`, `20`
data_t_1 = state_data %>% select(-state, `2`) %>% as.matrix() %>% t() # `2`, -`14`, -`16`, -`18`, -`20`

d_data_t = (data_t_1 - data_t)/2;

n_state = ncol(data_t)
n_time = nrow(data_t)

y = ham_dist %>% as.matrix() #%>% select(-state) 

# data_t = compare_data %>% filter(time != 20) %>% select(-time) %>% as.matrix()
# data_t_1 = compare_data %>% filter(time != 2) %>% select(-time) %>% as.matrix()

cmp_model = stan_model(file = "Research/RBC/Data Analysis/Mix_Gaussian_March2018/dirich_state.stan")

fit_Yule_func = function(cmp_model){
  
  num_chains = 4

  data_for_stan = list(n_state = n_state, n_time = n_time, x_t = data_t, y = y) 
  
  run_fit = sampling(cmp_model, iter = 100, warmup = 50, data = data_for_stan, chains = num_chains, thin = 1, cores = 14)
  
  run_fit
}

run_fit = fit_Yule_func(cmp_model)

traceplot(run_fit)

save(run_fit, file = "Desktop/Mar2018/fit_data_3.RData")
load(file = "Desktop/Mar2018/fit_data_1.RData")

x_t_1 = rstan::summary(run_fit, par = "x_t_1") 
x_t_1_new = x_t_1$summary
x_t_1_new$mean
tbl_summary = summary(run_fit)$summary %>% as.data.frame() %>% rownames_to_column() %>%  as_tibble()

lambda = tbl_summary$`50%`

lambda = lambda[7:46]

lambda_new = lambda %>% matrix(nrow =4)


#lambda = c(2.599334e-01, 3.2e-01, 0)
lambda1 = c(6.5e-01, 1.568821e-03, 5.224124e-04)

lambda1 = c(0.4, 0.07, 0.08)
lambda2 = c(3.0e-01, 4.0e-01, 0.025)
lambda3 = c(4.330774e+00, 5e-01, 4e-01 )

hill <- function(Kd=1e6, L, n =5){
  
  y = L^n / (Kd + L^n)
  
  vals = y
  
  return(vals)
  
}

qplot(0:3, hill(2, 0:3, 1)) + ylim(0,1)

lambda_func(lambda1, lambda2, lambda3, data_t)

lambda_func <- function(lambda1, lambda2, lambda3, data_t) {
  
prob_0 = data_t[1,]

time = seq(2, 20, 0.1)

prob_1 = prob_0[1]
prob_2 = prob_0[2]
prob_3 = prob_0[3]
prob_4 = prob_0[4]

prob1 <- c()
prob2 <- c()
prob3 <- c()
prob4 <- c()

ntime = length(time)

for(t in 1:ntime){
  
  tp = time[t]
  print(tp)
  
  if(tp==2) {
    
    prob1 = append(prob1, prob_1)
    prob2 = append(prob2, prob_2)
    prob3 = append(prob3, prob_3)
    prob4 = append(prob4, prob_4)
    
  }
  
  else {
    
  if(time[t] <11) lambda = lambda1
  if(time[t] > 11 & time[t] < 15) lambda = lambda2
  if(time[t] >= 15) {
    
    lambda = lambda3
    lambda[3] = lambda[3]*hill(2, time[t]-15, 1)
    
  }
  
  dt = time[t] - time[t-1]  
  
  prob_1_t = prob_1
  prob_2_t = prob_2
  prob_3_t = prob_3
  prob_4_t = prob_4
  
  prob_1 = ((1-lambda[1]*(dt))*prob_1_t) #+ lambda[4]
  prob_2 = (lambda[1]*prob_1_t*(dt) + (1-lambda[2]*(dt))*prob_2_t) #+ lambda[5]
  prob_3 = (lambda[2]*prob_2_t*(dt) + (1-lambda[3]*(dt))*prob_3_t) #+ lambda[6]
  
  #prob_4 = 1-(prob_1+prob_2+prob_3)
  prob_4 = (lambda[3]*(dt)*prob_3_t) + prob_4_t#+ lambda[7]
  
  prob1 = append(prob1, prob_1)
  prob2 = append(prob2, prob_2)
  prob3 = append(prob3, prob_3)
  prob4 = append(prob4, prob_4)
  }
  
}

compare_data = data.frame(time, State1 = prob1, State2 = prob2, State3 = prob3, State4 = prob4)
data_t_0 = state_data %>% select(-state) %>% as.matrix() %>% t() #, -`14`, -`16`, -`18`, -`20`
original_data = data_t_0 %>% as.tibble() %>% mutate(t=seq(2,20,2))

pp1 = compare_data %>% select(-time) %>% #rowSums()
  ggplot() +
  geom_line(data=compare_data, aes(x=time, y=State1)) +
  geom_point(data = original_data, aes(x=t, y=V1), colour="blue", size = 3.0) +
  geom_vline(xintercept = c(11, 15), col = "red", lty = 2, alpha= 0.7) +
  geom_line(data=compare_data_1, aes(x=time, y=State1), col = "blue", alpha = 0.5) +
  theme_classic() +
  theme(text = element_text(size=20)) 

pp2 = ggplot() +
  geom_line(data=compare_data, aes(x=time, y=State2)) +
  geom_point(data = original_data, aes(x=t, y=V2), colour="blue", size = 3.0) +
  geom_vline(xintercept = c(11, 15), col = "red", lty = 2, alpha= 0.7) +
  geom_line(data=compare_data_1, aes(x=time, y=State2), col = "blue", alpha = 0.5) +
  theme_classic() +
  theme(text = element_text(size=20)) 

pp3 = ggplot() +
  geom_line(data=compare_data, aes(x=time, y=State3)) +
  geom_point(data = original_data, aes(x=t, y=V3), colour="blue", size = 3.0) +
  geom_vline(xintercept = c(11, 15), col = "red", lty = 2, alpha= 0.7) +
  geom_line(data=compare_data_1, aes(x=time, y=State3), col = "blue", alpha = 0.5) +
  theme_classic() +
  theme(text = element_text(size=20)) 

pp4 = ggplot() +
  geom_line(data=compare_data, aes(x=time, y=State4)) +
  geom_point(data = original_data, aes(x=t, y=V4), colour="blue", size = 3.0) +
  geom_vline(xintercept = c(11, 15), col = "red", lty = 2, alpha= 0.7) +
  geom_line(data=compare_data_1, aes(x=time, y=State4), col = "blue", alpha = 0.5) +
  theme_classic() +
  theme(text = element_text(size=20)) 

gridExtra::grid.arrange(pp1, pp2, pp3, pp4)
}

load(file = "Desktop/Mar2018/comp_1.RData")

prob_0 = data_t[1,]

time = seq(2, 20, 0.1)

prob_1 = prob_0[1]
prob_2 = prob_0[2]
prob_3 = prob_0[3]
prob_4 = prob_0[4]

prob1 <- c()
prob2 <- c()
prob3 <- c()
prob4 <- c()

ntime = length(time)

for(t in 1:ntime){
  
  tp = time[t]
  print(tp)
  
  if(tp==2) {
    
    prob1 = append(prob1, prob_1)
    prob2 = append(prob2, prob_2)
    prob3 = append(prob3, prob_3)
    prob4 = append(prob4, prob_4)
    
  }
  
  else {
    
    # if(time[t] <11) lambda = lambda1
    # if(time[t] > 11 & time[t] < 15) lambda = lambda2
    # if(time[t] >= 15) {
    #   
    #   lambda = lambda3
    #   lambda[3] = lambda[3]*hill(2, time[t]-15, 1)
    #   
    # }
    # 
    dt = time[t] - time[t-1]  
    
    prob_1_t = prob_1
    prob_2_t = prob_2
    prob_3_t = prob_3
    prob_4_t = prob_4
    
    prob_1 = ((1-lambda[1]*(dt))*prob_1_t) #+ lambda[4]
    prob_2 = (lambda[1]*prob_1_t*(dt) + (1-lambda[2]*(dt))*prob_2_t) #+ lambda[5]
    prob_3 = (lambda[2]*prob_2_t*(dt) + (1-lambda[3]*(dt))*prob_3_t) #+ lambda[6]
    
    #prob_4 = 1-(prob_1+prob_2+prob_3)
    prob_4 = (lambda[3]*(dt)*prob_3_t) + prob_4_t#+ lambda[7]
    
    prob1 = append(prob1, prob_1)
    prob2 = append(prob2, prob_2)
    prob3 = append(prob3, prob_3)
    prob4 = append(prob4, prob_4)
  }}


# compare_data_1 = data.frame(time, State1 = prob1, State2 = prob2, State3 = prob3, State4 = prob4)
# save(compare_data_2, file = "Desktop/Mar2018/comp_2.RData")
# 
# compare_data_3 = data.frame(time, State1 = prob1, State2 = prob2, State3 = prob3, State4 = prob4)
# save(compare_data_3, file = "Desktop/Mar2018/comp_3.RData")





data_t_0 = state_data %>% select(-state) %>% as.matrix() %>% t() #, -`14`, -`16`, -`18`, -`20`
original_data = data_t_0 %>% as.tibble() %>% mutate(t=seq(2,20,2))

t = seq(2, 20, 2)

compare_data = data.frame(time = t, 
                      State1 = lambda_new[1,], 
                      State2 = lambda_new[2,], 
                      State3 = lambda_new[3,], 
                      State4 = lambda_new[4,]) %>% as.tibble()

compare_data %>% select(-time) %>% rowSums()

pp1 = ggplot() +
  geom_line(data=compare_data, aes(x=time, y=State1)) +
  geom_point(data = original_data, aes(x=t, y=V1), colour="red", size = 3.0) 

pp2 = ggplot() +
  geom_line(data=compare_data, aes(x=time, y=State2)) +
  geom_point(data = original_data, aes(x=t, y=V2), colour="red", size = 3.0)

pp3 = ggplot() +
  geom_line(data=compare_data, aes(x=time, y=State3)) +
  geom_point(data = original_data, aes(x=t, y=V3), colour="red", size = 3.0)

pp4 = ggplot() +
  geom_line(data=compare_data, aes(x=time, y=State4)) +
  geom_point(data = original_data, aes(x=t, y=V4), colour="red", size = 3.0)

gridExtra::grid.arrange(pp1, pp2, pp3, pp4)

xx = tbl_summary %>% mutate(parameter = y, state= z, type = p) %>% filter(type==0) 
xx$`50%` %>% matrix(30, 5)
#protein probability
xx = tbl_summary %>% mutate(parameter = y, state= z, type = p) %>% 
  filter(type==0) %>%  
  ggplot(aes(x=state, y = `50%`)) +
  geom_line() + geom_point(aes(colour=factor(state)), size = 2) + 
  geom_line(data=test_data_P_, aes(x=state, y=`50%`, colour="red"))+
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.25) + facet_wrap(~parameter, ncol = 5)

#cell cycle probability
tbl_summary %>% mutate(parameter = y, state= z, type = p) %>% filter(type==2) %>%  ggplot(aes(x=parameter, y = `50%`)) +
  geom_point(aes(colour="red"), size = 2) + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.25) + facet_wrap(~state, ncol = 3)


all_prot = names(gated_data)[-1*(1:2)]
pdf("Desktop/Jan2018/K1_code.pdf", 14, 9)

day_data = cell_code %>% select(obstime) 
day_data = day_data %>% mutate(obstime = obstime %>% as.character() %>% as.numeric()) 
day_data = day_data$obstime 

cell_data = cell_code %>% select(ATRX:TAL1) %>% as.tibble() %>% as.matrix()
   
dim(day_data)
nrow(cell_data)
   
protein_data = ncol(cell_data)
state_number = 5

day_index = cell_code %>% group_indices(obstime)
days = unique(day_data) %>% as.matrix() %>% sort(decreasing = FALSE)
day_value = length(unique(day_index))

cycle_data = cell_code %>% select(cycle_state)

cycle_code <- c()
  
for(i in 1:nrow(cycle_data)){
  
  t = cycle_data[i,]
  
  if(t==0) correct_vals = c(1, 0, 0, 0)
  if(t==1) correct_vals = c(0, 1, 0, 0)
  if(t==2) correct_vals = c(0, 0, 1, 0)
  if(t==3) correct_vals = c(0, 0, 0, 1)
  
  cycle_code = rbind(cycle_code, correct_vals)
  
}

cycle_code %>% as.matrix()
dim(cycle_code)

prot_comb_data = cell_data %>% unique.array() %>% apply(., 1, paste, collapse = "")
cell_comb_data = cell_data %>% apply(., 1, paste, collapse = "")
n_prot_comb = length(prot_comb_data)

prot_comb = matrix(0, nrow = nrow(cell_data), ncol = n_prot_comb)

for(i in 1:nrow(cell_data)){
  
  prot_comb_index = which(prot_comb_data == cell_comb_data[i])
  prot_comb[i,prot_comb_index] = 1
  
}

run_fit = fit_Yule_func(cell_data, day_data, protein_data, state_number, day_value, day_index, days, cycle_code, prot_comb, n_prot_comb, cmp_model)

cell_data[1,]

x = c("ATRX", "CD123", "CD235ab", "CD33", "CD34", "CD36", "CD38", "CD41", "CD44", "CD45RA", "CD49f", "CD71", "CD90", "CEBPa", "c_Jun", "c_Myc",
      "CXCR4", "FLI1", "GATA1", "GATA2", "HBA", "IKZF1", "KAT3B", "KLF1", "MAFG", "MEF2C", "NFE2p45", "PU1", "RUNX1", "TAL1")

c = c("G0", "G1", "S", "G2/M")

y = c(paste0(x, "_1"), paste0(x, "_2"), paste0(x, "_3"), paste0(x, "_4"), "lambda1", "lambda2", "lambda3", "lambda4", "lp_")

y = c(x,x,x,x,x,x, "lambda1", "lambda2", "lambda3", "lambda4","lambda5","lambda6",
      c,c,c,c,c,c, "lp_")

z = c(rep(1,30), rep(2,30), rep(3,30), rep(4,30), rep(5,30), rep(6,30), rep(7,6), 
      rep(1, 4), rep(2, 4), rep(3,4), rep(4,4), rep(5,4), rep(6,4), 14)

p = c(rep(0, 180), rep(1, 6), rep(2, 24), 3)

###5 states
y = c(rep(x, 5), "lambda1", "lambda2", "lambda3", "lambda4", "lambda5",
      rep(c, 5), rep("pi_", 55), "lp_")

z = c(rep(1,30), rep(2,30), rep(3,30), rep(4,30), rep(5,30), rep(7,5), 
      rep(1, 4), rep(2, 4), rep(3, 4), rep(4, 4), rep(5, 4), rep(6, 55), 14)

p = c(rep(0, 150), rep(1, 5), rep(2, 20), rep(3, 55), 4)

traceplot(run_fit)#, par="P_[5,3]")
save(run_fit, file = "Desktop/Jan2018/5states_simulate_j.RData")
load(file = "Desktop/Jan2018/5states_filter_nocycle.RData")

tbl_summary = summary(run_fit)$summary %>% as.data.frame() %>% rownames_to_column() %>%  as_tibble()
xx = tbl_summary %>% mutate(parameter = y, state= z, type = p) %>% filter(type==0) 
xx$`50%` %>% matrix(30, 5)
#protein probability
xx = tbl_summary %>% mutate(parameter = y, state= z, type = p) %>% 
  filter(type==0) %>%  
  ggplot(aes(x=state, y = `50%`)) +
  geom_line() + geom_point(aes(colour=factor(state)), size = 2) + 
  geom_line(data=test_data_P_, aes(x=state, y=`50%`, colour="red"))+
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.25) + facet_wrap(~parameter, ncol = 5)

#cell cycle probability
tbl_summary %>% mutate(parameter = y, state= z, type = p) %>% filter(type==2) %>%  ggplot(aes(x=parameter, y = `50%`)) +
  geom_point(aes(colour="red"), size = 2) + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.25) + facet_wrap(~state, ncol = 3)

#state transition rate
tbl_summary %>% mutate(parameter = y, state= z, type = p) %>% filter(type==3) %>%  ggplot(aes(x=parameter, y = `50%`)) +
  scale_y_log10() +
  geom_point(aes(colour=factor(state)), size = 2) + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.25) 

fit_out = tbl_summary %>% select(rowname,`2.5%`, `50%`, `97.5%`) #%>% mutate(type = c("mu1", "mu2", "phi1", "phi2", "theta", "lp__"))

if(any(tbl_summary$Rhat < 0.9 | tbl_summary$Rhat > 1.1)) print("convergence_problem")

pp_phi = bind_rows(pp_mu1, pp_mu2) %>% filter(type =="phi") %>% ggplot(aes(x = Day, y = `50%`, colour = data_set)) + geom_point(size = 2) +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%` ), width = 0.5) + facet_wrap(~type, scale = "free_y")

#write.csv(tbl_summary, file = direct0)

fit_out = tbl_summary %>% select(rowname,`2.5%`, `50%`, `97.5%`) 

day_vals = c("2", "8", "12", "14", "16", "20") #c("0", "2", "4", "6", "8", "10", "11", "12", "14", "16", "18", "20", "MNCs", "Jurkats")

for(d in day_vals){
  
  i = d
  
  kk = which(day_vals == d)
  print(kk)
  data_x = data_x1 %>% filter(Day == i) %>% select(prot_count)
  
  tbl_pred = tibble(x = c(0, 1:max(data_x)), 
                    y1 = c((fit_out$`50%`[15+kk-1])*(fit_out$`50%`[21+kk-1] + (1-fit_out$`50%`[21+kk-1])*dnbinom(0, mu = fit_out$`50%`[1], size = fit_out$`50%`[8])), 
                           (fit_out$`50%`[15+kk-1])*((1-fit_out$`50%`[21+kk-1])*dnbinom(x[-1], mu = fit_out$`50%`[1], size = fit_out$`50%`[8]))), 
                    y2 = (1-fit_out$`50%`[15+kk-1])*dnbinom(x, mu = fit_out$`50%`[2+kk-1], size = fit_out$`50%`[9+kk-1]))
  
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
#   ggplot(aes(x=prot_count)) + geom_bar(aes(y = ..count../sum(..count..))) + geom_point(data = tbl_pred, aes(x = x, y =y1), col = "red", alpha= 0.9) + 
#   geom_point(data = tbl_pred, aes(x = x, y =y2), col = "blue", alpha= 0.9) +
#   ggtitle(paste("Day", i)) + ylab("Proportion of counts") #+  xlim(-1, 50)
#  
#  print(paste0("Done", d))
# }

newpic = gridExtra::grid.arrange(grobs = pp_out[1:6], ncol = 3, top=grid::textGrob(protein, gp=grid::gpar(fontsize=15,font=8)))
print(paste0("Done", protein))
dev.off()
}
dev.off()


all_prot = names(gated_data_k1)[-1*(1:2)] 

all_prot = all_prot[-c(14, 37, 38, 43, 44)]

pdf("Desktop/Nov9/Mix_model/fit_mixed_K1.pdf", 14, 9)

## Loop through proteins
for(protein in all_prot){
  
  prot_data = gated_data_k1 %>% select(Day, !!!sym(protein)) %>% rename(prot_count = !!!sym(protein)) %>% mutate(prot_count = ceiling(prot_count)) %>%  
    group_by(Day) %>% sample_n(1000) %>% filter(prot_count<tukey_lim(prot_count, 3)) 
  
  day_index = prot_data %>% group_indices(.) 
  
  data_x1 = prot_data %>% ungroup() %>% select(Day, prot_count) %>% bind_cols(day_index = data.frame(day_index), .) %>% as.tibble()
  
  day_number = length(table(data_x1$day_index))
  
  direct1 = "Desktop/Nov9/Mix_model/"
  direct2 = "_K1.csv"
  direct3 = "_K1.RData"
  direct4 = "_K1.pdf"
  direct0 = paste(direct1, protein, direct2, sep="")
  direct = paste(direct1, protein, direct3, sep="")
  directp = paste(direct1, protein, direct4, sep="") 
  
  #pdf(directp, 14, 9)
  
  load(direct) # gated_data
 
  tbl_summary = summary(run_fit)$summary %>% as.data.frame() %>% rownames_to_column() %>%  as_tibble()
  
  if(any(tbl_summary$Rhat < 0.9 | tbl_summary$Rhat > 1.1)) print("convergence_problem")
  
  fit_out = tbl_summary %>% select(rowname,`2.5%`, `50%`, `97.5%`) 
  
  day_vals = c("0", "2", "4", "6", "8", "10", "11", "12", "14", "16", "18", "20", "MNCs", "Jurkats")
  
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

  newpic = gridExtra::grid.arrange(grobs = pp_out[1:14], ncol = 4, top=grid::textGrob(protein, gp=grid::gpar(fontsize=15,font=8)))

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

