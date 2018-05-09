library(rstan)
library(rlang)
library(tidyverse)
library(doParallel)
library(doMC)

doMC::registerDoMC(cores = 14)
rstan_options(auto_write = TRUE)

# options(mc.cores = parallel::detectCores())
tukey_lim = function(x, k = 2){
  q1_q3 = quantile(x, c(0.75, 0.9))
  upper_lim = q1_q3[2] + k*(q1_q3[2]-q1_q3[1])
  upper_lim
}


fit_neg_binom = function(x, cmp_model){
  num_chains = 4
  data_for_stan = list(N = length(x), x = as.integer(x))
  # run_fit = stan("model_poisson_gamma.stan", iter = 1000, data = data_for_stan, chains = num_chains, thin = 2, cores = 1)
  run_fit = sampling(cmp_model, iter = 1000, data = data_for_stan, chains = num_chains, thin = 2, cores = 10)
  tbl_summary = summary(run_fit)$summary %>% as_tibble() 
  if(any(tbl_summary$Rhat < 0.9 | tbl_summary$Rhat > 1.1)) print("convergence_problem")
  #tbl_summary %>% select(`2.5%`, `50%`, `97.5%`) %>% mutate(type = c("mu", "phi", "a", "lp__"))
  run_fit
}

cmp_model = stan_model(file = "Downloads/model_poisson_gamma.stan")

# x_true_1 = rnbinom(n = 1000, mu = 2000, size = 50)
# x1 = rpois(n = length(x_true_1), x_true_1) %>%  as.data.frame()
# 
# x_true_2 = rnbinom(n = 1000, mu = 20, size = 10)
# x2 = rpois(n = length(x_true_2), x_true_2) %>%  as.data.frame()
# 
# x_true_3 = rnbinom(n = 1000, mu = 10, size = 5)
# x3 = rpois(n = length(x_true_3), x_true_3) %>%  as.data.frame()
# 
# x_true_4 = rnbinom(n = 1000, mu = 10, size = 5)
# x4 = rpois(n = length(x_true_4), x_true_4) %>%  as.data.frame()
# 
# alpha1 = 0.1
# alpha2 = 0.01
# alpha3 = 0.01
# alpha4 = 0.02
# 
# y = rnbinom(n = 1000, mu = 10 + alpha1*x_true_1 + alpha2*x_true_2 + alpha3*x_true_3 + alpha4*x_true_4, size = 10) %>% as.data.frame()
# 
# x = bind_cols(x1, x2, x3, x4)
# 
# y %>% ggplot(aes(y)) + geom_histogram(binwidth = 0.5) 
# 
# 
# fit_out = fit_neg_binom(y, x, cmp_model)

load("Desktop/Dec06/Data1.rda")
load("Desktop/Dec06/Data2.rda")
#load("Desktop/Nov22/Data3.rda")

markers = c("089_CD41", "102Pd", "104Pd", "105Pd", "106Pd", "108Pd", 
            "110Pd", "120Sn", "127_IdU", "131Xe", "133Cs", "138Ba", "140Ce", "141_CD235ab", 
            "142_FlagTag", "143_CD45RA", "144_GlobinHBB", "145_C-EBPa", "146_MEF2C", "147_CD33", "148_p-HH3", "149_CD34", 
            "150_pRb", "151_CD123", "152_MAFG", "153_CD44", "154_NFE2p45", "155_CD36","156_GATA1", "158_RUNX1",
            "159_IKZF1", "160_ATRX", "161_CD90", "162_TAL1", "163_GATA2", "164_CD49F", "165_KLF1", "166_KAT3B",
            "167_PU1", "168_cyclinB1", "169_FLI1", "170Er", "171_CXCR4", "172_CD38", "173_c-JUN", "174_HistoneH3",
            "175_CD71", "176_c-Myc", "190_BCKG", "191Ir", "193Ir", "194Pt", "195Pt", "198Pt",          
            "208Pb")

#new_cell_data = cell_data %>% filter(`166_KAT3B` < tukey_lim(`166_KAT3B`, 20), `156_GATA1` < tukey_lim(`156_GATA1`, 20))

#new_cell_data = cell_data %>% filter(`175_CD71` < tukey_lim(`175_CD71`, 3), `176_c-Myc` < tukey_lim(`176_c-Myc`, 3))
library(gridExtra)

upplim <- function(count_value) {
  
  up_limit = tukey_lim(count_value, 4)
  
  if(up_limit >=10) return(up_limit)
  else return(10.)
  
}

pdf("Desktop/Dec06/hist_K1_K2.pdf", 12,7)
## Loop through proteins
pp_out = list()
i = 1
for(protein in markers){
  
  #protein = "150_pRb" 
    
  x1 = cell_data1 %>% rename(prot_count1 = !!!sym(protein)) %>% select(prot_count1) %>% 
    mutate(prot_count1 = ceiling(prot_count1)) %>% sample_n(1500)  %>% filter(prot_count1 <= upplim(prot_count1)) %>% as_vector()
  
  x2 = cell_data2 %>% rename(prot_count2 = !!!sym(protein)) %>% select(prot_count2) %>% 
    mutate(prot_count2 = ceiling(prot_count2)) %>% sample_n(1500) %>% filter(prot_count2 <= upplim(prot_count2)) %>% as_vector() 
  
  # x3 = cell_data3 %>% rename(prot_count3 = !!!sym(protein)) %>% select(prot_count3) %>% 
  #   mutate(prot_count3 = ceiling(prot_count3)) %>% sample_n(1500) %>% filter(prot_count3 <= upplim(prot_count3)) %>% as_vector()
  
  fit_out1 = fit_neg_binom(x1, cmp_model)
  fit_out2 = fit_neg_binom(x2, cmp_model)
  #fit_out3 = fit_neg_binom(x3, cmp_model)
  
  tbl_summary1 = summary(fit_out1)$summary %>% as.data.frame() %>% rownames_to_column() %>% as_tibble()
  tbl_summary2 = summary(fit_out2)$summary %>% as.data.frame() %>% rownames_to_column() %>% as_tibble()
  #tbl_summary3 = summary(fit_out3)$summary %>% as.data.frame() %>% rownames_to_column() %>% as_tibble()
  
  tbl_pred1 = tibble(x = 0:max(x1), y = dnbinom(x, mu = tbl_summary1$`50%`[1], size = tbl_summary1$`50%`[2]))
  tbl_pred2 = tibble(x = 0:max(x2), y = dnbinom(x, mu = tbl_summary2$`50%`[1], size = tbl_summary2$`50%`[2]))
  #tbl_pred3 = tibble(x = 0:max(x3), y = dnbinom(x, mu = tbl_summary3$`50%`[1], size = tbl_summary3$`50%`[2]))

  # bg_signal <- tbl_summary3$`50%`[1]
  # spill_signal <- tbl_summary2$`50%`[1] - tbl_summary3$`50%`[1]
  # stain_signal <- tbl_summary1$`50%`[1] - tbl_summary2$`50%`[1] - tbl_summary3$`50%`[1]
  # 
  # pp_out[[i]] <- data.frame(prot_name = protein, bg = bg_signal, spill = spill_signal, stain = stain_signal)
  # i = i +1
  pp_out[[1]] = as.tibble(x1) %>%
    ggplot(aes(x1)) + geom_bar(aes(y = ..count../sum(..count..))) + geom_point(data = tbl_pred1, aes(x = x, y = y), col = "red", alpha= 0.9) +
    ggtitle(paste0("K1 ", tbl_summary1$`50%`[1])) + ylab("Proportion of counts")

  pp_out[[2]] = as.tibble(x2) %>%
    ggplot(aes(x2)) + geom_bar(aes(y = ..count../sum(..count..))) + geom_point(data = tbl_pred2, aes(x = x, y =y), col = "red", alpha= 0.9) +
    ggtitle(paste0("K2 ", tbl_summary2$`50%`[1])) + ylab("Proportion of counts")

  # pp_out[[3]] = as.tibble(x3) %>%
  #   ggplot(aes(x3)) + geom_bar(aes(y = ..count../sum(..count..))) + geom_point(data = tbl_pred3, aes(x = x, y =y), col = "red", alpha= 0.9) +
  #   ggtitle(paste0("Unstained_cell ", tbl_summary3$`50%`[1])) + ylab("Proportion of counts")

  #grid.arrange(pp_out1, pp_out2, ncol = 1, main = protein)
  
  newpic = gridExtra::grid.arrange(grobs = pp_out, ncol = 1, top=grid::textGrob(protein, gp=grid::gpar(fontsize=15,font=8)))
  
  # plot(pp_out)
  # plot(pp_out2)
  # plot(pp_out3)
  print(paste0("Done", protein))
}

signal = do.call(rbind, pp_out)
save(signal, file = "Desktop/Nov22/signal.rda")

dev.off()

load("Desktop/Nov22/signal.rda")

signal = signal %>% t() 
colnames(signal) <- signal[1,]
signal = signal[-1,]
signal = signal %>% as.tibble() %>% select(`089_CD41`, `127_IdU`, `140Ce`:`190_BCKG`) 
signal0 = signal %>% data.frame() 
rownames(signal0) <- c("bg", "spill", "stain")

signal1 = signal0 %>% rownames_to_column(var = "signals") %>%  gather(Isotopes, Values, -signals) %>% as.tibble() %>% mutate(Values = as.numeric(Values))
                 
plot1 = signal1 %>% 
  ggplot(aes(x=Isotopes, y=asinh(Values), colour=signals)) + 
  geom_point(alpha = 0.9) +
  #scale_y_log10() +
  ggtitle("Channel Signals") + ylab("Number of counts") 

  geom_point(aes(y = stain), col = "green", alpha= 0.9)

  ggtitle(protein) + ylab("Proportion of counts") #+  xlim(-1, 50)
list(pp = pp_out, fit_out = fit_out)


all_prot = names(gated_data)[-1*(1:2)]
pdf("Neg_Bin_fit.pdf", 14, 9)
## Loop through proteins

for(protein in all_prot){
  
    protein = "142_FlagTag"
    
    prot_data = gated_data %>% rename(prot_count = !!!sym(protein)) %>% mutate(prot_count = ceiling(prot_count)) %>%  
      group_by(Sample) %>% sample_n(1000) %>% filter(prot_count<tukey_lim(prot_count, 3)) %>% ungroup() #%>% select(-Event_Length, -sum_prot, -mean_Ir, -mean_Pt, -Ir191, -Ir193, -Pt194, -Pt195)
    
    prot_data1 = prot_data %>% select(prot_count)
    
    prot_data2 = prot_data %>% select("140Ce", "141_CD235ab", "143_CD45RA", "144_GlobinHBB") %>% ceiling(.)
    
    #prot_data2 = prot_data %>% select(-starts_with(protein), -Cell_id, -prot_count) %>% select(-Day) %>% ceiling(.) %>% bind_cols(Day = prot_data$Day, .) %>% 
    #  group_by(Day) %>% ungroup()
    
    # prot_data1 = gated_data %>% select(Day, Cell_id, !!!sym(protein)) %>% rename(prot_count = !!!sym(protein)) %>% mutate(prot_count = ceiling(prot_count)) %>%  
    #group_by(Day) %>% sample_n(5000) %>% filter(prot_count<tukey_lim(prot_count, 3)) %>% ungroup()
    
    #prot_data2 = gated_data %>% select(-starts_with(protein)) %>% select(Cell_id = pr) select(-Day) %>% ceiling(.) %>% bind_cols(Day = gated_data$Day, .) %>% 
    #  group_by(Day) %>% sample_n(5000) %>% ungroup()
    
    ## Parallel loop through all days
    all_fits = foreach(i = unique(prot_data$Day))%dopar%{
      
    data_x = prot_data1 %>% as_vector()
    data_y = prot_data2 %>% as.matrix()
      
    length(data_x)
    ncol(data_y)
    colnames(data_y)
      
    fit_out = fit_neg_binom(data_x, data_y, cmp_model) #%>% mutate(Day = i)
      
    plot(fit_out, par = "a")
    
    tbl_summary = summary(fit_out)$summary %>% as.data.frame() %>% rownames_to_column() %>% as_tibble()
    
    write.csv(tbl_summary, file = direct0)
    
    save(tbl_summary, "Desktop/Nov22/fitting/tbl_.rda")
    
    load("Desktop/Nov22/fitting/tbl_.rda")
    
    lambda1 <- tbl_summary$rowname %>% stringr::str_extract(",1]")
    lambda2 <- tbl_summary$rowname %>% stringr::str_extract(",2]")
    lambda3 <- tbl_summary$rowname %>% stringr::str_extract(",3]")
    lambda4 <- tbl_summary$rowname %>% stringr::str_extract(",4]")
    
    lambda_id1 <- which(lambda1 == ",1]")
    lambda_id2 <- which(lambda2 == ",2]")
    lambda_id3 <- which(lambda3 == ",3]")
    lambda_id4 <- which(lambda4 == ",4]")
    
    lambda_01 <- tbl_summary[lambda_id1,] %>% select(`50%`) %>% as.data.frame()
    lambda_02 <- tbl_summary[lambda_id2,] %>% select(`50%`) %>% as.data.frame()
    lambda_03 <- tbl_summary[lambda_id3,] %>% select(`50%`) %>% as.data.frame()
    lambda_04 <- tbl_summary[lambda_id4,] %>% select(`50%`) %>% as.data.frame()
    
    selected_data = prot_data %>% select("144_GlobinHBB")
    
    pp_out1 = selected_data %>% ggplot(aes(x=`144_GlobinHBB`)) + geom_histogram(binwidth = 0.25) + xlab("GlobinHBB") + xlim(0,10) +ylim(0, 1000)
    pp_out2 = lambda_04 %>% ggplot(aes(x=`50%`)) + geom_histogram(binwidth = 0.5) + xlab("GlobinHBB_cor") + xlim(0,10)
                                      
    multiplot(pp_out1, pp_out2)
    
    lambda <- bind_cols(lambda_01, lambda_02, lambda_03, lambda_04) %>% as.matrix()
    
    a1 = tbl_summary[3,]$`50%` %>% as.data.frame()
    a2 = tbl_summary[4,]$`50%` %>% as.data.frame()
    a3 = tbl_summary[5,]$`50%` %>% as.data.frame()
    a4 = tbl_summary[6,]$`50%` %>% as.data.frame()
    
    a <- bind_rows(a1, a2, a3, a4) %>% as.matrix()
    
    xx = mean(lambda %*% a)
    
    fit_out1 = tbl_summary %>% select(rowname,`2.5%`, `50%`, `97.5%`) 
    
    tbl_pred = tibble(x = 0:max(data_x), y = dnbinom(x, mu = fit_out1$`50%`[1] + xx, size = fit_out1$`50%`[2]))
    
    pp_out = tibble(x = data_x) %>% 
    ggplot(aes(x)) + geom_bar(aes(y = ..count../sum(..count..))) + geom_point(data = tbl_pred, aes(x = x, y =y), col = "red", alpha= 0.9) +
    ggtitle(protein) + ylab("Proportion of counts") #+  xlim(-1, 50)
    list(pp = pp_out, fit_out = fit_out)
    }
    
    tbl_lambda = 1 #fit_out1 %>% select(rowname) %>% grep(',1]', ., value=TRUE)
    
    pp_mu = bind_rows(map(all_fits, 2)) %>% filter(type =="mu") %>% ggplot(aes(x = Day, y = `50%`)) + geom_point(size = 2) + 
      geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%` ), width = 0.5) + facet_wrap(~type, scale = "free_y")
    
    pp_phi = bind_rows(map(all_fits, 2)) %>% filter(type =="phi") %>% ggplot(aes(x = Day, y = `50%`)) + geom_point(size = 2) + 
      geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%` ), width = 0.5) + facet_wrap(~type, scale = "free_y")
    
    plots_list = map(all_fits, 1)
    
    gridExtra::grid.arrange(grobs = c(plots_list, list(pp_mu), list(pp_phi)),  top=grid::textGrob(protein, gp=grid::gpar(fontsize=15,font=8)))
    
    print(paste0("Done", protein))
    
}

dev.off()


# library(rstan)
# library(rlang)
# library(tidyverse)
# library(doParallel)
# library(doMC)
# 
# doMC::registerDoMC(cores = 14)
# rstan_options(auto_write = TRUE)
# 
# # options(mc.cores = parallel::detectCores())
# load("Data_K2.rda") # gated_data
# 
# tukey_lim = function(x, k = 2){
#   q1_q3 = quantile(x, c(0.75, 0.9))
#   upper_lim = q1_q3[2] + k*(q1_q3[2]-q1_q3[1])
#   upper_lim
# }
# 
# # prot_data %>% 
# #   ggplot(aes(x = prot_count))  + geom_bar(aes(y = ..count../sum(..count..))) + facet_wrap(~Day) #+ xlim(-1, 400) 
# # # + geom_histogram(stat = "density", bins = 100) 
# 
# fit_neg_binom = function(x, y, cmp_model){
#   num_chains = 4
#   data_for_stan = list(N = length(x), x = as.integer(x), num_proteins = ncol(y), y = y)
#   # run_fit = stan("model_poisson_gamma.stan", iter = 1000, data = data_for_stan, chains = num_chains, thin = 2, cores = 1)
#   run_fit = sampling(cmp_model, iter = 1000, data = data_for_stan, chains = num_chains, thin = 2, cores = 4)
#   tbl_summary = summary(run_fit)$summary %>% as_tibble() 
#   if(any(tbl_summary$Rhat < 0.9 | tbl_summary$Rhat > 1.1)) print("convergence_problem")
#   #tbl_summary %>% select(`2.5%`, `50%`, `97.5%`) %>% mutate(type = c("mu", "phi", "a", "lp__"))
#   run_fit
# }
# 
# cmp_model = stan_model(file = "Research/RBC/Data Analysis/Data_October/model_poisson_gamma.stan")
# 
# 
# all_prot = names(gated_data)[-1*(1:2)]
# pdf("Neg_Bin_fit.pdf", 14, 9)
# ## Loop through proteins
# 
# for(protein in all_prot){
#   
#   protein = "TAL1"
#   
#   prot_data = gated_data %>% rename(prot_count = !!!sym(protein)) %>% mutate(prot_count = ceiling(prot_count)) %>%  
#     group_by(Day) %>% sample_n(5000) %>% filter(prot_count<tukey_lim(prot_count, 3)) %>% ungroup() %>% select(-Event_Length, -sum_prot, -mean_Ir, -mean_Pt, -Ir191, -Ir193, -Pt194, -Pt195)
#   
#   prot_data1 = prot_data %>% select(Day, prot_count)
#   
#   prot_data2 = prot_data %>% select(RUNX1, IKZF1, ATRX, CD90, GATA2, CD49f,KLF1,KAT3B, Day) %>% select(-Day) %>% ceiling(.) %>% bind_cols(Day = prot_data$Day, .) %>% 
#     group_by(Day) %>% ungroup()
#   
#   #prot_data2 = prot_data %>% select(-starts_with(protein), -Cell_id, -prot_count) %>% select(-Day) %>% ceiling(.) %>% bind_cols(Day = prot_data$Day, .) %>% 
#   #  group_by(Day) %>% ungroup()
#   
#   # prot_data1 = gated_data %>% select(Day, Cell_id, !!!sym(protein)) %>% rename(prot_count = !!!sym(protein)) %>% mutate(prot_count = ceiling(prot_count)) %>%  
#   #group_by(Day) %>% sample_n(5000) %>% filter(prot_count<tukey_lim(prot_count, 3)) %>% ungroup()
#   
#   #prot_data2 = gated_data %>% select(-starts_with(protein)) %>% select(Cell_id = pr) select(-Day) %>% ceiling(.) %>% bind_cols(Day = gated_data$Day, .) %>% 
#   #  group_by(Day) %>% sample_n(5000) %>% ungroup()
#   
#   ## Parallel loop through all days
#   all_fits = foreach(i = unique(prot_data$Day))%dopar%{
#     
#     i = "2"
#     data_x = prot_data1 %>% filter(Day == i) %>% select(prot_count) %>% as_vector()
#     data_y = prot_data2 %>% filter(Day == i) %>% select(-Day) %>% as.matrix()
#     
#     length(data_x)
#     ncol(data_y)
#     colnames(data_y)
#     
#     fit_out = fit_neg_binom(data_x, data_y, cmp_model) #%>% mutate(Day = i)
#     
#     
#     tbl_pred = tibble(x = 0:max(data_x), y = dnbinom(x, mu = fit_out$`50%`[1], size = fit_out$`50%`[2]))
#     pp_out = tibble(x = data_x) %>% 
#       ggplot(aes(x)) + geom_bar(aes(y = ..count../sum(..count..))) + geom_point(data = tbl_pred, aes(x = x, y =y), col = "red", alpha= 0.9) +
#       ggtitle(paste("Day", i)) + ylab("Proportion of counts") #+  xlim(-1, 50)
#     list(pp = pp_out, fit_out = fit_out)
#   }
#   
#   pp_mu = bind_rows(map(all_fits, 2)) %>% filter(type =="mu") %>% ggplot(aes(x = Day, y = `50%`)) + geom_point(size = 2) + 
#     geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%` ), width = 0.5) + facet_wrap(~type, scale = "free_y")
#   
#   pp_phi = bind_rows(map(all_fits, 2)) %>% filter(type =="phi") %>% ggplot(aes(x = Day, y = `50%`)) + geom_point(size = 2) + 
#     geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%` ), width = 0.5) + facet_wrap(~type, scale = "free_y")
#   
#   plots_list = map(all_fits, 1)
#   
#   gridExtra::grid.arrange(grobs = c(plots_list, list(pp_mu), list(pp_phi)),  top=grid::textGrob(protein, gp=grid::gpar(fontsize=15,font=8)))
#   
#   print(paste0("Done", protein))
#   
# }
# 
# dev.off()

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
