#devtools::install_github('XycYuchenXu/eigTest', force = T, build_vignettes = F)
library(eigTest)
library(foreach)
library(doSNOW)
library(tidyverse)
library(tikzDevice)
library(ggpattern)


###### generate / load pvalues ######
# simulate pvalues from scratch
simu_pval = FALSE

d = 5
p = 2
k = d
samples = 500
SNRS = c(1000, 10)
n = sqrt(c(50, 250, 1000))

if (simu_pval) {
  set.seed(4000)
  means_c = generateMeans(d, p, k, snr = SNRS, control.g = TRUE)

  numCores = parallel::detectCores()
  cl <- makeCluster(numCores)
  registerDoSNOW(cl)
  simulated_c = simuSamples(means_c, n, samples, prl = T)

  totL = length(n)*samples*(length(SNRS)+1)
  pb <- txtProgressBar(max = totL, style = 3)
  progress <- function(n) {setTxtProgressBar(pb, n)}
  opts <- list(progress = progress)

  data_c = foreach(est_list = simulated_c, .inorder = F, .combine = bind_rows,
                   .options.snow = opts, .packages = c('eigTest', 'tidyverse')) %dopar% {
                     mu.bar = est_list$mu.bar; cov.bar = est_list$cov.bar
                     SNR = as.numeric(str_split(est_list$SNR, '=')[[1]][2])
                     CovRate = est_list$CovRate

                     data_temp = commutatorTest(mu.bar, cn = CovRate, cov.arr = cov.bar,
                                                param.out = T)[[1]] %>%
                       list2DF() %>% select(pvalue) %>%
                       mutate(SNR = paste0('SNR = ', ifelse(SNR==0, '$\\infty$', round(1/SNR))),
                              SampleSize = paste0('Sample size $n = ', round(CovRate^2), '$'))

                     return(data_temp)
                   }
  stopCluster(cl)
  orders = order(unique(data_c$SNR))
  orders = c(orders[1], rev(orders[2:length(orders)]))
  data_c$SNR = factor(data_c$SNR, levels = unique(data_c$SNR)[orders])
  data_c$SampleSize = factor(data_c$SampleSize, levels = paste0('Sample size $n = ', round(n^2), '$'))
  save(data_c, file = 'output/commutatorTest.RData')
} else {
  load('output/commutatorTest.RData')
}


###### plot histograms ######
binwidth = 0.05
breaks = seq(0, 1, 0.05)

tikz(file = "output/Plots/tikz/PvalueCommutator.tikz", standAlone=F,width = 6, height = 3)
ggplot(data_c) +
  geom_histogram_pattern(aes(x = pvalue, y = after_stat(density)*binwidth, fill = SNR),
                         pattern = rep(c(rep('none', 20),
                                         rep('stripe', 20*length(SNRS))), length(n)),
                         pattern_angle = rep(rep(seq(0, 180, length = length(SNRS)+2)[-(length(SNRS)+2)], each = 20),
                                             length(n)) - 90,
                         pattern_size = 0.1, pattern_colour = 'grey80', pattern_spacing = 0.015,
                         breaks = breaks,
                         position = position_dodge()) + theme_bw()+
  ggtitle('p-value histogram from commutator-based test') +
  facet_wrap(~SampleSize, ncol = length(n), dir = 'v') +scale_y_continuous(breaks = c(0.05, 0.25, 0.5, 0.75, 1), trans = 'sqrt')+
  theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(linewidth = 1), aspect.ratio = 1,
        panel.spacing.x = unit(4, 'mm'), legend.key.height = unit(0.5, 'cm'),
        legend.title = element_text(margin = margin(r = 3, unit = 'mm')),
        legend.text = element_text(margin = margin(r = 5, unit = 'mm')),
        strip.background =element_rect(fill="white", linewidth = 0)) +
  guides(fill=guide_legend(ncol=3,
                           override.aes = list(
                             pattern_size = 0.1, pattern_colour = 'grey80', pattern_spacing = 0.015,
                             pattern = c('none', rep('stripe', length(SNRS))),
                             pattern_angle = seq(0, 180, length = length(SNRS)+2)[-(length(SNRS)+2)] - 90
                           ))) +
  xlab('p-values') + ylab('Proportion') +
  labs(fill = 'SNR:')
dev.off()

p1 = ggplot(data_c) +
  geom_histogram_pattern(aes(x = pvalue, y = after_stat(density)*binwidth, fill = SNR),
                         pattern = rep(c(rep('none', 20),
                                         rep('stripe', 20*length(SNRS))), length(n)),
                         pattern_angle = rep(rep(seq(0, 180, length = length(SNRS)+2)[-(length(SNRS)+2)], each = 20),
                                             length(n)) - 90,
                         pattern_size = 0.1, pattern_colour = 'grey80', pattern_spacing = 0.015,
                         breaks = breaks,
                         position = position_dodge()) + theme_bw()+
  ggtitle('p-value histogram from commutator-based test') +
  facet_wrap(~SampleSize, ncol = length(n), dir = 'v', labeller = as_labeller(TeX, default = label_parsed)) +
  scale_y_continuous(breaks = c(0.05, 0.25, 0.5, 0.75, 1), trans = 'sqrt')+
  theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(linewidth = 1), aspect.ratio = 1,
        panel.spacing.x = unit(4, 'mm'), legend.key.height = unit(0.5, 'cm'),
        legend.title = element_text(margin = margin(r = 3, unit = 'mm')),
        legend.text = element_text(margin = margin(r = 5, unit = 'mm')),
        strip.background =element_rect(fill="white", linewidth = 0)) +
  guides(fill=guide_legend(ncol=3,
                           override.aes = list(
                             pattern_size = 0.1, pattern_colour = 'grey80', pattern_spacing = 0.015,
                             pattern = c('none', rep('stripe', length(SNRS))),
                             pattern_angle = seq(0, 180, length = length(SNRS)+2)[-(length(SNRS)+2)] - 90
                           ))) +
  xlab('p-values') + ylab('Proportion') +
  scale_fill_discrete(labels = TeX(levels(data_c$SNR))) +
  labs(fill = 'SNR:')
p1
ggsave(filename = 'output/Plots/png/PvalueCommutator.png', p1, width = 6, height = 3, units = 'in')


###### Type I/II errors ######
data_c %>% group_by(SNR, SampleSize) %>%
  summarise(RejRate = mean(pvalue <= 0.05)) %>% print(n = nrow(.))

