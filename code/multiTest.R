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

samples = 200
d = 8
p = 8
k = d
SNRS = c(10000, 100, 1)
n = sqrt(c(100, 1000, 10000, 100000))

if (simu_pval) {
  set.seed(2020)
  means_m = generateMeans(d, p, k, snr = SNRS, control.g = TRUE)
  v.t = eigen(means_m[1,1,,])$vectors

  numCores = parallel::detectCores()
  cl <- makeCluster(numCores)
  registerDoSNOW(cl)

  simulated_m = simuSamples(means_m, n, samples, prl = T)

  totL = length(n)*samples*(length(SNRS)+1)
  pb <- txtProgressBar(max = totL, style = 3)
  progress <- function(n) {setTxtProgressBar(pb, n)}
  opts <- list(progress = progress)

  data_m = foreach(est_list = simulated_m, .inorder = F, .combine = bind_rows,
                   .options.snow = opts, .packages = c('eigTest', 'tidyverse')) %dopar% {
                     mu.bar = est_list$mu.bar; cov.bar = est_list$cov.bar
                     SNR = as.numeric(str_split(est_list$SNR, '=')[[1]][2])
                     CovRate = est_list$CovRate

                     eigvPLG = JDTE(mu.bar)

                     data_eig = eigTest(mu.bar, cn = CovRate, cov.arr = cov.bar, V = eigvPLG,
                                        testType = 'chi', param.out = T)$chi %>%
                       list2DF() %>% select(pvalue) %>%
                       mutate(testType = 'Chi')
                     data_eig = eigTest(mu.bar, cn = CovRate, cov.arr = cov.bar, V = eigvPLG,
                                        testType = 'gam', param.out = T)$gam %>%
                       list2DF() %>% select(pvalue) %>%
                       mutate(testType = 'Gam') %>% bind_rows(data_eig)


                     if (SNR == 0) {
                       data_eig = eigTest(mu.bar, cn = CovRate, cov.arr = cov.bar,
                                          testType = 'chi', V = v.t, param.out = T)$chi %>%
                         list2DF() %>% select(pvalue) %>%
                         mutate(testType = 'Chi_0') %>% bind_rows(data_eig)
                     }

                     return(data_eig %>%
                              mutate(SNR = paste0('SNR = ', ifelse(SNR==0, '$\\infty$', round(1/SNR))),
                                     SampleSize = round(CovRate^2))
                     )
                   }

  stopCluster(cl)

  data_m$testType = factor(data_m$testType, levels = c('Chi_0', 'Chi', 'Gam'))
  orders = order(unique(data_m$SNR))
  orders = c(orders[1], rev(orders[2:length(orders)]))
  data_m$SNR = factor(data_m$SNR, levels = unique(data_m$SNR)[orders])
  data_m$SampleSize = paste0('Sample size $n = 10^', round(log10(data_m$SampleSize)), '$')
  save(data_m, file = 'output/multiTest.RData')
} else {
  load('output/multiTest.RData')
}


###### plot histograms ######
binwidth = 0.05
breaks = seq(0, 1, 0.05)

tikz(file = "output/Plots/tikz/PvalueMulti.tikz", standAlone=F,width = 7, height = 6)
ggplot(data_m) +
  geom_histogram_pattern(aes(x = pvalue, y = after_stat(density)*binwidth, fill = SNR),
                         pattern = c(rep('none', 20 * length(n)),
                                     rep(c(rep('none', 20),
                                           rep('stripe', 20*length(SNRS))),
                                         length(n) * 2)),
                         pattern_angle = c(rep(0, 20 * length(n)),
                                           rep(rep(seq(0, 180, length = length(SNRS)+2)[-(length(SNRS)+2)], each = 20),
                                               length(n) * 2)) - 90,
                         pattern_size = 0.1, pattern_colour = 'grey80', pattern_spacing = 0.015,
                         breaks = breaks, position = position_dodge()) + theme_bw()+
  ggtitle('p-value histogram from multi-sample test') +
  facet_grid(vars(testType), vars(SampleSize),# ncol = 4, dir = 'v',
             labeller = labeller(testType = c(`Chi_0` = 'Chi test with $V$',
                                              `Chi` = 'Chi test with $\\widehat{V}$',
                                              `Gam` = 'Gamma test with $\\widehat{V}$'))) +
  theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5),
        strip.background =element_rect(linewidth = 0),# fill="white"),
        panel.spacing.x = unit(4, 'mm'), legend.key.height = unit(0.5, 'cm'),
        legend.title = element_text(margin = margin(r = 3, unit = 'mm')),
        legend.text = element_text(margin = margin(r = 5, unit = 'mm')),
        panel.border = element_rect(linewidth = 1), strip.placement = 'inside') +
  guides(fill=guide_legend(ncol=4,
                           override.aes = list(
                             pattern_size = 0.1, pattern_colour = 'grey80', pattern_spacing = 0.015,
                             pattern = c('none', rep('stripe', length(SNRS))),
                             pattern_angle = seq(0, 180, length = length(SNRS)+2)[-(length(SNRS)+2)] - 90
                           ))) + xlab('p-values') + ylab('Proportion') +
  labs(fill = 'SNR:') +
  scale_fill_discrete()+
  scale_y_continuous(breaks = c(0.05, 0.25, 0.5, 0.75, 1), trans = 'sqrt')
dev.off()

p1 = ggplot(data_m %>%
              mutate(testType = recode_factor(testType,
                                              `Chi_0` = 'Chi test with $V$',
                                              `Chi` = 'Chi test with $\\widehat{V}$',
                                              `Gam` = 'Gamma test with $\\widehat{V}$'))) +
  geom_histogram_pattern(aes(x = pvalue, y = after_stat(density)*binwidth, fill = SNR),
                         pattern = c(rep('none', 20 * length(n)),
                                     rep(c(rep('none', 20),
                                           rep('stripe', 20*length(SNRS))),
                                         length(n) * 2)),
                         pattern_angle = c(rep(0, 20 * length(n)),
                                           rep(rep(seq(0, 180, length = length(SNRS)+2)[-(length(SNRS)+2)], each = 20),
                                               length(n) * 2)) - 90,
                         pattern_size = 0.1, pattern_colour = 'grey80', pattern_spacing = 0.015,
                         breaks = breaks, position = position_dodge()) + theme_bw()+
  ggtitle('p-value histogram from multi-sample test') +
  facet_grid(vars(testType), vars(SampleSize),# ncol = 4, dir = 'v',
             labeller = as_labeller(TeX, default = label_parsed)) +
  theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5),
        strip.background =element_rect(linewidth = 0),# fill="white"),
        panel.spacing.x = unit(4, 'mm'), legend.key.height = unit(0.5, 'cm'),
        legend.title = element_text(margin = margin(r = 3, unit = 'mm')),
        legend.text = element_text(margin = margin(r = 5, unit = 'mm')),
        panel.border = element_rect(linewidth = 1), strip.placement = 'inside') +
  guides(fill=guide_legend(ncol=4,
                           override.aes = list(
                             pattern_size = 0.1, pattern_colour = 'grey80', pattern_spacing = 0.015,
                             pattern = c('none', rep('stripe', length(SNRS))),
                             pattern_angle = seq(0, 180, length = length(SNRS)+2)[-(length(SNRS)+2)] - 90
                           ))) + xlab('p-values') + ylab('Proportion') +
  labs(fill = 'SNR:') +
  scale_fill_discrete(labels = TeX(levels(data_m$SNR))) +
  scale_y_continuous(breaks = c(0.05, 0.25, 0.5, 0.75, 1), trans = 'sqrt')
p1
ggsave(filename = 'output/Plots/png/PvalueMulti.png', p1, width = 7, height = 6, units = 'in')


###### Type I/II errors ######
data_m %>% group_by(testType, SNR, SampleSize) %>%
  summarise(RejRate = mean(pvalue <= 0.05)) %>% print(n = nrow(.))
