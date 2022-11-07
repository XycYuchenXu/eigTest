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
d = 4
p = 8
k = 2
SNRS = c(1000, 100, 10)
n = sqrt(c(100, 1000, 10000))

if (simu_pval) {
  set.seed(2020)
  means_p = generateMeans(d, p, k, snr = SNRS, control.g = TRUE)

  numCores = parallel::detectCores()
  cl <- makeCluster(numCores)
  registerDoSNOW(cl)

  simulated_p = simuSamples(means_p, n, samples, prl = T)

  totL = length(n)*samples*(length(SNRS)+1)
  pb <- txtProgressBar(max = totL, style = 3)
  progress <- function(n) {setTxtProgressBar(pb, n)}
  opts <- list(progress = progress)

  data_p = foreach(est_list = simulated_p, .inorder = F, .combine = bind_rows,
                   .options.snow = opts, .packages = c('eigTest', 'reshape2', 'tidyverse')) %dopar% {
                     mu.bar = est_list$mu.bar; cov.bar = est_list$cov.bar
                     SNR = as.numeric(str_split(est_list$SNR, '=')[[1]][2])
                     CovRate = est_list$CovRate

                     p_vector_partial = array(NA, dim = c(2, d-k+1))
                     for (kk in k:d) {
                       if (kk < d) {
                         Qk = expmPartSchur(mu.bar, kk, warmup = T)
                         B = array(0, c(p,kk,kk))
                         for (i in 1:p) {
                           Ai = tcrossprod(crossprod(Qk, mu.bar[i,,]), t(Qk))
                           B[i,,] = Ai[1:kk, 1:kk]
                         }
                         V = JDTE(B)
                       } else {
                         Qk = NULL; V = JDTE(mu.bar)
                       }
                       p_vector_partial[1, kk-k+1] = partialTest(mu.bar, cn = CovRate,
                                                                 cov.arr = cov.bar, Q = Qk, V = V,
                                                                 k = kk, testType = 'chi')
                       print(kk)
                       p_vector_partial[2, kk-k+1] = partialTest(mu.bar, cn = CovRate,
                                                                 cov.arr = cov.bar, Q = Qk, V = V,
                                                                 k = kk, testType = 'gam')
                     }
                     dimnames(p_vector_partial) = list(c('Chi', 'Gam'), paste('K =', k:d))

                     data_p = melt(p_vector_partial, value.name = 'pvalue', na.rm = T)
                     colnames(data_p)[1:2] = c('testType', 'K')

                     return(data_p %>%
                              mutate(SNR = paste0('SNR = ', ifelse(SNR==0, '$\\infty$', round(1/SNR))),
                                     SampleSize = round(CovRate^2))
                     )
                   }
  stopCluster(cl)

  data_p$testType = factor(data_p$testType, levels = c('Chi', 'Gam'))
  orders = order(unique(data_p$SNR))
  orders = c(orders[1], rev(orders[2:length(orders)]))
  data_p$SNR = factor(data_p$SNR, levels = unique(data_p$SNR)[orders])
  data_p$SampleSize = paste0('Sample size $n = 10^', round(log10(data_p$SampleSize)), '$')
  save(data_p, file = 'output/partTest.RData')
} else {
  load('output/partTest.RData')
}


###### plot histograms ######
breaks = seq(0,1,0.05); binwidth = 0.05

tikz(file = "output/Plots/PvaluePartial.tikz", standAlone=F, width = 6, height = 4.5)
ggplot(data_p %>% filter(K == 'K = 2')) +
  geom_histogram_pattern(aes(x = pvalue, y = after_stat(density)*binwidth, fill = SNR),
                         pattern = rep(c(rep('none', 20),
                                         rep('stripe', 20*length(SNRS))), length(n) * 2),
                         pattern_angle = rep(rep(seq(0, 180, length = length(SNRS)+2)[-(length(SNRS)+2)], each = 20),
                                             length(n) * 2) - 90,
                         pattern_size = 0.1, pattern_colour = 'grey80', pattern_spacing = 0.015,
                         breaks = breaks,
                         position = position_dodge()) + theme_bw()+
  ggtitle('p-value histogram from partial test') +
  facet_grid(vars(testType), vars(SampleSize), #ncol = length(n), dir = 'v',
             labeller = labeller(testType = c(`Chi` = 'Chi test with $\\widehat{V}$',
                                              `Gam` = 'Gamma test with $\\widehat{V}$'))) +
  theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5),
        strip.background =element_rect(linewidth = 0), aspect.ratio = 1,
        panel.spacing.x = unit(3, 'mm'), legend.key.height = unit(0.5, 'cm'),
        legend.title = element_text(margin = margin(r = 3, unit = 'mm')),
        legend.text = element_text(margin = margin(r = 5, unit = 'mm')),
        panel.border = element_rect(linewidth = 1), strip.placement = 'inside') +
  guides(fill=guide_legend(ncol=4,
                           override.aes = list(
                             pattern_size = 0.1, pattern_colour = 'grey80', pattern_spacing = 0.015,
                             pattern = c('none', rep('stripe', length(SNRS))),
                             pattern_angle = seq(0, 180, length = length(SNRS)+2)[-(length(SNRS)+2)] - 90
                           ))) +
  xlab('p-values') + ylab('Proportion') +
  labs(fill = 'SNR:') +
  scale_fill_discrete()+
  scale_y_continuous(breaks = c(0.05, 0.25, 0.5, 0.75, 1), trans = 'sqrt')
dev.off()

p1 = ggplot(data_p %>% filter(K == 'K = 2') %>%
              mutate(testType = recode_factor(testType,
                                              `Chi` = 'Chi test with $\\widehat{V}$',
                                              `Gam` = 'Gamma test with $\\widehat{V}$'))) +
  geom_histogram_pattern(aes(x = pvalue, y = after_stat(density)*binwidth, fill = SNR),
                         pattern = rep(c(rep('none', 20),
                                         rep('stripe', 20*length(SNRS))), length(n) * 2),
                         pattern_angle = rep(rep(seq(0, 180, length = length(SNRS)+2)[-(length(SNRS)+2)], each = 20),
                                             length(n) * 2) - 90,
                         pattern_size = 0.1, pattern_colour = 'grey80', pattern_spacing = 0.015,
                         breaks = breaks,
                         position = position_dodge()) + theme_bw()+
  ggtitle('p-value histogram from partial test') +
  facet_grid(vars(testType), vars(SampleSize), #ncol = length(n), dir = 'v',
             labeller = as_labeller(TeX, default = label_parsed)) +
  theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5),
        strip.background =element_rect(linewidth = 0), aspect.ratio = 1,
        panel.spacing.x = unit(3, 'mm'), legend.key.height = unit(0.5, 'cm'),
        legend.title = element_text(margin = margin(r = 3, unit = 'mm')),
        legend.text = element_text(margin = margin(r = 5, unit = 'mm')),
        panel.border = element_rect(linewidth = 1), strip.placement = 'inside') +
  guides(fill=guide_legend(ncol=4,
                           override.aes = list(
                             pattern_size = 0.1, pattern_colour = 'grey80', pattern_spacing = 0.015,
                             pattern = c('none', rep('stripe', length(SNRS))),
                             pattern_angle = seq(0, 180, length = length(SNRS)+2)[-(length(SNRS)+2)] - 90
                           ))) +
  xlab('p-values') + ylab('Proportion') +
  labs(fill = 'SNR:') +
  scale_fill_discrete(labels = TeX(levels(data_p$SNR))) +
  scale_y_continuous(breaks = c(0.05, 0.25, 0.5, 0.75, 1), trans = 'sqrt')
p1
ggsave(filename = 'output/Plots/PvaluePartial.png', p1, width = 6, height = 4.5, units = 'in')


###### Type I/II errors ######
data_p %>% group_by(testType, K, SNR, SampleSize) %>%
  summarise(RejRate = mean(pvalue <= 0.05)) %>% print(n = nrow(.))
