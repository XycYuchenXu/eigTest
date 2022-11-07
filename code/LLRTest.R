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
  means_l = generateMeans(d, p, k, snr = SNRS, control.g = TRUE)

  numCores = parallel::detectCores()
  cl <- makeCluster(numCores)
  registerDoSNOW(cl)
  simulated_l = simuSamples(means_l, n, samples, prl = T)

  eigvORC = JDTE(means_l[,1,,], iter = 1000)

  totL = length(n)*samples*(length(SNRS)+1)
  pb <- txtProgressBar(max = totL, style = 3)
  progress <- function(n) {setTxtProgressBar(pb, n)}
  opts <- list(progress = progress)

  data_l = foreach(est_list = simulated_l, .inorder = F, .combine = bind_rows,
                   .options.snow = opts, .packages = c('eigTest', 'tidyverse')) %dopar% {
                     mu.bar = est_list$mu.bar; cov.bar = est_list$cov.bar
                     SNR = as.numeric(str_split(est_list$SNR, '=')[[1]][2])
                     CovRate = est_list$CovRate

                     eigvPLG = JDTE(mu.bar)

                     data_LLR = projTest(mu.bar, cn = CovRate, param.out = T,
                                         CV = eigvPLG, poly.sp = F) %>%
                       list2DF() %>% select(pvalue) %>%
                       mutate(covType = 'Oracle', spaceType = 'Plugin eqref{eqn:eigvP}')
                     data_LLR = projTest(mu.bar, cn = CovRate, cov.arr = cov.bar,
                                         poly.sp = F, CV = eigvPLG, param.out = T) %>%
                       list2DF() %>% select(pvalue) %>%
                       mutate(covType = 'Plugin', spaceType = 'Plugin eqref{eqn:eigvP}') %>%
                       bind_rows(data_LLR)
                     data_LLR = projTest(mu.bar, cn = CovRate, param.out = T) %>% list2DF() %>%
                       select(pvalue) %>%
                       mutate(covType = 'Oracle', spaceType = 'Plugin eqref{eqn:polyP}') %>%
                       bind_rows(data_LLR)
                     data_LLR = projTest(mu.bar, cn = CovRate, cov.arr = cov.bar,
                                         param.out = T) %>%
                       list2DF() %>% select(pvalue) %>%
                       mutate(covType = 'Plugin', spaceType = 'Plugin eqref{eqn:polyP}') %>%
                       bind_rows(data_LLR)

                     if (SNR == 0) {
                       data_LLR = projTest(mu.bar, refMat = means_l[,1,,],
                                           cn = CovRate, param.out = T) %>%
                         list2DF() %>% select(pvalue) %>%
                         mutate(covType = 'Oracle', spaceType = 'Oracle eqref{eqn:polyP}') %>%
                         bind_rows(data_LLR)
                       data_LLR = projTest(mu.bar, cn = CovRate, param.out = T,
                                           CV = eigvORC, poly.sp = F) %>%
                         list2DF() %>% select(pvalue) %>%
                         mutate(covType = 'Oracle', spaceType = 'Oracle eqref{eqn:eigvP}') %>%
                         bind_rows(data_LLR)
                       data_LLR = projTest(mu.bar, cn = CovRate, refMat = means_l[,1,,],
                                           cov.arr = cov.bar, param.out = T) %>%
                         list2DF() %>% select(pvalue) %>%
                         mutate(covType = 'Plugin', spaceType = 'Oracle eqref{eqn:polyP}') %>%
                         bind_rows(data_LLR)
                       data_LLR = projTest(mu.bar, cn = CovRate, cov.arr = cov.bar,
                                           CV = eigvORC, poly.sp = F, param.out = T) %>%
                         list2DF() %>% select(pvalue) %>%
                         mutate(covType = 'Plugin', spaceType = 'Oracle eqref{eqn:eigvP}') %>%
                         bind_rows(data_LLR)
                     }
                     #                       })
                     return(data_LLR %>%
                              mutate(SNR = paste0('SNR = ', ifelse(SNR==0, '$\\infty$', round(1/SNR))),
                                     SampleSize = paste0('Sample size $n = ', round(CovRate^2), '$'))
                     )
                   }
  stopCluster(cl)
  orders = order(unique(data_l$SNR))
  orders = c(orders[1], rev(orders[2:length(orders)]))
  data_l$SNR = factor(data_l$SNR, levels = unique(data_l$SNR)[orders])
  data_l$SampleSize = factor(data_l$SampleSize, levels = paste0('Sample size $n = ', round(n^2), '$'))
  data_l$spaceType = factor(data_l$spaceType, levels = c('Oracle eqref{eqn:polyP}',
                                                         'Plugin eqref{eqn:polyP}',
                                                         'Oracle eqref{eqn:eigvP}',
                                                         'Plugin eqref{eqn:eigvP}'))
  save(data_l, file = 'output/LLRTest.RData')
} else {
  load('output/LLRTest.RData')
}


###### plot histograms ######
binwidth = 0.05
breaks = seq(0, 1, 0.05)

tikz(file = "output/Plots/tikz/PvalueLLR.tikz", standAlone=F, width = 6, height = 7.5)
ggplot(data_l %>% filter(covType == 'Plugin')) +
  geom_histogram_pattern(aes(x = pvalue, y = after_stat(density)*binwidth, fill = SNR),
                         pattern = rep(c(rep('none', 20 * length(n)),
                                         rep(c(rep('none', 20),
                                               rep('stripe', 20*length(SNRS))),
                                             length(n))), 2),
                         pattern_angle = rep(c(rep(0, 20 * length(n)),
                                               rep(rep(seq(0, 180, length = length(SNRS)+2)[-(length(SNRS)+2)], each = 20),
                                                   length(n))), 2) - 90,
                         pattern_size = 0.1, pattern_colour = 'grey80', pattern_spacing = 0.015,
                         breaks = breaks, position = position_dodge()) + theme_bw()+
  ggtitle('p-value histogram from LLR test') +
  facet_grid(vars(spaceType), vars(SampleSize)) +
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
  labs(fill = 'SNR:')
dev.off()

p1 = ggplot(data_l %>% filter(covType == 'Plugin') %>%
              mutate(spaceType = recode_factor(spaceType,
                                               `Oracle eqref{eqn:polyP}` = 'Oracle (B.5)',
                                               `Plugin eqref{eqn:polyP}` = 'Plugin (B.5)',
                                               `Oracle eqref{eqn:eigvP}` = 'Oracle (B.6)',
                                               `Plugin eqref{eqn:eigvP}` = 'Plugin (B.6)'))) +
  geom_histogram_pattern(aes(x = pvalue, y = after_stat(density)*binwidth, fill = SNR),
                         pattern = rep(c(rep('none', 20 * length(n)),
                                         rep(c(rep('none', 20),
                                               rep('stripe', 20*length(SNRS))),
                                             length(n))), 2),
                         pattern_angle = rep(c(rep(0, 20 * length(n)),
                                               rep(rep(seq(0, 180, length = length(SNRS)+2)[-(length(SNRS)+2)], each = 20),
                                                   length(n))), 2) - 90,
                         pattern_size = 0.1, pattern_colour = 'grey80', pattern_spacing = 0.015,
                         breaks = breaks, position = position_dodge()) + theme_bw()+
  ggtitle('p-value histogram from LLR test') +
  facet_grid(vars(spaceType), vars(SampleSize), labeller = as_labeller(TeX, default = label_parsed)) +
  scale_y_continuous(breaks = c(0.05, 0.25, 0.5, 0.75, 1), trans = 'sqrt')+
  scale_fill_discrete(labels = TeX(levels(data_l$SNR))) +
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
p1
ggsave(filename = 'output/Plots/png/PvalueLLR.png', p1, width = 6, height = 7.5, units = 'in')


###### Type I/II errors ######
data_l %>% filter(covType == 'Plugin') %>%
  group_by(spaceType, SNR, SampleSize) %>%
  summarise(RejRate = mean(pvalue <= 0.05)) %>% print(n = nrow(.))

