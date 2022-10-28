#devtools::install_github('XycYuchenXu/eigTest', force = T, build_vignettes = F)
library(eigTest)
library(foreach)
library(doSNOW)
library(tidyverse)
library(tikzDevice)


###### generate / load pvalues ######
# simulate pvalues from scratch
simu_pval = FALSE

samples = 200
d = 4
p = 8
k = d
SNRS = c(1000, 100, 10)
n = sqrt(c(100, 1000, 10000, 100000))

if (simu_pval) {
  set.seed(2020)
  means_m = generateMeans(d, p, k, snr = SNRS, control.g = TRUE)
  v.t = eigen(means_m[1,1,,])$vectors

  numCores = parallel::detectCores()/2
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
                     SNR = est_list$SNR; CovRate = est_list$CovRate

                     eigvPLG = JDTE(mu.bar)

                     data_eig = eigTest(mu.bar, cn = CovRate, cov.arr = cov.bar, V = eigvPLG,
                                        testType = 'chi', param.out = T)$chi %>%
                       list2DF() %>% select(pvalue) %>%
                       mutate(testType = 'Chi')
                     data_eig = eigTest(mu.bar, cn = CovRate, cov.arr = cov.bar, V = eigvPLG,
                                        testType = 'gam', param.out = T)$gam %>%
                       list2DF() %>% select(pvalue) %>%
                       mutate(testType = 'Gam') %>% bind_rows(data_eig)


                     if (SNR == '1/SNR=0') {
                       data_eig = eigTest(mu.bar, cn = CovRate, cov.arr = cov.bar,
                                          testType = 'chi', V = v.t, param.out = T)$chi %>%
                         list2DF() %>% select(pvalue) %>%
                         mutate(testType = 'Chi_0') %>% bind_rows(data_eig)
                     }

                     return(data_eig %>%
                              mutate(SNR = SNR, SampleSize = round(CovRate^2))
                     )
                   }

  stopCluster(cl)

  data_m$testType = factor(data_m$testType, levels = c('Chi_0', 'Chi', 'Gam'))
  data_m$SampleSize = paste0('Sample size $n = 10^', round(log10(data_m$SampleSize)), '$')
  save(data_m, file = 'output/multiTest.RData')
} else {
  data_m = load('output/multiTest.RData')
}


###### plot histograms ######
binwidth = 0.05
breaks = seq(0, 1, 0.05)

#tikz(file = "output/Plots/PvalueMulti.tikz", standAlone=F,width = 7, height = 6)
ggplot(data_m) +
  geom_histogram(aes(x = pvalue, y = ..density..*binwidth, fill = SNR),
                 breaks = breaks,
                 position = position_dodge()) + theme_bw()+
  ggtitle('P-value histogram from multi-sample test') +
  facet_grid(vars(testType), vars(SampleSize),# ncol = 4, dir = 'v',
             labeller = labeller(testType = c(`Chi_0` = 'Chi test with $V$',
                                              `Chi` = 'Chi test with $\\widehat{V}$',
                                              `Gam` = 'Gamma test with $\\widehat{V}$'))) +
  theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5),
        strip.background =element_rect(size = 0),# fill="white"),
        panel.border = element_rect(size = 1), strip.placement = 'inside') +
  guides(fill=guide_legend(ncol=4)) + xlab('P values') + ylab('Proportion') +
  labs(fill = 'SNR:') +
  scale_fill_discrete()+
  scale_y_continuous(breaks = c(0.05, 0.25, 0.5, 0.75, 1), trans = 'sqrt')
#dev.off()
