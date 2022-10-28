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
                     SNR = est_list$SNR; CovRate = est_list$CovRate

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
                     colnames(data_p)[1:2] = c('TestType', 'K')

                     return(data_p %>%
                              mutate(SNR = SNR, SampleSize = round(CovRate^2))
                     )
                   }
  stopCluster(cl)
  save(data_p, file = 'output/partTest.RData')
} else {
  data_p = load('output/partTest.RData')
}


###### plot histograms ######
breaks = seq(0,1,0.05); binwidth = 0.05

#tikz(file = "output/Plots/PvaluePartial.tikz", standAlone=F, width = 6, height = 4.5)
ggplot(data_p %>% filter(K == 'K = 2')) +
  geom_histogram(aes(x = pvalue, y = ..density..*binwidth, fill = SNR),
                 breaks = breaks,
                 position = position_dodge()) + theme_bw()+
  ggtitle('P-value histogram from partial test') +
  facet_grid(vars(TestType), vars(SampleSize), #ncol = length(n), dir = 'v',
             labeller = labeller(TestType = c(`Chi` = 'Chi test with $\\widehat{V}$',
                                              `Gam` = 'Gamma test with $\\widehat{V}$'))) +
  theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5),
        strip.background =element_rect(size = 0),
        panel.border = element_rect(size = 1), strip.placement = 'inside') +
  guides(fill=guide_legend(ncol=4)) + xlab('P values') + ylab('Proportion') +
  labs(fill = 'SNR:') +
  scale_fill_discrete()+
  scale_y_continuous(breaks = c(0.05, 0.25, 0.5, 0.75, 1), trans = 'sqrt')
#dev.off()

