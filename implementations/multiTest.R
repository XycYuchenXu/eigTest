#devtools::install_github('XycYuchenXu/eigTest', force = T, build_vignettes = F)
library(eigTest)
library(foreach)
library(doSNOW)
library(tidyverse)
#library(tikzDevice)

set.seed(2020)
samples = 200
d = 4
p = 8
k = d
SNRS = c(1000, 100, 10)
n = sqrt(c(100, 1000, 10000, 100000))
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

                   data_eig = eigTest(mu.bar, cn = CovRate, cov.arr = cov.bar,
                                      testType = 'chi', param.out = T) %>%
                     list2DF() %>% select(pvalue) %>%
                     mutate(testType = 'chi')
                   data_eig = eigTest(mu.bar, cn = CovRate, cov.arr = cov.bar,
                                      testType = 'gam', param.out = T) %>%
                     list2DF() %>% select(pvalue) %>%
                     mutate(testType = 'gam') %>% bind_rows(data_eig)


                   if (SNR == '1/SNR=0') {
                     data_eig = eigTest(mu.bar, cn = CovRate, cov.arr = cov.bar,
                                        testType = 'chi', V = v.t, param.out = T) %>%
                       list2DF() %>% select(pvalue) %>%
                       mutate(testType = 'chi_0') %>% bind_rows(data_LLR)
                   }

                   return(data_eig %>%
                            mutate(SNR = SNR, SampleSize = round(CovRate^2))
                   )
                 }

stopCluster(cl)

data_m$TestType = factor(dataP$TestType, levels = c('Chi_0', 'Chi', 'Gam'))
save(data_m, file = 'implementations/output/multiTest.RData')

binwidth = 0.05
breaks = seq(0, 1, 0.05)

#tikz(file = "implementations/Plots/PvalueMulti.tikz", standAlone=F,width = 7, height = 6)
ggplot(dataP) +
  geom_histogram(aes(x = P_value, y = ..density..*binwidth, fill = SNR),
                 breaks = breaks,
                 position = position_dodge()) + theme_bw()+
  ggtitle('P-value histogram from multi-sample test') +
  facet_grid(vars(TestType), vars(SampleSize),# ncol = 4, dir = 'v',
             labeller = labeller(TestType = c(`Chi_0` = 'Chi test with $V$',
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
