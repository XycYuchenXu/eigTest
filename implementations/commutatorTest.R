#devtools::install_github('XycYuchenXu/eigTest', force = T, build_vignettes = F)
library(eigTest)
library(foreach)
library(doSNOW)
library(tidyverse)
#library(tikzDevice)

set.seed(7202)
d = 5
p = 2
k = d
samples = 500
SNRS = c(50, 1)
n = sqrt(c(50, 100, 250))
means_c = generateMeans(d, p, k, snr = SNRS, control.g = TRUE)
simulated_c = simuSamples(means_c, n, samples)

numCores = parallel::detectCores()/2
totL = length(n)*samples*(length(SNRS)+1)
pb <- txtProgressBar(max = totL, style = 3)
progress <- function(n) {setTxtProgressBar(pb, n)}
opts <- list(progress = progress)

cl <- makeCluster(numCores)
registerDoSNOW(cl)
data_c = foreach(est_list = simulated_c, .inorder = F, .combine = bind_rows,
                 .options.snow = opts, .packages = c('eigTest', 'tidyverse')) %dopar% {
                   mu.bar = est_list$mu.bar; cov.bar = est_list$cov.bar
                   SNR = est_list$SNR; CovRate = est_list$CovRate

                   data_temp = commutatorTest(mu.bar, cn = CovRate, cov.arr = cov.bar,
                                              param.out = T)[[1]] %>%
                     list2DF() %>% select(pvalue) %>%
                     mutate(SNR = SNR, SampleSize = round(CovRate^2))

                   return(data_temp)
                 }
stopCluster(cl)

save(data_c, file = 'implementations/output/commutatorTest.RData')

binwidth = 0.05
breaks = seq(0, 1, 0.05)

#tikz(file = "implementations/Plots/PvalueCommutator.tikz", standAlone=F,width = 6, height = 3)
ggplot(data_c) +
  geom_histogram(aes(x = pvalue, y = ..density..*binwidth, fill = SNR),
                 breaks = breaks,
                 position = position_dodge()) + theme_bw()+
  ggtitle('P-value histogram from Commutator-based test') +
  facet_wrap(~SampleSize, ncol = length(n), dir = 'v') +scale_y_continuous(breaks = c(0.05, 0.25, 0.5, 0.75, 1), trans = 'sqrt')+
  theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(size = 1), aspect.ratio = 0.8,
        legend.spacing.y = unit(0.1, 'cm'), legend.key.height = unit(0.5, 'cm'),
        strip.background =element_rect(fill="white", size = 0)) +
  guides(fill=guide_legend(ncol=3)) + xlab('P values') + ylab('Proportion') +
  labs(fill = 'SNR:') #scale_fill_discrete(labels = c('Chi', 'Gam'))+
#dev.off()
