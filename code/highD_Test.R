library(eigTest)
library(foreach)
library(doSNOW)
library(tidyverse)
library(reshape2)
#library(tikzDevice)


###### generate / load pvalues ######
# simulate pvalues from scratch
simu_pval = FALSE

d = 2:20
p = 2
SNRS = c(4, 0.25)
samples = 200
n = sqrt(c(50, 100, 500))

if (simu_pval) {
  numCores = parallel::detectCores()
  totL = length(n)*samples*(length(SNRS)+1)
  pb <- txtProgressBar(max = totL, style = 3)
  progress <- function(n) {setTxtProgressBar(pb, n)}
  opts <- list(progress = progress)
  data_highD = c()

  cl <- makeCluster(numCores)
  registerDoSNOW(cl)

  for (i in 1:length(d)) {
    set.seed(7202)
    if (i == 1) {cat(paste('Dimension d =', d[i], ': \n'))}
    else {cat(paste('\nDimension d =', d[i], ': \n'))}

    means = generateMeans(d[i], p, snr = SNRS, control.g = T)
    simulated = simuSamples(means, n, samples, prl = T)
    gc()
    cat('\nTesting: \n')

    eigvORC = JDTE(means[,1,,], iter = 1000)

    data_i = foreach(est_list = simulated, m =icount(), .inorder = F,
                     .combine = bind_rows, .options.snow = opts,
                     .packages = c('eigTest', 'reshape2', 'tidyverse')) %dopar% {
                       if (m %% max(100, totL %/% 50) == 1) {gc()}
                       #                       print(m)
                       #                       profvis({
                       mu.bar = est_list$mu.bar; cov.bar = est_list$cov.bar
                       SNR = est_list$SNR; CovRate = est_list$CovRate

                       data_temp = commutatorTest(mu.bar, cn = CovRate, param.out = T)[[1]] %>%
                         list2DF() %>% select(df, pvalue) %>% mutate(covType = 'Oracle')
                       data_temp = commutatorTest(mu.bar, cn = CovRate, cov.arr = cov.bar,
                                                  param.out = T)[[1]] %>%
                         list2DF() %>% select(df, pvalue) %>% mutate(covType = 'Plugin') %>%
                         bind_rows(data_temp) %>%
                         mutate(SNR = SNR, Dimension = d[i], SampleSize = round(CovRate^2), testType = 'COM')

                       eigvPLG = JDTE(mu.bar)

                       data_LLR = projTest(mu.bar, cn = CovRate, param.out = T,
                                           CV = eigvPLG, poly.sp = F) %>%
                         list2DF() %>% select(df, pvalue) %>%
                         mutate(covType = 'Oracle', spaceType = 'eigv-PLG')
                       data_LLR = projTest(mu.bar, cn = CovRate, cov.arr = cov.bar,
                                           poly.sp = F, CV = eigvPLG, param.out = T) %>%
                         list2DF() %>% select(df, pvalue) %>%
                         mutate(covType = 'Plugin', spaceType = 'eigv-PLG') %>%
                         bind_rows(data_LLR)
                       data_LLR = projTest(mu.bar, cn = CovRate, param.out = T) %>% list2DF() %>%
                         select(df, pvalue) %>%
                         mutate(covType = 'Oracle', spaceType = 'poly-PLG') %>%
                         bind_rows(data_LLR)
                       data_LLR = projTest(mu.bar, cn = CovRate, cov.arr = cov.bar,
                                           param.out = T) %>%
                         list2DF() %>% select(df, pvalue) %>%
                         mutate(covType = 'Plugin', spaceType = 'poly-PLG') %>%
                         bind_rows(data_LLR)

                       if (SNR == '1/SNR=0') {
                         data_LLR = projTest(mu.bar, refMat = means[,1,,],
                                             cn = CovRate, param.out = T) %>%
                           list2DF() %>% select(df, pvalue) %>%
                           mutate(covType = 'Oracle', spaceType = 'poly-ORC') %>%
                           bind_rows(data_LLR)
                         data_LLR = projTest(mu.bar, cn = CovRate, param.out = T,
                                             CV = eigvORC, poly.sp = F) %>%
                           list2DF() %>% select(df, pvalue) %>%
                           mutate(covType = 'Oracle', spaceType = 'eigv-ORC') %>%
                           bind_rows(data_LLR)
                         data_LLR = projTest(mu.bar, cn = CovRate, refMat = means[,1,,],
                                             cov.arr = cov.bar, param.out = T) %>%
                           list2DF() %>% select(df, pvalue) %>%
                           mutate(covType = 'Plugin', spaceType = 'poly-ORC') %>%
                           bind_rows(data_LLR)
                         data_LLR = projTest(mu.bar, cn = CovRate, cov.arr = cov.bar,
                                             CV = eigvORC, poly.sp = F, param.out = T) %>%
                           list2DF() %>% select(df, pvalue) %>%
                           mutate(covType = 'Plugin', spaceType = 'eigv-ORC') %>%
                           bind_rows(data_LLR)
                       }
                       #                       })
                       return(data_LLR %>%
                                mutate(SNR = SNR, Dimension = d[i],
                                       SampleSize = round(CovRate^2),
                                       testType = 'LLR') %>%
                                bind_rows(data_temp)
                       )
                     }
    rm(simulated); gc()
    data_highD = rbind(data_highD, data_i)
  }
  stopCluster(cl)
  save(data_highD, file = 'output/highD_test.RData')
} else {
  data_highD = load('output/highD_Test.RData')
}


