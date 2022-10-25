library(eigTest)
library(tidyverse)
#library(tikzDevice)

data(hudsonDaily); data(hudsonWeekly)

hudsonPlot = rbind(cbind(hudsonDaily, reso = 'Daily discharge (2015 - 2020)'),
                   cbind(hudsonWeekly, reso = 'Weekly discharge (1979 - 2014)'))

p1 = ggplot(hudsonPlot, aes(x = datetime)) + facet_wrap(~reso, scales = 'free_x', ncol = 1) +
  geom_ribbon(aes(ymin = rep(0, length(datetime)), ymax = p25_va, fill = 'Drought'), alpha = 0.3) +
  geom_ribbon(aes(ymax = rep(Inf, length(datetime)), ymin = p75_va, fill = 'Flooding'), alpha = 0.3) +
  geom_ribbon(aes(ymin=p25_va, ymax=p75_va, fill = 'Normal'), alpha=0.5) +
  xlab('Date') +
  #ggtitle('Daily discharge (2015 - 2020)')+
  labs(color = 'State', fill = 'State') +
  geom_path(aes(y = Discharge), size = 0.3) +
  theme(legend.position = 'bottom', legend.title = element_blank(), strip.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5), strip.background = element_blank()) +
  scale_fill_discrete(breaks = c('Drought', 'Normal', 'Flooding'),
                      labels = c('Drought (0 - 25)', 'Normal (25 - 75)', 'Flooding (75 - 100)')) +
  geom_point(aes(y = as.double(Level)*15000 - 15000, color = as.character(Level)), size = 1, shape = 4) +
  scale_y_continuous('Discharge ($\\mbox{ft}^3/s$)',
                     sec.axis = sec_axis(~., breaks = (1:3)*15000 - 15000, labels = c('Drought', 'Normal', 'Flooding'))) +
  scale_color_discrete(breaks = c('Drought', 'Normal', 'Flooding'),
                       labels = c('Drought (0 - 25)', 'Normal (25 - 75)', 'Flooding (75 - 100)'))

#tikz(file = "implementations/Plots/Streamflow.tikz", standAlone=F, width = 6, height = 4.5)
p1
#dev.off()

labelTab = cbind(hudsonDaily$Level, hudsonWeekly$Level)

L = nrow(hudsonWeekly)

p = 2
matTran = array(0, dim = c(p, 3, 3))
matCov = array(0, dim = c(p, 9, 9))
dimnames(matTran) = list(c('Daily', 'Weekly'), c('Drought', 'Normal', 'Flooding'), c('Drought', 'Normal', 'Flooding'))
dimnames(matCov) = list(c('Daily', 'Weekly'), NULL, NULL)

i0 = as.integer(labelTab[1,1:p])
for (i in 2:L) {
  i1 = as.integer(labelTab[i,1:p])
  for (j in 1:p) {
    matTran[j, i0[j], i1[j]] = matTran[j, i0[j], i1[j]] + 1
  }
  i0 = i1
}

for (j in 1:p) {
  matTran[j,,] = matTran[j,,]/rowSums(matTran[j,,])
  Label = as.integer(unlist(labelTab[,j]))
  for (i in 1:3) {
    pi = sum(Label == i)/L
    Ei = matrix(0, ncol = 3, nrow = 3)
    Ei[i,i] = 1
    Qi = (diag(matTran[j,i,]) - matrix(matTran[j,i,], nrow = 3, ncol = 3) *
            matrix(matTran[j,i,], nrow = 3, ncol = 3, byrow = T))/pi
    matCov[j,,] = matCov[j,,] + kronecker(Qi, Ei)
  }
}

dailyTran = as_tibble(matTran[1,,]) %>%
  mutate(From = colnames(matTran[1,,])) %>%
  pivot_longer(cols = 1:3, names_to = 'To') %>%
  mutate(Resolution = 'Daily')

weeklyTran = as_tibble(matTran[2,,]) %>%
  mutate(From = colnames(matTran[2,,])) %>%
  pivot_longer(cols = 1:3, names_to = 'To') %>%
  mutate(Resolution = 'Weekly')

tranMats = rbind(dailyTran, weeklyTran)

#tikz(file = "implementations/Plots/tranMat.tikz", standAlone=F, width =7.5, height = 3)
ggplot(tranMats, aes(x = To, y = From, fill = value)) +
  geom_tile(color = 'white')+
  scale_fill_gradient2(low = "white", high = "black", limit = c(0,1), space = "Lab",
                       name="probability:") +
  theme_minimal()+ facet_wrap(~Resolution, ncol = 2) +
  geom_text(aes(label = sprintf("%0.3f", round(value, digits = 3)),
                color = value > 0.5), size = 4)+
  ggtitle('Transition probability matrices') +
  scale_color_manual(guide = 'none', values = c("black", "white")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(angle = 45, vjust = 3.5, hjust = 0.1),
        axis.text.x = element_text(),
        axis.title.y = element_text(vjust = -3.5),
        legend.position = 'bottom', strip.text = element_text(),
        #panel.grid.major = element_blank(),
        panel.background = element_blank(),
        legend.box.margin = margin(0,0,0,0), legend.margin = margin(0,0,0,0),
        panel.border = element_blank())+
  coord_fixed() +
  scale_y_discrete(limits = colnames(matTran[1,,])) +
  scale_x_discrete(limits = colnames(matTran[2,,])) +
  guides(fill = guide_colorbar(title.vjust = 0.8, barwidth = 9, barheight = 1))
#dev.off()

partialTest(matTran, cn = sqrt(L-1), cov.arr = matCov, nn = T, k = 1,
            testType = 'chi', eps = (L-1)^(-1/3), param.out = T)
partialTest(matTran, cn = sqrt(L-1), cov.arr = matCov, nn = T, k = 1,
            testType = 'gam', param.out = T)

statDist = nnPartSchur(matTran)[,1] %>%
  abs %>% proportions()
names(statDist) = colnames(matTran[1,,])
statDist
