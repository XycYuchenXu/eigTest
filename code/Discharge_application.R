library(eigTest)
library(tidyverse)
library(tikzDevice)
library(ggpattern)
library(latex2exp)

###### data summary & plot ######
data(hudsonDaily); data(hudsonWeekly)

hudsonPlot = rbind(cbind(hudsonDaily, reso = 'Daily discharge (2015 - 2020)'),
                   cbind(hudsonWeekly, reso = 'Weekly discharge (1979 - 2014)'))

tikz(file = "output/Plots/tikz/Streamflow.tikz", standAlone=F, width = 7, height = 4.5)
ggplot(hudsonPlot, aes(x = datetime)) + facet_wrap(~reso, scales = 'free_x', ncol = 1) +
  geom_ribbon_pattern(aes(ymin = rep(0, length(datetime)),
                          ymax = p25_va, fill = 'Drought'), alpha = 0.3,
                      pattern = 'stripe', pattern_angle = 45, pattern_alpha = 0.1,
                      pattern_size = 0.1, pattern_fill = 'black', pattern_spacing = 0.05) +
  geom_ribbon_pattern(aes(ymax = rep(Inf, length(datetime)),
                          ymin = p75_va, fill = 'Flooding'), alpha = 0.3,
                      pattern = 'stripe', pattern_angle = 135, pattern_alpha = 0.1,
                      pattern_size = 0.1, pattern_fill = 'black', pattern_spacing = 0.05) +
  geom_ribbon(aes(ymin=p25_va, ymax=p75_va, fill = 'Normal'), alpha=0.7) +
  xlab('Date') +
  labs(color = 'State', fill = 'State') +
  guides(fill=guide_legend(override.aes = list(
    alpha = c(0.3, 0.7, 0.3),
    pattern = c('stripe', 'none', 'stripe'),
    pattern_angle = c(45, 0, 135), pattern_alpha = 0.1,
    pattern_size = 0.1, pattern_fill = 'black', pattern_spacing = 0.015
  ))) + geom_path(aes(y = Discharge), linewidth = 0.3) +
  theme(legend.position = 'bottom', legend.title = element_blank(),
        strip.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(margin = margin(r = 5, unit = 'mm')),
        strip.background = element_blank()) +
  scale_fill_discrete(breaks = c('Drought', 'Normal', 'Flooding'),
                      labels = c('Drought (0 - 25)', 'Normal (25 - 75)', 'Flooding (75 - 100)')) +
  geom_point(aes(y = as.double(Level)*15000 - 15000, color = as.character(Level)), size = 1, shape = 4) +
  scale_y_continuous('Discharge ($\\mbox{ft}^3/s$)',
                     sec.axis = sec_axis(~., breaks = (1:3)*15000 - 15000, labels = c('Drought', 'Normal', 'Flooding'))) +
  scale_color_discrete(breaks = c('Drought', 'Normal', 'Flooding'),
                       labels = c('Drought (0 - 25)', 'Normal (25 - 75)', 'Flooding (75 - 100)'))
dev.off()

p1 = ggplot(hudsonPlot, aes(x = datetime)) + facet_wrap(~reso, scales = 'free_x', ncol = 1) +
  geom_ribbon_pattern(aes(ymin = rep(0, length(datetime)),
                          ymax = p25_va, fill = 'Drought'), alpha = 0.3,
                      pattern = 'stripe', pattern_angle = 45, pattern_alpha = 0.1,
                      pattern_size = 0.05, pattern_fill = 'black', pattern_spacing = 0.05) +
  geom_ribbon_pattern(aes(ymax = rep(Inf, length(datetime)),
                          ymin = p75_va, fill = 'Flooding'), alpha = 0.3,
                      pattern = 'stripe', pattern_angle = 135, pattern_alpha = 0.1,
                      pattern_size = 0.05, pattern_fill = 'black', pattern_spacing = 0.05) +
  geom_ribbon(aes(ymin=p25_va, ymax=p75_va, fill = 'Normal'), alpha=0.7) +
  xlab('Date') +
  labs(color = 'State', fill = 'State') +
  guides(fill=guide_legend(override.aes = list(
    alpha = c(0.3, 0.7, 0.3),
    pattern = c('stripe', 'none', 'stripe'),
    pattern_angle = c(45, 0, 135), pattern_alpha = 0.1,
    pattern_size = 0.1, pattern_fill = 'black', pattern_spacing = 0.015
  ))) + geom_path(aes(y = Discharge), linewidth = 0.3) +
  theme(legend.position = 'bottom', legend.title = element_blank(),
        strip.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(margin = margin(r = 5, unit = 'mm')),
        strip.background = element_blank()) +
  scale_fill_discrete(breaks = c('Drought', 'Normal', 'Flooding'),
                      labels = c('Drought (0 - 25)', 'Normal (25 - 75)', 'Flooding (75 - 100)')) +
  geom_point(aes(y = as.double(Level)*15000 - 15000, color = as.character(Level)), size = 1, shape = 4) +
  scale_y_continuous(TeX("Discharge ($ft^3/s$)"),
                     sec.axis = sec_axis(~., breaks = (1:3)*15000 - 15000, labels = c('Drought', 'Normal', 'Flooding'))) +
  scale_color_discrete(breaks = c('Drought', 'Normal', 'Flooding'),
                       labels = c('Drought (0 - 25)', 'Normal (25 - 75)', 'Flooding (75 - 100)'))
p1
ggsave(filename = 'output/Plots/png/Streamflow.png', p1, width = 7, height = 4.5, units = 'in')


###### transition probabilities ######
labelTab = cbind(hudsonDaily$Level, hudsonWeekly$Level)
L = nrow(hudsonWeekly)

p = 2; d = 3
matTran = array(0, dim = c(p, d, d))
matCov = array(0, dim = c(p, d^2, d^2))
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
  for (i in 1:d) {
    pi = sum(Label == i)/L
    Ei = matrix(0, ncol = d, nrow = d)
    Ei[i,i] = 1
    Qi = (diag(matTran[j,i,]) - matrix(matTran[j,i,], nrow = d, ncol = d) *
            matrix(matTran[j,i,], nrow = d, ncol = d, byrow = T))/pi
    matCov[j,,] = matCov[j,,] + kronecker(Qi, Ei)
  }
}

dailyTran = as_tibble(matTran[1,,]) %>%
  mutate(From = colnames(matTran[1,,])) %>%
  pivot_longer(cols = 1:d, names_to = 'To') %>%
  mutate(Resolution = 'Daily')

weeklyTran = as_tibble(matTran[2,,]) %>%
  mutate(From = colnames(matTran[2,,])) %>%
  pivot_longer(cols = 1:d, names_to = 'To') %>%
  mutate(Resolution = 'Weekly')

tranMats = rbind(dailyTran, weeklyTran)

tikz(file = "output/Plots/tikz/tranMat.tikz", standAlone=F, width =7.5, height = 4)
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
        legend.title = element_text(margin = margin(r = 3, unit = 'mm')),
        panel.background = element_blank(),
        legend.box.margin = margin(0,0,0,0), legend.margin = margin(0,0,0,0),
        panel.border = element_blank())+
  coord_fixed() +
  scale_y_discrete(limits = colnames(matTran[1,,])) +
  scale_x_discrete(limits = colnames(matTran[2,,])) +
  guides(fill = guide_colorbar(title.vjust = 0.8, barwidth = 9, barheight = 1))
dev.off()

p2 = ggplot(tranMats, aes(x = To, y = From, fill = value)) +
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
        legend.title = element_text(margin = margin(r = 3, unit = 'mm')),
        panel.background = element_blank(),
        legend.box.margin = margin(0,0,0,0), legend.margin = margin(0,0,0,0),
        panel.border = element_blank())+
  coord_fixed() +
  scale_y_discrete(limits = colnames(matTran[1,,])) +
  scale_x_discrete(limits = colnames(matTran[2,,])) +
  guides(fill = guide_colorbar(title.vjust = 0.8, barwidth = 9, barheight = 1))
p2
ggsave(filename = 'output/Plots/png/tranMat.png', p2, width = 7.5, height = 4, units = 'in')


###### static distribution estimate / test ######
# transpose for test setup
matTran_t = matTran
for (i in 1:p) {
  matTran_t[i,,] = t(matTran[i,,])
}
matCov_t = matCov[,
                  as.vector(matrix(1:d^2, ncol = d, byrow = T)),
                  as.vector(matrix(1:d^2, ncol = d, byrow = T))]

partialTest(matTran_t, cn = sqrt(L-1), cov.arr = matCov_t, nn = T, k = 1,
            testType = 'chi', eps = (L-1)^(-1/3), param.out = T)
partialTest(matTran_t, cn = sqrt(L-1), cov.arr = matCov_t, nn = T, k = 1,
            testType = 'gam', param.out = T)

statDist = nnPartSchur(matTran_t)[,1] %>%
  abs %>% proportions()
names(statDist) = colnames(matTran_t[1,,])
statDist
