#devtools::install_github('XycYuchenXu/eigTest', force = T, build_vignettes = F)
library(eigTest)
library(plotrix)
library(reshape2)
library(ggplot2)
library(tikzDevice)


###### time series plot ######
c2p = as_tibble(as.matrix(countryMacro[[1]])) %>% mutate(country = names(countryMacro)[1], Quarters = time(countryMacro[[1]]))
for (i in 2:length(countryMacro)) {
  c2p = as_tibble(as.matrix(countryMacro[[i]])) %>% mutate(country = names(countryMacro[i]), Quarters = time(countryMacro[[i]])) %>%
    bind_rows(c2p)
}
names(c2p)[1:3] = c('GDP', 'M2', 'REER')
c2p$GDP = c2p$GDP / 10^6
c2p$M2 = c2p$M2 / 10^7

labNames = c('GDP (USD, $\\times 10^{12}$)', 'M2 (USD, $\\times 10^{13}$)', 'REER ($2015=100$)')
names(labNames) = names(c2p)[1:3]

c2p = pivot_longer(c2p, cols = c(GDP, M2, REER), names_to = 'key', values_to = 'value')
c2p$country = factor(c2p$country, levels = rev(names(countryMacro)))

tikz('output/Plots/tikz/ts.tikz', standAlone = F, width = 6, height = 3)
ggplot(data = c2p) +
  geom_path(aes(x = Quarters, y = value, color = country, linetype = country), size = 1.5) +
  facet_wrap(~key, scales = 'free_y', labeller = labeller(key = labNames)) +
  labs(color = 'Country:', linetype = 'Country:') + theme_bw() + ylab('') +
  ggtitle('Quarterly macroeconomic indices') +
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_blank(), legend.position = 'bottom',
        legend.title = element_text(margin = margin(r = 3, unit = 'mm')),
        legend.key.width = unit(1.5, 'cm'), panel.border = element_rect(size = 1),
        legend.spacing.y = unit(0.1, 'cm'), legend.key.height = unit(0.5, 'cm'),
        strip.background =element_rect(fill="white", size = 0)) +
  guides(col=guide_legend(ncol=4))
dev.off()

p1 = ggplot(data = c2p) +
  geom_path(aes(x = Quarters, y = value, color = country, linetype = country), size = 1.5) +
  facet_wrap(~key, scales = 'free_y', labeller = labeller(key = TeX(labNames))) +
  labs(color = 'Country:', linetype = 'Country:') + theme_bw() + ylab('') +
  ggtitle('Quarterly macroeconomic indices') +
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_blank(), legend.position = 'bottom',
        legend.title = element_text(margin = margin(r = 3, unit = 'mm')),
        legend.key.width = unit(1.5, 'cm'), panel.border = element_rect(size = 1),
        legend.spacing.y = unit(0.1, 'cm'), legend.key.height = unit(0.5, 'cm'),
        strip.background =element_rect(fill="white", size = 0)) +
  guides(col=guide_legend(ncol=4))
p1
ggsave(filename = 'output/Plots/png/ts.png', p1, width = 6, height = 3, units = 'in')


###### data preprocessing ######
m = nrow(countryMacro[[1]]) - 1
d = ncol(countryMacro[[1]])
p = length(countryMacro)

# transpose for test setup
countryCoeff_t = countryCoeff
for (i in 1:p) {
  countryCoeff_t[i,,] = t(countryCoeff[i,,])
}
countryCovar_t = countryCovar[,
                              as.vector(matrix(1:d^2, ncol = d, byrow = T)),
                              as.vector(matrix(1:d^2, ncol = d, byrow = T))]

###### global test & estimation ######
# multi-sample test, full-rank, Corollary 4.1
eigTest(countryCoeff_t, cn = sqrt(m), cov.arr = countryCovar_t, testType = 'chi', param.out = T)
eigTest(countryCoeff_t, cn = sqrt(m), cov.arr = countryCovar_t, testType = 'gam', param.out = T)

# partial test, k = 1, Proposition 5.2 & Corollary 5.1
partialTest(countryCoeff_t, cn = sqrt(m), k = 1, cov.arr = countryCovar_t, testType = 'chi', param.out = T)
partialTest(countryCoeff_t, cn = sqrt(m), k = 1, cov.arr = countryCovar_t, testType = 'gam', param.out = T)

# partial test, k = 2, Corollary 5.1
partialTest(countryCoeff_t, cn = sqrt(m), k = 2, cov.arr = countryCovar_t, testType = 'gam', param.out = T)

# estimated (partially) common eigenvectors
V = JDTE(countryCoeff_t)
V

Q = expmPartSchur(countryCoeff_t, k = 2, warmup = T)
B = array(0, dim = c(p, d, d))
for (i in 1:p) {
  B[i,,] = tcrossprod(crossprod(Q, countryCoeff_t[i,,]), t(Q))
}
V = JDTE(B[,1:2, 1:2])
Q[,1:2] %*% V


###### continent-grouped test ######
# Asia, Corollary 4.1
eigTest(countryCoeff_t[1:3,,], cn = sqrt(m), cov.arr = countryCovar_t[1:3,,], testType = 'chi', param.out = T)
eigTest(countryCoeff_t[1:3,,], cn = sqrt(m), cov.arr = countryCovar_t[1:3,,], testType = 'gam', param.out = T)

# Europe, Corollary 4.1
eigTest(countryCoeff_t[4:6,,], cn = sqrt(m), cov.arr = countryCovar_t[4:6,,], testType = 'chi', param.out = T)
eigTest(countryCoeff_t[4:6,,], cn = sqrt(m), cov.arr = countryCovar_t[4:6,,], testType = 'gam', param.out = T)

# North America, Corollary 4.1
eigTest(countryCoeff_t[7:8,,], cn = sqrt(m), cov.arr = countryCovar_t[7:8,,], testType = 'chi', param.out = T)
eigTest(countryCoeff_t[7:8,,], cn = sqrt(m), cov.arr = countryCovar_t[7:8,,], testType = 'gam', param.out = T)


###### pairwise commutator test ######
countries = names(countryMacro)
comm.pair.test = 0.5*diag(length(countries))
for (i in 1:length(countries)) {
  for (j in i:length(countries)) {
    if (i == j) {next()}
    comm.pair.test[i,j] = commutatorTest(countryCoeff[c(i,j),,],
                                         cn = sqrt(m),
                                         cov.arr = countryCovar[c(i,j),,])
  }
}
colnames(comm.pair.test) = countries
rownames(comm.pair.test) = countries
comm.pair.test

comm.pair.complete = comm.pair.test + t(comm.pair.test)

pair.test2plot = melt(comm.pair.complete, na.rm = TRUE)
pair.test2plot$value = round(pair.test2plot$value, 3)

text2plot = pair.test2plot

groupMat = matrix(NA, nrow = length(countries), ncol = length(countries))
colnames(groupMat) = countries
rownames(groupMat) = countries
groupMat[1:2,1:2] = 'NA'
groupMat[3:5,3:5] = 'EU'
groupMat[6:8,6:8] = 'NA'
groupMat = melt(groupMat, na.rm = T)
groupMat$value = as.factor(groupMat$value)

tikz('output/Plots/tikz/pvalMat.tikz', standAlone = F, width = 5, height = 4.5)
ggplot(data = pair.test2plot, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "black",
                       limit = c(0,1), space = "Lab",
                       name="p-value:") +
  theme_minimal()+
  geom_text(data = text2plot,
            aes(Var2, Var1, label = sprintf("%0.3f", round(value, digits = 3)), color = value > 0.5),
            size = 3.5)+
  scale_color_manual(guide = 'none', values = c("black", "white")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        legend.position = 'right',#c(0.5, 0.75),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.justification = c(0, 0.5),
        legend.direction = 'vertical')+#"horizontal")+
  coord_fixed() + ggtitle('Simultaneous commutator test p-values') +
  scale_y_discrete(position = 'left', limits = rev(levels(pair.test2plot$Var2))) +
  guides(fill = guide_colorbar(barwidth = 1.5, barheight = 7, title.vjust = 1,
                               title.position = "top", title.hjust = 0.5,
                               title.theme = element_text(margin = margin(0,0,8,0))))+
  geom_tile(data = groupMat, aes(x = Var1, y = Var2), colour = "red", fill = NA, linewidth = 1)
dev.off()

p2 = ggplot(data = pair.test2plot, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "black",
                       limit = c(0,1), space = "Lab",
                       name="p-value:") +
  theme_minimal()+
  geom_text(data = text2plot, aes(Var2, Var1, label = sprintf("%0.3f", round(value, digits = 3)), color = value > 0.5), size = 3.5)+
  scale_color_manual(guide = 'none', values = c("black", "white")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        legend.position = 'right',#c(0.5, 0.75),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.justification = c(0, 0.5),
        legend.direction = 'vertical')+#"horizontal")+
  coord_fixed() + ggtitle('Simultaneous commutator test p-values') +
  scale_y_discrete(position = 'left') +scale_x_discrete(limits = rev(levels(pair.test2plot$Var2))) +
  guides(fill = guide_colorbar(barwidth = 1.5, barheight = 7, title.vjust = 1,
                               title.position = "top", title.hjust = 0.5,
                               title.theme = element_text(margin = margin(0,0,8,0))))+
  geom_tile(data = groupMat, aes(x = Var1, y = Var2), colour = "red", fill = NA, linewidth = 1)
p2
ggsave(filename = 'output/Plots/png/pvalMat.png', width = 5, height = 4.5, units = 'in')


