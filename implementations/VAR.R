library(eigTest)
library(plotrix)
library(reshape2)
library(ggplot2)
#library(tikzDevice)

m = nrow(countryMacro[[1]])
eigTest(countryCoeff, cn = sqrt(m), cov.arr = countryCovar, testType = 'chi')
eigTest(countryCoeff, cn = sqrt(m), cov.arr = countryCovar, testType = 'gam')

partialTest(countryCoeff, cn = sqrt(m), k = 2, cov.arr = countryCovar, testType = 'chi')
partialTest(countryCoeff, cn = sqrt(m), k = 2, cov.arr = countryCovar, testType = 'gam')

###  pairwise commutator test  ###
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
groupMat[1:3,1:3] = 'Asia'
groupMat[4:6,4:6] = 'EU'
groupMat[7:8,7:8] = 'NA'
groupMat = melt(groupMat, na.rm = T)
groupMat$value = as.factor(groupMat$value)

#tikz('implementations/Plots/pvalMat.tikz', standAlone = F, width = 5, height = 4.5)
ggplot(data = pair.test2plot, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "black",
                       limit = c(0,1), space = "Lab",
                       name="p-value:") +
  theme_minimal()+
  geom_text(data = text2plot, aes(Var2, Var1, label = sprintf("%0.3f", round(value, digits = 3)), color = value > 0.5), size = 4)+
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
  geom_tile(data = groupMat, aes(x = Var1, y = Var2), colour = "red", fill = NA, size = 1)
#dev.off()




