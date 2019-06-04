###################################################################
#Code to plot triads calculations, 2- vs. 3-host chains (Section S5)
#
#Date:
#07/05/2019 - First script
#
#Author:
#Gianluigi Rossi
###################################################################
#Packages
require(ggplot2)
require(scales)
###################################################################
rm(list = ls())

#load code for pair calculations
source("190424_KFE_pairs_calc_noPrints.R")
#load code for triad calculations
source("190510_KFE_triad_calc_noPrints.R")

###################################################################
load("KFE_SEI_triads_calculations.Rdata")


#transform the KFE data for 
pairs.list.fortri.final.edit1 <- pairs.list.fortri.final
names(pairs.list.fortri.final.edit1)[7:8] <- c("Prob", "Chain")
pairs.list.fortri.final.edit1[8] <- "A-U-B"

pairs.list.fortri.final.edit2 <- pairs.list.fortri.final
names(pairs.list.fortri.final.edit2)[7:8] <- c("Prob", "Chain")
pairs.list.fortri.final.edit2[7] <- pairs.list.fortri.final.edit2[8]
pairs.list.fortri.final.edit2[8] <- "A-B"

pairs.list.fortri.final.edit <- rbind(pairs.list.fortri.final.edit1, pairs.list.fortri.final.edit2)
pairs.list.fortri.final.edit$snps.d1 <- paste(pairs.list.fortri.final.edit$snps.d1, "(source)")
pairs.list.fortri.final.edit$snps.d2 <- paste(pairs.list.fortri.final.edit$snps.d2, "(infected)")

#quartz()

#pdf("Fig_S5.2.pdf", height = 9, width = 9)
ggplot(pairs.list.fortri.final.edit, aes(y = Prob, x = sample.t2/365, col = Chain)) +
  geom_line(size = 1.5) + 
  xlab(expression(paste(Delta,"t (years)",sep=""))) + ylab("Transmission probability")  +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  facet_grid(snps.d2 ~ snps.d1) +
  scale_x_continuous(breaks=c(0,2,4,6,8,10), labels=c(0,2,4,6,8,10))+
  theme_bw()
dev.off()

###################################################################




