#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

require(ape)
require(maps)
require(phytools)
require(readr)
require(ggplot2)

plot_history <- function(axis, input_path, title_str)
{
  #par(xpd=NA, oma=c(0,0,2,0), ask=F, cex=0.42, font.main = 1)
  tree_str = read_file(true_history_path)
  tree = read.simmap(text=tree_str)
  cols<-setNames(c("blue","red"), c(0,1))
  plotSimmap(tree,cols,ftype="i",fsize=2, pts=TRUE, add=TRUE)
  title(title_str)
}

# plot the true history
mp_history_path = args[1]
ml_history_path = args[2]
mapping1_path = args[3]
mapping2_path = args[4]
mapping3_path = args[5]
ouput_path = args[6]
pdf(file = output_path)
plot.new()
par(bg = "white")           # default is likely to be transparent
split.screen(c(2, 1))       # split display into two screens
split.screen(c(1, 2), screen = 1) # split the upper half into 2
split.screen(c(1, 3), screen = 2) # now split the bottom half into 3

#m <- rbind(c(1, 1, 1, 2, 2, 2), c(3, 3, 4, 4, 5, 5))
# plot mp history
plot_history(1, mp_history_path, "Maximum parsimony partition")
plot_history(2, ml_history_path, "Maximum parsimony partition")
plot_history(3, mapping1_path, "")
plot_history(4, mapping2_path, "Stochastic mappings")
plot_history(5, mapping3_path, "")
dev.off()