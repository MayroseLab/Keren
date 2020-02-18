#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

require(ape)
require(maps)
require(phytools)
require(readr)

if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
}

# function that plots a single stochastic mapping and writes it to a file
plot_sotchastic_mapping <- function(input_path, num) { ###CHANGED HERE###
  tree_str = read_file(input_path)
  tree = read.simmap(text=tree_str)
  cols<-setNames(c("blue","red"), c(0,1))
  plotSimmap(tree,cols,ftype="i",fsize=0.7, pts=TRUE, add=TRUE)
  if (num==0) {
    title("expected history")
  } else {
    title(paste0("history ",num)) 
  }
}

# plot the true history
true_history_path =args[1]
pdf(file = args[2])
par(xpd=NA, oma=c(0,0,2,0), ask=F, cex=0.42, font.main = 1)
plot.new()
tree_str = read_file(true_history_path)
tree = read.simmap(text=tree_str)
cols<-setNames(c("blue","red"), c(0,1))
plotSimmap(tree,cols,ftype="i",fsize=0.7, pts=TRUE, add=TRUE)
title("true history")
dev.off()
