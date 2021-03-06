#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly=TRUE)

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
#true_history_path =args[1]
#pdf(file = args[2])
par(xpd=NA, oma=c(0,0,2,0), ask=F, cex=0.42, font.main = 1)
plot.new()
#tree_str = read_file(true_history_path)
#tree_str = "((((((((S20:{0,0.0215273},S21:{0,0.0215273}):{0,0.10068},(S4:{0,0.0414948},S11:{0,0.0414948}):{0,0.080712}):{0,0.0459345},S5:{0,0.168141}):{0,0.00699474},((S28:{0,0.102636},S13:{0,0.102636}):{0,0.0718511},S31:{0,0.174487}):{0,0.000648973}):{0,0.0248507},((S14:{0,0.185971},(S19:{0,0.0977777},S8:{0,0.0977777}):{0,0.0881933}):{0,0.00385452},((((S18:{0,0.0874275},S25:{0,0.0874275}):{0,0.0544917},S24:{1,0.141919}):{0,0.0285361},S9:{0,0.170456}):{0,0.0140935},((S23:{0,0.0442562},S22:{0,0.0442562}):{0,0.0601431},S10:{0,0.104399}):{0,0.0801496}):{0,0.0052765}):{0,0.0101613}):{0,0.0198745},(S30:{0,0.0582782},S29:{0,0.0582782}):{0,0.161583}):{0,0.00101218},((((S12:{0,0.031},S17:{0,0.031}):{0,0.0541545},S2:{0,0.0851544}):{0,0.00286273},(S3:{0,0.000786646},S1:{0,0.000786646}):{0,0.0872305}):{0,0.122434},(S15:{0,0.0440658},S27:{0,0.0440658}):{0,0.166386}):{0,0.0104218}):{0,0.00230681},(((S16:{0,0.0636479},S7:{0,0.0636479}):{0,0.0530129},S26:{0,0.116661}):{0,0.0257351},(S6:{0,0.0602486},S32:{0,0.0602486}):{0,0.0821473}):{0,0.0807844});"
tree_str = "((((S16:{0,0.0636479},S7:{0,0.0636479}):{0,0.0530129},S26:{0,0.116661}):{0,0.0257351},(S6:{0,0.0602486},S32:{0,0.0602486}):{0,0.0821473}):{0,0.0807844},(((((((S20:{0,0.0215273},S21:{0,0.0215273}):{0,0.10068},(S4:{0,0.0414948},S11:{0,0.0414948}):{0,0.080712}):{0,0.0459345},S5:{0,0.168141}):{0,0.00699474},((S28:{0,0.102636},S13:{0,0.102636}):{0,0.0718511},S31:{0,0.174487}):{0,0.000648973}):{0,0.0248507},((S14:{0,0.185971},(S19:{0,0.0977777},S8:{0,0.0977777}):{0,0.0881933}):{0,0.00385452},((((S18:{0,0.0874275},S25:{0,0.0874275}):{0,0.0544917},S24:{1,0.0688384:0,0.0730806}):{0,0.0285361},S9:{0,0.170456}):{0,0.0140935},((S23:{0,0.0442562},S22:{0,0.0442562}):{0,0.0601431},S10:{0,0.104399}):{0,0.0801496}):{0,0.0052765}):{0,0.0101613}):{0,0.0198745},(S30:{0,0.0582782},S29:{0,0.0582782}):{0,0.161583}):{0,0.00101218},((((S12:{0,0.031},S17:{0,0.031}):{0,0.0541545},S2:{0,0.0851544}):{0,0.00286273},(S3:{0,0.000786646},S1:{0,0.000786646}):{0,0.0872305}):{0,0.122434},(S15:{0,0.0440658},S27:{0,0.0440658}):{0,0.166386}):{0,0.0104218}):{0,0.00230681});"
tree = read.simmap(text=tree_str)
cols<-setNames(c("blue","red"), c(0,1))
plotSimmap(tree,cols,ftype="i",fsize=2, pts=TRUE, add=TRUE)
#title(args[3])
dev.off()
