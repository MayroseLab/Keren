require(ape)
require(maps)
require(phytools)
require(readr)

# hard-coded details - please fill in dir
dir = # fill here your download directory
histories_dir = paste0(dir,"/simmap_format_histories/")
output_path = pastre0(dir+"/histories_visualization.pdf")

# function that plots a single stochastic mapping and writes it to a file
plot_sotchastic_mapping <- function(input_dir, fig_pos) {
  tree_str = read_file(input_path)
  tree = read.simmap(text=tree_str)
  cols<-setNames(c("blue","red"), c(0,1))
  plotSimmap(tree,cols,ftype="i",fsize=0.7, pts=TRUE, add=TRUE)
  
}
  
  
# plot sotichastic mappings in an input directory
file.names <- dir(histories_dir, pattern =".nwk")
pdf(onefile=TRUE, file = output_path) # failed to plot all the histories in the same file :(
par(mfrow=c(length(file.names)/2,length(file.names)/2), mfcol=c(length(file.names)/2,length(file.names)/2))
plot.new()
pos_lst = list(c(0,0), c(0,1), c(1,0), c(1,1))
for(i in 1:length(file.names)){
  input_path = paste0(histories_dir, file.names[i])
  res = plot_sotchastic_mapping(input_path, pos_lst[[i]])
}
dev.off()
