from ete3 import Tree, NodeStyle, TreeStyle

#output_dir = "C:/Users/ItayM5/Google Drive/MSc/posters and presentations/presentations/supplementary_materials/"
output_dir = "/groups/itay_mayrose/halabikeren/graphics/"
output_name = output_dir+ "tree.png"

tree_str = Tree('((S1:1,S2:1)N1:1,(S3:1,S4:1)N2:1;)')
tree = Tree(tree_str)

# Basic tree style
ts = TreeStyle()
ts.show_leaf_name = True
ts.show_scale = True
ts.rotation = 90
ts.branch_vertical_margin = 5

internals = NodeStyle()
internals["hz_line_type"] = 0
internals["vt_line_type"] = 0
internals["vt_line_width"] = 2
internals["hz_line_width"] = 2
internals["hz_line_color"] = "Silver" #454545" darker gray
internals["vt_line_color"] = "Silver" 
internals["shape"] = "circle"
internals["size"] = 3
internals["fgcolor"] = "Silver" 


clade = NodeStyle()
clade["hz_line_type"] = 0
clade["vt_line_type"] = 0
clade["vt_line_width"] = 2
clade["hz_line_width"] = 2
clade["hz_line_color"] = "Chocolate"
clade["vt_line_color"] = "Chocolate" #grey 
clade["shape"] = "circle"
clade["size"] = 3
clade["fgcolor"] = "Chocolate" 

cladeB = NodeStyle()
cladeB["hz_line_type"] = 0
cladeB["vt_line_type"] = 0
cladeB["vt_line_width"] = 2
cladeB["hz_line_width"] = 2
cladeB["hz_line_color"] = "Chocolate"
cladeB["vt_line_color"] = "Chocolate" #grey 
cladeB["shape"] = "circle"
cladeB["size"] = 8 
cladeB["fgcolor"] = "Chocolate" #olive green

# Applies the same static style to all nodes in the tree. Note that,
# if "nstyle" is modified, changes will affect to all nodes
'''
for n in t.traverse():
	n.set_style(internals)

group = t.get_common_ancestor('aaaaaaaaac', 'aaaaaaaaad')

for n in group.traverse():
	if n.is_leaf():
		n.set_style(cladeB)
	else:
		n.set_style(clade)
'''

'''
#t sim 35
for n in t.traverse():
	names = ['aaaaaaaaau', 'aaaaaaaaaw', 'aaaaaaaabi', 'aaaaaaaaaf', 
					'aaaaaaaaae', 'aaaaaaaaap', 'aaaaaaaaac']
	if n.name in names:
		n.set_style(cladeB)
	else:
		n.set_style(internals)
		n.set_style(internals)
t.prune(names)

for n in t.traverse():
	if n.is_leaf():
		n.set_style(cladeB)
	else:
		n.set_style(clade)		
'''

#t sim 18
for n in tree.traverse():
	names = ['S1','S2','S3','S4']
	if n.name in names:
		n.set_style(cladeB)
	else:
		n.set_style(internals)

tree.prune(names)

for n in tree.traverse():
	if n.is_leaf():
		n.set_style(cladeB)
	else:
		n.set_style(clade)

tree.show(tree_style=ts)
tree.render(output_name, w=300, units="mm", tree_style=ts)
