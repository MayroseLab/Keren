from ete3 import Tree
import re, argparse

def convert_history_to_simmap(tree_path, output_path):
    label_regex = re.compile("{(.*?)}")
    history = Tree(tree_path, format=1)
    #tree_str = history.write(format=8, outfile=None) # get a string of the tree newick only with names
    visited_nodes = [history.get_tree_root()]
    for node in history.traverse("postorder"):
        if not node in visited_nodes:
            node_label = label_regex.search(node.name).group(1)
            node_name = node.name.replace(label_regex.search(node.name).group(0),"")
            node_expression = ""
            if node.is_leaf():
                node_expression = node_expression + node_name
            node_expression = node_expression + ":{" + node_label + "," + str(node.dist)
            visited_nodes.append(node)
            if not node == history.get_tree_root():
                try:
                    curNode = node.up
                    while "mapping" in curNode.name:
                        node_label = label_regex.search(curNode.name).group(1)
                        node_expression = node_expression + ":" + node_label + "," + str(curNode.dist)
                        visited_nodes.append(curNode)
                        curNode = curNode.up
                except:
                    pass
            node_expression = node_expression + "}"
            node.name = node_expression
            #tree_str = tree_str.replace(node.name, node_expression)
    # remove mapping nodes
    for node in history.traverse("postorder"):
        if "mapping" in node.name:
            node.delete()
    # write tree on your own because ete3 writer has a bug
    tree_str = getTreeStr(history.get_tree_root())
    print("tree_str: ", tree_str)
    with open(output_path, "w") as output_file:
        output_file.write(tree_str + ";")


# writes tree with nodes names only in newick format
def getTreeStr(root):
    if root.is_leaf():
        return root.name

    tree_str = "("
    for child in root.get_children():
        tree_str = tree_str + getTreeStr(child) + ","
    tree_str = tree_str[:-1] + ")" + root.name
    return tree_str





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="converts history to simmap format")
    parser.add_argument("-i", "--input_tree_path", help="path to input tree, labeled by {-} format")
    parser.add_argument("-o", "--output_tree_path", help="path to output tree, labeled in simmap format")

    # input_tree_path = "C:/Users/ItayMNB7/Desktop/true_history.nwk"
    # output_tree_path = input_tree_path.replace(".nwk", "_simmap.nwk")

    args = parser.parse_args()
    input_tree_path = args.input_tree_path
    output_tree_path = args.output_tree_path

    convert_history_to_simmap(input_tree_path, output_tree_path)