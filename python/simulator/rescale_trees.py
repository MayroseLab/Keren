from ete3 import Tree
import re, argparse, os

def get_tree_size(tree):
    size = 0
    for node in tree.traverse():
        size+= node.dist
    return size


def rescale_tree(orig_tree_path, size, new_tree_path):
    tree = Tree(orig_tree_path)
    orig_size = get_tree_size(tree)
    factor = size / orig_size
    for node in tree.traverse():
        node.dist = node.dist * factor
    tree.write(outfile=new_tree_path)
    return 0



if __name__ == "__main__":

    # process input from cmd
    parser = argparse.ArgumentParser(description='rescales all trees in a given directory into a given size')
    parser.add_argument("--size", "-s", help="the size you wish the trees to be rescaled to", required=True)
    parser.add_argument("--orig_trees_path", "-i", help="directory which holds the trees to be rescaled", required=True)
    parser.add_argument("--new_trees_path", "-o", help="directory which holds the rescaled trees", required=True)
    args = parser.parse_args()
    size = int(args.size)
    orig_trees_path = args.orig_trees_path
    new_trees_path = args.new_trees_path
    if not os.path.exists(new_trees_path):
        res = os.system("mkdir " + new_trees_path)

    print("**** starting trees rescaling ****")
    for tree in os.listdir(orig_trees_path):
        orig_tree_path = orig_trees_path + tree
        new_tree_path = new_trees_path + tree
        res = rescale_tree(orig_tree_path, size, new_tree_path)

    print("**** process complete ****")


