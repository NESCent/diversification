from is_classes import Node 
from is_classes import Tree 
from is_classes import Level 

import levels

def Max(n1, n2):
    if (n1 > n2):
        return n1
    else:
        return n2

def concat(l1, l2):
    l = []
    l.extend(l1)
    l.extend(l2)
    return l


def NodeHeight(node):
    children_heights = [NodeHeight(child) for child in node.children]    
    # the starting value of -1 is provided to take care of leaves that have
    # no children.
    max_height = reduce(Max, children_heights, -1) + 1
    node.height = max_height
    return max_height


def TreeHeight(tree):
    tree.height = NodeHeight(tree.root)
    return tree.height

def NodeLeaves(node):
    if node.children == []:
        node.all_leaves = [node]
        node.num_leaves = 1
        return node.all_leaves
    children_leaves= [NodeLeaves(child) for child in node.children]
    node.all_leaves = reduce(concat, children_leaves)
    node.num_leaves = len(node.all_leaves)
    return node.all_leaves

def TreeLeaves(tree):
    tree.all_leaves = NodeLeaves(tree.root)
    tree.num_leaves = len(tree.all_leaves)
    return tree.all_leaves

#def MakeDelta(tree):
#    for j in range(0, tree.num_leaves):
#        tree.delta[j] = tree.all_leaves[j].state

def SubtreePartialOrder(subtree_root):
    if subtree_root.children == []:
        return [subtree_root]
    children_po = [SubtreePartialOrder(child) for child in subtree_root.children]
    subtree_root.partial_order = concat(reduce(concat, children_po), [subtree_root])
    return subtree_root.partial_order

def TreePartialOrder(tree):
    tree.partial_order = SubtreePartialOrder(tree.root)
    return tree.partial_order

def Average(l):
    return sum(l)/len(l)

def PrintLabel(node):
    print node.label
    return node.label

def PrintSelfParentLabel(node):
    if node.parent != None:
        print node.label + ':' + str(node.bifurcation_age) + '   ' + node.parent.label + ':' + str(node.parent_age)
    else:
        print node.label + ':' + str(node.bifurcation_age) + '   ' + str(node.parent)

def CalculateNodeAges(node):
    if node.children == []:
        node.bifurcation_age = 0.0
        node.parent_age = node.brlen
        print node.label, node.brlen, node.parent_age, node.bifurcation_age
        return
    child = node.children[0]
    node.bifurcation_age = child.parent_age
    node.parent_age = node.bifurcation_age + node.brlen
    print node.label, node.brlen, node.parent_age, node.bifurcation_age


def PostOrderTraverse(subtree_root, fn):
    """Traverse the subtree with the root subtree_root in post order, and
       perform fn on each node visited"""

    for child in subtree_root.children:
        PostOrderTraverse(child, fn)
    fn(subtree_root)

def PrintHistory(tree):
    for level in tree.levels:
        levels.PrintLevelHistory(level)

