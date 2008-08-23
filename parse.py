#from Python
import re
from string import Template

#user defined
import tree_traverse
import levels
from is_classes import Node 
from is_classes import Tree 

#module for parsing a newick tree. A newick tree is converted into a
#datastructure Tree composed of Node. See nodeAreas.py


def VerifyBrlen(token):
    """Verify that the given token conforms to a real floating point number.
       return True or False

       The definition follows python's definitition of a float, which is:
            floatnumber     ::= pointfloat | exponentfloat 
            pointfloat      ::= [intpart] fraction | intpart "." 
            exponentfloat   ::= (intpart | pointfloat) exponent 
            intpart         ::= digit+ 
            fraction        ::= "." digit+ 
            exponent        ::= ("e" | "E") ["+" | "-"] digit+ 
            """
    exponent = r'[eE][+-]\d+'
    fraction = r'\.\d+'
    intpart  = r'\d+'
    pf = Template('($ip)?($f)|($ip)\.')
    pointfloat = pf.substitute(ip=intpart, f=fraction)
    ef = Template('(($ip)|($pf))$e')
    exponentfloat = ef.substitute(ip=intpart, pf=pointfloat, e=exponent)
    fn = Template('($pf)|($ef)')
    floatnumber = fn.substitute(pf=pointfloat, ef=exponentfloat)
 
    re_floatnumber = re.compile(floatnumber)
    m = re_floatnumber.match(token)
    #assert (m != None), "Unrecognized format for brlen:"+token
    #assert (m.group == token), "unexpectedly long token:"+token
    if m != None:
        return True
    else:
        return False



def ExtractToken(line):
    """Extract the first proper newick tokens from the line;
       return a tuple (token_type, token). token_type is defined below

       The proper tokens are:
       Type 1   '('      

       Type 2   r',\s*' regular expression: comma followed by zero or more whitespace
                characters

       Type 3   ')[:floatingnumber]        - rightparen optionally followed by 
                                             : and a floatingpoint number denoting brlen

       Type 4   r'\S+'[:floatingnumber]    - any sequence of non-whitespace characters of 
                                             length atleast 1 optionally followed by : and a float

       Note that the above tokens are not fully refined: for example, (a) floatingnumber itself is an elaborate
       regular expression (b) '\S+' can be refined into r'(\w+)(#(\w+)(-\w+)?)?' with a taxon name possibly followed
       by island name which can possibly be followed by a zone name. The idea is that once the crudely defined
       tokens are identified, it can be verified that they conform to the proper syntax of taxon names, branch lengths
       etc.

       """
    
    #check if token is left par, and if so return its type (1)

    leftpar = r'\('
    re_leftpar = re.compile(leftpar)
    x = re_leftpar.match(line)
    if x != None:
        token = x.group()
        print 'Type 1: '+token
        return (1, token)

    # Check if Type 2: a
    # comma followed by zero or more whitespace characters

    not_names = r',\s*' 
    re_not_names = re.compile(not_names)
    x = re_not_names.match(line)
    if x != None:
        token = x.group()
        print 'Type 2: '+token
        return (2, token)

    # Check if Type 3: right par optionally followed by :float
    rightparfloat = r'\)(:[^,\)\s]+)?'
    re_rightparfloat = re.compile(rightparfloat)
    x = re_rightparfloat.match(line)
    if x != None:
        token = x.group()
        print 'Type 3: '+token
        l = token.split(':')
        if len(l) == 2:
            assert VerifyBrlen(l[1]), 'Possibly Illegal Branch Length:'+l[1]
        return (3, token)        
            
    # explaining the regular expression "names":
    # the token can take the following forms:
    # a. taxon[:branchlength]
    # b. taxon#area[:branchlength]
    # c. taxon#island-zone[:branchlength]
    # each of taxon, area, island and zone should themselves match \w+:
    #                                      (sequence of one or more alphanumeric characters)
    # branchlength should be a floatingpoint number as defined in python.
    # See function VerifyBrlen
    
    names = r'(\w+)(#(\w+)(-\w+)?)?(:[^,\)\s]+)?'
    re_names = re.compile(names)
    x = re_names.match(line)
    assert x != None, line+"Unrecognized token"
    token = x.group()
    print 'Type 4: '+token
    l = token.split(':')
    if len(l) == 2:
        assert VerifyBrlen(l[1]), 'Possibly Illegal Branch Length:'+l[0]+' '+l[1]
        
    return (4, token)        


def Read(tree_name, taxon_table, state_table, tree): 
    """Parse a newick string return the tree as a datastructure """    

    # variables for iterating through the string
    original_tree = tree
    tree_len = len(tree)
    parsed_len = 0

    # stack for parsing the string. Think push-down-automata
    stack = []

    # the logic for parsing is quite simple. It looks complicated only
    # because one has to catch syntax errors. The logic is as follows
    # (disregarding branch lengths):
    # Keeping pushing tokens on to a stack, until a ')' is encountered. At
    # that time, keep popping the items from the stack until a '(' is
    # encountered, and do any
    # operation with them (e.g., a node can be created that is the parent
    # of all the  popped nodes ) - essentially the operation should
    # result in the popped items being combined into one item. Now, there
    # should be a '(' at the top of the stack. Pop that, and push the
    # combined item to the stack. Essentially, replace, say, (x, y, z) with a single
    # item W. Now that doesn't take care of of expressions without the
    # outermost pair of parenthesis. The second outer while loop,
    # conditioned on len(stack) > 0, is to take care of those cases. Things
    # to keep in mind in this latter case: (a) only the outermost pair can
    # go missing, (b) this can be recognized only on reaching the end of
    # the whole string.


    while parsed_len < tree_len:
        (token_type, token) = ExtractToken(tree)
        token_len = len(token)
        parsed_len = parsed_len + token_len
        tree = original_tree[parsed_len:tree_len]
    
        # comma spaces token
        if token_type == 2:
            continue

        # right paren token
        if token_type == 3:
#           read in the potential branchlength
            z = token.split(':')
            prebrlen = z[0]
            brlen = 0.0
            if len(z) == 2:
                brlen = float(z[1])
            #else:
            #    print 'WARNING: UNSPECIFIED BRANCH LENGTH SET TO 0.0'

            node = Node()
            node.brlen = brlen
            try:
                top = stack.pop()
                while top.label != '(':
                    node.children.append(top)
                    top = stack.pop() 

                expr = ''
                for child in node.children:
                    child.parent = node
                    if expr == '':
                        expr = expr+child.label
                    else:
                        expr = expr+', '+child.label
                
                assert top.label == '('
                assert prebrlen == ')'
                node.label = '('+expr+')'
                stack.append(node)
                continue

            except IndexError:
                print """Malformed tree. A missing left paren?"""
        


        # left paren token
        if token_type == 1:
            node = Node()
            node.label = token
            stack.append(node)
            continue

        # type 4 token 
        z = token.split(':')
        name = z[0]
        brlen = 0.0
        if len(z) == 2:
            brlen = float(z[1])
        #else:
        #    print 'WARNING: UNSPECIFIED BRANCH LENGTH SET TO 0.0'
        node = Node()
        #node.label = name
        node.brlen = brlen

        taxon_state = name.split('#')
        taxon = taxon_state[0]
        node.label = taxon

        if len(taxon_table) == 0:
            node.taxon_name = taxon
        else:
            node.taxon_name = taxon_table[taxon]
            
        if len(taxon_state) > 1:
            state = taxon_state[1]
            if len(state_table) == 0:
                node.state = state
            else:
                node.state = state_table[taxon]
        else:
            node.state = 'NOWHERE'

        stack.append(node)
    
    tree = Tree()
    if len(stack) == 1:
        tree.root = stack.pop()
        tree.root.parent = None
        #return tree
    else:
        expr = ''
        node = Node()
        while len(stack) > 0:
            top = stack.pop()
            assert top.label != '(', """Malformed tree. A missing right paren?""" 
            node.children.append(top)
            top.parent = node
            if expr == '':
                expr = expr+top.label
            else:
                expr = expr+', '+top.label

        node.label = '(' + expr + ')'
        tree.root = node
        tree.root.parent = None

    tree.name = tree_name
    print 'Something'
    tree_traverse.PostOrderTraverse(tree.root, tree_traverse.CalculateNodeAges)
    tree_traverse.PrintLabel(tree.root)
    tree_traverse.PostOrderTraverse(tree.root, tree_traverse.PrintSelfParentLabel)
    tree_traverse.TreeLeaves(tree)
    levels.DemarcateLevels(tree)
    levels.PrintAllLevels(tree)
    temp_levels=tree.levels
    #k=len(tree.levels)
   # for i in range(1,2):
    #    tree.levels[i]=temp_levels[k-i]  
    return tree
