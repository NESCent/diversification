class Tree:
    """The class description of a tree.
    
    The class contains only one attribute: the __init__ function.
    But whenever an instance of Tree is created (i.e., during
    instantiation), __init__ is called, with its "self" parameter pointing
    to the instance object. The __init__ function, because of the
    assignments in it, cause the following data attributes to spring into
    existence in the instance object (but not in the class object).

    (all the attributes will get their proper value after appropriate
    functions are called during parsing.)

    all_leaves      a list of all leaves in the tree.
    
    partial_order   a list of all nodes in the tree such that a parent node
                    appears onnly after its children have appeared.
                        
    height          is the length of the longest path (measured in terms of
                    the number of branches) from the root of the tree to a
                    leaf.
    
    num_leaves      the total number of leaves in the tree.

    name            the name of the tree (like tree_13, tree_14 in
                    onetree.nex) as in the nexus file that contains the
                    tree in parenthesis notation.

    levels          A level is a portion of the tree between two
                    branching events (or the portion from the leaves to
                    the most recent branching event). 
                    
                    levels[] is a list that contains all the levels in the
                    tree. A level should be implemented as an instance of
                    the class Level (see below).
    """

    def __init__(self):
        self.all_leaves = []
        self.partial_order = []
        self.height = 0
        self.num_leaves = 0
        self.root = None
        self.name = None
        self.levels = []

class Level:
    """
    A level in the genealogical history of a sample is that portion 
    two successive branching events (or between the tips and the first
    branching event, as the case may be).

    The class contains only one attribute: the __init__ function.
    But whenever an instance of Level is created (i.e., during
    instantiation), __init__ is called, with its "self" parameter pointing
    to the instance object. The __init__ function, because of the
    assignments in it, cause the following data attributes to spring into
    existence in the instance object (but not in the class object).

    (all the attributes will get their proper value after appropriate
    functions are called during level demarcation.)

    
    begin_time      As mentioned earlier, a level is a portion of the
                    sample history between two successive branching events,
                    or between the tips and the most recent branching event
                    in the history of the sample. begin_time either 0.0
                    (for the most recent level) or is the age (date) of the
                    more recent of the two branching events that demarcate
                    it. Time is measured backwards, with the tips having
                    the age 0.0.

    begin_node      the parent node of the more recent branching event that
                    demarcates the level (a branching event involves one
                    parent node and two daughter nodes).

                    For the very first (i.e., most recent) level, the begin_node is one of the
                    leaves.

    end_time        The age of the less recent of the two branching events
                    that demarcate the level.

    end_node        the parent node of the less recent branching event that
                    demarcates the level 

    
    lineages        A list of lineages in a level. 
                    * Lineages and nodes are used interchangeably *, with a
                    lineage being identified with the more recent of its
                    two end points. For example, the following tree has
                    three nodes x, y and z; and three lineages 1, 2 & 3. 
                    We consider 1 = a, 2 = b and 3 = b.
                                
                                |
                              3 |
                                |z
                             1 / \ 2
                              /   \
                             x     y

    event_history   The (backwards) proposal algorithm will propose
                    migration events for each level. 
                    The character state assignment at the beginning (i.e.,
                    at its most recent) and the subsequent migration events
                    are recorded in event_history.

                    Suppose level = 4, and that there are 4 lineaages 1, 2, 3, & 4.
                    And let the character state assignments be {1:0, 2:1, 3:1,
                    4:1}. And let following events happen before the chain goes to
                    level 3:

                    a. lineage 2 goes from character state 1 to character state 0
                    b. lineage 3 goes from character state 1 to character state 0
                    c. lineages 2 & 3 coalesce to form a parent lineage 5; and this
                       ends level 4.

                    The event history for level 4 will be a list: [{1:0, 2:1, 3:1,
                    4:1}, 2, 3]. The first element is the character state at the
                    beginning. The other elements are the lineages that undergo
                    migration. Given that lineages only migrate between 2 states,
                    it is sufficient to list the lineages that migrated (their
                    states can be inferred).

                    For the next level, (i.e., level 3) the initial character state
                    assignment will be: {1:0, 4:1, 5:1} (since 2 & 3 were in state
                    1, their parent 5 also will be in state 1).
    """
    def __init__(self):
        self.begin_time = 0.0
        self.end_time = 0.0
        self.lineages = []
        self.begin_node = None
        self.end_node = None
        self.event_history = []                                     
    
    def __str__(self):
        return lineages

class Node:
    """The class description of a node in a tree

    The class contains only one attribute: the __init__ function.
    But whenever an instance of Node is created (i.e., during
    instantiation), __init__ is called, with its "self" parameter pointing
    to the instance object. The __init__ function, because of the
    assignments in it, cause the following data attributes to spring into
    existence in the instance object (but not in the class object).

    (all the attributes will get their proper value after appropriate
    functions are called during parsing.)

    children            a list containing the two children of a node. For leaves,
                        this list should be empty.

    all_leaves          a list containing all the leaves whose most recent 
                        common ancestor is the node.
                        

                        the list contains the *node itself* for nodes that are leaves.

    partial_order       each node is a root of a subtree. 
                        partial_order is a list of nodes of this subtree
                        such that a parent node appears in the list only after its
                        children have appeared.

    num_leaves          the number of leaves for which the node is the most
                        recent common ancestor.

    parent              the parent of the node in the tree. 
                        if the node is teh root of the whole tree, this
                        should be None.
    
    label               The label of a node is the subtree under the node
                        (i.e., one for which the node is the root)
                        described in parenthesis notation. 

                        This field is useful to see if parsing is being
                        done correctly, but otherwise is not used much.

    taxon_name          for leaves, this is the name of the taxon as in the
                        nexus file.

    state               this field is *deprecated* - which means it is not
                        really used anywhere, and in the future also it
                        shouldn't be used.

    
    brlen               the length of the branch (or equivalently, lineage) identified
                        with the node. Recall:
                        * Lineages and nodes are used interchangeably *, with a
                        lineage being identified with the more recent of its
                        two end points. For example, the following tree has
                        three nodes x, y and z; and three lineages 1, 2 & 3. 
                        We consider 1 = a, 2 = b and 3 = b.
                                    
                                    |
                                  3 |
                                    |z
                                 1 / \ 2
                                  /   \
                                 x     y

    bifurcation_age     The is really the age of the node. In the above
                        picture, z.bifurcation_age would have been 0.2 if
                        the branches (x, z) and (y, z) have the same
                        length, 0.2.

                        For leaves, this age is 0.0

    parent_age          The age of the the parent of the node; equals
                        bifurcation_age + brlen
    """
    def __init__(self):
        self.children = []
        self.all_leaves = []
        self.partial_order = []
        self.num_leaves = 0
        self.parent = None
        self.label = 'None'
        self.taxon_name = 'None'
        self.state = 'NOWHERE'
        self.brlen = 0.0
        self.bifurcation_age = 0.0
        self.parent_age = 0.0

    def __repr__(self):
        return self.label
        
# Tally & Perry: Ignore the following class definition.

# class definition begins.        
#class CondJumpChain:
#    """ The conditional jump chain class definition.
#
#        The class CondJumpChain itself is an object, and this object has only
#        function attributes, namely: 
#
#              fn. attrib.                          short desc.
#              ----------                           ----------
#
#            MakeTreeByLevels                makes a level-by-level representation
#                                            of the input tree G.
#            
#            MakeStateSpacesByLevels         constructs the state space of the
##                                            jump chain level by level.
#
#            MakeTransitionMatricesbyLevels  constructs the transition
#                                            matrix level-by-level.
#
#            SampleFromIS                    Runs the chain once and returns
#                                            the resulting event history
#                                            along with its density.
#
#            __init__                        This is called on
#                                            instantiation.
#
#        Except __init__, the other functions are actually defined outside
#        under a slightly different name. A function FunctionName is
#        will actually be a a local variable of the class; a function by the
#        slightly different name __FunctionName is defined outside and
#        assigned to the local variable FunctionName. 
#
#        Detailed descriptions of the functions are provided with their
#        definitions.
#
#        Data attributes of the instance objects:
#        ----------------------------------------
#
#        While the class CondJumpChain itself has only function attributes,
#        the instances of this class also have data attributes that spring
#        into existence upon instantiation. 
#
#        For example, the following statement creates an instance
#        cond_jump_chain_instance of the class CondJumpChain. 
#
#        cond_jump_chain_instance = CondJumpChain(G, delta, sigma)
#
#        At instantiation, the __init__ function of the class is automatically
#        called. Calling __init__ causes  the following variables
#        to spring into existence as data attributes of the
#        instance object cond_jump_chain_instance. 
#        
#        data attribute                      brief desc.
#        --------------                      ----------
#
#        tree_by_levels                      return value of call to MakeTreeByLevels
#        state_spaces_by_levels              return value of call to MakeStateSpacesByLevels
#        initial_state                       return value of call to GetInitialStateofCondJumpChain
#        transition_matrices                 return value of call to MakeTransitionMatricesbyLevels
#    """
#
#    # The following are functions defined outside the class.
#    MakeTreeByLevels = __MakeTreeByLevels
#    MakeStateSpacesByLevels = __MakeStateSpacesByLevels
#    GetInitialStateofCondJumpChain = __GetInitialStateofCondJumpChain
#    MakeTransitionMatricesbyLevels = __MakeTransitionMatricesbyLevels
#    SampleFromIS = __SampleFromIS
#
#    def __init__(G, delta, sigma):
#        
#        # When an instance object is initialized by, say, the statement
#        # cond_jump_chain_instance = CondJumpChain(G, delta, sigma), tree_by_levels, state_spaces_by_levels, 
#        # initial_state and transition_matrices, will spring into existence as the data
#        # attributes of the *instance* object cond_jump_chain_instance. 
#        # For each function attribute of the class object CondJumpChain,
#        # there will be a corresponding method attribute in the instance
#        # object.
#        self.tree_by_levels = self.MakeTreeByLevels(G)
#        self.state_spaces_by_levels = self.MakeStateSpacesByLevels(G, delta)
#        self.initial_state = self.GetInitialStateofCondJumpChain(G, delta)
#        self.transition_matrices = self.MakeTransitionMatricesbyLevels(sigma)
###### class definition ends.        
