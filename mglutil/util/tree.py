# ##################################################################################################
#  Disclaimer                                                                                      #
#  This file is a python3 translation of AutoDockTools (v.1.5.7)                                   #
#  Modifications made by Valdes-Tresanco MS (https://github.com/Valdes-Tresanco-MS)                #
#  Tested by Valdes-Tresanco-MS and Valdes-Tresanco ME                                             #
#  There is no guarantee that it works like the original distribution,                             #
#  but feel free to tell us if you get any difference to correct the code.                         #
#                                                                                                  #
#  Please use this cite the original reference.                                                    #
#  If you think my work helps you, just keep this note intact on your program.                     #
#                                                                                                  #
#  Modification date: 10/5/20 17:47                                                                #
#                                                                                                  #
# ##################################################################################################

#
# Last modified on Wed Oct 10 15:07:59 PDT 2001 by lindy


class TreeNode:
    """Base class of generic tree nodes.

    This class will work all by itself with data being
    whatever you want your nodes to be. It also could
    be used as the superclass of specific subclasses.
    """

    def __init__(self, parent=None, data=None):
        self.parent = parent
        self.children = []
        if parent:
            parent.add_child(self)
        if data:
            self.data = data

    def add_child(self, child):
        """Add given node to children of self.

        parent.add_child( node) adds node to children of parent.
        This is done in the __init__ if parent is given when node
        in instanciated (that's the prefered method).
        """
        self.children.append(child)

    def pre_traverse(self, f, *args, **kw):
        """Apply f to yourself and then your children recursively
        The function f must take the treeNode as it's first argument.
        """
        args = list(args)
        args[0] = self
        a = tuple([f] + args)

        f(*args, **kw)
        for c in self.children:
            c.pre_traverse(*a, **kw)

    def post_traverse(self, f, *args, **kw):
        """Traverse children with f and args, then visit parent.
        The function f must take the treeNode as it's first argument.
        """
        args = list(args)
        args[0] = self
        a = tuple([f] + args)
        for c in self.children:
            c.post_traverse(*a, **kw)
        f(*args, **kw)

    def get_iterator(self):
        """over-ride me to supply an appropriate subclass of TreeIterator"""
        raise NotImplementedError


class TreeIterator:
    """This iterator class is not finished yet.
    """

    def __init__(self, node):
        self.iterRoot = node
        self.currentNode = None  # set by subclass
        self.done = None

    def current(self):
        """return the currently-visited node"""
        return self.currentNode

    def done(self):
        """Returns false (None) until the traversal is finished"""
        return self.done

    def first(self):
        """Reset the currentNode to the initally specified node
        """
        self.currentNode = self.iterRoot

    def __next__(self):
        """Move currentNode on to the next item in the traversal.

        Over-ride this method to provide a specific type of traversal
        """
        # move on to next item
        raise NotImplementedError


class PostTreeIterator(TreeIterator):
    """This iterator class is not finished yet.
    """

    def __init__(self, node):
        TreeIterator.__init__(self)
        self.nodeList = []
        self.iterRoot.post_traverse(self.nodeList.append)
        self.currentIx = 0

    def __next__(self):
        self.currentNode = self.nodeList[self.currentIx]
        self.currentIx = self.currentIx + 1
        if self.currentIx == len(self.nodeList):
            self.done = 1


class PreTreeIterator(TreeIterator):
    """This iterator class is not finished yet.
    """

    def __init__(self, node):
        TreeIterator.__init__(self)
        self.nodeList = []
        self.iterRoot.pre_traverse(self.nodeList.append)
        self.currentIx = 0

    def __next__(self):
        self.currentNode = self.nodeList[self.currentIx]
        self.currentIx = self.currentIx + 1
        if self.currentIx == len(self.nodeList):
            self.done = 1
