#! /usr/bin/env python

import sys
import operator

oo = sys.maxint / 2 # oo is infinity

class Node:
    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.edges = {}
        self.link = 0

    def edge_length( self, position ):
        return min( self.end, position + 1 ) - self.start

    def get_edges( self ):
        return self.edges

    def set_start( self, start ):
        self.start = start

    def get_start( self ):
        return self.start

    def get_end( self ):
        return self.end

    def set_end( self, end ):
        self.end = end

    def put_edge( self, key, val ):
        self.edges[ key ] = val

    def set_link( self, node ):
        self.link = node

    def get_link( self ):
        return self.link

    def get_edge( self, key ):
        if key in self.edges:
            return self.edges[ key ]
        else:
            return None

    def str_val( self ):
        string = "start: " + str( self.start ) + " end: " + str( self.end ) + " edges: " + str( self.edges )
        return string

    def increment_start( self, val ):
        self.start += val


class SuffixTree:
    def __init__( self, t_len ):
        self.cur_node = 0
        self.nodes_added = 0
        self.nodes = [ Node( 0, oo ) ] * ( 2 * t_len + 2 )
        self.text = ""
        self.root = self.new_node( -1, -1 )
        self.pos = -1
        self.need_link = 0
        self.r = 0   # remainder - how many suffixes you haven't added in yet

        self.active_node = self.root
        self.active_edge = 0
        self.active_length = 0

    '''__author__ = agolotin '''
    ''' Performs recursive depth first search of the suffix tree  
        and returns the suffix array as a list '''
    def depthFirstSearch(self, main_edge = 1, traversed_so_far = 0):
        suffix_array = list()
        for edge_value, edge_index in sorted(self.nodes[main_edge].get_edges().iteritems(), key=operator.itemgetter(0)):
            if self.nodes[edge_index].get_end() == oo:
                suffix_array.append(self.nodes[edge_index].get_start() - traversed_so_far)
            else:
                traversed_so_far += self.nodes[edge_index].get_end() - self.nodes[edge_index].get_start()
                suffix_array += self.depthFirstSearch(edge_index, traversed_so_far)
                traversed_so_far -= self.nodes[edge_index].get_end() - self.nodes[edge_index].get_start()

        return suffix_array

    def add_char( self, c ):
        self.text += c
        self.pos += 1

        self.need_link = -1
        self.r += 1

        # while there are suffixes to be added
        while self.r > 0:
            if self.active_length == 0:
                self.active_edge = self.pos

            if self.get_active_edge() not in self.nodes[self.active_node].get_edges():
                leaf = self.new_node(self.pos, oo)
                self.nodes[self.active_node].put_edge(self.get_active_edge(), leaf)
                self.add_suffix_link(self.active_node)
            else:
                nex = self.nodes[self.active_node].get_edge(self.get_active_edge())#node that we are going to walk down
                if self.walk_down(nex):
                    continue
                if self.text[self.nodes[nex].get_start() + self.active_length] == c:
                    self.active_length += 1
                    self.add_suffix_link(self.active_node)
                    break

                split = self.new_node(self.nodes[nex].get_start(),
                                      self.nodes[nex].get_start() + self.active_length)
                self.nodes[self.active_node].put_edge(self.get_active_edge(), split)
                leaf = self.new_node(self.pos, oo)
                self.nodes[split].put_edge(c, leaf)
                self.nodes[nex].increment_start(self.active_length)
                self.nodes[split].put_edge(self.text[self.nodes[nex].get_start()], nex)
                self.add_suffix_link(split)

            self.r -= 1
            if self.active_node == self.root and self.active_length > 0:
                self.active_length -= 1
                self.active_edge = self.pos - self.r + 1
            else:
                self.active_node = self.nodes[self.active_node].get_link() if self.nodes[self.active_node].get_link() > 0 else self.root

    def walk_down( self, nex ):
        if self.active_length >= self.nodes[ nex ].edge_length( self.pos ):
            self.active_edge += self.nodes[ nex ].edge_length( self.pos )
            self.active_length -= self.nodes[ nex ].edge_length( self.pos )
            self.active_node = nex
            return True
        return False


    def get_active_edge( self ):
        return self.text[ self.active_edge ]

    def new_node( self, start, end ):
        self.cur_node += 1
        self.nodes[ self.cur_node ] = Node( start, end )
        self.nodes_added += 1
        return self.cur_node

    def add_suffix_link( self, node ):
        if self.need_link > 0:
            self.nodes[ self.need_link ].set_link( node )
        self.need_link = node

    ''' Modified to output the starting and ending indicies of an edge '''
    def edge_string( self, node ):
        start = self.nodes[ node ].get_start()
        end = min( self.pos + 1, self.nodes[ node ].get_end() )
        _str = self.text[start:end]
        _str += " " + str(self.nodes[node].get_start()) + ", "
        _str += "$" if self.nodes[node].get_end() == oo else str(self.nodes[node].get_end())
        return _str

    def print_tree(self ):
        self.print_edges( self.root )

    def print_graphviz_tree( self ):
        print "digraph {"
        print "\trankdir = LR;"
        print "\tedge [arrowsize=0.4,fontsize=10]"
        print "\tnode1 [label=\"\",style=filled,fillcolor=lightgrey,shape=circle,width=.1,height=.1];"
        print "//------leaves------"
        self.print_gv_leaves( self.root )
        print "//------internal nodes------"
        self.print_gv_internal_nodes( self.root )
        print "//------edges------"
        self.print_gv_edges( self.root )
        #print "//------suffix links------"
        #self.print_gv_suffix_links( self.root )
        print "}"

    def print_edges( self, x ):
        for key, child in self.nodes[ x ].get_edges().iteritems():
            print self.edge_string( child )
            self.print_edges( child )

    def print_gv_internal_nodes( self, x ):
        if x != self.root and len( self.nodes[ x ].get_edges() ) > 0:
            print "\tnode" + str( x ) + " [label=\"\",style=filled,fillcolor=lightgrey,shape=circle,width=.07,height=.07]"
        for key, child in self.nodes[ x ].get_edges().iteritems():
            self.print_gv_internal_nodes( child )

    def print_gv_suffix_links( self, x ):
        if self.nodes[ x ].get_link() > 0:
            print "\tnode" + str( x ) + " -> node" + str( self.nodes[ x ].get_link() ) + " [label=\"\",weight=1,style=dotted]"
        for key, child in self.nodes[ x ].get_edges().iteritems():
            self.print_gv_suffix_links( child )

    def print_gv_edges( self, x ):
        for key, child in self.nodes[ x ].get_edges().iteritems():
            print "\tnode" + str( x ) + " -> node" + str( child ) + " [label=\"" + self.edge_string( child ) + "\",weight=3]"
            self.print_gv_edges( child )

    def print_gv_leaves( self, x ):
        if len( self.nodes[ x ].get_edges() ) == 0:
            print "\tnode" + str( x ) + " [label=\"\",shape=point]"
        else:
            for key, child in self.nodes[ x ].get_edges().iteritems():
                self.print_gv_leaves( child )

if __name__ == "__main__":
    with open( sys.argv[ 1 ] ) as fh:
        text = fh.next().strip()

    tree = SuffixTree( len( text ) )
    for char in text:
        tree.add_char( char )

    if len( sys.argv ) > 2:
        if sys.argv[ 2 ] == "gv":
            tree.print_graphviz_tree()
    else:
        tree.print_tree()
