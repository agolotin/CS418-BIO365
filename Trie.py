#! /usr/bin/env python

class Trie(object):

    def __init__(self, words):

        # Function to create a node in a trie
        self.create_node = lambda p, d: {'parent':p, 'children':[], 'depth':d, 'end': False}

        self.nodes = {1:self.create_node(0,0)}
        self.edges = {}

        if type(words) is str:
            self.add_word(word)
        else:
            for w in words:
                self.add_word(w)
                
                
                
    def add_word(self, word):
       parent_node, insertion_substr = self.insertion_point(word, 1)

       for i in xrange(len(insertion_substr)):
           new_node_key = len(self.nodes) + 1

           self.nodes[new_node_key] = self.create_node(parent_node, self.nodes[parent_node]['depth']+1)
           self.nodes[parent_node]['children'].append(new_node_key)

           self.edges[parent_node, new_node_key] = insertion_substr[i]

           parent_node = new_node_key

       self.nodes[parent_node]['end'] = True


    def insertion_point(self, word_to_add, current_node):

        if word_to_add == '':
            return current_node, word_to_add

        for child_node in self.nodes[current_node]['children']:
            if self.edges[current_node, child_node] == word_to_add[0]:
                return self.insertion_point(word_to_add[1:], child_node)

        return current_node, word_to_add


