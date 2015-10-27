#! /usr/bin/env python

class Trie(object):

    def __init__(self, words):

        # Function to create a node in a trie
        self.create_node = lambda p, d: {'parent':p, 'children':[], 'depth':d, 'end': False}

        self.nodes = {1:self.create_node(0,0)}
        self.edges = {}

        # I just wanted to make it cool and allow for both a string and a list to be passed into the init argument
        if type(words) is str:
            self.add_word(word)
        else:
            for w in words:
                self.add_word(w)
                
                

                
    def add_word(self, word):
        # get the number of the parent node and the substring to insert
        # since part of the trie may already contain part of the word we want to insert
        #    we only want to get the substring that we have to add to the trie
       parent_node, insertion_substr = self.insertion_point(word, 1)

       for i in xrange(len(insertion_substr)):
           new_node_key = len(self.nodes) + 1

           self.nodes[new_node_key] = self.create_node(parent_node, self.nodes[parent_node]['depth']+1)
           self.nodes[parent_node]['children'].append(new_node_key)

           # Every edge in a dictionary is stroed as a key value pair, where key is a tuple comprised of the 
           #    parent node number and the node number of the value. Value is the actual nucleotide
           self.edges[parent_node, new_node_key] = insertion_substr[i]

           parent_node = new_node_key

       self.nodes[parent_node]['end'] = True




    # Traverses the trie to find where to insert the next new node
    # Function is used by add_word()
    def insertion_point(self, word_to_add, current_node):

        if word_to_add == '':
            return current_node, word_to_add

        for child_node in self.nodes[current_node]['children']:
            if self.edges[current_node, child_node] == word_to_add[0]:
                return self.insertion_point(word_to_add[1:], child_node)

        return current_node, word_to_add



    # Indentifies whether a prefix to a dna sequence is in the trie
    # Function is used only in trie matching problem(problem 26)
    # @param word is the dna sequence
    def prefix_in_trie(self, word, current_node):

        if self.nodes[current_node]['end'] == True:
            # If we reach the end of the trie we know that the suffix is in the trie
            return True
        elif word == '':
            return False;

        # Recursively traverse the trie and propogate the boolean up
        for child_node in self.nodes[current_node]['children']:
            if self.edges[current_node, child_node] == word[0]:
                return self.prefix_in_trie(word[1:], child_node)


        return False
