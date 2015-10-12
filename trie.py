#! /usr/bin/env python

class Trie:

    nuc = ""
    children = {'A':None, 'C':None, 'G':None, 'T':None}
    count = 0

    def __init__(self, nuc):
        self.nuc = nuc
        self.count = 0

    def getCount(self):
        return self.count

    def incrementCount(self):
        self.count += 1

    def addSeq(self, seq):
        if self.children[seq[0]] == None:
            self.children[seq[0]] = Trie(seq[0])

        if len(seq) == 1:
            self.children[seq].incrementCount()
            return True

        _seq = seq[1:len(seq)]
        return self.addSeq(_seq)

    def __getitem__(self, nuc):
        return self.children[nuc]

    def toString(self):
        main_str = set()
        temp_str = ""
        return self.toStringRec(main_str, temp_str)

    def toStringRec(self, main_str, temp_str):
        for nuc, trie in self.children.iteritems():
            if trie == None:
                continue
            temp_str += nuc
            if (trie[nuc].getCount() >= 1):
                main_str.add(temp_str)

            trie.toStringRec(main_str, temp_str)
            del(temp_str[-1])

        return main_str
