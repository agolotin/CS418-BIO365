import sys
from collections import defaultdict

with open(sys.argv[1]) as fd:
    text = fd.readline().strip()
    sa = fd.readline().strip()

#BWT from suffix array
_map = { i:text[i] for i in xrange(len(text)) }
print _map
text = ""
for num in sa.split(", "):
    if int(num) == 0:
        text += "$"
    else:
        text += _map[int(num)-1]
print text
