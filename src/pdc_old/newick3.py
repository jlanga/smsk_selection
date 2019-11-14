import string, sys
from shlex import shlex
from phylo3 import Node
import StringIO


class Tokenizer(shlex):
    """Provides tokens for parsing Newick-format trees"""
    def __init__(self, infile):
        shlex.__init__(self, infile)
        self.commenters = ''
        self.wordchars = self.wordchars+'-.@'
        self.quotes = "'"

    def parse_comment(self):
        while 1:
            token = self.get_token()
            if token == '':
                sys.stdout.write('EOF encountered mid-comment!\n')
                break
            elif token == ']':
                break
            elif token == '[':
                self.parse_comment()
            else:
                pass

def parse(input, ttable=None):
    """
    Parse a Newick-formatted tree description
    input is any file-like object that can be coerced into shlex,
    or a string (converted to StringIO)
    """
    if type(input) is str:
        input = StringIO.StringIO(input)
    
    start_pos = input.tell()
    tokens = Tokenizer(input)

    node = None; root = None
    lp=0; rp=0; rooted=1

    prev_tok = None

    while 1:
        token = tokens.get_token()
        #print token,
        if token == ';' or token == '':
            assert lp == rp, \
                   'unbalanced parentheses in tree description'
            break

        # internal node
        elif token == '(':
            lp = lp+1
            newnode = Node()
            newnode.istip = False
            if node:
                node.add_child(newnode)
            node = newnode

        elif token == ')':
            rp = rp+1
            node = node.parent
            
        elif token == ',':
##             if lp == rp:
            node = node.parent
            
        # branch length
        elif token == ':':
            token = tokens.get_token()

            if not (token == ''):
                try:
                    brlen = float(token)
                except ValueError:
                    raise "NewickError invalid literal for branch length, '%s'" % token
            else:
                raise "NewickError unexpected end-of-file (expecting branch length)"

            node.length = brlen
        # comment
        elif token == '[':
            tokens.parse_comment()

        # leaf node or internal node label
        else:
            if prev_tok != ')': # leaf node
                if ttable:
                    ttoken = ttable.get(token) or ttable.get(int(token))
                    if ttoken:
                        token = ttoken
                newnode = Node()
                newnode.label = token
                newnode.istip = True
                node.add_child(newnode)
                node = newnode
            else: # label
                # translation table for internal nodes labels?
                node.label = token

        prev_tok = token
        #print token, node
        

    input.seek(start_pos)

##     if rooted:
##         root = Fnode(isroot=1)
##         root.label = node.next.label; node.next.label = None
##         root.length = node.next.length; node.next.length = None
##         node.insert_fnode(root)

    #return root
    return node

def traverse(node):
    if node.istip: return node.back
    else: return node.next.back
        
def to_string(node, length_fmt=":%s"):
    if not node.istip:
        node_str = "(%s)%s" % \
                   (",".join([ to_string(child, length_fmt) \
                               for child in node.children ]),
                    node.label or ""
                    )
    else:
        node_str = "%s" % node.label

    if node.length is not None:
        length_str = length_fmt % node.length
        #length_str = ':%f' % node.length
        #length_str = str(length_str)
    else:
        length_str = ""

    s = "%s%s" % (node_str, length_str)
    return s

tostring = to_string
        
def parse_from_file(filename):
    if filename == '-':
        file = sys.stdin
    else:
        file = open(filename, 'r')
    content = string.strip(file.read())
    treedescs = string.split(content, ';')
    tree = parse(treedescs[0])
    file.close()
    return tree

if __name__ == "__main__":
    #import ascii
    s = "(a:3,(b:1e-05,c:1.3)int:5)root;"
    s = "(Nematostella_vectensis@146487:0.28570255971561625552,((Bargmannia@124461:0.11671371999896994198,(Nanomia@82407:0.02229037752330314398,(Frillagalma@103872:0.03352004060787890788,Agalma_elegans@38387:0.00723798292127085779):0.01418728724518653582):0.07522769713546764714):0.10356415483616157602,Clytia_hemisphaerica@165314:0.07583460439008869736):0.06930854841339750827,Hydra_magnipapillata@6480:0.08755522366843263016):0.0;"
    #s = "(a,b,c,d,e,f,g);"
    n = parse(s)
    print
    #print ascii.render(n)
    print(s)
    print(to_string(n))
    #print n.next.back.label
