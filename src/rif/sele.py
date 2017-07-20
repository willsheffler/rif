from parsimonious.grammar import Grammar, NodeVisitor

atom_name_grammar = Grammar(
    """
    expr = (anameORP ' ')* anameORP
    anameORP = anameOR / ('(' anameOR ')')
    anameOR = aname ((' OR ' / ' or ') aname)*
    aname    = ~'[a-zA-Z0-9]+'i
    """
)


def print_node(node):
    print(node)
    for m in 'expr_name match start text'.split():
        if hasattr(node, m):
            print('   ', m, getattr(node, m))


class AtomSeleVisitor(NodeVisitor):
    def __init__(self):
        self.atom_names = list()

    def visit_anameOR(self, node, visited_children):
        anames = [node.children[0].text.upper()]
        for c in node.children[1].children:
            anames.append(c.children[1].text.upper())
        self.atom_names.append(tuple(anames))

    def generic_visit(self, node, visited_children):
        pass


def parse_atom_names(sele):
    parse = atom_name_grammar.parse(sele)
    visitor = AtomSeleVisitor()
    visitor.visit(parse)
    return visitor.atom_names
