from parsimonious.grammar import Grammar, NodeVisitor

# todo: figure out how to properly use parsimonious
sele_grammar = Grammar(
    """
    expr = (obj ' ')* obj
    obj = ray / anameORP
    ray = anameORP '->' anameORP
    anameORP = anameOR / ('(' anameOR ')')
    anameOR = aname ((' OR ' / ' or ') aname)*
    aname    = ~'[a-zA-Z0-9]+'
    """
)


def print_node(node):
    print(node)
    for m in 'expr_name match start text'.split():
        if hasattr(node, m):
            print('   ', m, getattr(node, m))


def crappy_visit_anameOR(node):
    anames = [node.children[0].text.upper()]
    for c in node.children[1].children:
        anames.append(c.children[1].text.upper())
    return tuple(anames)


class AtomSeleVisitor(NodeVisitor):
    def __init__(self):
        self.atom_names = list()

    def visit_anameOR(self, node, vchild):
        anames = crappy_visit_anameOR(node)
        self.atom_names.append(anames)

    def visit_ray(self, node, vchild):
        raise ValueError('parsing ray with AtomSeleVisitor')

    def generic_visit(self, node, vchild):
        pass


def parse_atom_names(sele):
    parse = sele_grammar.parse(sele)
    visitor = AtomSeleVisitor()
    visitor.visit(parse)
    return visitor.atom_names


class RaySeleVisitor(NodeVisitor):
    def __init__(self):
        self.rays_atom_names = list()

    def visit_ray(self, node, vchild):
        ray = list()
        for c in (node.children[0], node.children[2]):
            ray.append(crappy_visit_anameOR(c.children[0]))
        self.rays_atom_names.append(tuple(ray))

    def generic_visit(self, node, vchild):
        pass


def parse_ray_atom_names(sele):
    parse = sele_grammar.parse(sele)
    visitor = RaySeleVisitor()
    visitor.visit(parse)
    return visitor.rays_atom_names
