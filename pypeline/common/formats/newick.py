import re
import collections

class NewickError(RuntimeError):
    pass

class NewickParseError(NewickError):
    pass

class Newick(tuple):
    def __new__(self, name = None, length = None, children = None):
        if (children is None):
            if not (name or length):
                raise NewickError("Leaf nodes MUST have either a name or a length")
        elif not children:
            raise NewickError("Empty internal nodes now allowed")
        return tuple.__new__(self, (name, length, children))

    @property
    def name(self):
        return self[0]

    @property
    def length(self):
        return self[1]

    @property
    def children(self):
        return self[2]

    @property
    def is_leaf(self):
        return self.children is None

    def get_leafs(self):
        if not self.is_leaf:
            for child in self.children:
                for leaf in child.get_leafs():
                    yield leaf
        else:
            yield self


    @classmethod
    def from_string(cls, string):
        tokens = _tokenize(string)
        if tokens and tokens[0] == "(":
            top_node = _parse_tokens(tokens)
        else:
            top_node = _parse_child(tokens)

        if tokens != [";"]:
            raise NewickParseError("Missing terminating semi-colon")

        return top_node

    def __repr__(self):
        return "%s;" % (self._to_str(),)

    def _to_str(self):
        fields = []
        if self.children:
            fields.append("(")
            for child in self.children:
                fields.append(child._to_str()) # pylint: disable=W0212
                fields.append(",")
            fields.pop()
            fields.append(")")
        if self.name is not None:
            fields.append(self.name)
        if self.length is not None:
            fields.append(":")
            fields.append(self.length)
        return "".join(fields)



################################################################################
################################################################################
## Functions related to NEWICK parsing

_TOKENIZER = re.compile("([():,;])")
_NODE_KEYS = frozenset(("name", "length", "children"))


def _tokenize(string):
    result = []
    for field in _TOKENIZER.split(string):
        field = field.strip()
        if field:
            result.append(field)
    return result


def _parse_tokens(tokens):
    assert tokens and tokens[0] == "("

    tokens.pop(0)
    child, children = None, []
    while tokens and (tokens[0] not in ");"):
        if tokens[0] == ",":
            children.append(child)
            tokens.pop(0)
        child = _parse_child(tokens)
    children.append(child)

    if any(child is None for child in children):
        raise NewickParseError("Implicit leaf nodes (no name OR length) are not allowed")
    elif not tokens or (tokens[0] != ")"):
        raise NewickParseError("Malformed Newick string, contains unbalanced parantheses")
    tokens.pop(0)

    return _parse_child(tokens, children = children)


def _parse_child(tokens, children = None):
    if tokens and tokens[0] == "(":
        return _parse_tokens(tokens)

    name, length = None, None
    while tokens and (tokens[0] not in ",);"):
        if (tokens[0] == ":"):
            if length is not None:
                raise NewickParseError("Node has multiple length values")
            tokens.pop(0)
            if tokens[0] in ",);":
                raise NewickParseError("Missing length value")
            length = tokens.pop(0).strip()
        else:
            name = tokens.pop(0).strip()

    return Newick(name     = name,
                  length   = length,
                  children = children)
