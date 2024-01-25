import sympy

def sympify(data):
    edges = []
    for i in range(data.c):
        edges.append([])
        for ed in data.edges[i]:
            a = sympy.Integer(ed[1][0])
            b = sympy.Integer(ed[1][1])
            edges[i].append((ed[0],(a,b)))
    data.edges = edges
    pos = []
    for i in range(data.c):
        a = sympy.Rational(data.pos[i][0]).limit_denominator(100000)
        b = sympy.Rational(data.pos[i][1]).limit_denominator(100000)
        pos.append((a,b))
    data.pos = pos
    return data

def is_undirected(data):
    for i in range(data.c):
        for ed in data.edges[i]:
            j = ed[0]
            ret = False
            for ed_ in data.edges[ed[0]]:
                if ed_[0]==i:
                    ok = True
                    for k in range(data.dim):
                        if ed[1][k]+ed_[1][k]!=0: ok = False
                    if ok:
                        ret = True
                        break
            if not ret: return False
    return True
