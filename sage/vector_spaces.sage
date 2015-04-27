# Stable implementation of vector spaces
########################################

class VectorSpace:
    def __init__(self, base, dim):
        self._base = base
        self._ambiant_dim = dim
        self._gens = [ ]
        self._indices = range(dim)

    def __repr__(self):
        s = "Basis:\n"
        gens = MatrixSpace(self._base, self.dimension(), self.ambiant_dimension())(self._gens)
        s += gens.str()
        return s

    def copy(self):
        other = VectorSpace(self._base, self._ambiant_dim)
        other._gens = list(self._gens)
        other._indices = list(self._indices)
        return other

    def gens(self):
        return self._gens

    def ambiant_dimension(self):
        return self._ambiant_dim

    def dimension(self):
        return len(self._gens)

    def orthogonal(self):
        K = self._base
        dim = self.ambiant_dimension()
        d = dim - self.dimension()
        VS = FreeModule(K, dim)

        indices = list(self._indices)
        indices.reverse()

        gens = [ ]
        for i in range(d):
            index = indices[i]
            coords = [ K(0) ] * dim
            coords[index] = K(1)
            for j in range(d,dim):
                coords[indices[j]] = -self._gens[dim-1-j][index]
            vector = VS(coords)
            gens.append(vector)

        orth = VectorSpace(K, dim)
        orth._indices = indices
        orth._gens = gens
        return orth

    def __add__(self,other):
        ans = self.copy()
        ans._add_in_place(other.gens())
        return ans

    def _add_in_place(self, vectors):
        K = self._base
        gens = self._gens
        indices = self._indices
        dim = len(gens)
        ambiant_dim = self._ambiant_dim
        if not isinstance(vectors, list):
           vectors = [ vectors ]
        for vector in vectors:
            for i in range(dim):
                vector -= vector[indices[i]] * gens[i]
            val = Infinity; ii = 0
            for i in range(dim, ambiant_dim):
                coeff = vector[indices[i]]
                if coeff == 0: continue
                v = coeff.valuation()
                if val > v:
                    val = v; ii = i
            if val is Infinity: continue
            index = indices[ii]
            indices[dim], indices[ii] = index, indices[dim]
            vector /= vector[index]
            for i in range(dim):
                gens[i] -= gens[i][index] * vector
            gens.append(vector)
            dim += 1

    def intersection(self,other):
        V = self.orthogonal() + other.orthogonal()
        return V.orthogonal()

    def meet(self,other):
        return self.intersection(other)

    def apply(self,f):
        vectors = [ ]
        for vector in self._gens:
            vectors = f(vector)
        ans = VectorSpace(self._base, self._ambiant_dim)
        ans._add_in_place(vectors)
        return ans



# Tests and statistics
######################

def val_vector(x):
    return min([ c.valuation() for c in x.list() ])
def prec_vector(x):
    return min([ c.precision_absolute() for c in x.list() ])

def stat_vectorspace(p=2, prec=100, iter=50, repeat=1000):
    K = Qp(p, prec=3*prec)
    Kint = K.integer_ring()
    V = K^3; Vint = Kint^3

    HOM = End(V)
    MS = MatrixSpace(Kint,3)

    theory = [ ]; theory2 = [ ]; practice = [ ]
    for no in range(repeat):
        L = VectorSpace(K,3)
        zero = K(0, 2*prec)
        one = K(1, 2*prec)
        delta = K(2^prec, 2*prec)
        L._add_in_place(V([one,zero,zero]))
        L1 = VectorSpace(K,3)
        L1._add_in_place(V([one,delta,zero]))
        L2 = VectorSpace(K,3)
        L2._add_in_place(V([one,zero,delta]))
        for i in range(iter):

            f1 = HOM(MS.random_element())
            f2 = HOM(MS.random_element())
            f3 = HOM(MS.random_element())
            f4 = HOM(MS.random_element())

            H1 = L.apply(f1) + L.apply(f2)
            H2 = L.apply(f3) + L.apply(f4)
            L = H1.intersection(H2)
            if L.dimension() != 1: break

            H1 = L1.apply(f1) + L1.apply(f2)
            H2 = L1.apply(f3) + L1.apply(f4)
            L1 = H1.intersection(H2)

            H1 = L2.apply(f1) + L2.apply(f2)
            H2 = L2.apply(f3) + L2.apply(f4)
            L2 = H1.intersection(H2)

        indices = L._indices
        v = L._gens[0]
        v1 = v - L1._gens[0]
        v2 = v - L2._gens[0]

        M = matrix(2,2,[v1[indices[1]], v1[indices[2]], v2[indices[1]], v2[indices[2]]])
        print [ x.valuation() for x in M.list() ]
        print RR(M.determinant().valuation()/2)
        theory.append( floor(M.determinant().valuation()/2) - prec )
        theory2.append( min(v1[indices[1]].valuation(), v2[indices[1]].valuation()) - prec )
        practice.append( prec_vector(v) - 2*prec )

        if no % 50 == 0: print no

    pts_theory = [ ]
    for i in range(min(theory), max(theory)+1):
        pts_theory.append((i, theory.count(i)))
    pts_theory2 = [ ]
    for i in range(min(theory2), max(theory2)+1):
        pts_theory2.append((i, theory2.count(i)))
    pts_practice = [ ]
    for i in range(min(practice), max(practice)+1):
        pts_practice.append((i, practice.count(i)))
    return pts_theory, pts_theory2, pts_practice
