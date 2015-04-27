from collections import defaultdict
import random

def comatrix(A):
    M = A.parent()
    d = A.det()
    if not d:
        return A.parent().zero()
    n = A.ncols()
    R = M.base_ring()
    K = R.fraction_field()
    B = A.change_ring(K).augment(M.one() * d)
    return B.echelon_form().change_ring(R).matrix_from_columns(range(n,2*n))

def prec_gens(A):
    """
    INPUT:

    - ``A`` -- a square matrix over Z

    OUTPUT:

    - the index of the precision increase inside Z^n
    - a sorted list of generators for the precision increase
    """
    R = A.parent().base_ring()
    Rx = PolynomialRing(R,'x')
    x = Rx.gen()
    n = A.ncols()
    B = x*identity_matrix(n) - A
    zero = R(0)
    C = comatrix(B)
    def padded(L):
        if len(L) < n:
            return L + [zero] * (n - len(L))
        else:
            return L[:n]
    gen_mat = matrix(R,[padded(c.list()) for c in C.list()]).transpose().change_ring(R)
    D, U, V = gen_mat.smith_form()
    rgen_mat = (~U)*D
    index = prod([D[i,i] for i in range(n)])
    return index, sorted([Rx(list(rgen_mat.column(c))) for c in range(n)])


class NonmaximalExamples(SageObject):
    def __init__(self):
        self.trials = defaultdict(int)
        self.examples = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        self.hodge_examples = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        self.hodge_stats = defaultdict(lambda: defaultdict(int))

    def add_examples(self, dim, p=2, k=10, num_examples=10, index_sought=None):
        """
        INPUT:

        - ``dim`` -- the dimension of the matrices to search for
        - ``p`` -- a prime
        - ``k`` -- ``p^k`` is the bound up to which integers are sampled uniformly in creating random matrices
        - ``num_examples`` -- the number of new examples to find
        - ``p`` -- a prime, given when searching for a particular index valuation
        """
        examples = self.examples[p]
        stats = self.hodge_stats[p]
        N = p^k
        for t in range(num_examples):
            A = random_matrix(ZZ,dim,x=N)
            if A.det() == 0: continue
            hodge_polygon = HP(A, p)
            hodge_vals = tuple([hodge_polygon[i] - hodge_polygon[i+1] for i in reversed(xrange(len(hodge_polygon)-1))])
            newton_polygon = NP(A,p)
            index, G = prec_gens(A)
            precision_polygon = CP(G,p)
            for i in range(len(precision_polygon)):
                try:
                    if precision_polygon[i] < hodge_polygon[i+1] or precision_polygon[i] > newton_polygon[i+1]:
                        raise RuntimeError("Matrix violates Hodge-Newton polygon conjecture!\n%s\n%s"%(A, G))
                except IndexError:
                    raise RuntimeError("Index Error!\nNP: %s\nHP: %s\nPP: %s\n%s\n%s"%(newton_polygon, hodge_polygon, precision_polygon, A, G))
            prec_hodge_diff = tuple(precision_polygon[i] - hodge_polygon[i+1] for i in range(len(precision_polygon)))
            prec_hodge_index = sum(prec_hodge_diff)
            lattice_prec_index = index.valuation(p) - sum(precision_polygon[i] for i in range(len(precision_polygon)))
            examples[hodge_vals][prec_hodge_index, lattice_prec_index].append((A,G))
            stats[prec_hodge_index, lattice_prec_index] += 1
    def add_examples_manydim(self, dims, p=2, k=10, iterations=500):
        for t in range(iterations):
            print "Iteration %s"%t
            for d in dims:
                self.add_examples(d, p, k)

    def add_smith_examples(self, dim, p, max_val=40, num_examples=10):
        trials = 0r
        found = 0r
        examples = self.examples[dim,0]
        while found < num_examples:
            vals = []
            while len(vals) < dim:
                new_val = ZZ.random_element()
                if 0 <= new_val <= max_val:
                    vals.append(p^new_val)
            D = diagonal_matrix(sorted(vals))
            U = random_matrix(ZZ,dim,algorithm="unimodular")
            V = random_matrix(ZZ,dim,algorithm="unimodular")
            A = U * D * V
            index, G = prec_gens(A)
            if index != 1r:
                examples[index].append((A,G))
                found += 1r
            trials += 1r
        self.trials[dim,0] += trials

    def add_hodge_examples(self, hodge_val_vector, p, num_examples=10, k=10, powers_of_p=False):
        """
        INPUTS:

        - ``hodge_val_vector`` -- a tuple with the list of valuations
                                  of diagonal entries in the Smith
                                  normal form of the matrices added.

        - ``p`` -- a prime

        - ``num_examples`` -- the number of examples to be added.

        - ``k`` -- an exponent of p modulo which we try to distribute
                   the entries of the matrix.

        - ``powers_of_p`` -- a boolean: whether to use powers of p or
                             to incorporate random, prime-to-p
                             factors.
        """
        hodge_val_vector = tuple(sorted(list(hodge_val_vector)))
        examples = self.hodge_examples[p][hodge_val_vector]
        dim = len(hodge_val_vector)
        det_valuation = sum(hodge_val_vector)
        hodge_polygon = [det_valuation]
        for v in reversed(hodge_val_vector):
            hodge_polygon.append(hodge_polygon[-1] - v)
        for example_num in range(num_examples):
            if not powers_of_p:
                # Since we're working over Z, to range over all matrices
                # with the given valuations in the diagonal part of the
                # smith decomposition, we need to allow p-adic units along
                # with just the powers of p.

                # This is biased a bit toward smaller elements, but that's okay.
                prime_to_p = ZZ.random_element(1,p^k).val_unit(p)[1]

                F = prime_to_p.factor()
                prime_to_p_vec = [1] * dim
                for q, e in F:
                    qexps = sorted([ZZ.random_element(e+1) for _ in range(dim)])
                    for i in range(dim):
                        prime_to_p_vec[i] *= q^qexps[i]
                D = diagonal_matrix([p^hodge_val_vector[i] * prime_to_p_vec[i] for i in range(dim)])
            else:
                D = diagonal_matrix([p^hodge_val_vector[i] for i in range(dim)])
            U = random_matrix(ZZ,dim,algorithm="unimodular")
            V = random_matrix(ZZ,dim,algorithm="unimodular")
            A = U*D*V
            newton_polygon = NP(A,p)
            index, G = prec_gens(A)
            precision_polygon = CP(G,p)
            for i in range(len(precision_polygon)):
                try:
                    if precision_polygon[i] < hodge_polygon[i+1] or precision_polygon[i] > newton_polygon[i+1]:
                        raise RuntimeError("Matrix violates Hodge-Newton polygon conjecture!\n%s\n%s"%(A, G))
                except IndexError:
                    raise RuntimeError("Index Error!\nNP: %s\nHP: %s\nPP: %s\n%s\n%s"%(newton_polygon, hodge_polygon, precision_polygon, A, G))
            prec_hodge_diff = tuple(precision_polygon[i] - hodge_polygon[i+1] for i in range(len(precision_polygon)))
            prec_hodge_index = sum(prec_hodge_diff)
            lattice_prec_index = index.valuation(p) - sum(precision_polygon[i] for i in range(len(precision_polygon)))
            examples[prec_hodge_index, lattice_prec_index].append((A,G))
            self.hodge_stats[p][prec_hodge_index, lattice_prec_index] += 1
    def add_hodge_examples_manyvec(self, dim, p, max_sum = None, iterations = 20):
        if max_sum is None: max_sum = dim + 4
        vecs = []
        for sum in range(max_sum+1):
            I = IntegerVectorsModPermutationGroup(SymmetricGroup(dim-1),sum=sum)
            for v in I:
                vecs.append(tuple([0] + list(reversed(v))))
        for t in range(iterations):
            if iterations > 1:
                print "Iteration %s"%t
            for i, v in enumerate(vecs):
                self.add_hodge_examples(v, p)
                if i%20 == 19:
                    print "    Hodge vector %s/%s"%(i+1,len(vecs))
    def add_hodge_examples_manydim(self, p, max_dim, max_sum = None, iterations = 20):
        for t in range(iterations):
            print "Iteration %s"%t
            for dim in range(3, max_dim+1):
                print "  Dimension %s"%dim
                self.add_hodge_examples_manyvec(dim, p, max_sum, 1)

    def _nontrivial(self, p, test_trivial, process_vec, filter_vec=None, hodge=False, show=True):
        if filter_vec is None:
            filter_vec = lambda vec: True
        nontrivial = defaultdict(int)
        trivial = defaultdict(int)
        if hodge:
            D = self.hodge_examples
        else:
            D = self.examples
        for vec, examples in D[p].iteritems():
            for indices, L in examples.iteritems():
                if filter_vec(vec):
                    if test_trivial(indices):
                        trivial[process_vec(vec)] += len(L)
                    else:
                        nontrivial[process_vec(vec)] += len(L)
        ratio = {}
        for measure in nontrivial.iterkeys():
            ratio[measure] = RR(nontrivial[measure]) / RR(nontrivial[measure] + trivial[measure])
        if show:
            for measure in sorted(ratio.keys()):
                print measure, ratio[measure]
            print "R value: %s"%(r_test(ratio.items()))
        else:
            return ratio

    def nontrivial_by_dimension(self, p, show=True):
        return self._nontrivial(p, lambda ind: ind[0] == 0 and ind[1] == 0, len, show=show)

    def nontrivial_by_sum(self, p, dim, show=True):
        return self._nontrivial(p, lambda ind: ind[0] == 0 and ind[1] == 0, sum, lambda vec: len(vec) == dim, show)

    def prec_above_hodge_by_dimension(self, p, show=True):
        return self._nontrivial(p, lambda ind: ind[0] == 0, len, show=show)

    def prec_above_hodge_by_sum(self, p, dim, show=True):
        return self._nontrivial(p, lambda ind: ind[0] == 0, sum, lambda vec: len(vec) == dim, show)

    def extra_lattice_prec_by_dimension(self, p, show=True):
        return self._nontrivial(p, lambda ind: ind[1] == 0, len, show=show)

    def extra_lattice_prec_by_sum(self, p, dim, show=True):
        return self._nontrivial(p, lambda ind: ind[1] == 0, sum, lambda vec: len(vec) == dim, show)

    def sample_Gs(self, p, n, hodge_val_vector=None, diagonal=None, hodge=False):
        Gs = []
        if hodge:
            H = self.hodge_examples[p]
        else:
            H = self.examples[p]
        hodge_vecs = H.keys()
        looping = 0
        while len(Gs) < n:
            if hodge_val_vector is None:
                vec = random.choice(hodge_vecs)
            else:
                vec = hodge_val_vector
            L = random.choice(H[vec].values())
            A, G = random.choice(L)
            if self._match_diag(diagonal,(A,G)):
                Gs.append(G)
            looping += 1
            if looping > 10*n:
                print "Ending search before %s examples found."%n
                break
        return Gs
    def all_exceeds(self,hodge=False):
        S = set([])
        if hodge:
            H = self.hodge_examples
        else:
            H = self.examples
        for examples in H.itervalues():
            for D in examples.itervalues():
                for exceed in D.iterkeys():
                    S.add(exceed)
        return sorted(list(S))
    def exceed_summary(self,p,dim=None,hodge=False):
        if hodge:
            H = self.hodge_examples[p]
        else:
            H = self.examples[p]
        to_print = []
        for vec in sorted(H.iterkeys()):
            if not self._match_dim(dim, len(vec)): continue
            D = H[vec]
            results = sorted([(len(v), k) for k, v in D.iteritems()])
            results = [("%s,%s"%(k[0],k[1]),str(v)) for v, k in reversed(results)]
            to_print.append([str(vec)] + results)
        col_width = [max(len(row[0]) for row in to_print)]
        col = 1
        while True:
            max_width_exceeds = -1
            max_width_counts = -1
            for row in to_print:
                if len(row) > col:
                    max_width_exceeds = max(len(row[col][0]), max_width_exceeds)
                    max_width_counts = max(len(row[col][1]), max_width_counts)
            if max_width_exceeds < 0: break
            col_width.append((max_width_exceeds, max_width_counts))
            col += 1
        for row in to_print:
            row_str = row[0] + " "*(1 + col_width[0] - len(row[0]))
            row_list = []
            for col, pair in enumerate(row):
                if col == 0: continue
                row_list.append(" " * (col_width[col][0] - len(pair[0])) + pair[0] + ":" + pair[1] + " " * (col_width[col][1] - len(pair[1])))
            print row_str + "-- " + " - ".join(row_list)

    @staticmethod
    def _match_dim(dim, d):
        return dim is None or d == dim
    @staticmethod
    def _match_N(p, N, strict):
        return N != 0 and (not strict or N.is_power_of(p))

    def index_stats(self, p, dim=None, strict=False, percentage=False):
        def match(d, N):
            return self._match_dim(dim, d) and self._match_N(p, N, strict)
        index_vals = defaultdict(int)
        total = 0r
        for key, examples in self.examples.iteritems():
            if match(*key):
                T = self.trials[key]
                total += T
                for index in examples.iterkeys():
                    index_vals[index.valuation(p)] += 1r
                    T -= 1
                # The unstored examples have index 1
                index_vals[0] += T
        if not index_vals: return []
        max_val = max(index_vals.keys())
        if percentage:
            total = RR(total)
            return [RR(index_vals[v]) / total for v in range(max_val+1)]
        else:
            return [index_vals[v] for v in range(max_val+1)]

    @staticmethod
    def _match_diag(need_diag, example):
        if need_diag is None: return True
        actually_diag = all(g.is_term() for g in example[1])
        return (need_diag and actually_diag) or (not need_diag and not actually_diag)
    @staticmethod
    def _match_index(p, index_required, actual_index):
        if index_required is None: return True
        return index_required.valuation(p) == actual_index.valuation(p)

    def find_example(self, p, dim, index=None, diagonal=None, all=False, k=10, max_extra_examples=1000):
        L = []
        for key, examples in self.examples.iteritems():
            if self._match_dim(dim, key[0]):
                for ind, exs in examples.iteritems():
                    if self._match_index(p, index, ind):
                        for ex in exs:
                            if self._match_diag(diagonal,ex):
                                if all:
                                    L.append(ex)
                                else:
                                    return ex
        if not all:
            if dim is not None:
                ex = self.add_examples(dim, p^k, max_extra_examples, p, index)
                if ex is not None:
                    return ex
            raise RuntimeError("No example found")
        return L

    def all_diagonal_examples(self):
        L = []
        for key, examples in self.examples.iteritems():
            for ind, exs in examples.iteritems():
                for ex in exs:
                    if self._match_diag(True, ex):
                        L.append(ex)
        return L

    def all_nondiagonal_examples(self):
        L = []
        for key, examples in self.examples.iteritems():
            for ind, exs in examples.iteritems():
                for ex in exs:
                    if self._match_diag(False, ex):
                        L.append(ex)
        return L

    def all_smith_examples(self):
        L = []
        for key, examples in self.examples.iteritems():
            if key[1] == 0:
                for ind, exs in examples.iteritems():
                    L.extend(exs)
        return L

def update(NME):
    NME2 = NonmaximalExamples()
    NME2.examples = NME.examples
    NME2.trials = NME.trials
    NME2.hodge_examples = NME.hodge_examples
    NME2.hodge_stats = NME.hodge_stats
    return NME2

def minval(A, p):
    return min([c.valuation(p) for c in A.list()])

def convexify(L):
    P = Polyhedron(vertices = list(enumerate(L)), rays=[(0,1)])
    verts = sorted(P.vertices_list())
    heights = [QQ(verts[0][1])]
    for i in range(len(verts)-1):
        cur = verts[i]
        next = verts[i+1]
        slope = (next[1] - cur[1])/(next[0] - cur[0])
        for x in range(cur[0],next[0]):
            heights.append(heights[-1] + slope)
    return heights

def NP(A, p):
    f = A.charpoly()
    return convexify([a.valuation(p) for a in f])

def newton_prec_diff(examples, p=2):
    diffs = set([])
    for D in examples.itervalues():
        for L in D.itervalues():
            for A, G in L:
                cp = CP(G, p)
                np = NP(A, p)
                diff = [np[i+1] - cp[i] for i in range(len(cp))]
                diffs.add(tuple(diff))
    return diffs

def HP(A, p):
    D, U, V = A.smith_form()
    n = min(D.nrows(), D.ncols())
    vals = sorted([D[i,i].valuation(p) for i in range(n)])
    S = sum(vals)
    HP = [S]
    for v in reversed(vals):
        HP.append(HP[-1] - v)
    return HP

def CP(G, p):
    vals = [min([g[i].valuation(p) for g in G if g[i]]) for i in range(len(G))]
    return convexify(vals)

def test_between(L, p):
    for A, G in L:
        np = NP(A, p)
        hp = HP(A, p)
        cp = CP(G, p)
        for i in range(len(cp)):
            if cp[i] < hp[i+1] or cp[i] > np[i+1]:
                print A
                print G
                raise RuntimeError

def HPexceed(L, p):
    exceed = defaultdict(int)
    for A, G in L:
        hp = HP(A, p)
        cp = CP(G, p)
        diff = tuple(cp[i] - hp[i+1] for i in range(len(cp)))
        exceed[diff] += 1
    return exceed

def HPexceed_display(L, p):
    exceed = HPexceed(L, p)
    excesses = sorted(exceed.keys(), cmp=lenlex)
    for key in excesses:
        print key, exceed[key]

def HPexceed_dictdisplay(examples, p):
    example_list = []
    for vecD in examples[p].itervalues():
        for L in vecD.itervalues():
            example_list.extend(L)
    HPexceed_display(example_list, p)

def lenlex(x, y):
    # For sorting difference tuples
    c = cmp(len(x),len(y))
    if c: return c
    return cmp(x, y)

def print_polygons(ex, p):
    A, G = ex
    print "NP: ", NP(A,p)
    print "HP: ", HP(A,p)
    print "CP:    ", CP(G,p)

def r_test(pairs):
    A = [a for a,b in pairs]
    B = [b for a,b in pairs]
    cross = [a*b for a,b in pairs]
    A2 = [a*a for a,b in pairs]
    B2 = [b*b for a,b in pairs]
    N = len(pairs)
    return RR(N*sum(cross) - sum(A)*sum(B)) / RR((N*sum(A2) - sum(A)^2)*(N*sum(B2) - sum(B)^2)).sqrt()
