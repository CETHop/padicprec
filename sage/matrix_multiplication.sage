p = 2
prec = 1000; d = 2
nb = 500
repeat = 100

R = Zp(p, prec=2*prec)
MS = MatrixSpace(R,d)

def prec_matrix(M):
    return min([ x.precision_absolute() for x in M.list() ])

Basis = [ ]
for i in range(d):
    for j in range(d):
        M = MS(0)
        M[i,j] = 1
        Basis.append(M)

theory = theory2 = valuation = practice = 0
for _ in range(repeat):
    Lattice = [ ]
    for i in range(d):
        for j in range(d):
            M = MS(0)
            M[i,j] = 1
            Lattice.append(M)

    M = MS(1)
    for c in range(nb):
        F = MS([ R.random_element().add_bigoh(prec) for _ in range(d^2) ])
        # dM * F + M * dF
        Lat = [ L*F for L in Lattice ] + [ M*L for L in Basis ]
        M = M * F
        #precM = prec_matrix(M)
        #M = MS([ x.add_bigoh(precM) for x in M.list() ])
        Lattice = [ ]
        for i in range(d):
            for j in range(d):
                val = Infinity; index = -1
                for ind in range(len(Lat)):
                    v = Lat[ind][i,j].valuation()
                    if v < val:
                        val = v; index = ind
                if index == -1: raise RuntimeError
                gen = Lat[index]
                Lattice.append(gen)
                del Lat[index]
                if c == nb-1:
                    theory += gen[i,j].valuation()
                for ind in range(len(Lat)):
                    scalar = Lat[ind][i,j] // gen[i,j]
                    scalar = scalar.lift_to_precision(prec)
                    Lat[ind][i,j] -= scalar * gen[i,j]

    #practice += min([ x.precision_absolute() for x in M.list() ]) - prec
    valuation += M[0,0].valuation()
    practice += M[0,0].precision_absolute() - prec
    #theory2 += Lattice[0][0,0].valuation()
    theory2 += min([ min([ x.valuation() for x in L.list() ]) for L in Lattice ])

theory /= d^2
valuation = 0
print "Theory:", RR((theory-valuation)/repeat)
print "Theory 2:", RR((theory2-valuation)/repeat)
print "Practice:", RR((practice-valuation)/repeat)
