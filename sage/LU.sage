#Computation of low(A) and up(A)#

def low(M):
    n = M.nrows()
    m = M.ncols()
    A = copy(M)
    for i in range(n):
        for j in range(i,m):
            A[i,j] = 0
    return A

def up(M):
    n = M.nrows()
    m = M.ncols()
    A = copy(M)
    for i in range(n):
        for j in range(i):
            A[i,j] = 0
    return A

def printval(M,p):
    n = M.nrows()
    m = M.ncols()
    Mval = MatrixSpace(QQ,n,m)(0)
    for i in range(n):
        for j in range(m):
            if M[i,j].valuation(p) ==+Infinity:
                Mval[i,j] = -42
            else:
                Mval[i,j] =M[i,j].valuation(p)
    print(Mval)
    return Mval            


def comparaison_LU(d,p,prec,nbiter,repeat):

    R = Zp(p, prec=2*prec)
    MS = MatrixSpace(R,d)

    MEX = MatrixSpace(QQ,d)

    #Basis for MEX

    Basis = [ ]
    listcouple = [ ]
    for ii in range(d):
        for jj in range(d):
            Mb = MEX(0)
            Mb[ii,jj] = 1
            Basis.append(Mb)
            listcouple.append([ii,jj])
    
    
    for qsd in range(repeat):
        listproj = [ ]
        listlatt = [ ]
        for wxc in range(nbiter):

            # Construction d'une matrice aleatoire #

            Mzp = MS.random_element()
            M = MEX(0)
            for i in range(d):
                for j in range(d):
                    M[i,j] = (Mzp[i,j]).lift()



            lu = M.LU(pivot ='nonzero')
            L = lu[1]
            U = lu[2]

            #Calcul de dM#

            MEX2 = MatrixSpace(QQ,d^2)

            DM = MEX2(0)


            for k in range(d^2):
                u = listcouple[k][0]
                v = listcouple[k][1]
                dM = Basis[k]
                meuv =L*low(L^(-1)*dM*U^(-1))+up(L^(-1)*dM*U^(-1))*U
                for kk in range(d^2):
                    uu = listcouple[kk][0]
                    vv = listcouple[kk][1]
                    DM[kk,k] = meuv[uu,vv]


            #Calcul de la valuation pour la projection#
            valproj = 0
            for l in range(d^2):
                valproj += min([DM[l,ll].valuation(p) for ll in range(d^2)])

            #Calcul de la perte de precision pour les lattices#
            vallatt = DM.det().valuation(p)


            #print("#############")
            #print("DM")
            #print(DM)
            #print("#############")
            #print("val de DM")
            #printval(DM,p)
            #print("#############")
            #print("valuation jagged")
            #print(valproj)
            #print("#############")
            #print("valuation lattice")
            #print(vallatt)
            listproj.append(valproj)
            listlatt.append(vallatt)
        print("################")
        print("################")
        print("################")
        print("################")
        print("moyenne proj")
        print(mean(listproj).numerical_approx())
        print("std proj")
        print(std(listproj).numerical_approx())
        print("################")
        print("moyenne lattice")
        print(mean(listlatt).numerical_approx())
        print("std latt")
        print(std(listlatt).numerical_approx())            

#choice of constants#

print("################")
print("################")
print("##### d = 2 ####")
print("################")
print("################")

p = 2
prec = 50; 
repeat = 2
d =2
nb = 2000

#print("nb de test")
#print(nb)
#comparaison_LU(d,p,prec,nb,2)

print("################")
print("################")
print("##### d = 3 ####")
print("################")
print("################")

d = 3
nb = 2000

print("nb de test")
print(nb)
#comparaison_LU(d,p,prec,nb,2)
#comparaison_LU(d,p,prec,2,2)

print("################")
print("################")
print("##### d = 4 ####")
print("################")
print("################")

d = 4
nb = 2000

print("nb de test")
print(nb)
comparaison_LU(d,p,prec,nb,2)
    
