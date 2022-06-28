import time
import sys
from random import randrange

from sage.geometry.polyhedron.parent import Polyhedra
from sage.geometry.polyhedron.backend_cdd import Polyhedron_QQ_cdd

from commun import dot, generateur, long_generateur, init_ineq_eq, generateur_compat

# dot(a, b) renvoie le produit scalaire des vecteurs a et b.
# generateur(v, m, u) engendre tous les v vérifiant les conditions v.m <= 2 et u.v >= 0. Modifie v en place.
# long_generateur(m, u) compte le nombre de vecteurs v vérifiant les conditions ci-dessus
# init_ineq_eq(n, angle_max, ineq, eq) initialise les listes ineq et eq avec les conditions que doivent vérifier les angles du polygone (compris entre 0 et angle_max, ordonnés et de somme (n-2))
# generateur_compat(v, beta) engendre tous les v tels que beta.v = 2. Modifie v en place

if len(sys.argv) != 3:
    print("Indiquer le nombre de côtés et l'angle maximal")
    exit()

n = int(sys.argv[1])
angle_max = int(sys.argv[2])

def is_good(X):    # Teste si X est bon
    ineq_test = [ tuple([0] + [v[i] for i in range(n)]) for v in X]
    B_test = Polyhedron_QQ_cdd(parent, None, [ineq_test, [tuple([0] + [1]*n)]]) & polytopes.hypercube(n)
    for som in B_test.vertices():
        for v in X:
            if dot(vector(som), v) > 0:
                return False
    return True

nb_appels = [0] ## appels recurse
t0 = time.time()

set_good = {} # set_good[profondeur] répertorie les bons ensembles pour une certaine profondeur
for profondeur in range(n):
    set_good[profondeur] = []

ineq = [] # conditions que doivent respecter les angles du polygone
eq = []
init_ineq_eq(n, angle_max, ineq, eq)

parent = Polyhedra(QQ, n, backend='cdd')

def recurse(X, profondeur, interdits):
    
    # Les premiers vecteurs de X forment une famille 'libre'
    # (il y en a autant que la profondeur)
    # les suivants ne sont que les vecteurs compatibles.

    nb_appels[0] += 1
    
    # Construire le polytope Bx
    eqx = [ tuple([-2] + [v[i] for i in range(n)]) for v in X]
    Bx = Polyhedron_QQ_cdd(parent, None, [ineq, eq+eqx])

    # S'assurer que le polytope est non vide
    if Bx.is_empty():
        return None

    # et qu'il contient des points non dégénérés
    for c in range(n):
        c_min = min([vector(som)[c] for som in Bx.vertices()])
        c_max = max([vector(som)[c] for som in Bx.vertices()])
        if c_min == c_max and c_min in [0, 1, 2]:
            return None
        
    # On cherche les vecteurs compatibles avec X
    # C'est un vecteur v tel que beta*v = 2 pour tout beta dans Bx
    beta = Bx.center()
    v = [0]*n
    for _ in generateur_compat(v, beta):
        if tuple(v) not in X:
            # tester si tous les points de Bx vérifient cette condition
            flag_appartient = True
            for som in Bx.vertices():
                if dot(vector(som), v) != 2:
                    flag_appartient = False
                    break
            if flag_appartient:
                X += [tuple(v)]
                if tuple(v) in interdits:
                    return None
        
    if is_good(X):      # bon ensemble trouvé
        set_good[profondeur].append(X)
        print('      ***', end = '')
        for i in range(n):
            print(len(set_good[i]), end = ' ')
        print("     = ", sum([len(set_good[i]) for i in range(n)]))

    # Calcul de m
    m = [ min([vector(s)[i] for s in Bx.vertices()]) for i in range(n) ]

    
    # Choix de u
        
    # Construction du polyèdre des possibilités pour u
    # u est de somme nulle,
    # de produit scalaire nul avec les vecteurs de X
    # on peut wlog le borner à [-1, 1]^n
    # et sur les composantes où m est nul, et u doit être < 0
    # (prendre un barycentre non dégénéré pour respecter le strict)
    
    eq_u = [ tuple([0] + [1]*n) ]
    eqx_u = [ tuple([0] + [v[i] for i in range(n)]) for v in X]
    ineq_u = [tuple([1] + [0]*i + [1] + [0]*(n-i-1)) for i in range(n)]
    ineq_u += [tuple([int(m[i] > 0)] + [0]*i + [-1] + [0]*(n-i-1)) for i in range(n)]
    Bu = Polyhedron_QQ_cdd(parent, None, [ineq_u, eq_u+eqx_u])
    vert = [vector(i) for i in Bu.vertices()]
    nb_som = len(vert)
    
    den = 10**8       # décimal aléatoire jusqu'à 8 chiffres après la virgule
    u_aleat = [0]*n
    poids = [0]*nb_som
    
    # nombre de tirages, en fonction du nombre de sommets du polyèdre, heuristique
    if nb_som <= 7:
        nb_times = [1, 5, 5, 10, 10, 20, 20][nb_som-1]
    elif nb_som <= 12:
        nb_times = 150
    else:
        nb_times = 1000

    wr = 10**6
    for rep in range(nb_times):
        for i in range(nb_som):
            poids[i] = randrange(1, den)
        for ind in range(n):
            s_ind = 0
            for k in range(nb_som):
                s_ind += poids[k] * vert[k][ind] / den
            u_aleat[ind] = s_ind
        if long_generateur(m, u_aleat) < wr:
            wr = long_generateur(m, u_aleat)
            u = tuple(u_aleat)
        

            
    # Backtrack 
    visites = set()
    v = [0]*n
    for _ in generateur(v, m, u):
        if tuple(v) not in X and tuple(v) not in interdits:
            
            if profondeur < 1:
                for _ in range(profondeur):
                    print('-  ', end = '')
                print(tuple(v))

            # éliminer les polyèdres vides et certains avec des angles dégénérés
            flag_negatif = False
            flag_positif = False
            count_ex = 0
            for som in Bx.vertices():
                val = dot(vector(som), v) - 2
                if val < 0:
                    flag_negatif = True
                elif val > 0:
                    flag_positif = True
                elif not (flag_negatif and flag_positif):
                    rem_som = vector(som)
                    count_ex += 1

            if not ( (count_ex >= 1) or (flag_positif and flag_negatif)):
                continue

            if ( not (flag_negatif and flag_positif) ) and count_ex == 1:
                # polytope réduit à rem_som, vérifier que les angles sont != 0, 180, 360
                if 0 in rem_som or 1 in rem_som or 2 in rem_som:
                    continue

            recurse([tuple(v)] + X, profondeur + 1, visites.union(interdits))
            visites.add(tuple(v))
            
    return None



recurse([], 0, set())
print("\n\n\n")
print(sum([len(set_good[i]) for i in range(n)]))
print("Time : ", time.time() - t0)

print("\nNombre d'appels récursifs : ", nb_appels[0])
for i in range(n):
    print(len(set_good[i]), end = ' ')
print()
print("\n\n\n")

for i in range(n):
    print(i)
    print(set_good[i])
    print()
print()
