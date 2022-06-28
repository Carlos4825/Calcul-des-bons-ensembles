
def dot(a, b):
    # renvoie le produit scalaire de a par b
    s = 0
    for i in range(len(a)):
        s += a[i]*b[i]
    return s


def generateur(v, m, u):
    # engendre tous les v vérifiant les conditions v.m <= 2 et u.v >= 0.
    # modifie v en place
    n = len(v)
    ind_cour = 0
    while True:
        if dot(v, m) > 2 or (dot(v, u) < 0 and m[ind_cour] == 0):
            # on ne pourra jamais revenir à dot(v, m) <= 2 ou dot(v, u) >= 0
            # (en utilisant le fait que u[i] < 0 pour ind_cour <= i <= n-1)
            v[ind_cour] = 0
            ind_cour -= 1
            if ind_cour == -1:
                break
            v[ind_cour] += 1
            continue
        if dot(v, u) >= 0:
            yield None
        ind_cour = n - 1
        v[ind_cour] += 1


def long_generateur(m, u):
    v = [0]*len(m)
    count = 0
    for _ in generateur(v, m, u):
        count += 1
    return count


def init_ineq_eq(n, angle_max, ineq, eq):
    
    # les angles sont positifs
    # eg (0, 1, 0, 0, 0, 0) / ... / (0, 0, 0, 0, 0, 1)
    ineq += [ tuple([0]*i + [1] + [0]*(n-i)) for i in range(1, n+1) ]
    
    # les angles valent au plus angle_max
    # eg (angle_max, -180, 0, 0, 0, 0) / ... / (angle_max, 0, 0, 0, 0, -180)
    ineq += [ tuple([angle_max] + [0]*(i-1) + [-180] + [0]*(n-i)) for i in range(1, n+1)]
    
    # les angles sont par ordre décroissant
    # eg (0, 1, -1, 0, 0, 0) / ... / (0, 0, 0, 0, 1, -1)
    ineq += [ tuple([0]*i + [1, -1] + [0]*(n-i-1)) for i in range(1, n) ]

    # les angles ont (n-2) pour somme
    # eg (-3, 1, 1, 1, 1, 1)
    eq += [ tuple([-n+2] + [1]*n) ]


def generateur_compat(v, beta):
    # engendre tous les v tels que beta.v = 2
    # modifie v en place
    n = len(v)
    ind_cour = n - 1
    while True:
        val_beta_v = dot(beta, v)
        if val_beta_v > 2:
            v[ind_cour] = 0
            ind_cour -= 1
            if ind_cour == -1:
                break
            v[ind_cour] += 1
            continue
        if val_beta_v == 2:
            yield None
        ind_cour = n - 1
        v[ind_cour] += 1
    
