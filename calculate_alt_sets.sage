def weyl_alternation_set(lam, mu=(0,0,0)):
    """
    Finds the Weyl alternation set for lam = mw_1 + nw_2 + kw_3
    and mu = c_1w_1 + c_2w_2 + c_3w+3
    """
    (m_, n_, k_) = lam
    (c1_, c2_, c3_) = mu

    alt_set = set()
    sigma_results_dict = get_sigma_results_dict(lam,mu)

    for s in W:
        if all( [ sigma_results_dict[s][j] >= 0 for j in range(0,3) ] ):
            alt_set.add(s)
    #print('alternation set')
    #print(alt_set)
    return alt_set

def weyl_alternation_set_17(lam, mu=(0,0,0)):
    """
    Finds the Weyl alternation set for lam = mw_1 + nw_2 + kw_3
    and mu = c_1w_1 + c_2w_2 + c_3w+3
    doesn't look at sigmas we know are never in the alternation set to make program more efficient
    """
    (m_, n_, k_) = lam
    (c1_, c2_, c3_) = mu

    alt_set = set()
    sigma_results_dict = get_sigma_results_dict(lam,mu)

    for s in W_17: 
        if all( [ sigma_results_dict[s][j] >= 0 for j in range(0,3) ] ):
            alt_set.add(s)
    #print('alternation set')
    #print(alt_set)
    return alt_set

'''DEFINITIONS'''

# Takes ambient space vector and converts to column matrix over SR
def ambient_to_list(v):
    return [v[i] for i in range(0,3)]

# Takes vector and converts to a list
def vector_to_list(v):
    return [v[i] for i in range(0,3)]

# Take R3 vector and gives coordinates in terms of alpha's
# v = c1 * a1 + c2 * a2 + c3 * a3
def ambient_to_alpha_coords(v):
    return root_matrix.solve_right(ambient_to_vector(v))

def vector_to_alpha_coords(v):
    return root_matrix.solve_right(v)

def ambient_to_vector(v):
    return matrix([v[i] for i in range(0,3)],ring=SR).transpose()

var('m n k c1 c2 c3')
W = WeylGroup(['C', 3], prefix='s')
a = W.domain().simple_roots()
P = W.domain().positive_roots()
[s1, s2, s3] = W.simple_reflections()
e=s1*s1
W_17 = [e, s1, s2, s3, s1*s2, s2*s1, s2*s3, s3*s1, s3*s2, s1*s2*s1, s2*s3*s1, s2*s3*s2, s3*s2*s1, s3*s1*s2, s3*s2*s3, s3*s1*s2*s1, s3*s2*s3*s2] 

# Simple root vectors for the Lie algebra of type C_3
a1 = ambient_to_vector(a[1]); #print('a1'); print(a1)
a2 = ambient_to_vector(a[2]); #print('a1'); print(a2)
a3 = ambient_to_vector(a[3]); #print('a1'); print(a3)

root_matrix = matrix([ambient_to_list(a[i]) for i in range(1,4)],ring=SR).transpose()

# Fundamental weight vectors for the Lie algebra of type C_3
w1 = (1/2) * (2*a1 + 2*a2 + 1*a3)
w2 = (1*a1 + 2*a2 + 1*a3)
w3 = (1/2) * (2*a1 + 4*a2 + 3*a3)

lam = m*w1 + n*w2 + k*w3
mu = c1*w1 + c2*w2 + c3*w3
rho = ambient_to_vector((1/2)*sum(P))

''' OTHER FUNCTIONS WE NEED '''

def get_sigma_results_dict(lam, mu):
    # lam + mu = (m,n,k) + (c1,c2,c3) = (m,n,k,c1,c2,c3)
    #     ^ This is NOT vector addition! We're concatenating the tuples
    #       to pass into the callable.
    return dict(
    [ (s,[component(*(lam+mu)) for component in weyl_action_callable_dict[s]] ) for s in W
    ]
    )

def weyl_actions():
    """
    Returns a dict with key,value pairs: (s, (b_1, b_2, b_3)) where s is an
    element of the Weyl group and b_1, b_2, b_3 are the coefficients
    of a_1, a_2, and a_3, respectively, of s(lambda + rho) - (rho + mu)
    """
    return [(s,vector_to_alpha_coords(s.matrix() * (lam + rho) - (rho + mu))) for s in W]

weyl_action_callable_dict = dict([ (p[0],[fast_callable(p[1][i][0], vars=[m,n,k,c1,c2,c3]) for i in range(0,3)]) for p in weyl_actions()])

def alt_sets_with_mus(lamnum, munum):
    '''
    iterates through lam=(0,0,0) to (lamnum, lamnum, lamnum) and mu=(0,0,0) to (munum, munum, munum)
    prints alternation sets which show up with a lambda and mu that induce it
    '''
    alt_sets = []
    alt_set_exp = []
    
    #for lam = (m,n,k), mu=(x,y,z)
    for m in range(lamnum+1):
        for n in range(lamnum+1):
            for k in range(lamnum+1):
                for x in range(munum+1):
                    for y in range(munum+1):
                        for z in range(munum+1):
                            alt_set = weyl_alternation_set_17(lam=(m,n,k), mu=(x,y,z))
                            if alt_set not in alt_sets: #if a new set
                                alt_sets.append(alt_set) #record alt set
                                alt_set_exp.append((alt_set,(m,n,k),(x,y,z))) #record how we got it
                                print((alt_set,(m,n,k),(x,y,z)))
    return alt_set_exp


def alt_sets_new_lattice(lamnum, munum):
    '''
    iterates through lam=(0,0,0) to (lamnum, lamnum, lamnum) and mu=(0,0,0) to (munum, munum, munum)
    For m+k+x+z divisible by 2,
    prints alternation sets which show up with a lambda and mu that induce it
    '''
    alt_sets = []
    alt_set_exp = []
    
    #for lam = (m,n,k), mu=(x,y,z)
    for m in range(lamnum+1):
        for n in range(lamnum+1):
            for k in range(lamnum+1):
                for x in range(munum+1):
                    for y in range(munum+1):
                        for z in range(munum+1):
                            if (m+k)%2 == 0 and (x+z)%2 == 0: #if its in our lattice
                                alt_set = weyl_alternation_set_17(lam=(m,n,k), mu=(x,y,z))
                                if alt_set not in alt_sets: #if a new set
                                    alt_sets.append(alt_set) #record alt set
                                    alt_set_exp.append((alt_set,(m,n,k),(x,y,z))) #record how we got it
                                    print((alt_set,(m,n,k),(x,y,z)))
    return alt_set_exp

#gives alternation set for a lambda and a mu
#weyl_alternation_set(lam=(2,1,0), mu=(0,0,0))
#weyl_alternation_set_17(lam=(2,1,0), mu=(0,0,0))

#gives all alt sets for lam going up to lamnum and mu going up to munum
#alt_sets_with_mus(lamnum=0,munum=50)

#gives all alt sets for lam going up to lamnum and mu going up to munum for m+k+x+z divisible by 2
alt_sets_new_lattice(lamnum=10, munum=10)
