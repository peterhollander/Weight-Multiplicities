'''
Authors: David Melendez, Lucy Martinez
MSRI-UP, June - July 2019

- Wherever not explicitly stated, we are working in the Lie algebra of type A_3.
- This program implements closed forms for the weight (q)-multiplicity of a weight
in the Lie algebra of type A_3.
- Throughout the program we denote by a_1, a_2, a_3 the simple roots of the Lie
algebra of type A_3, and we denote by w_1, w_2, w_3 the fundamental weights.
'''

def kostant_partition_function(xi, q_analog=True):
    """
    Evaluates the q-analog of Kostant's partition function for the Lie algebra
    of type A_3 on the weight xi.
    Parameters:
    - xi : 3-tuple (m,n,k) representing the weight m*a_1 + n*a_2 + k*a_3
    """
    (m,n,k) = xi
    result=0
    q=var('q')

    # We use a for loop to count partitions of xi.
    for f in range(0,min(m,n,k)+1):
        ''' f represents the coefficient of our highest weight al_1+al_2+al_3
            in our partition '''
        for d in range(0,min(m-f,n-f)+1):
            ''' d represents the coefficient of one of our next highest weights
                al_1+al_2 in our partition. '''
            for e in range(0, min(n-f-d,k-f)+1):
                ''' e represents the coefficient of the other of our next
                    highest weights al_2+al_3 in our partition.
                    A partition of xi with f (al_1+al_2+al_3)'s,
                    d (al_1+al_2)'s, and e (al_2+al_3)'s has
                    m+n+k-2f-e-d parts. Thus we add the following term to our
                    ourput polynomial: '''
                if(q_analog):
                    result += q^(m+n+k-2*f-e-d)
                else:
                    result += 1
    return result


def weight_multiplicity_A3(lam, mu=(0,0,0), q_analog=True):
    """
    Computes the weight q-multiplicity of the weight mu in the irreducible
    representation of sl_4(C) with highest weight lam
    Parameters:
    - lam : 3-tuple (m,n,k) representing the weight m*w_1 + n*w_2 * n*w_3
    - mu : 3-tuple (c1,c2,c3) representing the weight c1*w_1+c2*w_2+c3*w_3
    - q_analog : Boolean. If true, returns the weight q-multiplicity. If false,
                 returns the weight multiplicity.
    """
    sigma_results_dict = get_sigma_results_dict(lam,mu)

    abs_q_multiplicity = 0
    for s in W:
        vec = sigma_results_dict[s]
        sign = (-1)^(s.length())
        if(any([j < 0 or floor(j) != j for j in vec])):
            continue
        abs_q_multiplicity += sign * kostant_partition_function(vec)
    return sign*abs_q_multiplicity

###############################################################################
###############################################################################
###############################################################################
###############################################################################

def get_sigma_results_dict(lam, mu):
    return dict(
    [ (s,[component(*(lam+mu)) for component in weyl_action_callable_dict[s]] ) for s in W
    ]
    )

# Takes ambient space vector and converts to column matrix over SR
def ambient_to_list(v):
    return [v[i] for i in range(0,4)]

# Takes vector and converts to a list
def vector_to_list(v):
    return [v[i] for i in range(0,4)]


# Take R4 vector and gives coordinates in terms of alpha's
# v = c1 * a1 + c2 * a2 + c3 * a3
def ambient_to_alpha_coords(v):
    return root_matrix.solve_right(ambient_to_vector(v))

def vector_to_alpha_coords(v):
    return root_matrix.solve_right(v)

var('m n k c1 c2 c3')
W = WeylGroup(['A', 3], prefix='s')
a = W.domain().simple_roots()
P = W.domain().positive_roots()
[s1, s2, s3] = W.simple_reflections()

# Simple root vectors for the Lie algebra of type A_3
a1 = ambient_to_vector(a[1])
a2 = ambient_to_vector(a[2])
a3 = ambient_to_vector(a[3])

root_matrix = matrix([ambient_to_list(a[i]) for i in range(1,4)],ring=SR).transpose()

# Fundamental weight vectors for the Lie algebra of type A_3
w1 = (1/4) * (3*a1 + 2*a2 + 1*a3)
w2 = (1/2) * (1*a1 + 2*a2 + 1*a3)
w3 = (1/4) * (1*a1 + 2*a2 + 3*a3)

lam = m*w1 + n*w2 + k*w3
mu = c1*w1 + c2*w2 + c3*w3
rho = ambient_to_vector((1/2)*sum(P))

def weyl_actions():
    """
    Returns a dict with key,value pairs: (s, (b_1, b_2, b_3)) where s is an
    element of the Weyl group and b_1, b_2, b_3 are the coefficients
    of a_1, a_2, and a_3, respectively, of s(lambda + rho) - (rho + mu)
    """
    return [(s,vector_to_alpha_coords(s.matrix() * (lam + rho) - (rho + mu))) for s in W]

weyl_action_callable_dict = dict([ (p[0],[fast_callable(p[1][i][0], vars=[m,n,k,c1,c2,c3]) for i in range(0,3)]) for p in weyl_actions()])


# Reload file in interactive Sage environment
def rl():
    load('partitions.sage')
