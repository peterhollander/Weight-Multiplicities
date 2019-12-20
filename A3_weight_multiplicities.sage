'''
Authors: David Melendez, Lucy Martinez
MSRI-UP, June - July 2019

- Wherever not explicitly stated, we are working in the Lie algebra of type A_3.
- This program implements closed forms for the weight (q)-multiplicity of a weight
in the Lie algebra of type A_3.
- Throughout the program we denote by a_1, a_2, a_3 the simple roots of the Lie
algebra of type A_3, and we denote by w_1, w_2, w_3 the fundamental weights.
'''

def kostant_partition_function(xi, q_analog=True, triple_sum=False):
    """
    Evaluates the q-analog of Kostant's partition function for the Lie algebra
    of type A_3 on the weight xi.
    Parameters:
    - xi : 3-tuple (m,n,k) representing the weight m*a_1 + n*a_2 + k*a_3
    - q_analog: Use q-analog of Kostant partition function?
    - triple_sum: Use (slower) triple summation formula for partition function?
    """
    (m,n,k) = xi
    result=0
    q=var('q')

    if (triple_sum):
         # Triple summation to count partitions of xi
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
    else:
        # Closed formulas deduced in the paper
        if (m >= n and k >= n):
            for j in range(m+k-n, m+n+k+1):
                t = m + n + k - j
                L = min(floor(t/2), n-ceil(t/2))
                result += (L+1) * (L + (t%2) + 1) * q^j
        elif (m >= n >= k):
            for j in range(m, m+n+k+1):
                t = m + n + k - j
                F2 = min(floor(t/2), k)
                J = k - F2 - (t%2)
                F1 = max(t-n, 0)
                L = F2-F1
                if (J < 0):
                    c = (L+1)*(2*k - 2*F2 + L + 2)/2
                elif (0 <= J <= L):
                    c = ((L+1)*(L+2) - F2*(F2+1) - k*(k-1))/2 + F2*k + k*L - F2*L + (k-F2)*(t%2)
                elif (J > L):
                    c = (L+1)*(L+(t%2)+1)
                result += c*q^j
        elif (k >= n >= m):
            for j in range(k, k+n+m+1):
                t = k + n + m - j
                F2 = min(floor(t/2), m)
                J = m - F2 - (t%2)
                F1 = max(t-n, 0)
                L = F2-F1
                if (J < 0):
                    c = (L+1)*(2*m - 2*F2 + L + 2)/2
                elif (0 <= J <= L):
                    c = ((L+1)*(L+2) - F2*(F2+1) - m*(m-1))/2 + F2*m + m*L - F2*L + (m-F2)*(t%2)
                elif (J > L):
                    c = (L+1)*(L+(t%2)+1)
                result += c*q^j
        elif (n >= m >= k):
            for j in range(n, m+n+k+1):
                t = m + n + k - j
                F1 = max(0, t - min(m+k,n))
                F2 = min(floor(t/2), k)
                if (t-m < 0):
                    S1 = (t-F2) * (F2+1)
                    #S1 = t*(F2 - F1 + 1) - F2*(F2 + 1)
                elif (0 <= t-m <= F2):
                    S1 = ((t-m) * (t-m+1) + F1*(F1-1))/2 - F2*(F2+1) + m*(t-m-F1+1) + t*(F2-t+m)
                if (0 <= t-k <= F2):
                    S2 = (t-k)*(t-k-F1+1) - ((t-k)*(t-k+1) - F1*(F1-1))/2
                elif (t-k > F2):
                    S2 = (t-k)*(F2-F1+1) - (F2*(F2+1) - F1*(F1-1))/2
                    #S2 = F2*(F2-F1+1) - (F2*(F2+1) - F1*(F1-1))/2
                elif (t-k < 0):
                    S2 = 0
                c = S1 - S2 + F2 - F1 + 1
                result += c * q^j
        elif (n >= k >= m):
            for j in range(n, m+n+k+1):
                t = m + n + k - j
                F1 = max(0, t - min(m+k,n))
                F2 = min(floor(t/2), m)
                if (t-k < 0):
                    S1 = (t-F2) * (F2+1)
                elif (0 <= t-k <= F2):
                    S1 = ((t-k) * (t-k+1) + F1*(F1-1) - 2*F2*(F2+1) + 2*k*(t-k-F1+1) + 2*t*(F2-(t-k)))/2
                if (0 <= t-m <= F2):
                    S2 = (t-m)*(t-m-F1+1) - ((t-m)*(t-m+1) - F1*(F1-1))/2
                elif (t-m > F2):
                    S2 = (t-m)*(F2-F1+1) - (F2*(F2+1) - F1*(F1-1))/2
                elif (t-m < 0):
                    S2 = 0
                c = S1 - S2 + F2 - F1 + 1
                result += c * q^j
        
    return result

def kostant_weight_multiplicity(lam, mu=(0,0,0), q_analog=True):
    """
    Computes the weight q-multiplicity of the weight mu in the irreducible
    representation of sl_4(C) with highest weight lam
    Parameters:
    - lam : 3-tuple (m,n,k) representing the weight m*w_1 + n*w_2 * n*w_3
    - mu : 3-tuple (c1,c2,c3) representing the weight c1*w_1+c2*w_2+c3*w_3
    - q_analog : Boolean. If true, returns the weight q-multiplicity. If false,
                 returns the weight multiplicity.
    """

    # Dictionary mapping Weyl group elements s to the coordinates of 
    # s(lam+rho) - (mu + rho) with respect to the simple roots
    sigma_results_dict = get_sigma_results_dict(lam,mu)

    alternation_set = weyl_alternation_set(lam,mu)

    multiplicity = 0

    # Note that we only need to sum over the Weyl alternation set
    for s in alternation_set:
        # vec: Coefficients of the simple roots in s(lam+rho) - (rho + mu)
        vec = sigma_results_dict[s]
        sign = (-1)^(s.length())
        multiplicity += sign * kostant_partition_function(vec, q_analog=q_analog)
    return multiplicity

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
    return alt_set


###############################################################################
###############################################################################
###############################################################################
###############################################################################

def get_sigma_results_dict(lam, mu):
    # lam + mu = (m,n,k) + (c1,c2,c3) = (m,n,k,c1,c2,c3)
    #     ^ This is NOT vector addition! We're concatenating the tuples
    #       to pass into the callable.
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

def ambient_to_vector(v):
    return matrix([v[i] for i in range(0,4)],ring=SR).transpose()

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
