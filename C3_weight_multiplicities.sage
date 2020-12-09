'''
Authors: David Melendez, Lucy Martinez
MSRI-UP, June - July 2019

- Wherever not explicitly stated, we are working in the Lie algebra of type A_3. (now C_3)
- This program implements closed forms for the weight (q)-multiplicity of a weight
in the Lie algebra of type A_3. (now C_3)
- Throughout the program we denote by a_1, a_2, a_3 the simple roots of the Lie
algebra of type A_3 (now C_3), and we denote by w_1, w_2, w_3 the fundamental weights.
'''

def kostant_partition_function(xi, q_analog=True, triple_sum=True):
    """
    Evaluates the q-analog of Kostant's partition function for the Lie algebra
    of type C_3 on the weight xi.
    Parameters:
    - xi : 3-tuple (m,n,k) representing the weight m*a_1 + n*a_2 + k*a_3
    - q_analog: Use q-analog of Kostant partition function?
    - triple_sum: Use (slower) triple summation formula for partition function? -- this is the only formula we have at the moment
    """
    (m,n,k) = xi
    result=0
    q=var('q')

    if (triple_sum):
         # Triple summation to count partitions of xi
        for h in range(0,min(floor(m/2),floor(n/2),k)+1):
            ''' h represents the coefficient of 2al_1+2al_2+al_3
                in our partition '''
            for g in range(0,min(m-2*h,floor((n-2*h)/2),k-h)+1):
                ''' g represents the coefficient of al_1+2al_2+al_3
                    in our partition. '''
                for f in range(0, min(m-2*h-g,n-2*h-2*g,k-h-g)+1):
                    ''' f represents the coefficient of al_1+al_2+al_3
                        in our partition.'''
                    for i in range(0, min(floor((n-2*h-2*g-f)/2),k-h-g-f)+!):
                        '''i represents the coefficient of 2al_2+al_3 
                        in our partition'''
                        for d in range(0, min(m-2*h-g-f,n-2*h-2*g-f-2*i)+1):
                            '''d represents the coefficient of al_1+al_2 
                            in our partition'''
                            for e in range(0, min(n-2*h-2*g-f-2*i-d,k-h-g-f-i)+1):
                                '''e represents the coefficient of al_2+al_3
                                in our partition'''
                                '''A partition of xi with these coefficients has
                                m+n+k-d-e-2f-3g-4h-2i parts. Thus we add the following term to our
                                output polynomial: '''
                                if(q_analog):
                                    result += q^(m+n+k-d-e-2*f-3*g-4*h-2*i)
                                else:
                                    result += 1
    ''' Will edit when we have closed formulas
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
       '''

    if not q_analog:
        result = result.subs({q:1})
        
    return result




def kostant_weight_multiplicity(lam, mu=(0,0,0), q_analog=True):
    """
    Computes the weight q-multiplicity of the weight mu in the irreducible
    representation of sp_6(C) with highest weight lam
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



###vvv###
def plot_empty_region(mu, dist=15, color='red', size=25):
    """
    Plots the empty region (points lambda for which no Weyl group elements contribute)
    for a given mu = (m,n,k) correspondingn to m*w_1 + n*w_2 + k*w_3
    Parameters:
    - dist: Radius of region to check & plot points
    - color: Color of points
    - size: Radius of points (which are spheres)
    """
    return show(point_plot_reversed(dist, mu, W, color, size), frame=False)

def plot_alternation_diagram(sigmas, mu=(0,0,0), restricted=False, color='red', dist=20, size=15):
    '''
    Plots the Weyl alternation diagram for a given weight mu=m*w_1 + n*w+2 + k*w_3
    Parameters:
    restricted: If true, plot weights lambda such that A(lambda, mu) is exactly sigmas.
                If false, plot weights lambda such that A(lambda, mu) contains sigmas.
    '''
    the_plot = alt_set_plot(dist=dist, mu=mu, sigmas=sigmas, 
        color=color, size=size, restricted=restricted)
    return the_plot


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
W = WeylGroup(['C', 3], prefix='s')
a = W.domain().simple_roots()
P = W.domain().positive_roots()
[s1, s2, s3] = W.simple_reflections()
e=s1*s1

# Simple root vectors for the Lie algebra of type C_# (?)
a1 = ambient_to_vector(a[1])
a2 = ambient_to_vector(a[2])
a3 = ambient_to_vector(a[3])

root_matrix = matrix([ambient_to_list(a[i]) for i in range(1,4)],ring=SR).transpose()

# Fundamental weight vectors for the Lie algebra of type C_3
w1 = (1/2) * (2*a1 + 2*a2 + 1*a3)
w2 = (1*a1 + 2*a2 + 1*a3)
w3 = (1/2) * (2*a1 + 4*a2 + 3*a3)

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

def weyl_actions_sub(a,b,c):
    return [(s,vector_to_alpha_coords(s.matrix() * (lam + rho) - (rho + mu)).substitute([m==a,n==b,k==c])) for s in W]


weyl_action_callable_dict = dict([ (p[0],[fast_callable(p[1][i][0], vars=[m,n,k,c1,c2,c3]) for i in range(0,3)]) for p in weyl_actions()])

###############################################################################
###############################################################################
###############################################################################
###############################################################################
# Diagram code

# Create alpha & omega column matrices over Q

w1_ = vector(w1.transpose())
w2_ = vector(w2.transpose())
w3_ = vector(w3.transpose())

a1_ = vector(a1.transpose())
a2_ = vector(a2.transpose())
a3_ = vector(a3.transpose())

# Matrix whose rows are w1, w2, w3
omegas = matrix([w1_, w2_, w3_], ring=QQ)

# Matrix whose rows are a1, a2, a3
alphas = matrix([a1_, a2_, a3_], ring=QQ)


def gram_schmidt_symb(M):
    return [v.normalized() for v in M.gram_schmidt()[0].rows()]

# Substitutions of m,n,k in terms of x,y,z,c1,c2,c3 to ensure that alpha coefficients
# in partition function input are integers for integral x,y,z
var('x y z')
sub_1 = matrix([[3,2,1],[1,2,1],[1,2,3]],ring=SR).solve_right(vector([4*x+3*c1+2*c2+c3,2*y+c1+2*c2+c3,4*z+c1+2*c2+3*c3]))

# sigma(lam+rho) - (rho + mu) with integrality substitutions for m,n,k
# applied
sub_1_result = weyl_actions_sub(*sub_1);
sub_1_dict = dict(sub_1_result)

sub_1_callable = [[fast_callable(p[1][i][0], vars=[x,y,z,c1,c2,c3]) for i in range(0,3)] for p in sub_1_result]
sub_1_callable_dict = dict([(p[0], [fast_callable(p[1][i][0], vars=[x,y,z,c1,c2,c3]) for i in range(0,3)]) for p in sub_1_result])

xyz_to_mnk = [fast_callable(X, vars=[x,y,z,c1,c2,c3]) for X in sub_1]

alpha_projection_rows = gram_schmidt_symb(alphas)
alpha_projection_cols = matrix(alpha_projection_rows).transpose()

wp1 = alpha_projection_cols.solve_right(w1_)
wp2 = alpha_projection_cols.solve_right(w2_)
wp3 = alpha_projection_cols.solve_right(w3_)

ap1 = alpha_projection_cols.solve_right(a1_)
ap2 = alpha_projection_cols.solve_right(a2_)
ap3 = alpha_projection_cols.solve_right(a3_)


def point_plot_count_alpha(mu, dist=10):
    # Get the xyz coordinates
    sigmas = W
    A = matrix([[3/4,1/2,1/4],[1/2,1,1/2],[1/4,1/2,3/4]])
    [c1_, c2_, c3_] = A * vector(mu)
    coords_xyz = [
            (x_, y_, z_)
            for x_ in range(-dist, dist) for y_ in range(-dist, dist)
            for z_ in range(-dist, dist)
            if all( [
                any([
                sub_1_callable[list(sub_1_dict.keys()).index(s)][j](x_,y_,z_,c1_,c2_,c3_) < 0
                for j in range(0,3)
                ])
                for s in sigmas])
        ]
    return len(coords_xyz)

def point_plot_reversed(dist, mu, sigmas, color, size=10):
    if(isinstance(color,basestring)):
        col = tuple(colors[color])
    else:
        col = color
    # Get the xyz coordinates
    c1_ = mu[0]
    c2_ = mu[1]
    c3_ = mu[2]
    coords_xyz = [
            (x_, y_, z_)
            for x_ in range(-dist, dist) for y_ in range(-dist, dist)
            for z_ in range(-dist, dist)
            if all( [
                any([
                sub_1_callable[list(sub_1_dict.keys()).index(s)][j](x_,y_,z_,c1_,c2_,c3_) < 0
                for j in range(0,3)
                ])
                for s in sigmas])
        ]

    # Substitute in m,n,k
    coords_mnk = [tuple( [xyz_to_mnk[j](x_,y_,z_,c1_,c2_,c3_) for j in [0,1,2]] )
        for (x_,y_,z_) in coords_xyz]

    # Transform into omega coordinates
    coords_mnk_omega = [m_*wp1 + n_*wp2 + k_*wp3 for (m_,n_,k_) in coords_mnk]
    max_z = max([pt[2] for pt in coords_mnk_omega]) + 1
    min_z = min([pt[2] for pt in coords_mnk_omega]) - 1
    range_z = max_z-min_z
    points = [point3d(pt, size, color=tuple(min(1,col[j] * max(0.3,(max_z-pt[2]+0.3)/(range_z))) for j in range(0,3)), opacity=1) for pt in coords_mnk_omega]
    return show(sum(points), frame=False)

def alt_set_plot(dist, mu, sigmas, restricted=False, color='red', size=15):
    others = set(W).difference(sigmas)
    # Get the xyz coordinates
    if(isinstance(color,basestring)):
        col = tuple(colors[color])
    else:
        col = color
    c1_ = mu[0]
    c2_ = mu[1]
    c3_ = mu[2]
    sub1dict = dict(sub_1_result)
    if restricted:
        coords_xyz = [
                (x_, y_, z_)
                for x_ in range(-dist, dist) for y_ in range(-dist, dist)
                for z_ in range(-dist, dist)
                if all( [
                    sub_1_callable_dict[s][j](x_,y_,z_,c1_,c2_,c3_) >= 0
                    for j in range(0,3) for s in sigmas])
                and all( [ # For every other sigma
                    any([sub_1_callable_dict[s][j](x_,y_,z_,c1_,c2_,c3_) < 0 for j in range(0,3)]) # this sigma does not contribute
                    for s in others])
            ]
    else:
        coords_xyz = [
                (x_, y_, z_)
                for x_ in range(-dist, dist) for y_ in range(-dist, dist)
                for z_ in range(-dist, dist)
                if all( [
                    sub_1_callable_dict[s][j](x_,y_,z_,c1_,c2_,c3_) >= 0
                    for j in range(0,3) for s in sigmas])
            ]

    # A * [alpha coeffs vector] = [omega coeffs vector]
    A = matrix([[3/4, 1/2, 1/4], [1/2, 1, 1/2], [1/4, 1/2, 3/4]])
    #print(A.solve_right(vector([0,1,2])))

    # Substitute in m,n,k
    coeffs_mnk_omega = [tuple( [xyz_to_mnk[j](x_,y_,z_,c1_,c2_,c3_) for j in [0,1,2]] )
        for (x_,y_,z_) in coords_xyz]
    coeffs_mnk_omega_nn = [vec for vec in coeffs_mnk_omega 
            if all([vec[j] >= 0 for j in [0,1,2]])]
    coeffs_mnk_alpha = [A*vector(omega_coeffs) 
            for omega_coeffs in coeffs_mnk_omega]
    coeffs_mnk_alpha_nn = [vec for vec in coeffs_mnk_alpha 
            if all([vec[j] >= 0 for j in [0,1,2]])]

    # Transform into omega coordinates
    #coords_mnk_omega = [m_*wp1 + n_*wp2 + k_*wp3 for (m_,n_,k_) in coords_mnk]
    coords_mnk_alpha = [m_*ap1 + n_*ap2 + k_*ap3 for (m_,n_,k_) in coeffs_mnk_alpha]
    coords_mnk_omega = [m_*wp1 + n_*wp2 + k_*wp3 for (m_,n_,k_) in coeffs_mnk_omega]
    coords_mnk_alpha_nn = [m_*ap1 + n_*ap2 + k_*ap3 
        for (m_,n_,k_) in coeffs_mnk_alpha_nn]
    coords_mnk_omega_nn = [m_*wp1 + n_*wp2 + k_*wp3 
        for (m_,n_,k_) in coeffs_mnk_omega_nn]
    pts = coeffs_mnk_omega
    pts_plot = coords_mnk_omega
    if (not pts):
        print('empty')
        return 0;
    #coords_mnk_omega = [m_*wp1 + n_*wp2 + k_*wp3 for (m_,n_,k_) in coords_mnk if m_>=0 and n_>=0 and k_>=0]
    max_z = max([pt[2] for pt in pts_plot])+1
    min_z = min([pt[2] for pt in pts_plot])-1
    range_z = max_z-min_z
    if (range_z == 0):
        range_z = 1
    points = [point3d(pt, size, color=tuple(col[j] * (max_z-pt[2]+0.1)/range_z for j in range(0,3)), opacity=.75) for pt in pts_plot]
    return show(sum(points) + fundamental_weight_plot(range_z), frame=False)


def fundamental_weight_plot(scale=1):
    pO1 = tuple(ap1)
    pO2 = tuple(ap2)
    pO3 = tuple(ap3)

    width = scale
    txt = 400
    O1 = arrow3d(start=(0,0,0), end=tuple(wp1*scale), color='red', width=width)
    O2 = arrow3d(start=(0,0,0), end=tuple(wp2*scale), color='green', width=width)
    O3 = arrow3d(start=(0,0,0), end=tuple(wp3*scale), color='blue', width=width)
    O1_label = text3d("w_1", tuple((scale+0.3)*wp1), color='black', 
            fontsize='xx-large', fontweight='black', background_color='white')
    O2_label = text3d("w_2", tuple((scale+0.3)*wp2), color='black', 
            fontsize='xx-large', fontweight='black', background_color='white')
    O3_label = text3d("w_3", tuple((scale+0.3)*wp3), color='black', 
            fontsize='xx-large', fontweight='black', background_color='white')

    return O1 + O2 + O3 + O1_label + O2_label + O3_label


def simple_root_plot(scale=1):
    pO1 = tuple(ap1)
    pO2 = tuple(ap2)
    pO3 = tuple(ap3)

    width = scale
    txt = 400
    O1 = arrow3d(start=(0,0,0), end=tuple(ap1*scale), color='red', width=width)
    O2 = arrow3d(start=(0,0,0), end=tuple(ap2*scale), color='green', width=width)
    O3 = arrow3d(start=(0,0,0), end=tuple(ap3*scale), color='blue', width=width)
    O1_label = text3d("a_1", tuple((scale+0.3)*ap1), color='black', 
            fontsize='xx-large', fontweight='black', background_color='white')
    O2_label = text3d("a_2", tuple((scale+0.3)*ap2), color='black', 
            fontsize='xx-large', fontweight='black', background_color='white')
    O3_label = text3d("a_3", tuple((scale+0.3)*ap3), color='black', 
            fontsize='xx-large', fontweight='black', background_color='white')

    return O1 + O2 + O3 + O1_label + O2_label + O3_label


# Reload file in interactive Sage environment
def rl():
    load('partitions.sage')


rl()
kostant_partition_function((1,2,3), q_analog=True, triple_sum=True)
#weyl_alternation_set(lam=(2,3,4), mu=(0,0,0))
#kostant_weight_multiplicity(lam=(2,3,4), mu=(0,0,0), q_analog=True)

#plot_empty_region((0,0,0), dist=15, color='red', size=25)
