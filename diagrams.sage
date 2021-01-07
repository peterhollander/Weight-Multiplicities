from sage.plot.plot3d.shapes2 import Line

load('weightmult.sage')
init()

weyl_group = list(W) #W is from weightmult, it's the weyl group


def alternation_diagram(sigmas, color): #doesn't need changing? i think?
    return (omega_plot(dist=15,color='black')
        + point_plot(dist=10, mu=(0,0,0), sigmas=sigmas, color=color))


def cool_pic(): #changed
    return (omega_plot(dist=10,color='black')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[e], color='red')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s1], color='blue')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s2], color='green')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s3], color='yellow')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s1*s2], color='brown')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s2*s1], color='purple')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s3*s1], color='blueviolet')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s2*s3], color='cadetblue')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s3*s2], color='darkmagenta')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s1*s2*s3], color='darksalmon')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s1*s2*s1], color='firebrick')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s3*s2*s1], color='fuchsia')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s3*s1*s2], color='gold')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s2*s3*s1], color='greenyellow')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s2*s3*s2], color='darkturquoise')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s1*s2*s3*s1], color='lavender')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s1*s2*s3*s2], color='lightcoral')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s2*s3*s1*s2], color='navy')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s2*s3*s2*s1], color='olive')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s3*s1*s2*s1], color='plum')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s1*s2*s3*s2*s1], color='tan')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s1*s2*s3*s1*s2], color='tomato')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s2*s3*s1*s2*s1], color='yellowgreen')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s1*s2*s3*s1*s2*s1], color='violet') #--
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s3*s2*s3], color='lightyellow')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s3*s2*s3*s1], color='')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s3*s2*s3*s2], color='')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s3*s1*s2*s3], color='')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s2*s3*s1*s2*s3], color='')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s3*s1*s2*s3*s1], color='')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s3*s1*s2*s3*s2], color='')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s3*s2*s3*s2*s1], color='')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s3*s2*s3*s1*s2], color='')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s2*s3*s1*s2*s3*s1], color='')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s2*s3*s1*s2*s3*s2], color='')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s3*s1*s2*s3*s1*s2], color='')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s3*s1*s2*s3*s2*s1], color='')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s3*s2*s3*s1*s2*s1], color='')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s3*s2*s3*s1*s2*s3], color='')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s2*s3*s1*s2*s3*s2*s1], color='')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s2*s3*s1*s2*s3*s1*s2], color='')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s3*s2*s3*s1*s2*s3*s2], color='')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s3*s1*s2*s3*s1*s2*s1], color='')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s3*s2*s3*s1*s2*s3*s1], color='')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s2*s3*s1*s2*s3*s1*s2*s1], color='')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s3*s2*s3*s1*s2*s3*s2*s1], color='')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s3*s2*s3*s1*s2*s3*s1*s2], color='')
        + point_plot(dist=10, mu=(0,0,0), sigmas=[s3*s2*s3*s1*s2*s3*s1*s2*s1], color=''))

def omega_plot(dist, color): #are these the fancy w's omegas? whatever, then they might not need changing -- what is wp1 2 and 3, where are they defined? --
    #this draws the wp lines, their negatives and their labels, or at least returns them
    O1 = Line([(0,0,0), tuple(dist*wp1)], arrow_head=True, color=color)
    O2 = Line([(0,0,0), tuple(dist*wp2)], arrow_head=True, color=color)
    O3 = Line([(0,0,0), tuple(dist*wp3)], arrow_head=True, color=color)
    O1_label = text3d("w_1", tuple((dist+1.3)*wp1), color=color)
    O2_label = text3d("w_2", tuple((dist+1.3)*wp2), color=color)
    O3_label = text3d("w_3", tuple((dist+1.3)*wp3), color=color)
    O1_negative = Line([(0,0,0), tuple(-dist*wp1)], arrow_head=True, color=color)
    O2_negative = Line([(0,0,0), tuple(-dist*wp2)], arrow_head=True, color=color)
    O3_negative = Line([(0,0,0), tuple(-dist*wp3)], arrow_head=True, color=color)
    return O1 + O2 + O3 + O1_negative + O2_negative + O3_negative + O1_label + O2_label + O3_label

def operation_1_plot(): #what are these vectors? what  are ap 1 2 and 3,where are they defined? are they the alphas? (ap1; ap1+ap2; ap1, ap1+ap2)
    vec = Line([(0,0,0), tuple(ap1+ap2)], color='blue', arrow_head=True)
    vec_1 = Line([(0,0,0), tuple(ap1)], color='gray', arrow_head=True)
    vec_2 = Line([tuple(ap1), tuple(ap1+ap2)], color='gray', arrow_head=True)
    return vec + vec_1 + vec_2

def operation_2_plot(): #(ap2; ap2+ap3; ap2, ap2+ap3)
    vec = Line([(0,0,0), tuple(ap2+ap3)], color='red', arrow_head=True)
    vec_1 = Line([(0,0,0), tuple(ap2)], color='gray', arrow_head=True)
    vec_2 = Line([tuple(ap2), tuple(ap2+ap3)], color='gray', arrow_head=True)
    return vec + vec_1 + vec_2

def operation_3_plot(): #(ap1; ap1+ap2+ap3; ap1, ap1+ap2; ap1+ap2, ap1+ap2+ap3)
    vec = Line([(0,0,0), tuple(ap1+ap2+ap3)], color='green', arrow_head=True)
    vec_1 = Line([(0,0,0), tuple(ap1)], color='gray', arrow_head=True)
    vec_2 = Line([tuple(ap1), tuple(ap1+ap2)], color='gray', arrow_head=True)
    vec_3 = Line([tuple(ap1+ap2), tuple(ap1+ap2+ap3)], color='gray', arrow_head=True)
    return vec + vec_1 + vec_2 + vec_3

def partition_plot(): #plots partitions of mu = 2al_1 + 2al_2 + 3al_3 in sl_4 -- would have to change it to all partitions in sp_6 -- this seems not that useful??
    #draws all possible lines to mu = 2al_1 + 2al_2 + 3al_3 through its partitions i think (like vector addition)
    vec = Line([(0,0,0), tuple(2*ap1+2*ap2+2*ap3)], color='violet', arrow_head=True)

    # First partition: al_1 + al_2 + al_3 + (al_1+al_2+al_3)
    #+al_1
    vec_1_2 = Line([(0,0,0), tuple(ap1)], color='gray', arrow_head=True)
    #+(a_1+al_2+al_3)
    vec_1_1 = Line([tuple(ap1), tuple(2*ap1+ap2+ap3)], color='green', arrow_head=True)
    #+al_2
    vec_1_3 = Line([tuple(2*ap1+ap2+ap3), tuple(2*ap1+2*ap2+ap3)], color='gray', arrow_head=True)
    #+al_3
    vec_1_4 = Line([tuple(2*ap1+2*ap2+ap3), tuple(2*ap1+2*ap2+2*ap3)], color='gray', arrow_head=True)

    # Second partition: al_1 + al_3 + (al_1+al_2) + (al_2+al_3)
    #+al_1
    vec_2_1 = Line([(0,0,0), tuple(ap1)], color='gray', arrow_head=True)
    #+al_3
    vec_2_2 = Line([tuple(ap1), tuple(ap1+ap3)], color='gray', arrow_head=True)
    #+(al_1+al_2)
    vec_2_3 = Line([tuple(ap1+ap3), tuple(2*ap1+ap2+ap3)], color='blue', arrow_head=True)
    #+(al_2+al_3)
    vec_2_4 = Line([tuple(2*ap1+ap2+ap3), tuple(2*ap1+2*ap2+2*ap3)], color='red', arrow_head=True)

    #return vec + vec_1_1 + vec_1_2 + vec_1_3 + vec_1_4
    return vec + vec_2_1 + vec_2_2 + vec_2_3 + vec_2_4

def positive_root_plot(dist, color, color1, color2, color3, thickness): #needs to be updated to reflect roots and positive roots of sp_6 -- might need Harris's help for this one
    t = thickness

    pO1 = tuple(ap1)
    pO_1 = tuple(-ap1)
    pO2 = tuple(ap2)
    pO_2 = tuple(-ap2)
    pO3 = tuple(ap3)
    pO_3 = tuple(-ap3)

    pO12 = tuple(ap1+ap2)
    pO1_2 = tuple(ap1-ap2)
    pO_12 = tuple(-ap1+ap2)
    pO_1_2 = tuple(-ap1-ap2)

    pO23 = tuple(ap2+ap3)
    pO2_3 = tuple(ap2-ap3)
    pO_23 = tuple(-ap2+ap3)
    pO_2_3 = tuple(-ap2-ap3)

    pO123 = tuple(ap1+ap2+ap3)
    pO_123 = tuple(-ap1+ap2+ap3)
    pO1_23 = tuple(ap1-ap2+ap3)
    pO12_3 = tuple(ap1+ap2-ap3)
    pO_1_23 = tuple(-ap1-ap2+ap3)
    pO1_2_3 = tuple(ap1-ap2-ap3)
    pO_12_3 = tuple(-ap1+ap2-ap3)
    pO_1_2_3 = tuple(-ap1-ap2-ap3)

    label_A1 = text3d("a_1", tuple((norm(ap1))*ap1), color=color1)
    label_A2 = text3d("a_2", tuple((norm(ap2))*ap2), color=color1)
    label_A3 = text3d("a_3", tuple((norm(ap3))*ap3), color=color1)
    label_A12 = text3d("a_1+a_2", tuple((norm(ap1+ap2))*(ap1+ap2)), color=color2)
    label_A23 = text3d("a_2+a_3", tuple((norm(ap2+ap3))*(ap2+ap3)), color=color2)
    label_A123 = text3d("a_1+a_2+a_3", tuple((norm(ap1+ap2+ap3))*(ap1+ap2+ap3)), color=color3)
    labels = label_A1 + label_A2 + label_A3 + label_A12 + label_A23 + label_A123

    #vertices, has all the points
    verts = [pO1, pO_1, pO2, pO_2, pO3, pO_3, pO12, pO1_2, pO_12, pO_1_2,
    pO23, pO2_3, pO_23, pO_2_3, pO123, pO_123, pO1_23, pO12_3, pO_1_23,
    pO1_2_3, pO_12_3, pO_1_2_3]
    verts = [tuple(RDF(j) for j in v) for v in verts]

    #a line from 0 to each point
    O1 = Line([(0,0,0), tuple(ap1)], arrow_head=True, color=color1, thickness=t)
    O_1 = Line([(0,0,0), tuple(-ap1)], arrow_head=True, color=color1, thickness=t)
    O2 = Line([(0,0,0), tuple(ap2)], arrow_head=True, color=color1, thickness=t)
    O_2 = Line([(0,0,0), tuple(-ap2)], arrow_head=True, color=color1, thickness=t)
    O3 = Line([(0,0,0), tuple(ap3)], arrow_head=True, color=color1, thickness=t)
    O_3 = Line([(0,0,0), tuple(-ap3)], arrow_head=True, color=color1, thickness=t)

    O12 = Line([(0,0,0), tuple(ap1+ap2)], arrow_head=True, color=color2, thickness=t)
    O1_2 = Line([(0,0,0), tuple(ap1-ap2)], arrow_head=True, color=color, thickness=t)
    O_12 = Line([(0,0,0), tuple(-ap1+ap2)], arrow_head=True, color=color, thickness=t)
    O_1_2 = Line([(0,0,0), tuple(-ap1-ap2)], arrow_head=True, color=color2, thickness=t)


    O23 = Line([(0,0,0), tuple(ap2+ap3)], arrow_head=True, color=color2, thickness=t)
    O2_3 = Line([(0,0,0), tuple(ap2-ap3)], arrow_head=True, color=color, thickness=t)
    O_23 = Line([(0,0,0), tuple(-ap2+ap3)], arrow_head=True, color=color, thickness=t)
    O_2_3 = Line([(0,0,0), tuple(-ap2-ap3)], arrow_head=True, color=color2, thickness=t)

    O123 = Line([(0,0,0), tuple(ap1+ap2+ap3)], arrow_head=True, color=color3, thickness=t)
    O_123 = Line([(0,0,0), tuple(-ap1+ap2+ap3)], arrow_head=True, color=color, thickness=t)
    O1_23 = Line([(0,0,0), tuple(ap1-ap2+ap3)], arrow_head=True, color=color, thickness=t)
    O12_3 = Line([(0,0,0), tuple(ap1+ap2-ap3)], arrow_head=True, color=color, thickness=t)
    O_1_23 = Line([(0,0,0), tuple(-ap1-ap2+ap3)], arrow_head=True, color=color, thickness=t)
    O1_2_3 = Line([(0,0,0), tuple(ap1-ap2-ap3)], arrow_head=True, color=color, thickness=t)
    O_12_3 = Line([(0,0,0), tuple(-ap1+ap2-ap3)], arrow_head=True, color=color, thickness=t)
    O_1_2_3 = Line([(0,0,0), tuple(-ap1-ap2-ap3)], arrow_head=True, color=color3, thickness=t)


    #O1_label = text3d("a_1", tuple((dist+1.3)*ap1), color=color)
    #O2_label = text3d("a_2", tuple((dist+1.3)*ap2), color=color)
    #O3_label = text3d("a_3", tuple((dist+1.3)*ap3), color=color)
    #P = show(Polyhedron(vertices=verts))
    #P = point3d(verts, size=20)
    return (O1 + O_1 + O2 + O_2 + O3 + O_3 +
    O12 + O1_2 + O_12 + O_1_2 +
    O23 + O2_3 + O_23 + O_2_3 +
    O123 + O_123 + O1_23 + O12_3 +
    O_1_23 + O1_2_3 + O_12_3 + O_1_2_3
    + labels
    ) #returns lines to point and their labels

def simple_root_plot(scale=1): #doesn't need to be changed since we have the same simple roots
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

    return O1+O2+O3

def normalPlane(vec, dist): #plots a vector?? I think?
    tup = tuple(vec)
    plane = x*vec[0] + y*vec[1] + z*vec[2] == 0 #boolean, true if they're all 0, false otherwise
    return implicit_plot3d(plane, (x, -dist, dist), (y,-dist,dist), (z,-dist,dist),opacity=.5)

def point_plot(dist, mu, sigmas, color, size=15): #chack the change into omega coordinates part
    if(isinstance(color,basestring)):
        color = tuple(colors[color])
    else:
        color = color
    # Get the xyz coordinates
    c1_ = mu[0]
    c2_ = mu[1]
    c3_ = mu[2]
    sub1dict = dict(sub_1_result)
    coords_xyz = [
            (x_, y_, z_)
            for x_ in range(-dist, dist) for y_ in range(-dist, dist)
            for z_ in range(-dist, dist)
            if all( [
                sub_1_callable[list(sub1dict.keys()).index(s)][j](x_,y_,z_,c1_,c2_,c3_) >= 0
                for j in range(0,3) for s in sigmas])
        ]

    # Substitute in m,n,k
    coords_mnk = [tuple( [xyz_to_mnk[j](x_,y_,z_,c1_,c2_,c3_) for j in [0,1,2]] )
        for (x_,y_,z_) in coords_xyz]

    # Transform into omega coordinates
    coords_mnk_omega = [m_*wp1 + n_*wp2 + k_*wp3 for (m_,n_,k_) in coords_mnk]
    points = point3d(coords_mnk_omega, size, color=color, opacity=.5)
    return points

'''
Finds the Weyl alternation set for lam = mw_1 + nw_2 + kw_3
and mu = c_1w_1 + c_2w_2 + c_3w+3
'''
def alternation_set(lam, mu=(0,0,0)): #may need changing, but idk
    (x_, y_, z_) = lam
    (c1_, c2_, c3_) = mu
    alt_set = set()
    sub1dict = dict(sub_1_result) #defined in weightmult, sub_1_result = weyl_actions_sub(*sub_1), but  idk what that means
    for s in W: #here if it's not 0 we want to add it to the alt set
        if all( [
            sub_1_callable_dict[s][j](x_,y_,z_,c1_,c2_,c3_) >= 0
            for j in range(0,3)]): #sub_1_callable_dict = dict([(p[0], [fast_callable(p[1][i][0], vars=[x,y,z,c1,c2,c3]) for i in range(0,3)]) for p in sub_1_result])
            print((s,[sub_1_callable[list(sub1dict.keys()).index(s)][j](x_,y_,z_,c1_,c2_,c3_) for j in [0,1,2]]))
            alt_set.add(s)
    return alt_set

''' rewritten later on?
def point_plot_fade_only(dist, mu, sigmas, color, size=15):
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
    coords_xyz = [
            (x_, y_, z_)
            for x_ in range(-dist, dist) for y_ in range(-dist, dist)
            for z_ in range(-dist, dist)
            if all( [
                sub_1_callable[list(sub1dict.keys()).index(s)][j](x_,y_,z_,c1_,c2_,c3_) >= 0
                for j in range(0,3) for s in sigmas])
            and all( [ # For every other sigma
                any([sub_1_callable[list(sub1dict.keys()).index(s)][j](x_,y_,z_,c1_,c2_,c3_) < 0 for j in range(0,3)]) # this sigma does not contribute
                for s in others])
        ]

    # A * [alpha coeffs vector] = [omega coeffs vector]
    A = matrix([[3/4, 1/2, 1/4], [1/2, 1, 1/2], [1/4, 1/2, 3/4]])

    # Substitute in m,n,k
    coeffs_mnk_omega = [tuple( [xyz_to_mnk[j](x_,y_,z_,c1_,c2_,c3_) for j in [0,1,2]] )
        for (x_,y_,z_) in coords_xyz]
    coeffs_mnk_alpha = [A.solve_right(vector(omega_coeffs))
            for omega_coeffs in coeffs_mnk_omega]

    # Transform into omega coordinates
    #coords_mnk_omega = [m_*wp1 + n_*wp2 + k_*wp3 for (m_,n_,k_) in coords_mnk]
    coords_mnk_alpha = [m_*ap1 + n_*ap2 + k_*ap3 for (m_,n_,k_) in coeffs_mnk_alpha]
    print(coords_mnk_alpha)
    if (not coords_mnk_alpha):
        print('empty')
        return 0;
    #coords_mnk_omega = [m_*wp1 + n_*wp2 + k_*wp3 for (m_,n_,k_) in coords_mnk if m_>=0 and n_>=0 and k_>=0]
    max_z = max([pt[2] for pt in coords_mnk_alpha])
    min_z = min([pt[2] for pt in coords_mnk_alpha])
    range_z = max_z-min_z
    if (range_z == 0):
        range_z = 1
    points = [point3d(pt, size, color=tuple(col[j] * (max_z-pt[2]+0.1)/range_z for j in range(0,3)), opacity=.5) for pt in coords_mnk_alpha]
    return sum(points)
'''

color_dict = { #what is this a set of???
        frozenset([e]) : (0,0,0),
        frozenset([e,s2]) : (1,0,0),
        frozenset([e,s2,s3]) : (0,1,0),
        frozenset([e,s1,s2]) : (0,0,1),
        frozenset([e,s2,s3,s2*s3]) : (1,1,0),
        frozenset([e,s1,s3,s3*s1]) : (0,1,1),
        frozenset([e,s1,s2,s2*s1]) : (.25,1,1),
        frozenset([e,s1,s2,s3,s3*s1]) : (1,1,.25),
        frozenset([e,s1,s2,s3,s2*s3,s3*s1]) : (1,.25,.25),
        frozenset([e,s1,s2,s3,s2*s1,s3*s1]) : (.25,.5,1),
        frozenset([e,s2,s3,s2*s3,s3*s2,s2*s3*s2]) : (.25,1,.25),
        frozenset([e,s1,s2,s1*s2,s2*s1,s1*s2*s1]) : (1,.25,.75),
        }
pos_lambda_alts = color_dict.keys()

def plot_only(sigmas, size=15, dist=20): #plots the sigmas?
    sigmas_set = frozenset(sigmas)
    print(sigmas)
    point_plot_fade_only(dist=dist, mu=(0,0,0), sigmas=sigmas,
            color=color_dict[sigmas_set], size=size)

def print_only(sigmas, dist=20): #check substitute mnk and omega coeff
    mu = (0,0,0)
    others = set(W).difference(sigmas)
    c1_ = mu[0]
    c2_ = mu[1]
    c3_ = mu[2]
    sub1dict = dict(sub_1_result)
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
    #coeffs_mnk_alpha = [(m_, n_, k_) for (m_,n_,k_) in coeffs_mnk_alpha if m_ >= 0 and n_>=0 and k_>=0]

    # Transform into omega coordinates
    #coords_mnk_omega = [m_*wp1 + n_*wp2 + k_*wp3 for (m_,n_,k_) in coords_mnk]
    coords_mnk_alpha = [m_*ap1 + n_*ap2 + k_*ap3 for (m_,n_,k_) in coeffs_mnk_alpha]
    coords_mnk_omega = [m_*wp1 + n_*wp2 + k_*wp3 for (m_,n_,k_) in coeffs_mnk_omega]
    coords_mnk_alpha_nn = [m_*ap1 + n_*ap2 + k_*ap3
        for (m_,n_,k_) in coeffs_mnk_alpha_nn]
    coords_mnk_omega_nn = [m_*wp1 + n_*wp2 + k_*wp3
        for (m_,n_,k_) in coeffs_mnk_omega_nn]
    pts = coeffs_mnk_alpha_nn
    print(pts)

def point_plot_fade_only(dist, mu, sigmas, color, size=15): #check turn to mnk and omega coeff
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
    #coeffs_mnk_alpha = [(m_, n_, k_) for (m_,n_,k_) in coeffs_mnk_alpha if m_ >= 0 and n_>=0 and k_>=0]

    # Transform into omega coordinates
    #coords_mnk_omega = [m_*wp1 + n_*wp2 + k_*wp3 for (m_,n_,k_) in coords_mnk]
    coords_mnk_alpha = [m_*ap1 + n_*ap2 + k_*ap3 for (m_,n_,k_) in coeffs_mnk_alpha]
    coords_mnk_omega = [m_*wp1 + n_*wp2 + k_*wp3 for (m_,n_,k_) in coeffs_mnk_omega]
    coords_mnk_alpha_nn = [m_*ap1 + n_*ap2 + k_*ap3
        for (m_,n_,k_) in coeffs_mnk_alpha_nn]
    coords_mnk_omega_nn = [m_*wp1 + n_*wp2 + k_*wp3
        for (m_,n_,k_) in coeffs_mnk_omega_nn]
    pts = coeffs_mnk_alpha_nn
    print(len(pts))
    pts_plot = coords_mnk_alpha_nn
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
    return show(sum(points) + simple_root_plot(range_z), frame=False)

def point_plot_reversed_old(dist, mu, sigmas, color, size=10): #check mnk and omega coeff
    # Get the xyz coordinates
    c1_ = mu[0]
    c2_ = mu[1]
    c3_ = mu[2]
    sub1dict = dict(sub_1_result)
    coords_xyz = [
            (x_, y_, z_)
            for x_ in range(-dist, dist) for y_ in range(-dist, dist)
            for z_ in range(-dist, dist)
            if all( [
                sub_1_callable[list(sub1dict.keys()).index(s)][j](x_,y_,z_,c1_,c2_,c3_) < 0
                for j in range(0,3) for s in sigmas])
        ]

    # Substitute in m,n,k
    coords_mnk = [tuple( [xyz_to_mnk[j](x_,y_,z_,c1_,c2_,c3_) for j in [0,1,2]] )
        for (x_,y_,z_) in coords_xyz]

    # Transform into omega coordinates
    coords_mnk_omega = [m_*wp1 + n_*wp2 + k_*wp3 for (m_,n_,k_) in coords_mnk]
    points = point3d(coords_mnk_omega, size, color=color, opacity=.5)
    return points

def point_plot_reversed_nofade(dist, mu, sigmas, color, size=10):
    # Get the xyz coordinates
    c1_ = mu[0]
    c2_ = mu[1]
    c3_ = mu[2]
    sub1dict = dict(sub_1_result)
    coords_xyz = [
            (x_, y_, z_)
            for x_ in range(-dist, dist) for y_ in range(-dist, dist)
            for z_ in range(-dist, dist)
            if all( [
                any([
                sub_1_callable[list(sub1dict.keys()).index(s)][j](x_,y_,z_,c1_,c2_,c3_) < 0
                for j in range(0,3)
                ])
                for s in sigmas])
        ]

    # Substitute in m,n,k
    coords_mnk = [tuple( [xyz_to_mnk[j](x_,y_,z_,c1_,c2_,c3_) for j in [0,1,2]] )
        for (x_,y_,z_) in coords_xyz]

    # Transform into omega coordinates
    coords_mnk_omega = [m_*wp1 + n_*wp2 + k_*wp3 for (m_,n_,k_) in coords_mnk]
    points = point3d(coords_mnk_omega, size, color=color, opacity=.5)
    return points

def point_plot_reversed(dist, mu, sigmas, color, size=10): #check mnk and omega coeff
    if(isinstance(color,basestring)):
        col = tuple(colors[color])
    else:
        col = color
    # Get the xyz coordinates
    c1_ = mu[0]
    c2_ = mu[1]
    c3_ = mu[2]
    sub1dict = dict(sub_1_result)
    coords_xyz = [
            (x_, y_, z_)
            for x_ in range(-dist, dist) for y_ in range(-dist, dist)
            for z_ in range(-dist, dist)
            if all( [
                any([
                sub_1_callable[list(sub1dict.keys()).index(s)][j](x_,y_,z_,c1_,c2_,c3_) < 0
                for j in range(0,3)
                ])
                for s in sigmas])
        ]

    # Substitute in m,n,k
    coords_mnk = [tuple( [xyz_to_mnk[j](x_,y_,z_,c1_,c2_,c3_) for j in [0,1,2]] )
        for (x_,y_,z_) in coords_xyz]

    # Transform into omega coordinates
    coords_mnk_omega = [m_*wp1 + n_*wp2 + k_*wp3 for (m_,n_,k_) in coords_mnk]
    max_z = max([pt[2] for pt in coords_mnk_omega])
    min_z = min([pt[2] for pt in coords_mnk_omega])
    range_z = max_z-min_z
    print(col)
    points = [point3d(pt, size, color=tuple(min(1,col[j] * max(0.3,(max_z-pt[2]+0.3)/(range_z))) for j in range(0,3)), opacity=1) for pt in coords_mnk_omega]
    return sum(points)

def point_plot_reversed_polytope(dist, mu, sigmas, color, size=10): #check mnk and omega coeff
    # Get the xyz coordinates
    c1_ = mu[0]
    c2_ = mu[1]
    c3_ = mu[2]
    sub1dict = dict(sub_1_result)
    coords_xyz = [
            (x_, y_, z_)
            for x_ in range(-dist, dist) for y_ in range(-dist, dist)
            for z_ in range(-dist, dist)
            if any( [
                sub_1_callable[list(sub1dict.keys()).index(s)][j](x_,y_,z_,c1_,c2_,c3_) == 0
                for j in range(0,3) for s in sigmas])
        ]

    # Substitute in m,n,k
    coords_mnk = [tuple( [xyz_to_mnk[j](x_,y_,z_,c1_,c2_,c3_) for j in [0,1,2]] )
        for (x_,y_,z_) in coords_xyz]

    # Transform into omega coordinates
    coords_mnk_omega = [m_*wp1 + n_*wp2 + k_*wp3 for (m_,n_,k_) in coords_mnk]
    coords_mnk_omega_rdf = [tuple(RDF(j) for j in i) for i in coords_mnk_omega]
    #points = point3d(coords_mnk_omega, size, color=color, opacity=.5)
    polytope = Polyhedron(vertices=coords_mnk_omega_rdf)
    return plot(polytope)

def center_plot(dist, mu, color, size=25):
    return point_plot_reversed(dist, mu, weyl_group, color, size)

def center_polytope_pts_plot(dist, mu, color, size=25):
    return point_plot_reversed(dist, mu, weyl_group, color, size)

#def center_polytope_plot(mu):
#    pts = getPolytopePts3d(mu)
#    pts_ = [[RDF(j) for j in v] for v in verts]
#    print(pts_)
#    return plot(Polyhedron(pts_))

def center_polytope_plot(mu, poly=False):
    #mu given as linear combo of w1, w2, w3
    #convert to actual coords
    polytope = getPolytope3d(mu)
    verts = polytope.vertices()
    verts = [tuple(v[0] * wp1 + v[1] * wp2 + v[2] * wp3) for v in verts]
    verts = [[RDF(j) for j in v] for v in verts]
    new_polytope = Polyhedron(verts)
    if(poly):
        return new_polytope
    return plot(new_polytope)

#def center_polytope_verts(mu, dist=10):
#    return (omega_plot(dist=dist, color='black')
#        + )

def lattice_plot(dist):
    lines_xy = [
        Line( [ tuple(x_*wp1 + y_*wp2 - dist*wp3),
                tuple(x_*wp1 + y_*wp2 + dist*wp3) ],
              color=(.7,.7,.7),
              opacity=.1,
              thickness=.2)
        for x_ in range(-dist, dist)
        for y_ in range(-dist, dist)
    ]
    lines_xz = [
        Line( [ tuple(x_*wp1 - dist*wp2 + z_*wp3),
                tuple(x_*wp1 + dist*wp2 + z_*wp3) ],
              color=(.7,.7,.7),
              opacity=.1,
              thickness=.2)
        for x_ in range(-dist, dist)
        for z_ in range(-dist, dist)
    ]
    lines_yz = [
        Line( [ tuple(-dist*wp1 + y_*wp2 + z_*wp3),
                tuple(dist*wp1 + y_*wp2 + z_*wp3) ],
              color=(.7,.7,.7),
              opacity=.1,
              thickness=.2)
        for y_ in range(-dist, dist)
        for z_ in range(-dist, dist)
    ]
    return sum(lines_xy) + sum(lines_xz)

def init_diagrams():
    # Create omega column matrices over Q
    global w1_, w2_, w3_
    w1_ = vector(w1.transpose())
    w2_ = vector(w2.transpose())
    w3_ = vector(w3.transpose())

    global a1_, a2_, a3_
    a1_ = vector(a1.transpose())
    a2_ = vector(a2.transpose())
    a3_ = vector(a3.transpose())

    # Matrix whose rows are w1, w2, w3
    global omegas
    omegas = matrix([w1_, w2_, w3_], ring=QQ)

    # Matrix whose rows are a1, a2, a3
    global alphas
    alphas = matrix([a1_, a2_, a3_], ring=QQ)

    # bb is projection basis for omegas
    global bb, bbMat
    bb = gram_schmidt_symb(omegas)
    bbMat = matrix(bb).transpose()

    # alpha_projection is projection basis for alphas
    global alpha_projection_rows, alpha_projection_cols
    alpha_projection_rows = gram_schmidt_symb(alphas)
    alpha_projection_cols = matrix(alpha_projection_rows).transpose()

    global wp1, wp2, wp3
    wp1 = alpha_projection_cols.solve_right(w1_)
    wp2 = alpha_projection_cols.solve_right(w2_)
    wp3 = alpha_projection_cols.solve_right(w3_)

    global ap1, ap2, ap3
    ap1 = alpha_projection_cols.solve_right(a1_)
    ap2 = alpha_projection_cols.solve_right(a2_)
    ap3 = alpha_projection_cols.solve_right(a3_)

def gram_schmidt_symb(M):
    return [v.normalized() for v in M.gram_schmidt()[0].rows()]

def rl2():
    print('Reloading diagrams.sage...')
    load('diagrams.sage')

init_diagrams()
