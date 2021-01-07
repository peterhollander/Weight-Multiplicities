load('diagrams.sage')

#do this with all 48 elements
'''
sigmas added:
# s3s2s3 (after o)
# s3s2s3s1, s3s2s3s2, s3s1s2s3, (after t)
# s2s3s1s2s3, s3s1s2s3s1, s3s1s2s3s2, s3s2s3s2s1, s3s2s3s1s2, (after w)
# s2s3s1s2s3s1, s2s3s1s2s3s2,
# s3s1s2s3s1s2, s3s1s2s3s2s1, s3s2s3s1s2s1, s3s2s3s1s2s3,
# s2s3s1s2s3s2s1, s2s3s1s2s3s1s2, s3s2s3s1s2s3s2, s3s1s2s3s1s2s1, s3s2s3s1s2s3s1,
# s2s3s1s2s3s1s2s1, s3s2s3s1s2s3s2s1, s3s2s3s1s2s3s1s2,
# s3s2s3s1s2s3s1s2s1 (after x)
'''


#THIS IS WHERE THE 24 ELEMENTS START  #changed
# a=sigmas=[e]
def diagram_e():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[e], color='red'), frame=False)

#b=sigmas=[s1]
def diagram_s1():
    return show(omega_plot(dist=20,color='black')
        + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s1], color='blue'), frame=False)

#c=sigmas=[s2]
def diagram_s2():
    return show(omega_plot(dist=20,color='black')
        + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s2], color='green'), frame=False
        )

#d=sigmas=[s3]
def diagram_s3():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s3], color='yellow'), frame=False
    )

#e=sigmas=[s1*s2]
def diagram_s1s2():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s1*s2], color='orange'), frame=False
    )

#f=sigmas=[s2*s1]
def diagram_s2s1():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s2*s1], color='purple'), frame=False
    )

#g=sigmas=[s3*s1]
def diagram_s3s1():
    return show(omega_plot(dist=20,color='black')
        + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s3*s1], color='blueviolet'), frame=False
        )

#h=sigmas=[s2*s3]
def diagram_s2s3():
    return show(omega_plot(dist=20,color='black')
        + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s2*s3], color='cadetblue'), frame=False
        )

#i=sigmas=[s3*s2] --
def diagram_s3s2():
    return show(omega_plot(dist=20,color='black')
        + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s3*s2], color='cyan'),frame=False
        )

#j=sigmas=[s1*s2*s3] --
def diagram_s1s2s3():
    return show(omega_plot(dist=20,color='black')
        + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s1*s2*s3], color='pink'), frame=False
        )

#k=sigmas=[s1*s2*s1] --
def diagram_s1s2s1():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s1*s2*s1], color='limegreen'), frame=False
    )

#l=sigmas=[s3*s2*s1] --
def diagram_s3s2s1():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s3*s2*s1], color='fuchsia'), frame=False
    )

#m=sigmas=[s3*s1*s2] --
def diagram_s3s1s2():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s3*s1*s2], color='gold'), frame=False
    )

#n=sigmas=[s2*s3*s1] --
def diagram_s2s3s1():
    return show(omega_plot(dist=10,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s2*s3*s1], color='greenyellow'), frame=False
    )

#o=sigmas=[s2*s3*s2] --
def diagram_s2s3s2():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s2*s3*s2], color='papayawhip'), frame=False
    )

#p=sigmas=[s1*s2*s3*s1] --
def diagram_s1s2s3s1():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s1*s2*s3*s1], color='lavender'), frame=False
    )

#q=sigmas=[s1*s2*s3*s2] --
def diagram_s1s2s3s2():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s1*s2*s3*s2], color='lightcoral'), frame=False
    )

#r=sigmas=[s2*s3*s1*s2] --
def diagram_s2s3s1s2():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s2*s3*s1*s2], color='navy'), frame=False
    )

#s=sigmas=[s2*s3*s2*s1] --
def diagram_s2s3s2s1():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s2*s3*s2*s1], color='olive'), frame=False
    )

#t=sigmas=[s3*s1*s2*s1] --
def diagram_s3s1s2s1():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s3*s1*s2*s1], color='plum'), frame=False
    )

#u=sigmas=[s1*s2*s3*s2*s1] --
def diagram_s1s2s3s2s1():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s1*s2*s3*s2*s1], color='tan'), frame=False
    )

#v=sigmas=[s1*s2*s3*s1*s2] --
def diagram_s1s2s3s1s2():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s1*s2*s3*s1*s2], color='tomato'), frame=False
    )

#w=sigmas=[s2*s3*s1*s2*s1] --
def diagram_s2s3s1s2s1():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s2*s3*s1*s2*s1], color='dimgrey'), frame=False
    )

#x=sigmas=[s1*s2*s3*s1*s2*s1] --
def diagram_s1s2s3s1s2s1():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s1*s2*s3*s1*s2*s1], color='lightgrey'), frame=False
    )

#THE REST OF THE 48  (FIND COLORS FOR THEMMMM)

#y=sigmas=[s3*s2*s3] (after o)
def diagram_s3s2s3():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s3*s2*s3], color= 'lightyellow'), frame=False
    )

#z=sigmas=[s3*s2*s3*s1] (after t)
def diagram_s3s2s3s1():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s3*s2*s3*s1], color=''), frame=False
    )

#aa=sigmas=[s3*s2*s3*s2] (after t)
def diagram_s3s2s3s2():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s3*s2*s3*s2], color=''), frame=False
    )

#bb=sigmas=[s3*s1*s2*s3] (after t)
def diagram_s3s1s2s3():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s3*s1*s2*s3], color=''), frame=False
    )

#cc=sigmas=[s2*s3*s1*s2*s3] (after w)
def diagram_s2s3s1s2s3():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s2*s3*s1*s2*s3], color=''), frame=False
    )

#dd=sigmas=[s3*s1*s2*s3*s1] (after w)
def diagram_s3s1s2s3s1():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s3*s1*s2*s3*s1], color=''), frame=False
    )

#ee=sigmas=[s3*s1*s2*s3*s2] (after w)
def diagram_s3s1s2s3s2():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s3*s1*s2*s3*s2], color=''), frame=False
    )

#ff=sigmas=[s3*s2*s3*s2*s1] (after w)
def diagram_s3s2s3s2s1():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s3*s2*s3*s2*s1], color=''), frame=False
    )

#gg=sigmas=[s3*s2*s3*s1*s2] (after w)
def diagram_s3s2s3s1s2():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s3*s2*s3*s1*s2], color=''), frame=False
    )

#hh=sigmas=[s2*s3*s1*s2*s3*s1] (after x)
def diagram_s2s3s1s2s3s1():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s2*s3*s1*s2*s3*s1], color=''), frame=False
    )

#ii=sigmas=[s2*s3*s1*s2*s3*s2] (after x)
def diagram_s2s3s1s2s3s2():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s2*s3*s1*s2*s3*s2], color=''), frame=False
    )

#jj=sigmas=[s3*s1*s2*s3*s1*s2] (after x)
def diagram_s3s1s2s3s1s2():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s3*s1*s2*s3*s1*s2], color=''), frame=False
    )

#kk=sigmas=[s3*s1*s2*s3*s2*s1] (after x)
def diagram_s3s1s2s3s2s1():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s3*s1*s2*s3*s2*s1], color=''), frame=False
    )

#ll=sigmas=[s3*s2*s3*s1*s2*s1] (after x)
def diagram_s3s2s3s1s2s1():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s3*s2*s3*s1*s2*s1], color=''), frame=False
    )

#mm=sigmas=[s3*s2*s3*s1*s2*s3] (after x)
def diagram_s3s2s3s1s2s3():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s3*s2*s3*s1*s2*s3], color=''), frame=False
    )

#nn=sigmas=[s2*s3*s1*s2*s3*s2*s1] (after x)
def diagram_s2s3s1s2s3s2s1():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s2*s3*s1*s2*s3*s2*s1], color=''), frame=False
    )

#oo=sigmas=[s2*s3*s1*s2*s3*s1*s2] (after x)
def diagram_s2s3s1s2s3s1s2():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s2*s3*s1*s2*s3*s1*s2], color=''), frame=False
    )

#pp=sigmas=[s3*s2*s3*s1*s2*s3*s2] (after x)
def diagram_s3s2s3s1s2s3s2():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s3*s2*s3*s1*s2*s3*s2], color=''), frame=False
    )

#qq=sigmas=[s3*s1*s2*s3*s1*s2*s1] (after x)
def diagram_s3s1s2s3s1s2s1():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s3*s1*s2*s3*s1*s2*s1], color=''), frame=False
    )

#rr=sigmas=[s3*s2*s3*s1*s2*s3*s1] (after x)
def diagram_s3s2s3s1s2s3s1():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s3*s2*s3*s1*s2*s3*s1], color=''), frame=False
    )

#ss=sigmas=[s2*s3*s1*s2*s3*s1*s2*s1] (after x)
def diagram_s2s3s1s2s3s1s2s1():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s2*s3*s1*s2*s3*s1*s2*s1], color=''), frame=False
    )

#tt=sigmas=[s3*s2*s3*s1*s2*s3*s2*s1] (after x)
def diagram_s3s2s3s1s2s3s2s1():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s3*s2*s3*s1*s2*s3*s2*s1], color=''), frame=False
    )

#uu=sigmas=[s3*s2*s3*s1*s2*s3*s1*s2] (after x)
def diagram_s3s2s3s1s2s3s1s2():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s3*s2*s3*s1*s2*s3*s1*s2], color=''), frame=False
    )

#vv=sigmas=[s3*s2*s3*s1*s2*s3*s1*s2*s1] (after x)
def diagram_s3s2s3s1s2s3s1s2s1():
    return show(omega_plot(dist=20,color='black')
    + point_plot_fade(dist=10, mu=(0,0,0), sigmas=[s3*s2*s3*s1*s2*s3*s1*s2*s1], color=''), frame=False
    )

#the following give the outer shape of all the 24 weyl elements
mu_= (6,-3,0)
def cool_pic2(): #changed
    return show(omega_plot(dist=10,color='black')
        + point_plot(dist=10, mu=mu_, sigmas=[e], color='red')
        + point_plot(dist=10, mu=mu_, sigmas=[s1], color='blue')
        + point_plot(dist=10, mu=mu_, sigmas=[s2], color='green')
        + point_plot(dist=10, mu=mu_, sigmas=[s3], color='yellow')
        + point_plot(dist=10, mu=mu_, sigmas=[s1*s2], color='orange')
        + point_plot(dist=10, mu=mu_, sigmas=[s2*s1], color='purple')
        + point_plot(dist=10, mu=mu_, sigmas=[s3*s1], color='blueviolet')
        + point_plot(dist=10, mu=mu_, sigmas=[s2*s3], color='cadetblue')
        + point_plot(dist=10, mu=mu_, sigmas=[s3*s2], color='cyan')
        + point_plot(dist=10, mu=mu_, sigmas=[s1*s2*s3], color='pink')
        + point_plot(dist=10, mu=mu_, sigmas=[s1*s2*s1], color='limegreen')
        + point_plot(dist=10, mu=mu_, sigmas=[s3*s2*s1], color='fuchsia')
        + point_plot(dist=10, mu=mu_, sigmas=[s3*s1*s2], color='gold')
        + point_plot(dist=10, mu=mu_, sigmas=[s2*s3*s1], color='greenyellow')
        + point_plot(dist=10, mu=mu_, sigmas=[s2*s3*s2], color='papayawhip')
        + point_plot(dist=10, mu=mu_, sigmas=[s1*s2*s3*s1], color='lavender')
        + point_plot(dist=10, mu=mu_, sigmas=[s1*s2*s3*s2], color='lightcoral')
        + point_plot(dist=10, mu=mu_, sigmas=[s2*s3*s1*s2], color='navy')
        + point_plot(dist=10, mu=mu_, sigmas=[s2*s3*s2*s1], color='olive')
        + point_plot(dist=10, mu=mu_, sigmas=[s3*s1*s2*s1], color='plum')
        + point_plot(dist=10, mu=mu_, sigmas=[s1*s2*s3*s2*s1], color='tan')
        + point_plot(dist=10, mu=mu_, sigmas=[s1*s2*s3*s1*s2], color='tomato')
        + point_plot(dist=10, mu=mu_, sigmas=[s2*s3*s1*s2*s1], color='dimgrey')
        + point_plot(dist=10, mu=mu_, sigmas=[s1*s2*s3*s1*s2*s1], color='lightgrey') #--
        + point_plot(dist=10, mu=mu_, sigmas=[s3*s2*s3], color='lightyellow')
        + point_plot(dist=10, mu=mu_, sigmas=[s3*s2*s3*s1], color='')
        + point_plot(dist=10, mu=mu_, sigmas=[s3*s2*s3*s2], color='')
        + point_plot(dist=10, mu=mu_, sigmas=[s3*s1*s2*s3], color='')
        + point_plot(dist=10, mu=mu_, sigmas=[s2*s3*s1*s2*s3], color='')
        + point_plot(dist=10, mu=mu_, sigmas=[s3*s1*s2*s3*s1], color='')
        + point_plot(dist=10, mu=mu_, sigmas=[s3*s1*s2*s3*s2], color='')
        + point_plot(dist=10, mu=mu_, sigmas=[s3*s2*s3*s2*s1], color='')
        + point_plot(dist=10, mu=mu_, sigmas=[s3*s2*s3*s1*s2], color='')
        + point_plot(dist=10, mu=mu_, sigmas=[s2*s3*s1*s2*s3*s1], color='')
        + point_plot(dist=10, mu=mu_, sigmas=[s2*s3*s1*s2*s3*s2], color='')
        + point_plot(dist=10, mu=mu_, sigmas=[s3*s1*s2*s3*s1*s2], color='')
        + point_plot(dist=10, mu=mu_, sigmas=[s3*s1*s2*s3*s2*s1], color='')
        + point_plot(dist=10, mu=mu_, sigmas=[s3*s2*s3*s1*s2*s1], color='')
        + point_plot(dist=10, mu=mu_, sigmas=[s3*s2*s3*s1*s2*s3], color='')
        + point_plot(dist=10, mu=mu_, sigmas=[s2*s3*s1*s2*s3*s2*s1], color='')
        + point_plot(dist=10, mu=mu_, sigmas=[s2*s3*s1*s2*s3*s1*s2], color='')
        + point_plot(dist=10, mu=mu_, sigmas=[s3*s2*s3*s1*s2*s3*s2], color='')
        + point_plot(dist=10, mu=mu_, sigmas=[s3*s1*s2*s3*s1*s2*s1], color='')
        + point_plot(dist=10, mu=mu_, sigmas=[s3*s2*s3*s1*s2*s3*s1], color='')
        + point_plot(dist=10, mu=mu_, sigmas=[s2*s3*s1*s2*s3*s1*s2*s1], color='')
        + point_plot(dist=10, mu=mu_, sigmas=[s3*s2*s3*s1*s2*s3*s2*s1], color='')
        + point_plot(dist=10, mu=mu_, sigmas=[s3*s2*s3*s1*s2*s3*s1*s2], color='')
        + point_plot(dist=10, mu=mu_, sigmas=[s3*s2*s3*s1*s2*s3*s1*s2*s1], color=''), frame=False
        )


#the following gives the shape of the empty space of all the 24 weyl elements
def cool_pic3(mu_a, dist=20, color='red', size=25): #may need changing
    return (omega_plot(dist=dist,color='black')
        + center_plot(dist=dist,
        mu=(2*mu_a[0]-mu_a[1], -mu_a[0] + 2*mu_a[1] - mu_a[2],
        -mu_a[1] + 2*mu_a[2]), color=color, size=size)) #is this our vectors? the positive root ones?

def center_polytope_pts(mu_a, dist=20, color='red', size=25): #may need changing
    return show(omega_plot(dist=dist,color='black')
        + center_polytope_pts_plot(dist=dist,
        mu=(2*mu_a[0]-mu_a[1], -mu_a[0] + 2*mu_a[1] - mu_a[2],
        -mu_a[1] + 2*mu_a[2]), color=color, size=size), frame=False) #same vectors as above

def center_polytope_pts_2(mu, dist=20, color='red', size=25): #this shouldn't need changing
    return (omega_plot(dist=dist,color='black')
        + center_polytope_pts_plot(dist=dist,
        mu=mu, color=color, size=size))


def center_polytope(mu_a, dist=20): #may need changing
    new_mu = mu=(2*mu_a[0]-mu_a[1],
    -mu_a[0] + 2*mu_a[1] - mu_a[2],
    -mu_a[1] + 2*mu_a[2]) #same vectors as above I think
    print(new_mu)
    return show(omega_plot(dist=dist,color='black')
        + center_polytope_plot(new_mu), frame=False
        )
