#P=[p1,p3,p4]
#Q=[q1,q4,q5,q6]
#R=[r1,r3,r4]
[s1,s2,s3]=WeylGroup(['A',3],prefix='s').simple_reflections()


def find_subset(P,Q,R):
    subset=set()
    dict={}
    if(P[0] and Q[0] and R[0]):
        subset.add(1)
        dict.update({1:'P[0],Q[0],R[0]'})
    if(P[1] and Q[1] and R[0]):
        subset.add(s1*s2*s1)
        dict.update({s1*s2*s1:'P[1],Q[1],R[0]'})
    if(P[0] and Q[0] and R[2]):
        subset.add(s3)
        dict.update({s3:'P[0],Q[0],R[2]'})
    if(P[1] and Q[3] and R[1]):
        subset.add(s3*s1*s2)
        dict.update({s3*s1*s2:'P[1],Q[3],R[1]'})
    if(P[1] and Q[3] and R[0]):
        subset.add(s1*s2)
        dict.update({s1*s2:'P[1],Q[3],R[0]'})
    if(P[0] and Q[2] and R[2]):
        subset.add(s2*s3)
        dict.update({s2*s3:'P[0],Q[2],R[2]'})
    if(P[2] and Q[0] and R[0]):
        subset.add(s1)
        dict.update({s1:'P[2],Q[0],R[0]'})
    if(P[0] and Q[2] and R[1]):
        subset.add(s2*s3*s2)
        dict.update({s2*s3*s2:'P[0],Q[2],R[1]'})
    if(P[2] and Q[1] and R[0]):
        subset.add(s2*s1)
        dict.update({s2*s1:'P[2],Q[1],R[0]'})
    if(P[2] and Q[0] and R[2]):
        subset.add(s3*s1)
        dict.update({s3*s1:'P[2],Q[0],R[2]'})
    if(P[0] and Q[3] and R[1]):
        subset.add(s3*s2)
        dict.update({s3*s2:'P[0],Q[3],R[1]'})
    if(P[0] and Q[3] and R[0]):
        subset.add(s2)
        dict.update({s2:'P[0],Q[3],R[0]'})
    return subset


def all_subsets():
    P=[False,False,False]
    Q=[False,False,False,False]
    R=[False,False,False]

    P_ind=[0,1,2]
    Q_ind=[0,1,2,3]
    R_ind=[0,1,2]

    P_s=list(powerset(P_ind))
    Q_s=list(powerset(Q_ind))
    R_s=list(powerset(R_ind))

    #S=set()
    S=[]
    for p in P_s:
        for q in Q_s:
            for r in R_s:
                for i in p:
                    P[i]=True
                for j in q:
                    Q[j]=True
                for k in r:
                    R[k]=True
                #S.add(find_subset(P,Q,R))
                S.append( (find_subset(P,Q,R),P,Q,R) )
                R=[False,False,False]
            Q=[False,False,False,False]
        P=[False,False,False]
    return S


def rl():
    load('coef.sage')
