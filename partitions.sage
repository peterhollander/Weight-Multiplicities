
#lists partitions
def listPartitions(m,n,k): #needs to be updated for sp_6
    #each element of v has the coefficients in this order : (a1, a2, a3, a1+a2, a2+a3, a1+a2+a3, a1+2a2+a3, 2a1+2a2+a3, 2a2+a3)
    # m-..., d, e, f, g, h, i
    v=[]
    for h in range(0,min(floor(m/2),floor(n/2),k)+1):
        ''' h represents the coefficient of 2al_1+2al_2+al_3
            in our partition '''
        for g in range(0,min(m-2*h,floor((n-2*h)/2),k-h)+1):
            ''' g represents the coefficient of al_1+2al_2+al_3
                in our partition. '''
            for f in range(0, min(m-2*h-g,n-2*h-2*g,k-h-g)+1):
                ''' f represents the coefficient of al_1+al_2+al_3
                    in our partition.'''
                for i in range(0, min(floor((n-2*h-2*g-f)/2),k-h-g-f)+1):
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
                            v.append((m-d-f-g-2*h, n-e-f-2*g-2*h-2*i, k-e-f-g-h-i, d, e, f, g, h, i))
                            #coefficients are in this order: (a1, a2, a3, a1+a2, a2+a3, a1+a2+a3, a1+2a2+a3, 2a1+2a2+a3, 2a2+a3)
                            #i'm pretty sure this is right


    '''#these are the sl_4 partitions:
    for f in range(0,min(m,n,k)+1):
        for d in range(0,min(m-f,n-f)+1):
            for e in range(0, min(n-f-d,k-f)+1):
                v.append((m-d-f,n-d-e-f,k-e-f,d,e,f))
    '''
    return v

#returns kostant's q-analog polynomial
def qPoly(m,n,k): #updated for sp_6
    v=0
    q=var('q')
    for h in range(0,min(floor(m/2),floor(n/2),k)+1):
        ''' h represents the coefficient of 2al_1+2al_2+al_3
            in our partition '''
        for g in range(0,min(m-2*h,floor((n-2*h)/2),k-h)+1):
            ''' g represents the coefficient of al_1+2al_2+al_3
                in our partition. '''
            for f in range(0, min(m-2*h-g,n-2*h-2*g,k-h-g)+1):
                ''' f represents the coefficient of al_1+al_2+al_3
                    in our partition.'''
                for i in range(0, min(floor((n-2*h-2*g-f)/2),k-h-g-f)+1):
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
                            v += q^(m+n+k-d-e-2*f-3*g-4*h-2*i)
    ''' #this is the A_3 polynomial:
    for f in range(0,min(m,n,k)+1):
        for d in range(0,min(m-f,n-f)+1):
            for e in range(0, min(n-f-d,k-f)+1):
                v+=q^(m+n+k-2*f-e-d)
    '''
    return v



#we just dont have cases for sp_6, sadly
'''
# m,k => n
def caseOne(m,n,k,t):
    L=min(floor(t/2), n-ceil(t/2))
    return sum([2*j+t%2+1 for j in range(0,L+1)])


# m => n => k
def caseTwo(m,n,k,t):
    B= min(floor(t/2),k)
    A= max(t-n,0)
    L=B-A
    return sum([min((2*j+(t%2)), k-B+j)+1 for j in range(0,L+1)])
'''


# Reload file in interactive Sage environment
def rl():
    load('partitions.sage')

rl()
