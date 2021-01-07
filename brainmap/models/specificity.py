"""
@name: specificity
@description:
Specficity and basal connectivity models

@author: Christopher Brittin
@email: "cabrittin"+ <at>+ "gmail"+ "."+ "com"
@date: 2020-March
"""

import numpy as np
from scipy.special import comb

class Model:
    def __init__(self,n=5):
       self.Z = np.zeros(n)
       self.space = 0

    def error(self,data):
        return abs(data - self.Z).sum() 

    def min_error(self,error,derror=1):
        idx = np.argsort(error[:,-1])[:derror]
        return error[idx,:]
   
    def get_error(self,D):
        i = 0
        error = np.zeros(self.dim)
        for p in self.loop(D):
            error[i,:] = p
            i += 1
        return error

    def normalize(self):
        norm = 1 - self.Z[0]
        if norm > 0: self.Z[1:] /= norm
 

class Model_P1(Model):
    def __init__(self):
        Model.__init__(self)
        self.dim = (101,2) 
    
    def loop(self,data):
        for self.p in np.linspace(0,1,101):
            self.compute() 
            yield [self.p,self.error(data)]
    
    def set(self,params):
        self.p = params
        self.compute() 
        
    def compute(self):
        p = self.p
        mp = 1 - p
        self.Z[0] = mp**4
        self.Z[1] = 4*(p*(mp**3))
        self.Z[2] = 6*((p**2)*(mp**2))
        self.Z[3] = 4*((p**3)*mp)
        self.Z[4] = p**4

        self.normalize()

class Model_P2(Model):
    def __init__(self):
        Model.__init__(self)
        self.dim = (101**2,3)

    def loop(self,data):
        for self.p in np.linspace(0,1,101):
            for self.q in np.linspace(0,1,101):
                self.compute() 
                yield [self.p,self.q,self.error(data)]
    
    def set(self,params):
        self.p = params[0]
        self.q = params[1]
        self.compute() 
        
    def compute(self):
        p,q = self.p,self.q
        mp = 1-p
        mq = 1-q
        
        self.Z[0] = (mp**4)*(mq**4)
        self.Z[1] = 4*(mp**3)*(mq**3)*(p + q*mp)
        self.Z[2] = 6*(mp**2)*(mq**2)*(p**2 + 2*p*q*mp + (q*mp)**2)
        self.Z[3] = 4*mp*mq*(p**3 + 3*(p**2)*q*mp + 3*p*((q*mp)**2) + (q*mp)**3) 
        self.Z[4] = p**4 + (q**4)*mp**4 

        self.normalize()

class Model_P3_Basal(Model):
    def __init__(self):
        Model.__init__(self)
    
    def error(self,data):
        return abs(data[1:] - self.Z[1:]).sum() 

    def min_error(self,error,derror=1):
        idx = np.argsort(error[:,-1])[:derror]
        return error[idx,:]
   
    def loop(self,data):
        for self.p in np.linspace(0,1,101):
            for self.q in np.linspace(0,1,101):
                for self.f in np.linspace(0,1,101): 
                    self.compute() 
                    yield [self.p,self.q,self.f,self.error(data)]
    
    def set(self,params):
        self.p = params[0]
        self.q = params[1]
        self.f = params[2]
        self.compute() 
    

    def compute(self):
        p,q,f = self.p,self.q,self.f
        mp = 1-p
        mq = 1-q
        mf = 1-f

        self.Z[0] = f*(mp**4)*(mq**4) + mf*(mq**4)
        self.Z[1] = 4*(f*(mp**3)*(mq**3)*(p + q*mp) + mf*q*(mq**3))
        self.Z[2] = 6*(f*(mp**2)*(mq**2)*(p**2 + 2*p*mp*q + (q**2)*(mp**2)) \
                    + mf*(q**2)*(mq**2))
        self.Z[3] = 4*(f*mp*mq*(p**3 + 3*(p**2)*mp*q + 3*p*(mp**2)*(q**2) + (mp**3)*(q**3))\
                    + mf*(q**3)*mq)
        self.Z[4] = f*((p**4) + 4*(p**3)*mp*q + 6*(p**2)*(mp**2)*(q**2) \
                    + 4*p*(mp**3)*(q**3) + (mp**4)*(q**4)) \
                    + mf*(q**4)
        
        self.normalize()

class Model_P3_Avoidance(Model_P3_Basal):
    def __init__(self):
        Model_P3_Basal.__init__(self)

    def compute(self):
        p,q,f = self.p,self.q,self.f
        mp = 1-p
        mq = 1-q
        mf = 1-f
        self.Z[0] = f*(mp**4) + mf*(mq**4)
        self.Z[1] = 4*(f*p*(mp**3) + mf*q*(mq**3))
        self.Z[2] = 6*(f*(p**2)*(mp**2) + mf*(q**2)*(mq**2))
        self.Z[3] = 4*(f*(p**3)*mp + mf*(q**3)*mq)
        self.Z[4] = f*(p**4) + mf*(q**4)
    
        self.normalize()

class Model_P3_Avoidance_Gen(Model_P3_Basal):
    def __init__(self,n=5):
        Model_P3_Basal.__init__(self)
        self.Z = np.zeros(n)
        
    def compute(self):
        p,q,f = self.p,self.q,self.f
        mp = 1-p
        mq = 1-q
        mf = 1-f

        n = len(self.Z) - 1
        for i in range(n+1):
            c = comb(n,i,exact=True)
            prob = c*(f*(p**i)*(mp**(n-i)) + mf*(q**i)*(mq**(n-i)))
            self.Z[i] = prob
        
        self.normalize()

 
