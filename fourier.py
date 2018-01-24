import math

def shift_x(x_min,x_max,x):
    # shifts x into range -pi to pi, linear transformation
    s = (x-x_min) / (x_max-x_min) #gives how far through the range we are
    r = (2*s - 1) * math.pi # maps onto desired range
    return r

class fourier_series:
    def __init__(self):
        self.k = 0.0
        self.c = []
        self.s = []
        self.x_min = 0.0
        self.x_max = 0.0
    def __call__(self,x):
        nx = shift_x(x_min,x_max,x)
        r = self.k / 2.0
        for i in range(len(self.c)):
            r += self.c[i] * math.cos((i+1)*nx)
        for i in range(len(self.s)):
            r += self.s[i] * math.sin((i+1)*nx)
        return r

class fourier_block:
    def __init__(self,value,x_left,x_right):
        self.value = value
        self.x_left = x_left
        self.x_right = x_right

def fourier(x_min,x_max,fbl,cos_count,sin_count):
    # fbl is a list of fourier_block, order not important
    # cos_count and sin_count are integers
    # setup
    ou = fourier_series()
    ou.x_min = x_min
    ou.x_max = x_max
    # for faster computation, we want to shift all x values in fbl first
    nfbl = []
    for p in fbl:
        nfbl.append( fourier_block( p.value , shift_x(x_min,x_max,p.left) , shift_x(x_min,x_max,p.right) ) )
    # we can now use values in nfbl without shifting them
    # now find k (a_0)
    a = 0.0
    for p in nfbl:
        a += p.value * (p.right - p.left)
    ou.k = a / math.pi
    # now do the cosines
    for n in range(1,cos_count+1):
        a = 0.0
        for p in nfbl:
            a += p.value * (math.sin(n*p.right) - math.sin(n*p.left)) / n
        ou.c.append(a)
    # now do the sines
    for n in range(1,sin_count+1):
        a = 0.0
        for p in nfbl:
            a += p.value * (math.cos(n*p.left) - math.cos(n*p.right)) / n
        ou.s.append(a)
    # and we are done
    return ou
