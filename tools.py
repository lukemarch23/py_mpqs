

def modular_sqrt(a, p):
    if legendre_symbol(a, p) != 1:return 0
    elif a == 0:return 0
    elif p == 2:return p
    elif p % 4 == 3:return pow(a, (p + 1) / 4, p)
    s = p - 1
    e = 0
    while s % 2 == 0:
        s /= 2
        e += 1
    n = 2
    while legendre_symbol(n, p) != -1:
        n += 1
    x = pow(a, (s + 1) / 2, p)
    b = pow(a, s, p)
    g = pow(n, s, p)
    r = e
    while True:
        t = b
        m = 0
        for m in xrange(r):
            if t == 1:break
            t = pow(t, 2, p)
        if m == 0:
            return x
        gs = pow(g, 2 ** (r - m - 1), p)
        g = (gs * gs) % p
        x = (x * gs) % p
        b = (b * g) % p
        r = m
def mod_inv(a, m):
  a = int(a%m)
  x, u = 0, 1
  while a:
    x, u = u, x - (m/a)*u
    m, a = a, m%a
  return x
def legendre_symbol(a, p):
    ls = pow(a, (p - 1) / 2, p)
    return -1 if ls == p - 1 else ls
'''
Input: 
    M: 2d binary matrix represented as a 1d array of integers
    H: 2d binary matrix represented as a 1d array of integers
    columnCount: number of columns
Return:
    None (in place updates)
'''
def reduceRowEchelonForm(M,H,columnCount):
    if not M: return
    lead = 0
    rowCount = len(M)
    for r in xrange(rowCount):
        if lead >= columnCount:
            return
        i = r
        while M[i] & (1<<lead)==0:
            i += 1
            if i == rowCount:
                i = r
                lead += 1
                if columnCount == lead:
                    return
        M[i],M[r] = M[r],M[i]
        H[i],H[r] = H[r],H[i]
        for i in xrange(rowCount):
            if i != r and M[i] & (1<<lead):
                M[i]^=M[r]
                H[i]^=H[r]
        lead += 1
 
def isqrt(n):
  c = n*4/3
  d = c.bit_length()
  a = d>>1
  if d&1:
    x = 1 << a
    y = (x + (n >> a)) >> 1
  else:
    x = (3 << a) >> 2
    y = (x + (c >> a)) >> 1
  if x != y:
    x = y
    y = (x + n/x) >> 1
    while y < x:
      x = y
      y = (x + n/x) >> 1
  return x
 
def prime_sieve(n):
    from math import sqrt;from array import array
    ar = array('B')
    [ar.append(0) for x in xrange(n/16+1)]
    primes = [2,3]
    p=1;sq=int(sqrt(n)+1);i=2
    def setv(n): ar[(n>>3)] |= (1<<(n&7))
    while i<=((n-1)>>1):
        x=(i<<1)+1
        if x<sq and not ar[(i>>3)] & (1<<(i&7)): 
                for k in xrange(x,n/x+1,2):setv((k*x-1)>>1)
        if not ( ar[(i>>3)] & (1<<(i&7)) ) :primes.append((i<<1)+1)
        i+=p;p = 2 if p==1 else 1
    return primes

def miller_rabin_pass(a, s, d, n):
  a_to_power = pow(a, d, n)
  if a_to_power == 1:
    return True
  for i in range(s-1):
    if a_to_power == n - 1:
      return True
    a_to_power = (a_to_power * a_to_power) % n
  return a_to_power == n - 1
 
def miller_rabin(n,steps=20):
  from random import randrange
  d = n - 1
  s = 0
  while d % 2 == 0:
    d >>= 1
    s += 1
  for repeat in range(steps):
    a = 0
    while a == 0:
      a = randrange(n)
    if not miller_rabin_pass(a, s, d, n):
      return False
  return True
'''
Not 100% prime, just very likely which may cause problems
with a very small probability
'''
primes = prime_sieve(1000)
def next_prime(n):
    if n==2:return 3
    n+=2
    while 1:
        q=1
        for p in primes:
            if p*p>n:break
            if n%p==0:q= 0;break
        if q==0:n+=2;continue
        else : break
        if miller_rabin(n,10):break
        n+=2
    return n