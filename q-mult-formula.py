if __name__ == "__main__":
  
  lam = input("Enter Lambda (format: (#,#,#)): ")
  mu = input("Enter Mu (format: (#,#,#)): ")
  
  m = lam[1]
  n = lam[3]
  k = lam[5]
  
  x = mu[1]
  y = mu[3]
  z = mu[5]
  
  a = m + n + k - x - y - z
  b = n + k - x - y - z - 1
  c = k - x - y - z - 2
  d = m + 2*n + 2*k - x - 2*y - 2*z
  e = m + n + 2*k - x - 2*y - 2*z - 1
  f = n + 2*k - x - 2*y - 2*z - 2
  g = m + n - x - 2*y - 2*z - 3
  h = n - x - 2*y - 2*z - 4
  i = m - x - 2*y - 2*z - 4
  j = m/2 + n + 3*k/2 - x/2 - y - 3*z/2
  l = m/2 + n + k/2 - x/2 - y - 3*z/2 - 1
  o = m/2 + k/2 - x/2 - y - 3*z/2 - 2
  p = -m/2 + k/2 - x/2 - y - 3*z/2 - 3
  r = m/2 - k/2 - x/2 - y - 3*z/2 - 3
  
  vals = [a,b,c,d,e,f,g,h,i,j,l,o,p,r]
  N = set()
  for val in vals:
    if (int(val) == val) and (val >= 0):
      N = N.union(val)
    
  
  mult = ""
  if N == set(a,b,c,d,e,f,g,h,i,l,r,o,j):
    mult = "A - B - C - D + E + F + G + H + I - J - K - L - N - O + Q"
  elif N == set():
    mult = ""
  elif N == set(a,b,c,d,e,f,g,i,l,r,o,j):
    mult = "A - B - C - D + E + F + G + H + I - J - L - N - O + Q"
  elif N == set(a,b,c,d,e,f,g,i,l,o,p,j):
    mult = "A - B - C - D + E + F + G + H + I - J - L - M - N + P"
  elif N == set(a,b,c,d,e,f,g,h,l,o,p,j):
    mult = "A - B - C - D + E + F + G + H + I - J - K - M - N + P"
  elif N == set(a,b,c,d,e,f,g,l,o,p,j):
    mult = "A - B - C - D + E + F + G + H + I - J - M - N + P"
  elif N == set(a,b,c,d,e,f,g,h,i,l,o,j):
    mult = "A - B - C - D + E + F + G + H + I - J - K - L - N"
  elif N == set(a,b,c,d,e,f,g,i,l,o,j):
    mult = "A - B - C - D + E + F + G + H + I - J - L - N"
  elif N == set(a,b,c,d,e,f,l,o,p,j):
    mult = "A - B - C - D + E + F + H + I - J - M - N + P"
  elif N == set(a,b,c,d,e,f,g,h,l,o,j):
    mult = "A - B - C - D + E + F + G + H + I - J - K - N"
  elif N == set(a,b,d,e,f,g,h,i,l,r,o,j):
    mult = "A - B - C - D + F + G + H + I - K - L - O + Q"
  elif N == set(a,b,c,d,e,f,g,l,o,j):
    mult = "A - B - C - D + E + F + G + H + I - J - N"
  elif N == set(a,b,d,e,f,g,i,l,r,o,j):
    mult = "A - B - C - D + F + G + H + I - L - O + Q"
  elif N == set(a,b,c,d,e,f,g,h,l,j):
    mult = "A - B - C - D + E + F + G + H - J - K"
  elif N == set(a,b,c,d,e,f,l,o,j):
    mult = "A - B - C - D + E + F + H + I - J - N"
  elif N == set(a,b,d,e,g,i,l,r,o,j):
    mult = "A - B - C - D + G + H + I - L - O + Q"
  elif N == set(a,b,d,e,f,g,h,i,l,o,j):
    mult = "A - B - C - D + F + G + H + I - K - L"
  elif N == set(a,b,c,d,e,f,g,l,j):
    mult = "A - B - C - D + E + F + G + H - J"
  elif N == set(a,b,d,e,f,g,i,l,o,j):
    mult = "A - B - C - D + F + G + H + I - L"
  elif N == set(a,b,d,e,f,g,h,l,o,j):
    mult = "A - B - C - D + F + G + H + I - K"
  elif N == set(a,b,c,d,e,f,l,j):
    mult = "A - B - C - D + E + F + H - J"
  elif N == set(a,b,d,e,f,g,h,l,j):
    mult = "A - B - C - D + F + G + H - K"
  elif N == set(a,b,d,e,f,g,l,o,j):
    mult = "A - B - C - D + F + G + H + I"
  elif N == set(a,b,d,e,g,l,o,j):
    mult = "A - B - C - D + G + H + I - L"
  elif N == set(a,d,e,g,i,l,r,o,j):
    mult = "A - C - D + G + I - L - O + Q"
  elif N == set(a,b,d,e,f,g,l,j):
    mult = "A - B - C - D + F + G + H"
  elif N == set(a,b,d,e,f,l,o,j):
    mult = "A - B - C - D + F + H + I"
  elif N == set(a,b,d,e,g,l,o,j):
    mult = "A - B - C - D + G + H + I"
  elif N == set(a,b,c,d,e,f,j):
    mult = "A - B - C + E + F - J"
  elif N == set(a,b,d,e,f,l,j):
    mult = "A - B - C - D + F + H"
  elif N == set(a,b,d,e,g,l,j):
    mult = "A - B - C - D + G + H"
  elif N == set(a,d,e,g,i,l,o,j):
    mult = "A - C - D + G + I - L"
  elif N == set(a,b,d,e,l,o,j):
    mult = "A - B - C - D + H + I"
  elif N == set(a,b,d,e,l,j):
    mult = "A - B - C - D + H"
  elif N == set(a,d,e,g,l,o,j):
    mult = "A - C - D + G + I"
  elif N == set(a,b,d,e,f,j):
    mult = "A - B - C + F"
  elif N == set(a,b,d,l,j):
    mult = "A - B - D + H"
  elif N == set(a,d,e,g,l,j):
    mult = "A - C - D + G"
  elif N == set(a,d,e,l,o,j):
    mult = "A - C - D + I"
  elif N == set(a,b,d,e,j):
    mult = "A - B - C"
  elif N == set(a,d,e,l,j):
    mult = "A - C - D"
  elif N == set(a,b,d,j):
    mult = "A - B"
  elif N == set(a,d,e,j):
    mult = "A - C"
  elif N == set(a,d,l,j):
    mult = "A - D"
  elif N == set(a,d,j):
    mult = "A"
  else:
    mult = "0"
  
print("m_q({},{}) = {}".format(Lambda, Mu, mult))
