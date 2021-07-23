if __name__ == "__main__":
  
  lam = input("Enter Lambda (format: (#,#,#)): ")
  mu = input("Enter Mu (format: (#,#,#)): ")
  
  m = int(lam[1])
  n = int(lam[3])
  k = int(lam[5])
  
  x = int(mu[1])
  y = int(mu[3])
  z = int(mu[5])
  
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
  N = {a}
  for val in vals:
    if (int(val) == val) and (val >= 0):
      N.add(val)
    
  
  mult = ""
  if N == {a,b,c,d,e,f,g,h,i,l,r,o,j}:
    mult = "A - B - C - D + E + F + G + H + I - J - K - L - N - O + Q"
  elif N == {}:
    mult = ""
  elif N == {a,b,c,d,e,f,g,i,l,r,o,j}:
    mult = "A - B - C - D + E + F + G + H + I - J - L - N - O + Q"
  elif N == {a,b,c,d,e,f,g,i,l,o,p,j}:
    mult = "A - B - C - D + E + F + G + H + I - J - L - M - N + P"
  elif N == {a,b,c,d,e,f,g,h,l,o,p,j}:
    mult = "A - B - C - D + E + F + G + H + I - J - K - M - N + P"
  elif N == {a,b,c,d,e,f,g,l,o,p,j}:
    mult = "A - B - C - D + E + F + G + H + I - J - M - N + P"
  elif N == {a,b,c,d,e,f,g,h,i,l,o,j}:
    mult = "A - B - C - D + E + F + G + H + I - J - K - L - N"
  elif N == {a,b,c,d,e,f,g,i,l,o,j}:
    mult = "A - B - C - D + E + F + G + H + I - J - L - N"
  elif N == {a,b,c,d,e,f,l,o,p,j}:
    mult = "A - B - C - D + E + F + H + I - J - M - N + P"
  elif N == {a,b,c,d,e,f,g,h,l,o,j}:
    mult = "A - B - C - D + E + F + G + H + I - J - K - N"
  elif N == {a,b,d,e,f,g,h,i,l,r,o,j}:
    mult = "A - B - C - D + F + G + H + I - K - L - O + Q"
  elif N == {a,b,c,d,e,f,g,l,o,j}:
    mult = "A - B - C - D + E + F + G + H + I - J - N"
  elif N == {a,b,d,e,f,g,i,l,r,o,j}:
    mult = "A - B - C - D + F + G + H + I - L - O + Q"
  elif N == {a,b,c,d,e,f,g,h,l,j}:
    mult = "A - B - C - D + E + F + G + H - J - K"
  elif N == {a,b,c,d,e,f,l,o,j}:
    mult = "A - B - C - D + E + F + H + I - J - N"
  elif N == {a,b,d,e,g,i,l,r,o,j}:
    mult = "A - B - C - D + G + H + I - L - O + Q"
  elif N == {a,b,d,e,f,g,h,i,l,o,j}:
    mult = "A - B - C - D + F + G + H + I - K - L"
  elif N == {a,b,c,d,e,f,g,l,j}:
    mult = "A - B - C - D + E + F + G + H - J"
  elif N == {a,b,d,e,f,g,i,l,o,j}:
    mult = "A - B - C - D + F + G + H + I - L"
  elif N == {a,b,d,e,f,g,h,l,o,j}:
    mult = "A - B - C - D + F + G + H + I - K"
  elif N == {a,b,c,d,e,f,l,j}:
    mult = "A - B - C - D + E + F + H - J"
  elif N == {a,b,d,e,f,g,h,l,j}:
    mult = "A - B - C - D + F + G + H - K"
  elif N == {a,b,d,e,f,g,l,o,j}:
    mult = "A - B - C - D + F + G + H + I"
  elif N == {a,b,d,e,g,l,o,j}:
    mult = "A - B - C - D + G + H + I - L"
  elif N == {a,d,e,g,i,l,r,o,j}:
    mult = "A - C - D + G + I - L - O + Q"
  elif N == {a,b,d,e,f,g,l,j}:
    mult = "A - B - C - D + F + G + H"
  elif N == {a,b,d,e,f,l,o,j}:
    mult = "A - B - C - D + F + H + I"
  elif N == {a,b,d,e,g,l,o,j}:
    mult = "A - B - C - D + G + H + I"
  elif N == {a,b,c,d,e,f,j}:
    mult = "A - B - C + E + F - J"
  elif N == {a,b,d,e,f,l,j}:
    mult = "A - B - C - D + F + H"
  elif N == {a,b,d,e,g,l,j}:
    mult = "A - B - C - D + G + H"
  elif N == {a,d,e,g,i,l,o,j}:
    mult = "A - C - D + G + I - L"
  elif N == {a,b,d,e,l,o,j}:
    mult = "A - B - C - D + H + I"
  elif N == {a,b,d,e,l,j}:
    mult = "A - B - C - D + H"
  elif N == {a,d,e,g,l,o,j}:
    mult = "A - C - D + G + I"
  elif N == {a,b,d,e,f,j}:
    mult = "A - B - C + F"
  elif N == {a,b,d,l,j}:
    mult = "A - B - D + H"
  elif N == {a,d,e,g,l,j}:
    mult = "A - C - D + G"
  elif N == {a,d,e,l,o,j}:
    mult = "A - C - D + I"
  elif N == {a,b,d,e,j}:
    mult = "A - B - C"
  elif N == {a,d,e,l,j}:
    mult = "A - C - D"
  elif N == {a,b,d,j}:
    mult = "A - B"
  elif N == {a,d,e,j}:
    mult = "A - C"
  elif N == {a,d,l,j}:
    mult = "A - D"
  elif N == {a,d,j}:
    mult = "A"
  else:
    mult = "0"
  
print("m_q({},{}) = {}".format(Lambda, Mu, mult))
