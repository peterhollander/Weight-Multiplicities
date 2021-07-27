# (c) Pamela E. Harris, Peter Hollander, Maria Rodriguez-Hertz, Daniel Qin
# A program which calculates a Kostant's Weight q-multiplicity formula for a given lambda-mu pair
# To run, type python3 q-mult-formula.py
# You will be prompted to enter lambda, then mu -- do so in the format requested
# The output will be a formula using capital letters A through Q

if __name__ == "__main__":
  
  lam = input("Enter Lambda (format: (#,#,#)): ")[1:-1].split(",") #remove parens, then break into components
  mu = input("Enter Mu (format: (#,#,#)): ")[1:-1].split(",") #remove parens, then break into components
  
  m = int(lam[0])
  n = int(lam[1])
  k = int(lam[2])
  
  x = int(mu[0])
  y = int(mu[1])
  z = int(mu[2])
  
  a = m + n + k - x - y - z
  b = n + k - x - y - z - 1
  c = k - x - y - z - 2
  d = m + 2*n + 2*k - x - 2*y - 2*z
  e = m + n + 2*k - x - 2*y - 2*z - 1
  f = n + 2*k - x - 2*y - 2*z - 2
  g = m + n - x - 2*y - 2*z - 3
  h = n - x - 2*y - 2*z - 4
  i = m - x - 2*y - 2*z - 4
  j = 0.5*m + n + 1.5*k- 0.5*x - y - 1.5*z
  l = 0.5*m + n + 0.5*k - 0.5*x - y - 1.5*z - 1
  o = 0.5*m + 0.5*k - 0.5*x - y - 1.5*z - 2
  p = -0.5*m + 0.5*k - 0.5*x - y - 1.5*z - 3
  r = 0.5*m - 0.5*k - 0.5*x - y - 1.5*z - 3
  
  vals = [a,b,c,d,e,f,g,h,i,j,l,o,p,r]
  N = ""
  #create a binary string representing which variables are nonnegative integers (1) and which are not (0)
  for val in vals:
    if (int(val) == val) and (val >= 0):
      N += "1"
    else:
      N += "0"
    
  #based on the string N, print a formula
  mult = ""
  if N == "11111111111101": #1
    mult = "A - B - C - D + E + F + G + H + I - J - K - L - N - O + Q"
  elif N == "11111111111110": #2
    mult = "A - B - C - D + E + F + G + H + I - J - K - L - M - N - P"
  elif N == "11111110111101": #3
    mult = "A - B - C - D + E + F + G + H + I - J - L - N - O + Q"
  elif N == "11111110111110": #4
    mult = "A - B - C - D + E + F + G + H + I - J - L - M - N + P"
  elif N == "11111111011110": #5
    mult = "A - B - C - D + E + F + G + H + I - J - K - M - N + P"
  elif N == "11111110011110": #6
    mult = "A - B - C - D + E + F + G + H + I - J - M - N + P"
  elif N == "11111111111100": #7
    mult = "A - B - C - D + E + F + G + H + I - J - K - L - N"
  elif N == "11111110111100": #8
    mult = "A - B - C - D + E + F + G + H + I - J - L - N"
  elif N == "11111100011110": #9
    mult = "A - B - C - D + E + F + H + I - J - M - N + P"
  elif N == "11111111011100": #10
    mult = "A - B - C - D + E + F + G + H + I - J - K - N"
  elif N == "11011111111101": #11
    mult = "A - B - C - D + F + G + H + I - K - L - O + Q"
  elif N == "11111110011100": #12
    mult = "A - B - C - D + E + F + G + H + I - J - N"
  elif N == "11011110111101": #13
    mult = "A - B - C - D + F + G + H + I - L - O + Q"
  elif N == "11111111011000": #14
    mult = "A - B - C - D + E + F + G + H - J - K"
  elif N == "11111100011100": #15
    mult = "A - B - C - D + E + F + H + I - J - N"
  elif N == "11011010111101": #16
    mult = "A - B - C - D + G + H + I - L - O + Q"
  elif N == "11011111111100": #17
    mult = "A - B - C - D + F + G + H + I - K - L"
  elif N == "11111110011000": #18
    mult = "A - B - C - D + E + F + G + H - J"
  elif N == "11011110111100": #19
    mult = "A - B - C - D + F + G + H + I - L"
  elif N == "11011111011100": #20
    mult = "A - B - C - D + F + G + H + I - K"
  elif N == "11111100011000": #21
    mult = "A - B - C - D + E + F + H - J"
  elif N == "11011111011000": #22
    mult = "A - B - C - D + F + G + H - K"
  elif N == "11011110011100": #23
    mult = "A - B - C - D + F + G + H + I"
  elif N == "11011010011100": #24
    mult = "A - B - C - D + G + H + I - L"
  elif N == "10011010111101": #25
    mult = "A - C - D + G + I - L - O + Q"
  elif N == "11011110011000": #26
    mult = "A - B - C - D + F + G + H"
  elif N == "11011100011100": #27
    mult = "A - B - C - D + F + H + I"
  elif N == "11011010011100": #28
    mult = "A - B - C - D + G + H + I"
  elif N == "11111100010000": #29
    mult = "A - B - C + E + F - J"
  elif N == "11011100011000": #30
    mult = "A - B - C - D + F + H"
  elif N == "11011010011000": #31
    mult = "A - B - C - D + G + H"
  elif N == "10011010111100": #32
    mult = "A - C - D + G + I - L"
  elif N == "11011000011100": #33
    mult = "A - B - C - D + H + I"
  elif N == "11011000011000": #34
    mult = "A - B - C - D + H"
  elif N == "10011010011100": #35
    mult = "A - C - D + G + I"
  elif N == "11011100010000": #36
    mult = "A - B - C + F"
  elif N == "11010000011000": #37
    mult = "A - B - D + H"
  elif N == "10011010011000": #38
    mult = "A - C - D + G"
  elif N == "10011000011100": #39
    mult = "A - C - D + I"
  elif N == "11011000010000": #40
    mult = "A - B - C"
  elif N == "10011000011000": #41
    mult = "A - C - D"
  elif N == "11010000010000": #42
    mult = "A - B"
  elif N == "10011000010000": #43
    mult = "A - C"
  elif N == "10010000011000": #44
    mult = "A - D"
  elif N == "10010000010000": #45
    mult = "A"
  else: #46
    mult = "0"
  
print("m_q({},{}) = {}".format(lam, mu, mult))
