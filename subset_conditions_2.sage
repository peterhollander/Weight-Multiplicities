W = WeylGroup(['C', 3], prefix='s')
[s1, s2, s3] = W.simple_reflections()
e=s1*s1

var('P Q R')
p = var('P1 P2 P3 P4')
q = var('Q1 Q2 Q3 Q4 Q5 Q6')
r = var('R1 R2 R3 R4')

[P1, P2, P3, P4] = p
[Q1, Q2, Q3, Q4, Q5, Q6] = q
[R1, R2, R3, R4] = r

pqr = p + q + r

alts = load('alts_mu_positive')
weight_mult_alts = load('weight_mult_alts')

positive_conds= {
  e:{ P:1, Q:1, R:1 },
  s1*s2*s1:{ P:3, Q:4, R:1 },
  s3:{ P:1, Q:1, R:4 },
  s1*s2:{ P:3, Q:6, R:1 },
  s2*s3:{ P:1, Q:5, R:4 },
  s1:{ P:4, Q:1, R:1 },
  s2*s3*s2:{ P:1, Q:5, R:3 },
  s2*s1:{ P:4, Q:4, R:1 },
  s3*s1:{ P:4, Q:1, R:4 },
  s3*s2:{ P:1, Q:6, R:3 },
  s2:{ P:1, Q:6, R:1 },
  s3*s1*s2:{ P:3, Q:6, R:3 }
}

all_conds = {
  e:{ P:1, Q:1, R:1 },
  s1*s2*s3*s2*s1:{ P:2, Q:3, R:2 },
  s1*s2*s3*s1:{ P:2, Q:3, R:4 },
  s1*s2*s3*s2:{ P:2, Q:5, R:3 },
  s1*s2*s1:{ P:3, Q:4, R:1 },
  s3*s1*s2*s1:{ P:3, Q:4, R:2 },
  s3:{ P:1, Q:1, R:4 },
  s3*s1*s2:{ P:3, Q:6, R:3 },
  s1*s2:{ P:3, Q:6, R:1 },
  s2*s3*s1*s2*s1:{ P:3, Q:2, R:2 },
  s2*s3:{ P:1, Q:5, R:4 },
  s2*s3*s1*s2:{ P:3, Q:2, R:3 },
  s1:{ P:4, Q:1, R:1 },
  s2*s3*s2*s1:{ P:4, Q:3, R:2 },
  s2*s3*s1:{ P:4, Q:3, R:4 },
  s2*s3*s2:{ P:1, Q:5, R:3 },
  s2*s1:{ P:4, Q:4, R:1 },
  s3*s2*s1:{ P:4, Q:4, R:2 },
  s3*s1:{ P:4, Q:1, R:4 },
  s3*s2:{ P:1, Q:6, R:3 },
  s2:{ P:1, Q:6, R:1 },
  s1*s2*s3*s1*s2*s1:{ P:2, Q:2, R:2 },
  s1*s2*s3:{ P:2, Q:5, R:4 },
  s1*s2*s3*s1*s2:{ P:2, Q:2, R:3 }
}

def find_conditions(alt, conds):
    (true_, false_) = get_general_formula(alt, conds)

    # First pass: Remove redundancies with nonnegative assumptinos
    for var in true_:
        for s in false_:
            if var in s:
                s.remove(var)

    # Get max length of false set
    max_length = max([len(s) for s in false_])

    '''
    Remove duplicates until 'removing duplicates' does nothing,
    i.e. until there are no more duplicates
    '''
    prev_false = []
    while(false_ != prev_false):
        prev_false = false_[:]
        for size in range(1,max_length+1):
            for f in false_:
                if len(f) == size:
                    for suspect in false_:
                        if(f != suspect and f.issubset(suspect)):
                            false_.remove(suspect)
                            #false_[j] = false_[j].difference(f)

    return (true_,remove_duplicates(false_))


def get_general_formula(alt, conds):
    true_ = set()
    false_ = []
    for sigma in conds.keys():
        term = {p[conds[sigma][P]-1] , q[conds[sigma][Q]-1] , r[conds[sigma][R]-1]}
        if(sigma in alt):
            for t in term:
                true_.add(t)
        else:
            false_.append(term)
    return (true_,false_)


def remove_duplicates(L):
    result = []
    for l in L:
        if l not in result:
            result.append(l)
    return result

# Get the big formula for everything

def get_formula(alt):
    f = propcalc.formula('x')
    for sigma in conds.keys():
        sigma_conds = ''
        if(sigma not in alt):
            sigma_conds += '~'
        sigma_conds += '(p%d&q%d&r%d)' % (conds[sigma][P], conds[sigma][Q], conds[sigma][R])
        sigma_statement = propcalc.formula(sigma_conds)
        f = f.add_statement(sigma_statement,'&')
    return f


alt_conds = dict((alt, find_conditions(alt, all_conds)) for alt in alts)
alt_conds_pos = dict((alt, find_conditions(alt, positive_conds)) for alt in weight_mult_alts)

def pretty_print_conds(conds_dict, sublist, num):
    subset = frozenset(sublist)
    (pos, neg) = conds_dict[subset]
    latex = ''

    # String for Weyl group subset
    subset_str = '\\item \\altcirc{%s} $\A(\lam,\mu)=\{%s\}$' % (str(num), (''.join(str(s)+',' for s in sublist))[:-1])

    # String for nonnegative conditions
    #pos_str = '%s\\in\\N' % (''.join(str(cond)+',' for cond in pos))[:-1]
    pos_str = '%s' % (''.join("$%s$, " % str(cond) for cond in pos))[:-1]
    # Remove last ", "
    #pos_str = pos_str[:-2]

    latex += '%s if %s' % (subset_str, pos_str)

    neg_singletons = [l for l in neg if len(l) == 1]
    neg_multiples = [l for l in neg if len(l) != 1]

    for s in neg_singletons:
        for cond in s:
            latex += "$\\neg %s$, " % str(cond)
    for s in neg_multiples:
        for cond in s:
            latex += "$\\neg %s$ $\\vee$ " % str(cond)
        # Get rid of last " $\vee$ "
        # Add comma
        latex = latex[:-8]
        latex += ", "
    #Get rid of last ", "
    #latex = latex[:-2]

    # String for single negative conditions
    #neg_ones = set.union(*[s for s in neg if len(s) == 1])
    #neg_others = [s for s in neg if len(s) > 1]

    #latex += '\\\\'
    latex = latex.replace('P1', 'K_1')
    latex = latex.replace('P2', 'K_2')
    latex = latex.replace('P3', 'K_3')
    latex = latex.replace('P4', 'K_4')
    latex = latex.replace('Q1', 'K_5')
    latex = latex.replace('Q2', 'K_6')
    latex = latex.replace('Q3', 'K_7')
    latex = latex.replace('Q4', 'K_8')
    latex = latex.replace('Q5', 'K_9')
    latex = latex.replace('Q6', 'K_{10}')
    latex = latex.replace('R1', 'K_{11}')
    latex = latex.replace('R2', 'K_{12}')
    latex = latex.replace('R3', 'K_{13}')
    latex = latex.replace('R4', 'K_{14}')
    latex = latex.replace('s1', 's_1')
    latex = latex.replace('s2', 's_2')
    latex = latex.replace('s3', 's_3')
    latex = latex.replace('*', '')

    #latex = latex.replace('P', 'P_')
    #latex = latex.replace('Q', 'Q_')
    #latex = latex.replace('R', 'R_')
    #latex = latex.replace('s', 's_')
    #latex = latex.replace('*', '')
    #latex = '%s & %s \\text{and}\\ %s\\ \\text{and}\\  %s' % (subset_str,
    #pos_str, neg_ones_str, neg_others_str)
    return latex

def all_conds_latex(fname):
    f = open(fname, 'w+')
    alt_conds_sorted = load('alts_sorted.sobj')
    j = 1
    for sublist in alt_conds_sorted:
        subset = frozenset(sublist)
        if(subset != frozenset()):
            f.write('%s\r\n' % pretty_print_conds(alt_conds, sublist, j))
            j += 1
    f.close()
