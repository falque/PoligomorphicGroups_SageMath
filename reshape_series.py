
###### Auxiliary, conveniency functions ######

def identify_cyclotomic_index(poly, higher_bound) :
    z = poly.parent().gen()
    for n in range(1, higher_bound+1) :
        if poly == cyclotomic_polynomial(n, z) :
            return n
    raise ValueError("Argument (=%s) should be a cyclotomic polynomial"%poly)

def _multiplicities_of_factors(denom, bound_on_degree) :
    # * Input:  -``denom``: a product of cyclotomic polynomials
    # * Output: a list describing a "factorization" (multiplicities in ZZ)
    #           involving factors of the form (1 - z^i)
    # This uses the Moebius inversion formula, applied to cyclotomic polynomials
    factors = list(factor(denom))
    indices = []
    for f in factors :
        n = identify_cyclotomic_index(f[0], bound_on_degree)
        indices.append(n)
    mult = []  # creation and initialization to 0 of the result list
    for i in range(max(indices)+1) : 
        # max(indices) will be the highest degree in the new factorization
        mult.append(0)
    for i in range(len(factors)) :
        n = indices[i]  # the factor is the n-th cyclotomic poly
        for d in divisors(n) :
            mult[d] += moebius(n/d)*factors[i][1]  # values in ZZ
            # factors[i][1] is the multiplicity of the factor
    return mult

def _pos_of_first_non_zero(L) :
    for i in range(len(L)) :
        if L[i] <> 0 :
            return i
    return -1

def _build_factorization(list_of_mult, variable = 'z') :
    # build a polynomial factorization of (1-z^i) with multiplicities 
    # specified by the list
    from sage.rings.integer_ring import ZZ
    if isinstance(variable, str):
        R = PolynomialRing(ZZ, Integer(1), order = 'neglex', names=(variable,))
        (variable,) = R._first_ngens(1)
    z = variable
    facto = []
    for i in range(1, len(list_of_mult)) :
        facto.append(((1-z**i), list_of_mult[i]))
    return Factorization(facto)

def _facto_to_product(facto) :
    # build a polynomial product of (1-z^i) according to the factorization
    res = 1
    for f in facto :
        res *= f[0]**f[1]
    return res

##### MAIN functions #####

def test_denom(mult, series) : 
    # Test if a denominator that is a product of (1-z^i) with the 
    # specified multiplicities makes all coefficients of the numerator >= 0
    # mult should define (via _build_product) a denominator that is a 
    # multiple of the current one in series
    num = series.numerator()
    z = series.parent().gen()
    new_denom = _facto_to_product(_build_factorization(mult, z))
    new_num = ZZ[z](num * (new_denom/series.denominator()).numerator()) # should be a polynomial
    # .numerator() just ensures we get a polynomial type, ZZ[z] makes it considered univariate
    # Are the coefficients all positive ? (or negative)
    all_positive = True
    all_negative = True
    for coeff in new_num.coefficients(sparse=False) :
        if coeff < 0 :
            all_positive = False
        elif coeff > 0 :
            all_negative = False
    return all_positive or all_negative

def find_denom(denom_as_list, pos, m, size, series) :
    # Recursive search (backtracking type)
    # - size : should be len(denom_as_list) (to avoid systematic recomputing)
    # Test all denominators with factor in position pos pushed at least 
    # to pos*m and preceding factors untouched
    # Research on list, but heuristically changes one (1-z^i) in the
    # denominator defined by the list into (1-z^(i*m)) (by multiplying
    # numerator and denominator by the appropriate cyclotomic polynomials)
    # and tests if it makes the numerator positive; 
    # if not, calls or returns
    new_pos = pos*m
    # safety test
    if new_pos > size :
        return []
    # make the safe move
    new_denom = deepcopy(denom_as_list)
    new_denom[pos] -= 1
    new_denom[new_pos] += 1
    # test of success
    if test_denom(new_denom, series) :
        return new_denom
    # recursive calls
    # 1) try to push the factor further, recursion on m
    return_value = find_denom(denom_as_list, pos, m+1, length, series)
    if return_value <> [] :
        return return_value
    next_pos = pos_first_non_zero(new_denom)
    # last hope; if fail, backtrack by returning to the calling function
    return find_denom(new_denom, next_pos, 2, length, series)

def nice_factorization(series, bound_on_degree, print_facto=False):
    # series must not be symbolic
    denom = series.denominator()
    mult = _multiplicities_of_factors(denom, bound_on_degree)
    size = len(mult)
    first_pos = _pos_of_first_non_zero(mult)
    new_mult = find_denom(mult, first_pos, 1, size, series)
    if new_mult <> [] :
        z = series.parent().gen()
        facto = _build_factorization(new_mult, z)
        new_denom = _facto_to_product(facto)
        correction_factor = new_denom/denom
        new_num = expand(series.numerator()*correction_factor)
        if print_facto:
            print repr(new_num) + " / " + repr(facto)
        new_num = Factorization([(new_num, 1)])
        return new_num.__mul__(facto.__invert__())
    else :
        raise Error("Argument series (=%s) could not be reshaped"%series)

#############################


# Now useless (been added as a method in PoligomorphicGroup)
def nice_factorization_of_group_series(group, variable='z', print_facto=False) :
    bound_on_degree = group.finite_group.degree()
    # bound_on_degree : higher bound for the degree of generators in the
    #                   orbit algebra (will be the degree of the finite
    #                   diagonal action)
    from sage.rings.integer_ring import ZZ
    if isinstance(variable, str):
        R = PolynomialRing(ZZ, Integer(1), order = 'neglex', names=(variable,))
        (variable,) = R._first_ngens(1)
        z = variable
        # in theory equivalent to R.<z> = PolynomialRing(ZZ, 1, order = 'neglex')
        # but preparsed so it is accepted in a function
        # (something like z = PolynomialRing(ZZ, order='neglex') just won't work
        # in the sense that the printing of the series will not follow the specified order)
    else:
        R = PolynomialRing(ZZ, Integer(1), order = 'neglex', names=('z',))
        (z,) = R._first_ngens(1)
    # un peu brutal mais bon... on a besoin que la serie ne soit pas une expression symbolique
    # sans quoi factor ne renvoie pas une factorisation et c'est un pb
    series = group.profile_series(z)
    return nice_factorization(series, bound_on_degree, print_facto)


