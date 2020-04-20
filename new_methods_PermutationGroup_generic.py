
def polynomial_of_profile(perm_group):
    # méthode à rajouter à la classe PermutationGroup_generic façon monkey patching
    
    r"""
    Return the generating polynomial of the (finite) profile of the group.
    
    OUTPUT: 

    * A symbolic polynomial in z with nonnegative integer coefficients

    EXAMPLES::
        
        sage: C5 = FinitePoligomorphicGroup(CyclicPermutationGroup(5))
        sage: C5.polynomial_of_profile()
        z^5 + z^4 + 2*z^3 + 2*z^2 + z + 1
    
    """

    if not (isinstance(perm_group, PermutationGroup_generic) and perm_group.is_finite()):
        raise TypeError("Argument must be a finite permutation group.")
    var('z')
    poly = 0
    CCl = libgap.ConjugacyClasses(perm_group)
    for gbar in CCl:
        g = libgap.Representative(gbar)
        W_g = 1
        for k in libgap.CycleLengths(g, list(perm_group.domain())):
            W_g *= (1+z**Integer(k))   # choice on each cycle of g: either take (z**1) or not (z**0)
        poly += W_g * Integer(libgap.Size(gbar))
    return expand(poly / perm_group.cardinality())

def profile_series(gp, variable='z'):
    r"""
    Return the (finite) generating series of the (finite) profile of the group.

    The profile of a permutation group G is the counting function that
    maps each nonnegative integer n onto the number of orbits of the
    action induced by G on the n-subsets of its domain.
    If f is the profile of G, f(n) is thus the number of orbits of
    n-subsets of G.

    INPUT:

    - ``variable`` -- a variable, or variable name as a string (default: `'z'`)

    OUTPUT:

    - A polynomial in ``variable`` with nonnegative integer coefficients.
      By default, a polynomial in z over ZZ.

    EXAMPLES::

        sage: C8 = CyclicPermutationGroup(8)
        sage: C8.profile_series()
        z^8 + z^7 + 4*z^6 + 7*z^5 + 10*z^4 + 7*z^3 + 4*z^2 + z + 1
        sage: D8 = DihedralGroup(8)
        sage: poly_D8 = D8.profile_series()
        sage: poly_D8
        z^8 + z^7 + 4*z^6 + 5*z^5 + 8*z^4 + 5*z^3 + 4*z^2 + z + 1
        sage: poly_D8.parent()
        Univariate Polynomial Ring in z over Rational Field
        sage: D8.profile_series(variable='y')
        y^8 + y^7 + 4*y^6 + 5*y^5 + 8*y^4 + 5*y^3 + 4*y^2 + y + 1
        sage: u = var('u')
        sage: D8.profile_series(u).parent()
        Symbolic Ring

    """
    from sage.rings.integer_ring import ZZ

    if isinstance(variable, str):
        variable = ZZ[variable].gen()
    cycle_poly = gp.cycle_index()
    return cycle_poly.expand(2).subs(x0 = 1, x1 = variable)

profile_polynomial = profile_series

def profile(gp, n, using_polya=True):
    r"""
    Return the value in ``n`` of the profile of the group ``self``.

    Optional argument ``using_polya`` allows to change the default method.

    INPUT:

    - ``n`` -- a nonnegative integer

    - ``using_polya`` (optional) -- a boolean: if ``True`` (default), the computation
      uses P\'{o}ly\`{a} enumeration (and all values of the profile are cached, so this
      should be the method used in case several of them are needed);
      if ``False``, uses the GAP interface to compute the orbit.

    OUTPUT:

    - A nonnegative integer that is the number of orbits of ``n``-subsets
      under the action induced by ``self`` on the subsets of its domain
      (i.e. the value of the profile of ``self`` in ``n``)

    EXAMPLES::

        sage: C6 = CyclicPermutationGroup(6)
        sage: C6.profile(2)
        3
        sage: C6.profile(3)
        4
        sage: D8 = DihedralGroup(8)
        sage: D8.profile(4, using_polya=False)
        8

    """

    if using_polya:
        return gp.profile_polynomial()[n]
    else:
        subs_n = libgap.Combinations(list(gp.domain()), n)
        return len(libgap.Orbits(gp, subs_n, libgap.OnSets))

def age_explicit(perm_group, homogeneous_component='unspecified'):
    if not (isinstance(perm_group, PermutationGroup_generic) and perm_group.is_finite()):
        raise TypeError("Argument must be a finite permutation group.")
    # remove when properly added to the class
    if homogeneous_component == 'unspecified':
        subs = libgap.Combinations(list(perm_group.domain()))
        age = libgap.Orbits(perm_group, subs, libgap.OnSets)
        return age.sage()
    n = homogeneous_component
    if not (isinstance(n, Integer) and n >= 0):
        raise TypeError("Optional argument 2 homogeneous_component must be a positive integer.")
    subs_n = list(libgap.Combinations(list(perm_group.domain()), n))
    # gp_n = libap.Action(perm_group, subs_n, libgap.OnSets)
    age_n = libgap.Orbits(perm_group, subs_n, libgap.OnSets).sage()
    for i in range(len(age_n)):
        age_n[i] = sorted(age_n[i])
    return age_n

PermutationGroup_generic.polynomial_of_profile = polynomial_of_profile
PermutationGroup_generic.profile_series = profile_series
PermutationGroup_generic.profile = profile
PermutationGroup_generic.age_explicit = age_explicit


def age_representatives(perm_group):
    if not (isinstance(perm_group, PermutationGroup_generic) and perm_group.is_finite()):
        raise TypeError("Argument must be a finite permutation group.")
    # remove when properly added to the class
    subs = libgap.Combinations(list(finite_group.domain()))
    finite_age = list(libgap.Orbits(finite_group,subs,libgap.OnSets))
    #finite_age.pop(0)   # removes the orbit of the empty set, which is the unit of the algebra
    age_representatives = []
    for finite_orbit in finite_age:
        repres = finite_orbit[0]
        age_representatives.append(repres)
    return age_representatives

