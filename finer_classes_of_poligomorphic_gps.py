
class DirectProductOfPoligomorphicGroups(PoligomorphicGroup):
    # knowledge of a direct product structure may ease the computation of the profile
    def __init__(self, groups):   
        # groups = list of P-oligomorphic groups and maybe finite permutation groups
        if not (isinstance(groups, list) and all([isinstance(g, PoligomorphicGroup) or isinstance(g, PermutationGroup_generic) for g in groups])):
            raise TypeError("Argument should be a list of P-oligomorphic groups or finite permutation groups.")
        self.groups = groups

    def __repr__(self):
        to_print = "DirectProduct([ "
        first = True
        for gp in self.groups:
            if not first:
                to_print += ", "
            to_print += repr(gp)
            first = False
        to_print += " ])"
        return to_print

    def __str__(self):
        to_print = "Direct product of P-oligomorphic groups: [ "
        first = True
        for gp in self.groups:
            if not first:
                to_print += ", "
            to_print += repr(gp)
            first = False
        to_print += " ]"
        return to_print

    def kernel(self):
        ker = []
        for group in groups:
            ker += group.kernel()
        return ker

    def restriction_to_kernel(self):
        list_gp_ker = []
        for group in groups:
            list_gp_ker.append(group.restriction_to_kernel())
        return direct_product_permgroups(list_gp_ker)
        
    def profile_series(self, variable='z'):
        r"""
        Return the generating series of the profile of the group.
        
        OUTPUT:

        * A symbolic series in z with nonnegative integer coefficients
    
        EXAMPLES::
        
            sage: C5 = CyclicPermutationGroup(5)
            sage: SwrC5 = WreathProductInfiniteBlocks(SymInfinity(), C5)
            sage: S3 = SymmetricGroup(3)
            sage: G = DirectProductOfPoligomorphicGroups([SwrC5, S3, RevQQ()])
            sage: G.profile_series()
            (z^4 - 3*z^3 + 5*z^2 - 3*z + 1)*(z^2 + 1)*(z + 1)/((z^4 + z^3 + z^2 + z + 1)*(z - 1)^6)
        
        """
        from sage.rings.integer_ring import ZZ
        if isinstance(variable, str):
            variable = ZZ[variable].gen()
        z = variable
        series = 1
        for group in self.groups:
            series *= group.profile_series(z)
        return simplify(series)

    # TODO: conversion method to PoligomorphicGroup_generic
        



class PermutingSuperblocks(PoligomorphicGroup):   
    # most general class for now (other than PoligomorphicGroup_generic)
    def __init__(self, single_superblock, finite_group):
        if not isinstance(single_superblock, SingleSuperblock):
            raise TypeError("Argument single_superblock (=%s) must be a P-oligomorphic group with a single superblock."%single_superblock)
        if not (isinstance(finite_group, PermutationGroup_generic) and finite_group.is_finite()):
            raise TypeError("Argument finite_group (=%s) must be a finite permutation group."%finite_group)
        self.group_on_superblocks = finite_group
        self.one_superblock = single_superblock

    def kernel(self):
        return []

    def restriction_to_kernel(self):
        return PermutationGroup([])

    def action_on_superblocks(self):
        return self.group_on_superblocks

    def number_of_superblocks(self):
        return self.action_on_superblocks().degree()

    def restriction_to_superblock(self, block=[]): # block is actually useless here
        return self.one_superblock

    def diagonal_action(self): # TODO
        superblock = self.restriction_to_superblock()
        H0 = superblock._H0()
        dir_prod = PermutationGroup([])
        basis = list(range(self.number_of_superblocks()))
        for i in basis:
            dir_prod = dir_prod.direct_product(H0, False)
        
        return NotImplemented

    def restriction_to_finite_block(self, block=[]): # block is actually useless here
        return self.restriction_to_superblock().restriction_to_finite_fblock()

    def size_of_max_blocks(self):
        return self.restriction_to_finite_block().degree()

    def action_on_finite_blocks(self):  # or not ?
    # could be implemented via WreathProductInfiniteBlocks, the finite blocks being considered singletons
        pass

    def _highly_homogeneous_group(self):  # should be SymInfinity() if non trivial finite blocks
        return self.restriction_to_superblock().action_on_finite_blocks()

    def is_highly_homogeneous(self):
        return self.number_of_superblocks() == 1 and self.restriction_to_superblock().is_highly_homogeneous()

    # TODO: methods for profile and series


class SingleSuperblock(PermutingSuperblocks):
    r"""
    A P-oligomorphic group with only one superblock.

    """

    def __init__(self, restriction_to_block, wreath_base=None, hhomogeneous_group=SymInfinity()):
        if not (isinstance(hhomogeneous_group, HighlyHomogeneousGroup)):
            raise TypeError("Argument hhomogeneous_group (=%s) must be a closed highly homogeneous group."%hhomogeneous_group)
        if not isinstance(restriction_to_block, PermutationGroup_generic):
            raise TypeError("First argument restriction_to_block (=%s) must be a finite permutation group."%restriction_to_block)
        if wreath_base == None:
            wreath_base = restriction_to_block # TODO: check math vocabulary
        if not isinstance(wreath_base, PermutationGroup_generic):  
            raise TypeError("Optional argument wreath_base (=%s) must be a finite permutation group."%wreath_base)
        if restriction_to_block.degree() != wreath_base.degree():
            raise ValueError("Arguments restriction_to_block (=%s) and wreath_base (optional, =%s) must be finite permutation groups of same degree."%(restriction_to_block, wreath_base))
        H = libgap.NormalClosure(restriction_to_block, wreath_base)
        self.H = PermutationGroup(gap_group = H)
        if (restriction_to_block.degree() != 1) and (hhomogeneous_group != SymInfinity()):
            raise ValueError("Action on blocks must be symmetric or the blocks should be trivial for the group to be P-oligomorphic.")
            # TODO: this case could be accepted though, and trigger an error if the user
            # tries to use the conversion method towards a generic P-oligomorphic group,
            # but it would mostly be interesting if the method for computing the profile
            # is consequently adapted/enlarged (this should be done at the level of wreath
            # products with finite blocks: one does not count subsets any more, but rather
            # words with certain symmetries, which hands an exponential profile...).
        H0 = libgap.ClosureGroup(restriction_to_block, wreath_base)
        # closure in case wreath_base is not a subgroup of restriction_to_block
        self.H0 = PermutationGroup(gap_group = H0)
        self.group_on_blocks = hhomogeneous_group

    def __eq__(self, other):
        if not isinstance(other, SingleSuperblock):
            return False
        data1 = (self.group_on_blocks, self.H0, self.H)
        data2 = (other.group_on_blocks(), other.restriction_to_block(), other.wreath_base())
        return data1 == data2

    def action_on_superblocks(self):
        return PermutationGroup_generic([[1]])

    def restriction_to_superblock(self):
        return self

    def action_on_finite_blocks(self):
        return self.group_on_blocks

    def restriction_to_finite_block(self):
        return self.H0

    diagonal_action = restriction_to_finite_block
    _H0 = restriction_to_finite_block

    def wreath_base(self):
        return self.H

    _H1 = wreath_base

    @cached_method
    def _synchro(self):
        return self.restriction_to_finite_block().quotient(self.wreath_base())

    def size_of_max_blocks(self):
        return self.wreath_base().degree()

    def is_highly_homogeneous(self):
        return self.wreath_base().degree() == 1

    def is_wreath_product(self):
        return self.wreath_base() == self.restriction_to_finite_block()

    def underlying_wreath_product(self):
        return WreathProductFiniteBlocks(self.wreath_base(), self.action_on_finite_blocks())

    def number_of_superblocks(self):
        return 1

    def has_independent_superblocks(self):
        return True
    
    def __repr__(self):
        if self.is_highly_homogeneous():
            return repr(self.action_on_finite_blocks())
        elif self.is_wreath_product():
            return repr(self._H0()) + "  wreath  " + repr(self.action_on_finite_blocks())
        else:
            return repr(self._H1()) + "  wreath  " + repr(self.action_on_finite_blocks()) + "  with diagonal action of  " + repr(self._H0())

    def __str__(self):
        r"""
        Return a user readable description of the object.

        """
        if self.is_highly_homogeneous():
            return repr(self.action_on_finite_blocks())
        elif self.is_wreath_product():
            return "The P-oligomorphic wreath product  " + repr(self._H0()) + "  wreath  " + repr(self.action_on_finite_blocks())
        else:
            return "The P-oligomorphic group  " + repr(self._H1()) + "  wreath  " + repr(self.action_on_finite_blocks()) + "  with diagonal action of  " + repr(self._H0()) + "  within blocks altogether  "

    def print_tower(self):
        # TODO: find a way to access libgap.StructureDescription (to be called on gap(self._Hi()))
        str_H0 = "H0 "
        str_H1 = "H  "
        print(str_H0, str_H1, str_H1, str_H1, str_H1, str_H1, " ...\nwith H0 =", self._H0(), "\n     H =", self._H1())

    def profile_series(self, variable='z'):
        from sage.rings.integer_ring import ZZ

        if isinstance(variable, str):       # creation or recuperation of the variable
            variable = ZZ[variable].gen()
        z = variable
        M = self.size_of_max_blocks()
        dom = [None]
        homom = [None]
        for d in range(1, M+1):      # will enable to consider the actions on orbits of degree d
            age_d = self.wreath_base().age_explicit(Integer(d))
            dom.append(list(range(1, len(age_d)+1)))
            homom.append(libgap.ActionHomomorphism(self._H0(), age_d, libgap.OnSetsSets))
        CCl = libgap.ConjugacyClasses(self._H0())
        polya_sum = 0
        for gbar in CCl:
            g = libgap.Representative(gbar)          # computation of the weight of g
            W_g = 1
            for d in range(1, M+1):           # run through the orbital degrees of the orbits of the H1's
                g_d = libgap.Image(homom[d], g)      # g acting on the orbits of degree d
                CT_d = libgap.CycleLengths(g_d, dom[d])
                for k in CT_d:
                    W_g /= (1 - z**(d*Integer(k)))
            polya_sum += W_g*len(libgap.List(gbar))
        result = simplify(polya_sum / self._H0().order())
        return result


    def as_PoligomorphicGroup(self):
        H0 = self._H0()
        repres_of_superblock = H0.domain()
        return PoligomorphicGroup_generic(H0, block_system=[repres_of_superblock], wreath_bases=[[repres_of_superblock, self._H1()]], hhomogeneous_groups=[[repres_of_superblock, self.action_on_finite_blocks()]])


class WreathProductFiniteBlocks(SingleSuperblock):
    r"""
    An (unrestricted) wreath product as a P-oligomorphic group with finite maximal blocks.

    """

    def __init__(self, finite_group, action_on_blocks=SymInfinity()):
        if not isinstance(action_on_blocks, HighlyHomogeneousGroup):
            raise TypeError("Second (optional) argument must be a closed highly homogeneous group.")
        if not (isinstance(finite_group, PermutationGroup_generic) and finite_group.is_finite()):
            raise TypeError("First argument must be a finite permutation group.")
        if not (finite_group.degree() == 1 or action_on_blocks == SymInfinity()):
            raise ValueError("Action on blocks must be symmetric or the blocks should be trivial for the group to be P-oligomorphic.")
        self.group_on_blocks = action_on_blocks
        self.H = finite_group
        # doit y avoir moyen de reutiliser l'init de SingleSuperblock

    def action_on_finite_blocks(self):  # does the same as the inherited method action_on_blocks
        return self.group_on_blocks

    def restriction_to_finite_block(self):
        return self.H

    diagonal_action = restriction_to_finite_block
    _H0 = restriction_to_finite_block
    _H1 = restriction_to_finite_block

    def __repr__(self):
        if self.is_highly_homogeneous():
            return repr(self.action_on_finite_blocks())
        return repr(self._H0()) + "  wreath  " + repr(self.action_on_finite_blocks())

    def as_PoligomorphicGroup(self):
        H0 = self._H0()
        repres_of_superblock = list(H0.domain())
        return PoligomorphicGroup_generic(H0, block_system=[repres_of_superblock], hhomogeneous_groups=[[repres_of_superblock, self.action_on_finite_blocks()]])


class WreathProductInfiniteBlocks(PermutingSuperblocks):
    r"""
    An (unrestricted) wreath product as a P-oligomorphic group with infinite blocks.

    """
    def __init__(self, finite_group, hhomogeneous_group=SymInfinity()):
        if not (isinstance(finite_group, PermutationGroup_generic) and finite_group.is_finite()):
            raise TypeError("First argument must be a finite permutation group.")
        if not isinstance(hhomogeneous_group, HighlyHomogeneousGroup):
            raise TypeError("Second argument must be a closed highly homogeneous group.")
        self.group_on_blocks = finite_group
        self.hhomogeneous_group = hhomogeneous_group

    def action_on_infinite_blocks(self):
        return self.group_on_blocks

    diagonal_action = action_on_infinite_blocks
    action_on_superblocks = action_on_infinite_blocks

    def restriction_to_infinite_block(self):
        return self.hhomogeneous_group

    def restriction_to_superblock(self):
        return self.hhomogeneous_group

    def is_highly_homogeneous(self):
        return self.action_on_finite_blocks().degree() == 1

    def __repr__(self):
        if self.is_highly_homogeneous():
            return repr(self.restriction_to_infinite_block())
        return repr(self.restriction_to_infinite_block()) + "  wreath  " + repr(self.action_on_infinite_blocks())

    def size_of_max_blocks(self):
        return Infinity

    def number_of_superblocks(self):  
        return self.action_on_superblocks().degree()
        
    def profile_series(self, variable='z'):
        r"""
        Return the generating series of the profile of the group.

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

            sage: C5 = CyclicPermutationGroup(5)
            sage: SwrC5 = WreathProductInfiniteBlocks(SymInfinity(), C5)
            sage: SwrC5.profile_series()
            -(z^4 - 3*z^3 + 5*z^2 - 3*z + 1)/((z^4 + z^3 + z^2 + z + 1)*(z - 1)^5)
        
        """
        from sage.rings.integer_ring import ZZ

        if isinstance(variable, str):
            variable = ZZ[variable].gen()
        z = variable
        CCl = libgap.ConjugacyClasses(self.group_on_blocks)
        polya_sum = 0
        for gbar in CCl:
            g = libgap.Representative(gbar)
            W_g = 1
            for k in libgap.CycleLengths(g, list(self.group_on_blocks.domain())):
                W_g /= (1 - z**Integer(k))   # choice of image of "h" on each cycle of g (how many elements in "this" block)
            polya_sum += W_g * Integer(libgap.Size(gbar))
        return simplify(polya_sum / self.group_on_blocks.cardinality())

    def as_PoligomorphicGroup(self):
        repres_of_blocks = []  # trivial block system
        hhomogeneous_extensions = []
        not_sym = self.hhomogeneous_group != SymInfinity()
        for a in self.group_on_blocks.domain():
            repres_of_blocks.append([a])
            if not_sym:
                hhomogeneous_extensions.append([[a], self.hhomogeneous_group])
        return PoligomorphicGroup_generic(self.group_on_blocks, block_system=repres_of_blocks, hhomogeneous_groups=hhomogeneous_extensions)
    


def restriction(gp, subset):
    stab = libgap.Stabilizer(gp, subset, libgap.OnSets)
    return libgap.Action(stab, subset)

def blocks_intransitive(gp):
    if gp.is_transitive():
        return list(libgap.Blocks(gp, gp.domain()))
    block_syst = []
    orbits = gp.orbits()
    for orbit in orbits:
        restr = restriction(gp, orbit)
        block_syst += list(libgap.Blocks(restr, orbit))
    return block_syst

