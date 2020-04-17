

################################### the generic class ########################################


class PoligomorphicGroup_generic(PoligomorphicGroup):
    def is_block(block, system):
    # TODO: find something more appropriate
        for b in system:
            if Set(b) == Set(block):
                return True
        return False

    def __init__(self, finite_group, block_system=None, block_seed=None, wreath_bases=[], hhomogeneous_groups=[], decorated_block_system=None):
    # wreath_bases = liste de 2-listes, 1er elt = le block dont l'orbite est concernee
    # et 2e elt = le sg normal (ou juste sg dont on prend la cloture normale)
    # meme fonctionnement pour les hhomogeneous_groups
    # par defaut tout est wreath avec sg_\infty (H0=H1)
    # decorated_block_system = liste de 3-listes, 2e moyen de creer une instance (NotImplemented)
    # attribut = decorations (dictionnaire dont les valeurs sont des listes, et les clés des blocs (tuples))
    # ou bien 3 dictionnaires attributs ? (peut-être un peu plus lourd pour la machine mais plus de sémantique)
        if not (isinstance(finite_group, PermutationGroup_generic) and finite_group.is_finite()):
            raise TypeError("Argument finite_group (=%s) must be a finite permutation group."%finite_group)
        # First attribute
        self.finite_group = finite_group
        if decorated_block_system == None:
            if block_system == None:
                if block_seed != None:
                    if not finite_group.is_transitive(): 
                    # TODO: or use homemade function (see blocks_intransitive)
                        raise ValueError("The computation of a block system via argument block_seed (=%s) requires argument finite_group (=%s) to be transitive. Use argument block_system for non transitive groups."%(block_seed, finite_group))
                    block_system = list(libgap.Blocks(finite_group, finite_group.domain(), list(block_seed)))
                else:                                           # default
                    if not (finite_group.is_transitive()):
                        block_system = finite_group.orbits()
                    else:
                        dom = list(range(1,finite_group.degree()+1))
                        block_system = list(libgap.Blocks(finite_group, dom))
            # check argument hhomogeneous_groups
            # check type
            if not all([isinstance(hhg, list) 
                            and sorted(hhg[0]) in block_system 
                            and (isinstance(hhg[1], HighlyHomogeneousGroup)
                                 or hhg[1] == PermutationGroup([])) for hhg in hhomogeneous_groups]):
                raise TypeError("Argument hhomogeneous_groups (=%s) must be a list of lists of first element a block of finite_group and second element a highly homogeneous group (or the trivial group)."%hhomogeneous_groups)
            # check value
            if not all([len(hhg[0]) == 1 
                        or hhg[1] == SymInfinity() 
                        or hhg[1] == PermutationGroup([]) for hhg in hhomogeneous_groups]):
                raise ValueError("Non trivial blocks must be endowed with the highly homogeneous group SymInfinity() (or the trivial group).")
            # creation of a dictionary block -> index of orbit, in order to test the validity 
            # of arguments (values, later)
            self.orbit_indices = {}  # TODO: keep as attribute or clear afterwards?
            orbits_of_blocks = libgap.Orbits(finite_group, block_system, libgap.OnSets)
            i = 0
            for orbit_of_b in orbits_of_blocks:
                for block in orbit_of_b:
                    block_as_tuple = tuple(sorted(block.sage()))
                    self.orbit_indices[block_as_tuple] = i
                i += 1
            # check that the choices (values) in wreath_bases and then in hhomogeneous_groups are compatible
            # ie no different choices for blocks in a same orbit
            l = len(wreath_bases)
            for i in range(l-1):
                for j in range(i+1, l):
                # i,j run through pairs of blocks mentionned in argument wreath_bases
                    block1 = tuple(sorted(wreath_bases[i][0]))
                    block2 = tuple(sorted(wreath_bases[j][0]))
                    wb1 = wreath_bases[i][1]
                    wb2 = wreath_bases[j][1]
                    # rem: if some blocks of a given orbit are specified while others are not
                    # the specification will win over the default choice (see fourth attribute below)
                    if self.orbit_indices[block1] == self.orbit_indices[block2] and wb1 != wb2:
                        raise ValueError("Choice must be the same for blocks in the same orbit in argument list wreath_base"%(wreath_bases[i][0], wreath_bases[j][0]))
            # second compatibility check (in hhomogeneous_groups), plus identification of kernel
            l = len(hhomogeneous_groups)
            kernel = Set([])
            for i in range(l):
            # even if there can be no pair for the last one, we need to get its hhg for the kernel
                block1 = tuple(sorted(hhomogeneous_groups[i][0]))
                hhg1 = hhomogeneous_groups[i][1]
                for j in range(i+1, l):
                # i,j run through pairs of blocks mentionned in argument hhomogeneous_groups
                    block2 = tuple(sorted(hhomogeneous_groups[j][0]))
                    hhg2 = hhomogeneous_groups[j][1]
                    if self.orbit_indices[block1] == self.orbit_indices[block2] and hhg1 != hhg2:
                        raise ValueError("Choice must be the same for blocks in the same orbit (%s and %s) in argument list hhomogeneous_groups."%(hhomogeneous_groups[i][0], hhomogeneous_groups[j][0]))
                if hhg1 == PermutationGroup([]):
                    kernel = kernel.union(Set(block1))
            keys = list(self.orbit_indices.keys())
            # remove the blocks that are part of the kernel
            i = 0
            while i < len(keys):
                if Set(keys[i]).issubset(kernel):
                    self.orbit_indices.pop(keys[i]) # no need any more
                    keys.pop(i)
                    i -= 1
                i += 1
            # Second attribute
            self.ker = list(kernel)
            # Third attribute
            self.repres_from_superblocks = keys
            # TODO: réintégrer les autres tests ici (temporairement stockés l'autre fichier)
            # dictionary of "decorations" used to convert finite blocks of finite_group into superblocks of self
            # (un peu lourd mais pas catastrophique pour ce qu'on veut en faire)
            # par contre, obligation de convertir les listes en tuples pour pouvoir servir de cles...
            # Fourth attribute
            decorations = {}
            # 1st decoration: the restriction H0 to whichever finite block of the superblock
            #                 (not part of the classification since the information is
            #                  technically already included in finite_group; worth keeping ?)
            for block in self.repres_from_superblocks: 
                decorations[tuple(block)] = []
                block_stab = libgap.Stabilizer(self.finite_group, block, libgap.OnSets)
                H0 = libgap.Action(block_stab, block)
                H0 = PermutationGroup(gap_group = H0)
                decorations[tuple(block)].append(H0)
            # 2nd decoration: the wreath base H1 (restriction to whichever f.b. of the sb. after fixing another one)
            done = []       # already managed indices
            for wb in wreath_bases:                     # if H1 specified
                block = tuple(sorted(wb[0]))
                if block in self.repres_from_superblocks: # the kernel is already handled
                    index = self.orbit_indices[block]
                    if not index in done:
                        H0 = decorations[block][0] # already added in decorations (above)
                        H1 = wb[1]
                        if not libgap.IsSubgroup(H0, H1):
                            raise ValueError("The wreath base H1 (=%s) must be a subgroup of the restriction to the block B (=%s), which is H0 = %s"%(H1, wb[0], H0))
                        H1 = libgap.NormalClosure(H0, H1)   # H1 ajusted to normal closure of wb[1] in H0
                        H1 = PermutationGroup(gap_group = H1)
                        for other_block in self.orbit_indices:
                        # if the wreath base is not specified for all blocks in the orbit, 
                        # automatically uses the specified one(s) for the others   
                            if self.orbit_indices[other_block] == index:
                                decorations[other_block].append(H1)
                        done.append(index)
            for block in self.repres_from_superblocks:
            # default value for H1, if not specified at all (then H1 = H0)
                deco = decorations[tuple(block)]       # list (equal pointers)
                if len(deco) < 2:
                    deco.append(deco[0])
            # 3rd decoration: the highly homogeneous group permuting the finite blocks of the superblock
            done = []       # already managed indices
            for hhg in hhomogeneous_groups:       
            # if some hhg is specified, it is passed to all blocks in the same orbit
                block = tuple(sorted(hhg[0]))
                if block in self.repres_from_superblocks: # the kernel is already handled
                    index = self.orbit_indices[block]
                    if not index in done:
                        for other_block in self.orbit_indices:
                            if self.orbit_indices[other_block] == index:
                                decorations[other_block].append(hhg[1])
            for block in self.orbit_indices:
            # default value for hhg, if not specified
                deco = decorations[block]
                if len(deco) < 3:
                    deco.append(SymInfinity())
            # eventuellement rajouter une 4e deco avec l'indice d'orbite ? ...
            # third attribute
            self.decorations = decorations
            
        else:
            return NotImplemented
            # TODO: to implement
            dbs = decorated_block_system
            if not (isinstance(dbs, list) and all([isinstance(decor_block, list) for decor_block in dbs])):
                raise TypeError("not good")
            # self.decorations = ...

    ####
    # Methods for accessing attributes and data
    ####

    def diagonal_action(self):
        return self.finite_group

    def kernel(self):
        return self.ker

    def has_kernel(self):
        return self.kernel() != []
    
    def restriction_to_kernel(self):
        kernel_stab = libgap.Stabilizer(self.diagonal_action(), self.kernel(), libgap.OnSets)
        gp = libgap.Action(kernel_stab, self.kernel())
        return PermutationGroup(gap_group = gp)
        
    def blocks_from_superblocks(self):
        # return a list of sorted tuples (can be used as keys for self.decorations)
        return self.repres_from_superblocks

    def blocks_from_superblocks_as_lists(self):
        # return a list of lists
        result = []
        for block in self.repres_from_superblocks:
            result.append(list(block))
        return res

    def restriction_to_finite_block(self, block):  # block is an element from repres_from_superblocks
        return self.decorations[tuple(sorted(block))][0]

    _H0 = restriction_to_finite_block

    def wreath_base_on_block(self, block):         # block is an element from repres_from_superblocks
        return self.decorations[tuple(sorted(block))][1]

    _H1 = wreath_base_on_block

    def highly_homogeneous_extension(self, block): # block is an element from repres_from_superblocks
        return self.decorations[tuple(sorted(block))][2]

    _hhg = highly_homogeneous_extension

    ####
    # Methods for obtaining and building more usable data
    ####

    def list_of_superblocks(self):
        result = []
        for block in self.repres_from_superblocks:
            H0 = self._H0(block)
            H1 = self._H1(block)
            hhg = self._hhg(block)
            superblock = SingleSuperblock(H0, H1, hhg)
            result.append(superblock)
        return result

    def superblock_of_block(self, block):     # block is an element from repres_from_superblocks
        H0 = self._H0(block)
        H1 = self._H1(block)
        hhg = self._hhg(block)
        return SingleSuperblock(H0, H1, hhg)

    def subwreath_of_superblock(self, block): # block is an element from repres_from_superblocks
        H1 = self._H1(block)
        hhg = self._hhg(block)
        return WreathProductFiniteBlocks(H1, hhg)
                
    def normal_sg_finite_index(self): # ignores the kernel
        restrictions_of_normal_sg = []
        for block in self.blocks_from_superblocks():
            restrictions_of_normal_sg.append(self.subwreath_of_superblock(block))
        return DirectProductOfPoligomorphicGroups(restrictions_of_normal_sg)

    ####
    # Testing methods and __repr__
    ####

    def is_highly_homogeneous(self):
        return (not self.has_kernel()) and self.diagonal_action().degree() == 1
    
    def has_single_superblock(self):
        return len(self.blocks_from_superblocks()) == 1

    def is_wreath_product(self):
        # in case self has a single superblock, test wether H0 = H1
        return (not self.has_kernel()) and self.has_single_superblock() and self.decorations[self.blocks_from_superblocks()[0]][0] == self.decorations[self.blocks_from_superblocks()[0]][1]

    def __repr__(self):
        first_block = self.blocks_from_superblocks()[0]
        result = ""
        if self.is_highly_homogeneous():
            return repr(self.decorations[first_block][2])
        elif self.is_wreath_product():
            return repr(self.decorations[first_block][0]) + "  wreath  " + repr(self.decorations[first_block][2])
        elif self.has_single_superblock():
            if self.list_of_superblocks()[0].is_wreath_product():
                result = repr(self.decorations[first_block][1]) + "  wreath  " + repr(self.decorations[first_block][2])
            else:
                result = repr(self.decorations[first_block][1]) + "  wreath  " + repr(self.decorations[first_block][2]) + "  with diagonal action of  " + repr(self.decorations[first_block][0])
        else:
            result = "P-oligomorphic group of diagonal group (" + repr(self.diagonal_action()) + ") and normal subgroup of finite index " + repr(self.normal_sg_finite_index())
        if self.has_kernel():
            result += " with kernel domain " + repr(self.kernel())
        return result

    ####
    # Utilitary and computational methods 
    ####
    
    def _max_size_of_finite_block(self):
        max_size = 0
        for block in self.repres_from_superblocks:
            if max_size < len(block):
                max_size = len(block)
        return max_size

    def _union_finite_ages_per_degree(self):
        ages_per_deg = []
        max_deg = self._max_size_of_finite_block()
        ages_per_deg.append([[]])             # orbit of empty set
        for i in range(1, max_deg+1):
            ages_per_deg.append([])
        shift = 0
        for block in self.blocks_from_superblocks():
            for d in range(1, max_deg+1):     # for every orbital degree
                H1 = self._H1(block)
                age_d = age_explicit(H1, Integer(d))
                for orbit in age_d : 
                # shifts the domain of H1 so it is disjoint from the former
                    for subset in orbit :
                        for i in range(d) :   # d = orbital degree = size of the subsets
                            subset[i] += shift
                ages_per_deg[d] += age_d
            shift += len(block)
        return ages_per_deg
        

    def profile_series(self, variable='z'):
        r"""
        Return the generating series of the profile of the P-oligomorphic group.

        INPUT:

        - ``variable`` -- a variable, or variable name as a string (default: `'z'`)
        
        OUTPUT:
        
        - A series in ``variable`` with nonnegative integer coefficients.
          By default, a series in z over ZZ.
        
        EXAMPLES::
            
            sage: var('u')
            sage: S4 = SymmetricGroup(4)
            sage: G = PoligomorphicGroup_generic(S4, block_system = [[1],[2],[3],[4]], hhomogeneous_groups = [[[1], AutQQ()]])
            sage: G.profile_series(u)
            -1/4/(u^4 - 1) + 1/8/(u^2 - 1)^2 + 1/3/((u^3 - 1)*(u - 1)) - 1/4/((u^2 - 1)*(u - 1)^2) + 1/24/(u - 1)^4
            sage: G.profile_series(u) == 1/((1-u)*(1-u**2)*(1-u**3)*(1-u**4))
            True


        """
        from sage.rings.integer_ring import ZZ

        if isinstance(variable, str):       # creation or recuperation of the variable
            variable = ZZ[variable].gen()
        z = variable
        max_deg = self._max_size_of_finite_block()           # upper bound of computation
        finite_orbits = self._union_finite_ages_per_degree() # finite orbits except in kernel
        homom =     [None]
        dom =       [None]
        homom_ker = [None]
        dom_ker =   [None]
        for d in range(1, max_deg+1): # will enable to consider the actions on orbits of degree d(>1)
            # actions on orbits outside the kernel
            homom.append(libgap.ActionHomomorphism(self.finite_group, finite_orbits[d], libgap.OnSetsSets))
            dom.append(list(range(1, len(finite_orbits[d])+1))) # to compute cycle lengths below
            # same with the kernel, if relevant
            if self.has_kernel() :
                subsets_d = list(Combinations(self.kernel(), d))
                homom_ker.append(libgap.ActionHomomorphism(self.finite_group, subsets_d, libgap.OnSets)) 
                # in the kernel, each set is alone in its K-orbit,
                # so here (above) sets are identified to their orbits,
                # hence OnSets rather then OnSetsSets
                dom_ker.append(list(range(1, len(subsets_d)+1)))
        CCl = libgap.ConjugacyClasses(self.finite_group)
        polya_sum = 0
        for gbar in CCl:
            g = libgap.Representative(gbar)      
            # computation of the contribution of g
            W_g = 1
            for d in range(1, max_deg+1):        
            # run through the orbital degrees of the orbits of the H1's
                g_d = libgap.Image(homom[d], g)  # g acting on the orbits of degree d
                CT_d = libgap.CycleLengths(g_d, dom[d])
                for k in CT_d:
                    W_g /= (1 - z**(d*Integer(k)))
                # same with the kernel, only each subset may only be selected once
                if self.has_kernel() and d <= len(self.kernel()):
                    g_ker_d = libgap.Image(homom_ker[d], g)
                    CT_d = libgap.CycleLengths(g_ker_d, dom_ker[d])
                    for k in CT_d:
                        W_g *= 1 + z**(d*Integer(k))
            polya_sum += W_g*len(libgap.List(gbar))
        result = polya_sum / self.finite_group.order()
        return result
                


        
    
