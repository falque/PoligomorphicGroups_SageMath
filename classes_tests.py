
r"""
Test file

Classes to be tested:

class PoligomorphicGroup(Parent)
    class HighlyHomogeneousGroup(PoligomorphicGroup)
        class AutQQ(HighlyHomogeneousGroup)
        class RevQQ(HighlyHomogeneousGroup)
        class AutQQCircle(HighlyHomogeneousGroup)
        class RevQQCircle(HighlyHomogeneousGroup)
        class SymInfinity(HighlyHomogeneousGroup)
    class DirectProductOfPoligomorphicGroups(PoligomorphicGroup)
    class PermutingSuperblocks(PoligomorphicGroup)
        class WreathProductInfiniteBlocks(PermutingSuperblocks)
        class SingleSuperblock(PermutingSuperblocks)
            class WreathProductFiniteBlocks(SingleSuperblock)
    class PoligomorphicGroup_generic(PoligomorphicGroup)
"""

C = [PermutationGroup([])]
S = [PermutationGroup([])]
D = [PermutationGroup([])]
A = [PermutationGroup([])]

for i in range(1,16):
  C.append(CyclicPermutationGroup(i))
  S.append(SymmetricGroup(i))
  D.append(DihedralGroup(i))
  A.append(AlternatingGroup(i))

# highly homogeneous groups

AutQ = AutQQ()
RevQ = RevQQ()
AutQc = AutQQCircle()
RevQc = RevQQCircle()
S_infty = SymInfinity()
assert S_infty.is_highly_homogeneous()

var('u')
assert S_infty.profile_series(u) == 1 / (1-u)

# wreath products with infinite blocks

AQwS5 = WreathProductInfiniteBlocks(S[5], AutQ)
assert AQwS5.action_on_infinite_blocks() == S[5]
assert AQwS5.restriction_to_infinite_block() == AutQQ()
hilbsym5 = 1/((1-u)*(1-u**2)*(1-u**3)*(1-u**4)*(1-u**5))
assert AQwS5.profile_series(u) == hilbsym5
assert AQwS5.profile_first_values(20) == hilbsym5.series(u, 21).coefficients(u, sparse=False)

# wreath products with finite blocks

S_inf = WreathProductFiniteBlocks(PermutationGroup([]))
assert S_inf.is_highly_homogeneous()
C7wS = WreathProductFiniteBlocks(C[7])
assert C7wS.size_of_max_blocks() == 7
assert C7wS._synchro() == PermutationGroup([])
hilbC7 = 1/((1-u)*(1-u**2)**3*(1-u**3)**5*(1-u**4)**5*(1-u**5)**3*(1-u**6)*(1-u**7))
assert hilbC7 == C7wS.profile_series(u)
assert hilbC7.series(u, 21).coefficients(u, sparse=False) == C7wS.profile_first_values(20)
assert C7wS.as_PoligomorphicGroup().diagonal_action() == C[7]
assert C7wS.as_PoligomorphicGroup().profile_series(u) == C7wS.profile_series(u)

# single superblock and direct product

C7wS_asSingleSb = SingleSuperblock(C[7])
assert C7wS.is_wreath_product()
assert C7wS_asSingleSb.profile_series(u) == C7wS.profile_series(u)
G1 = SingleSuperblock(S[7], A[7], S_infty)
G2 = WreathProductInfiniteBlocks(S[5], RevQ)
h1 = G1.profile_series()
h2 = G2.profile_series()
G = DirectProductOfPoligomorphicGroups([G1, G2])
assert G.profile_series() == h1*h2

# general class

G = PoligomorphicGroup_generic(C[7])
assert G.profile_series(u) == C7wS.profile_series(u)
G = PoligomorphicGroup_generic(S[4], block_system = [[1],[2],[3],[4]], hhomogeneous_groups = [[[1], AutQQ()]])
assert G.profile_series(u) == 1/((1-u)*(1-u**2)*(1-u**3)*(1-u**4)) # hilbert series of symmetric polynomials


















