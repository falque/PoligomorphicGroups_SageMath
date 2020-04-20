# sage-poligomorphic-groups

## Context
Let G be a permutation group of a denumerable set E. The **profile** of G is the function f which counts, for each integer n, the (possibly infinite) number f(n) of orbits of G acting on the n-subsets of E. 
When this function takes only finite values and is in addition bounded by a polynomial, G is said to be **P-oligomorphic**. It was a conjecture by Peter Cameron that the profile of these groups was actually asymptotically equivalent to a polynomial.
In a recent [paper](https://www.lri.fr/~falque/paper.pdf) yet to be submitted (the ["short version"](https://www.mat.univie.ac.at/~slc/wpapers/FPSAC2018/83-Falque-Thiery.html) of which was published by the international conference FPSAC 2018), the conjecture was proved and a *classification* of all P-oligomorphic groups was described. 

The SageMath code hosted on this repository leans on this classification to offer an implementation of these groups. It allows the user to build and manipulate P-oligomorphic groups, and to compute their profile, using some variants of the *P'olya enumeration*.
A more extensive description can be found in Appendix A.2 of the [PhD thesis](https://www.lri.fr/~falque/manuscrit.pdf) of the author.

## Classes implemented
```py=
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
```


## How to use 
To use, launch SageMath and type commands:

```py=
load("new_methods_PermutationGroup_generic.py")
load("PoligomorphicGroup.py")
load("finer_classes_of_poligomorphic_gps.py")
load("PoligomorphicGroup_generic.py")
load("reshape_series.py")
load("classes_tests.py")
```
See then the Example notebook for examples of standard uses.

