#! @Chapter Relalgs
#! @Section Free groups
#! @Description
#! Decides whether $u$ and $v$ are automorphically equivalent in $F$.
#! @Arguments F, u, v
DeclareOperation( "AreAutomorphicEquivalent", [ IsFreeGroup, IsAssocWordWithInverse, IsAssocWordWithInverse ] );

#! @Chapter Automorphisms
#! @Section Automorphisms

#! @Description
#! Returns an automorphism such that the images of the basis elements are the given words $u$ and $v$.
#! If the automorphism does not exist, it returns <C>false</C>.
#! @Arguments F, u, v
DeclareOperation( "AutomorphismOfF2ByImages", [ IsFreeGroup, IsAssocWordWithInverse, IsAssocWordWithInverse ] );
#! @Description
#! Returns an automorphism whose image in $\mathrm{GL}_2(\mathbb{Z})$ is the given matrix.
#! @Arguments F, M
DeclareOperation( "AutomorphismOfF2ByMatrix", [ IsFreeGroup, IsList ] );