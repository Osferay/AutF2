#! @Description
#! Decides whether u and v are automorphic equivalent in F.
#! @Arguments F, u, v
DeclareOperation( "AreAutomorphicEquivalent", [ IsFreeGroup, IsAssocWordWithInverse, IsAssocWordWithInverse ] );
#! @Description
#! Returns an automorphism such that the images of the basis are the given words u and v. If the automorphism does not exists, returns false.
#! @Arguments F, u, v
DeclareOperation( "AutomorphismOfF2ByImages", [ IsFreeGroup, IsAssocWordWithInverse, IsAssocWordWithInverse ] );
#! @Description
#! Returns an automorphism that has the given matrix as image in $GL_2(\mathbb{Z})$. If the automorphism does not exists, returns false.
#! @Arguments F, M.
DeclareOperation( "AutomorphismOfF2ByMatrix", [ IsFreeGroup, IsList ] );