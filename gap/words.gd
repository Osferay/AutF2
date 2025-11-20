#! @Description
#! Returns whether the word u is cyclicaly reduced or not.
#! @Arguments u.
DeclareProperty( "IsCyclicalyReducedWord", IsAssocWordWithInverse );
#! @Description
#! Returns whether the word v is a subword of u or not.
#! @Arguments u,v.
DeclareOperation( "IsSubword", [ IsAssocWordWithInverse, IsAssocWordWithInverse ] );
#! @Description
#! Returns whether the word v is a prefix of u or not.
#! @Arguments u,v.
DeclareOperation( "IsPrefix", [ IsAssocWordWithInverse, IsAssocWordWithInverse ] );
#! @Description
#! Returns whether the word v is a suffix of u or not.
#! @Arguments u,v.
DeclareOperation( "IsSuffix", [ IsAssocWordWithInverse, IsAssocWordWithInverse ] );