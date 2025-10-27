#
# AutF2: Computations for the automorphisms group of F2
#
#! @Chapter Introduction
#!
#! AutF2 is a package which does some
#! interesting and cool things
#!
#! @Chapter Automorphisms
#! @ChapterLabel autos
#! @ChapterTitle Automorphisms

#! @Description
#! Object that represents an automorphism.
DeclareCategory( "IsAutomorphismOfF2", IsObject);

DeclareCategoryCollections( "IsAutomorphismOfF2" );
DeclareRepresentation( "RepAutomorphismOfF2", 
                        IsAttributeStoringRep, 
                        ["functions", "freeGroup", "images", "word"] );
DeclareAttribute( "AutomorphismOfF2Family", IsFamily );

#! @Description
#! Constructor of the automorphism object given a list of functions. 
#! @Arguments F, list
DeclareOperation("AutomorphismOfF2", [ IsFreeGroup, IsList ] );
#! @Description
#! Returns whether or not an automorphism is trivial.
#! @Arguments aut
DeclareAttribute( "IsIdentityAutomorphismOfF2", IsAutomorphismOfF2 );
#! @Description
#! Returns the representation of an automorphism in $GL_2(\mathbb{Z})$.
#! @Arguments aut
DeclareAttribute( "MatrixRepresentationOfAutomorphismOfF2", IsAutomorphismOfF2);
#! @Description
#! Returns whether or not an automorphism is special, that is its representation in $GL_2(\mathbb{Z})$ has determinant 1.
#! @Arguments aut
DeclareAttribute( "IsSpecialAutomorphismOfF2", IsAutomorphismOfF2 );