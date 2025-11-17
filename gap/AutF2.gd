#
# AutF2: Computations for the automorphisms group of F2
#
#! @Chapter Introduction
#!
#! AutF2 is a package to compute with automorphisms of F2.
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
#! Constructor of the automorphism object given a word of automorphisms. 
#! @Arguments F, list
DeclareOperation("AutomorphismOfF2", [ IsFreeGroup, IsList ] );
#! @Description
#! Returns the word of the given automorphism.
#! @Arguments aut
DeclareAttribute( "WordOfAutomorphismOfF2", IsAutomorphismOfF2 );
#! @Description
#! Returns the images of the generators of $F_2$ by the given automorphism.
#! @Arguments aut
DeclareAttribute( "ImagesAutomorphismOfF2", IsAutomorphismOfF2 );
#! @Description
#! Returns the image of the given word by the given automorphism.
#! @Arguments aut, w
DeclareOperation( "ImageAutomorphismOfF2", [ IsAutomorphismOfF2, IsAssocWordWithInverse ] );
#! @Description
#! Returns whether or not an automorphism is trivial.
#! @Arguments aut
DeclareProperty( "IsIdentityAutomorphismOfF2", IsAutomorphismOfF2 );
#! @Description
#! Returns the representation of an automorphism in $GL_2(\mathbb{Z})$.
#! @Arguments aut
DeclareAttribute( "MatrixRepresentationOfAutomorphismOfF2", IsAutomorphismOfF2);
#! @Description
#! Returns whether or not an automorphism is special, that is its representation in $GL_2(\mathbb{Z})$ has determinant 1.
#! @Arguments aut
DeclareAttribute( "IsSpecialAutomorphismOfF2", IsAutomorphismOfF2 );
#! @Description
#! Returns whether the given automorphism is determined by conjugating an element
#! @Arguments aut
DeclareProperty( "IsConjugacyAutomorphismOfF2", IsAutomorphismOfF2 );
#! @Description
#! Getter of the element that determines a conjugacy automorphism
#! @Arguments aut
DeclareAttribute( "ConjugacyElementConjugacyAutomorphismOfF2", IsAutomorphismOfF2 );
#! @Description
#! Returns the automorphism determinated by the conjugacy of the given word in $F_2$
#! @Arguments F, w
DeclareOperation( "ConjugacyOfAutomorphismOfF2", [ IsFreeGroup, IsAssocWordWithInverse ] );
#! @Description
#! Computes the left canonical form of the given automorphism
#! @Arguments aut
DeclareAttribute( "LeftCanonicalFormAutomorphismOfF2", IsAutomorphismOfF2 );