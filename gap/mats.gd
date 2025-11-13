#! @Chapter Matrix groups
#!
#! @Chapter Matrix groups
#! @ChapterLabel mats
#! @ChapterTitle Matrix groups

#! @Description
#! Decides whether a 2x2 matrix is reduced.
#! @Arguments A
DeclareGlobalFunction( "IsReducedMatrix22" );
#! @Description
#! Reduces a 2x2 matrix.
#! @Arguments A
DeclareGlobalFunction( "ReduceMatrix22" );
#! @Description
#! Decides whether a pair of 2x2 matrices are conjugate. If so returns the element C such that $A=C^{-1}BC$.
#! @Arguments A,B
DeclareGlobalFunction( "ConjugacyMatrix22" );
#! @Description
#! Computes a generating set of the centralizer of a 2x2 matrix.
#! @Arguments A,B
DeclareGlobalFunction( "CentralizerMatrix22" );