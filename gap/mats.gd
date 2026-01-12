#! @Chapter Matrix groups
#!
#! @Chapter Matrix groups
#! @ChapterLabel mats
#! @ChapterTitle Matrix groups

#! @Description
#! Decides whether a pair of 2x2 matrices in $\operatorname{GL}_2(\mathbb{Z})$ are conjugate. If so returns the element $C$ such that $A=C^{-1}BC$.
#! @Arguments A,B
DeclareGlobalFunction( "ConjugacyGL2Z" );
#! @Description
#! Given a matrix $A$, computes a matrix $C$ such that $\langle -C, C \rangle = C_{\operatorname{GL}_2(\mathbb{Z})}(A)$ and an integer $n$ such that $C^n=A$.
#! @Arguments A,B
DeclareGlobalFunction( "CentralizerGL2Z" );