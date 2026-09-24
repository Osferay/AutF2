#! @Chapter Relalgs
#! @Section Matrix groups

#! @Description
#! Decides whether two $2 \times 2$ matrices in $\mathrm{GL}_2(\mathbb{Z})$ are conjugate.
#! If so, returns the element $C$ such that $A = C^{-1} B C$.
#! @Arguments A,B
DeclareGlobalFunction( "ConjugacyGL2Z" );
#! @Description
#! Given a matrix $A$, computes a matrix $C$ and an integer $n$ such that
#! $\operatorname{Cent}_{\mathrm{GL}_2(\mathbb{Z})}(A) = \langle -C, C \rangle$ and $C^n = A$.
#! @Arguments A,B
DeclareGlobalFunction( "CentralizerGL2Z" );
#! @Description
#! Given a matrix $A$, decides whether $A \in \langle \mathrm{gens} \rangle$ and,
#! if so, computes a word $w$ in the generators representing $A$.
#! @Arguments gens,A
DeclareGlobalFunction( "MembershipSubgroupSL2Z" );