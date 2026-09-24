
#! @Chapter Preface
#! @ChapterLabel preface
#! @ChapterTitle Preface
#!
#! In this package we provide methods to compute with automorphisms of the free group of rank $2$.
#! The details of the algorithms are available in <Cite Key="autf2" />.
#! We use the following notation:
#! $$ \sigma=\begin{cases}
#! x \mapsto y,\\
#! y \mapsto x,
#! \end{cases}\quad
#! \phi_1=\begin{cases}
#! x \mapsto xy^{-1},\\
#! y \mapsto y,
#! \end{cases}\quad
#! \phi_2=\begin{cases}
#! x \mapsto x,\\
#! y \mapsto yx,
#! \end{cases}\quad
#! \phi_3=\begin{cases}
#! x \mapsto y^{-1}x,\\
#! y \mapsto y.
#! \end{cases} $$
#! Then $\operatorname{Aut}(F_2)$ has the following presentation.
#! $$ \begin{aligned}
#! \operatorname{Aut}(F_2) = \langle \sigma,\phi_1,\phi_2,\phi_3 \mid & \phi_1\phi_3=\phi_3\phi_1, \quad \phi_2\phi_3\phi_2 = \phi_3\phi_2\phi_3, \quad \phi_2\phi_1\phi_2 = \phi_1\phi_2\phi_1, \\
#!    & (\phi_1\phi_2\phi_3)^4=1, \quad \sigma^2=1, \quad \phi_1^{\sigma} = \phi_2^{-1}, \\ 
#!    & \phi_3\sigma = \sigma\phi_3^{-1}\phi_2^{-1}\phi_1^{-1}\phi_2\phi_3, \quad \phi_3^{-1}\sigma= \sigma\phi_3^{-1}\phi_2^{-1}\phi_1\phi_2\phi_3\rangle.
#! \end{aligned} $$
#! Denote the commutator subgroup of a group $G$ by $G^{\prime}$. 
#! For $\psi\in \operatorname{Aut}(F_2)$, we define the map $\psi^\ast$ by $\psi^\ast(x F_2^{\prime}) = \psi(x) F_2^{\prime}$.
#! If $\overline{\cdot}$ denotes the isomorphism between $\operatorname{Aut}(F_2/F_2^\prime)$ and $\operatorname{GL}_2(\mathbb{Z})$, then the map
#! \begin{align*}
#!    \Psi \colon & \operatorname{Aut}(F_2) \to \operatorname{GL}_2(\mathbb{Z}) \\
#!    & \psi \mapsto \overline{(\psi^\ast)}
#! \end{align*}
#! is a group homomorphism and the kernel of this map is $\operatorname{Inn}(F_2)$.
#! We denote by $\operatorname{SA}_2$ the subgroup of $\operatorname{Aut}(F_2)$ consisting of those automorphisms that are mapped into $\operatorname{SL}_2(\mathbb{Z})$ under $\Psi$.
#! Let $B_4$ denote the braid group on $4$ strands.
#! Then $\operatorname{SA}_2 \cong B_4 / Z(B_4)$.
#! Using this, we have that each automorphism can be written uniquely as a word of the form
#! $$ \alpha = \sigma^{\varepsilon_1}\Delta^{\varepsilon_2}\alpha_1\alpha_2\cdots\alpha_k $$
#! where $\varepsilon_1,\varepsilon_2 \in \{0,1\}$, $\Delta = \phi_1(\phi_2\phi_1)(\phi_3\phi_2\phi_1)$, and $\alpha_1,\dots,\alpha_k\in \operatorname{SA}_2$ whose images under the isomorphism between $\operatorname{SA}_2$ and $B_4 / Z(B_4)$ are simple braids.
#! We call this expression the left canonical form of $\alpha$.

#! @Chapter Automorphisms
#! @ChapterLabel autos
#! @ChapterTitle Automorphisms

#! @Section Automorphisms

#! @Description
DeclareCategory( "IsAutomorphismOfF2", IsObject);

DeclareCategoryCollections( "IsAutomorphismOfF2" );
DeclareRepresentation( "RepAutomorphismOfF2", 
                        IsAttributeStoringRep, 
                        ["freeGroup", "lcf"] );
DeclareAttribute( "AutomorphismOfF2Family", IsFamily );

#! @Description
#! Constructor of the automorphism object given a word of automorphisms. 
#! @Arguments F, list
DeclareOperation("AutomorphismOfF2", [ IsFreeGroup, IsList ] );
#! @Description
#! Returns the left canonical form of the given automorphism.
#! @Arguments aut
DeclareAttribute( "WordOfAutomorphismOfF2", IsAutomorphismOfF2 );
#! @Description
#! Returns the images of the generators of $F_2$ under  the given automorphism.
#! @Arguments aut
DeclareAttribute( "ImagesAutomorphismOfF2", IsAutomorphismOfF2 );
#! @Description
#! Returns the image of the given word under the given automorphism.
#! @Arguments aut, w
DeclareOperation( "ImageByAutomorphismOfF2", [ IsAutomorphismOfF2, IsAssocWordWithInverse ] );
#! @Description
#! @Arguments aut
DeclareProperty( "IsIdentityAutomorphismOfF2", IsAutomorphismOfF2 );
#! @Description
#! Returns the image of an automorphism in $\mathrm{GL}_2(\mathbb{Z})$ under the map $\Psi$.
#! @Arguments aut
DeclareAttribute( "MatrixRepresentationOfAutomorphismOfF2", IsAutomorphismOfF2);
#! @Description
#! @Arguments aut
DeclareProperty( "IsSpecialAutomorphismOfF2", IsAutomorphismOfF2 );
#! @Description
#! @Arguments aut
DeclareProperty( "IsConjugacyAutomorphismOfF2", IsAutomorphismOfF2 );
#! @Description
#! Getter of the element that determines a conjugacy automorphism.
#! @Arguments aut
DeclareAttribute( "ConjugacyElementConjugacyAutomorphismOfF2", IsAutomorphismOfF2 );
#! @Description
#! Returns the automorphism determinated by conjugacy of the given word in $F_2$
#! @Arguments F, w
DeclareOperation( "ConjugacyAutomorphismOfF2", [ IsFreeGroup, IsAssocWordWithInverse ] );