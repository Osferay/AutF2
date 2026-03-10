gap> F := FreeGroup( 2 );;
gap> w := F.1^3*F.2^5*F.1^-5*F.2^-3;;
gap> WordEqualToInverseBySigma( F, 2, w );
rec( sym := <identity ...>, w1 := f1^3*f2^5 )
gap> w := F.1^3*F.2^2*F.1^3;;
gap> WordEqualToInverseBySigma( F, 3, w );
rec( sym := <identity ...>, w1 := f1^3*f2 )
gap> w := F.2*F.1^7*F.2^-1*F.1^7*F.2;;
gap> WordEqualToInverseBySigma( F, 3, w );
rec( sym := f2, w1 := f2*f1^7*f2^-1 )
gap> w := F.2*F.1^-1*F.2^-4*F.1*F.2;;
gap> WordEqualToInverseBySigma( F, 4, w );
rec( sym := <identity ...>, w1 := f2*f1^-1*f2^-2 )
gap> w := F.1^-5*F.2^-3*F.1^5;;
gap> WordEqualToInverseBySigma( F, 4, w );
rec( sym := f2, w1 := f1^-5*f2^-2 )
gap> w := F.1^-2*F.2^-3*F.1^3*F.2^2*F.1*F.2^2*F.1^-2*F.2^-1;;
gap> WordConjugateToInverseBySigma( F, 2, w );
rec( sym := [ <identity ...>, <identity ...> ], u := f2*f1^2*f2^-2*f1^-1, 
  u1 := f2*f1^2, w1 := (f1^-2*f2^-1)^2*f2^-2, 
  w2 := (f1^-2*f2^-1)^2*f2^-2*f1^2*(f1*f2^2)^2 )
gap> h := SolveQuestion2Order2( F, 2, w );
rec( h := f2*f1^2*f2^3*f1^2, u := f1*f2^2*f1^-2*f2^-1 )
gap> w := F.1^-1*F.2^3*F.1;;
gap> WordConjugateToInverseBySigma( F, 3, w );
rec( sym := [ f2, <identity ...> ], u := f1^-2, u1 := f1^-1, w1 := f2, 
  w2 := f2^3 )
gap> w := F.1*F.2*F.1^-5*F.2*F.1*F.2^2*F.1^7*F.2^2;;
gap> WordConjugateToInverseBySigma( F, 3, w );
rec( sym := [ f1, f1 ], u := f2^-2*f1^-7*f2^-2, u1 := f2^-2*f1^-4, 
  w1 := f1^4*f2^2*f1*f2*f1^-3, w2 := f1^4*f2^2*f1*f2*f1^-5*f2*f1*f2^2*f1^3 )
gap> h := SolveQuestion2Order2( F, 3, w );
rec( h := f2^-1*(f2^-1*f1^-1)^2, u := f2^2*f1^7*f2^2 )
gap> w := F.1^-3*F.2*F.1^3*F.2^-4;; 
gap> WordConjugateToInverseBySigma( F, 4, w );
rec( sym := [ f2, <identity ...> ], u := f2^4, u1 := f2^2, w1 := f2^-2*f1^-3, 
  w2 := f2^-2*f1^-3*f2*f1^3*f2^-2 )
gap> w := F.1*F.2^-5*F.1^-1*F.2^2*F.1^-1*F.2^3*F.1*F.2^2;;
gap> WordConjugateToInverseBySigma( F, 4, w );
rec( sym := [ f2, f2 ], u := f2^-2*f1^-1*f2^-3*f1*f2^-2, 
  u1 := f2^-2*f1^-1*f2^-2, w1 := (f2^2*f1)^2*f2^-3, 
  w2 := (f2^2*f1)^2*f2^-5*f1^-1*f2^2*f1^-1*f2 )
gap> h := SolveQuestion2Order2( F, 4, w );
rec( h := f2^-2*f1^-1*f2*f1^-1, u := f2^2*f1^-1*f2^3*f1*f2^2 )
gap> a := AutomorphismOfF2( F, [ "s", 2 ] );;
gap> b := AutomorphismOfF2( F, [ "s", -1, 2, -1, 3, 2, -1, 3 ] );;
gap> AreConjugateAutomorphismsOfF2( a, b );
Automorphism of F2 with word [ -1, 2, 3 ]
gap> a := AutomorphismOfF2( F, [ "s", 3, -2 ] );;
gap> b := AutomorphismOfF2( F, [ "s", "d", 2, 3, 2, 2, 3, 2, -1, -1 ] );;
gap> AreConjugateAutomorphismsOfF2( a, b );
Automorphism of F2 with word [ "s", -1 ]