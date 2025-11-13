gap> F := FreeGroup( 2 );;
gap> phi1 := AutomorphismOfF2( F, [1] );;
gap> phi2 := AutomorphismOfF2( F, [2] );;
gap> phi3 := AutomorphismOfF2( F, [3] );;
gap> sigma := AutomorphismOfF2( F, ["s"] );;
gap> IsIdentityAutomorphismOfF2( Inverse( phi1 )*Inverse( phi3 )*phi1*phi3 );
true
gap> IsIdentityAutomorphismOfF2( Inverse( phi3 )*Inverse( phi2 )*Inverse(phi3)*phi2*phi3*phi2 );
true
gap> IsIdentityAutomorphismOfF2( phi2*Inverse( phi1 )*phi2*phi1*Inverse( phi2 )*phi1 );
true
gap> IsIdentityAutomorphismOfF2( (Inverse( phi1 )*phi2*phi3 )^4 );
true
gap> IsIdentityAutomorphismOfF2( phi1^sigma*Inverse( phi2 ) );
true
gap> IsIdentityAutomorphismOfF2( phi3^sigma*Inverse( phi3 )*Inverse( phi2 )*Inverse( phi1 )*phi2*phi3 );
true
gap> Order( sigma );
2
gap> aut := AutomorphismOfF2( F, [ 1, 3, 2, -3, -2, 1, 2, 3 ] );;
gap> ConjugacyElementConjugacyAutomorphismOfF2( aut );
f2*f1