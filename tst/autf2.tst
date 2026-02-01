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
gap> lcf := LeftCanonicalFormAutomorphismOfF2( aut );;
gap> AutomorphismOfF2( F, lcf );
f1 -> f1^-1*f2^-1*f1*f2*f1
f2 -> f1^-1*f2*f1
gap> phi1*phi3 = phi3*phi1;
true
gap> phi2*phi3*phi2 = phi3*phi2*phi3;
true
gap> phi1^-1*phi2*phi1^-1 = phi2*phi1^-1*phi2;
true
gap> phi1^sigma = phi2;
true
gap> phi3^sigma = phi3^-1*phi2^-1*phi1*phi2*phi3;
true
gap> phi4 := phi3^-1*phi2^-1*phi1*phi2*phi3;;
gap> cent := CentralizerAutomorphismOfF2( phi4 );;
gap> ForAll( cent, x -> phi4*x = x*phi4 );
true
gap> ImageAutomorphismOfF2( phi1*phi3*phi4, F.1*F.2 );
f2^-1*f1*f2*f1^-1*f2
gap> imgs := ImagesAutomorphismOfF2( aut );;
gap> AutomorphismOfF2ByImages( F, imgs[1], imgs[2] );
f1 -> f1^-1*f2^-1*f1*f2*f1
f2 -> f1^-1*f2*f1
gap> B := [ [ 1, 1 ], [ -3, -2 ] ];;
gap> AutomorphismOfF2ByMatrix( F, B );
f1 -> f1*f2
f2 -> (f2^-1*f1^-1)^2*f1^-1