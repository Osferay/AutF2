gap> F := FreeGroup( 2 );;
gap> a := AutomorphismOfF2( F, [ "s", 2 ] );;
gap> b := AutomorphismOfF2( F, [ "s", -1, 2, -1, 3, 2, -1, 3 ] );;
gap> AreConjugateAutomorphismsOfF2( a, b );
Automorphism of F2 with word [ -1, 2, 3 ]