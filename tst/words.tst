gap> F := FreeGroup( 2 );;
gap> V := [ F.1*F.2^2, F.1*F.2^-1, F.2^2 ];;
gap> NielsenReducedSetBacktrack( V, [[1],[2],[3]] );
[ [ f1, f2 ], [ [ 1, -3, -3, -2, 1, -3, -3, -2, 1 ], [ -1, 2, 3, 3, -1, 2, 3, 3, -1, 1, -3, -3, -2, 1 ] ] ]
gap> CosetRepresentativeReducedNielsenSet( NielsenReducedSet(V), F.1*F.2 );
[ f1*f2, <identity ...> ]