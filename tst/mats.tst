gap> A := [ [ 1, 0 ], [ 0, -1 ] ];;
gap> B := [ [ 1, 0 ], [ 6, -1 ] ];;
gap> ConjugacyGL2Z( A, B );
[ [ 1, 0 ], [ -3, 1 ] ]
gap> A := [ [ 1, 1 ], [ 0, -1 ] ];;
gap> B := [ [ 4, -1 ], [ 15, -4 ] ];;
gap> ConjugacyGL2Z( A, B );
[ [ 4, -1 ], [ -3, 1 ] ]
gap> A := [ [ -4, 1 ], [ -13, 3 ] ];;
gap> B := [ [ 1, 1 ], [ -3, -2 ] ];;
gap> ConjugacyGL2Z( A, B );
[ [ 1, -1 ], [ -2, 3 ] ]
gap> A := [ [ -3, -1 ], [ 10, 3 ] ];;
gap> B := [ [ 2, -1 ], [ 5, -2 ] ];;
gap> C := ConjugacyGL2Z( A, B );
[ [ 1, 0 ], [ -5, 1 ] ]
gap> A := [ [ -3, -1 ], [ 13, 4 ] ];;
gap> B := [ [ 2, -1 ], [ 3, -1 ] ];;
gap> ConjugacyGL2Z( A, B );
[ [ 1, -1 ], [ -2, 3 ] ]
gap> A := [ [ 5, 2 ], [ -8, -3 ] ];;
gap> B := [ [ -9, -2 ], [ 50, 11 ] ];;
gap> ConjugacyGL2Z( A, B );
[ [ -1, 0 ], [ 7, 1 ] ]
gap> CentralizerGL2Z(B);
rec( exponent := 2, gen := [ [ -4, -1 ], [ 25, 6 ] ] )
gap> A := [ [ 14, 3 ], [ -75, -16 ] ];;
gap> B := [ [ -1, 0 ], [ -3, -1 ] ];;
gap> ConjugacyGL2Z( A, B );
[ [ 0, -1 ], [ 1, 5 ] ]
gap> CentralizerGL2Z(A);
rec( exponent := 3, gen := [ [ 4, 1 ], [ -25, -6 ] ] )
gap> A := [ [ -5, -1 ], [ 29, 6 ] ];;
gap> B := [ [ -2, -5 ], [ 1, 3 ] ];;
gap> ConjugacyGL2Z( A, B );
[ [ -1, -2 ], [ 5, 11 ] ]
gap> A := [ [ 5, 14 ], [ -1, -3 ] ];;
gap> B := [ [ -1, -2 ], [ 1, 3 ] ];;
gap> C := ConjugacyGL2Z( A, B );
[ [ 3, 4 ], [ -1, -1 ] ]
gap> A := [ [ -17, -54 ], [ 11, 35 ] ];;
gap> B := [ [ 25, 29 ], [ -6, -7 ] ];;
gap> ConjugacyGL2Z( A, B );
[ [ -2, -5 ], [ 1, 2 ] ]
gap> CentralizerGL2Z(B);
rec( exponent := 1, gen := [ [ 25, 29 ], [ -6, -7 ] ] )
gap> M := [ [ 97, 216 ], [ 22, 49 ] ];;
gap> MembershipCommutatorSL2Z(M);
[ 1, 1, 2, 1, 2, 2, 1 ]
gap> M := [ [ 216, 97 ], [ 49, 22 ] ];;
gap> w := WordGL2ZinST( M );
rec( det := -1, wS := [ 1, 1, 1, 1, 1, 1, 1, 1, 0 ], 
  wT := [ 1, 2, 2, 2, 4, 2, 6, 2, 1 ] )
gap> MatrixGL2ZbyWordST( w );
[ [ 216, 97 ], [ 49, 22 ] ]
gap> gens := [ [ [ -1, 1 ], [ -1, 0 ] ], [ [ -1, 0 ], [ -1, -1 ] ], [ [ 0, 1 ], [ -1, -1 ] ] ];;