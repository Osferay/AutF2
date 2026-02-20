ConjugacySA2 := function( a, b )
	local	b1, b2, j, conj, v, d, a1, f, z0, h0;

	b1 := WordOfAutomorphismOfF2( a );
	b2 := WordOfAutomorphismOfF2( b );

	if IsSpecialAutomorphismOfF2( a ) and IsSpecialAutomorphismOfF2( b ) then
		b1 := WordOfSpecialAutomorphismOfF2ToBraidWord( b1 );
		b2 := WordOfSpecialAutomorphismOfF2ToBraidWord( b2 );
		j  := rec( word1 := b1, word2 := b2 );

		AutF2WriteJSON( j );
		AutF2CallCpp( "conj" );
		conj := AutF2ReadJSON();

		if IsBool( conj[1] ) then
			return false;
		else
			conj := BraidWordToWordOfSpecialAutomorphismOfF2( conj );
			return AutomorphismOfF2( a!.freeGroup, conj );
		fi;
	fi;
end;

ReduceToQuestion1 := function( a, b )

	local	F, A, B, C, c, v;

	F := a!.freeGroup;
	A := MatrixRepresentationOfAutomorphismOfF2( a );
	B := MatrixRepresentationOfAutomorphismOfF2( b );

	C := ConjugacyGL2Z( A, B );

	if IsBool( C ) then
		return false;
	fi;

	c := AutomorphismOfF2ByMatrix( F, C );

	v := b^-1*a^c;
	if IsConjugacyAutomorphismOfF2( v ) then
		return [v,c];
	else
		Error( "v should be a conjugacy automorphism.");
	fi;

end;

ReduceToQuestion2 := function( b, v, b1 )

	local	z, a1, t, cent, gens, M, w, t0, i, f, z0;

	z  := b1^-1*b^-1*b1*b*v;
	a1 := b^b1;
	t  := ConjugacySA2( a1^2, (a1*z)^2 );

	if IsBool(t) then
		return false;
	fi;

	cent := CentralizerAutomorphismOfF2( a1^2 );
	gens := List( cent, MatrixRepresentationOfAutomorphismOfF2 );
	gens := ReduceParallelGeneratorSetSubgroupSL2Z( gens, cent );
	M := MatrixRepresentationOfAutomorphismOfF2( t^-1 );
	w := MembershipSubgroupSL2Z( gens[1], M );

	if IsBool(w) then
		return false;
	fi;

	t0 := ProductAutomorphismsOfF2ByWord( gens[2], w );
	f  := t0*t;

	if not IsConjugacyAutomorphismOfF2( f ) then
		Error( "f should be a conjugacy automorphism.");
	fi;

	z0  := f*z^-1*a1^-1*f^-1*a1;

	if IsConjugacyAutomorphismOfF2( z0 ) then
		z0 := ConjugacyElementConjugacyAutomorphismOfF2( z0 );
		return [a1, z0, f];
	else
		Error( "z0 has to be a conjugacy automorphism");
	fi;

end;

SolveQuestion2 := function( a, z )

	local A, F, t, d, r, C, p, b, w, a0, A0, gens, c, c1;

	if z = z^0 then
		return z;
	fi;

	if ImageByAutomorphismOfF2( a, z ) <> z^-1 then
		return false;
	fi;

	A := MatrixRepresentationOfAutomorphismOfF2( a );
	F := a!.freeGroup;
	t := Trace( A^2 );
	d := DeterminantIntMat( A^2 );
	r := 1;

	if AbsInt(t) = 2 and d = 1 then
		C := ConjugacyClassParabolicMatrix( A ).rep;
		t := C[1][2]; 
		p := AutomorphismOfF2( F, [ -3, -2, -1, 2, 3 ] );
		b := ConjugacySA2( a^2, p^t );

		if not IsBool( b ) then
			r := 2;
		fi;
	fi;

	if r = 1 then
		w := RootFreeGroup( z );
		r := w[2];
		w := w[1];

		if IsEvenInt( r ) then	
			return w^(-r/2);
		else
			return false;
		fi;

	else
		a0 := a^b;
		A0 := MatrixRepresentationOfAutomorphismOfF2( a0 );
		if A0[1] = 1 then 
			return false;
		fi;

		w    := RootFreeGroup( z );
		r    := w[2];
		w    := w[1];
		gens := GeneratorsOfGroup( F );
		C    := [ gens[1], gens[1]^2, gens[2]^-1*gens[1]*gens[2], gens[2]^-1*gens[1]^-1*gens[2] ];
		
		for c in C do
			c1   := RepresentativeAction( F, c, w );
			if not IsBool( c1 ) then
				if IsEvenInt( r ) then
					return w^(r/2);
				else
					return false;
				fi;
			fi;
		od;

		if IsEvenInt( r ) then
			return w^(-r/2);
		fi;
	fi;
end;	

InstallGlobalFunction( AreConjugateAutomorphismsOfF2, function( a, b )
	local	b1, v, d, F, s3, B, C, i, j, a1, f, z0, h0;


	v  := ReduceToQuestion1( a, b );

	if IsBool(v) then
		return false;
	fi;

	d  := v[2];
	v  := v[1];

	F  := a!.freeGroup;
	s3 := AutomorphismOfF2( F, [-1,2,-1,3,2,3] );
	B  := MatrixRepresentationOfAutomorphismOfF2( b );
	C  := CentralizerGL2Z( B );
	
	for i in [0,1] do
		for j in [0..(C.exponent-1)] do	
			b1 := s3^i*AutomorphismOfF2ByMatrix( F, C.gen )^j;
			a1 := ReduceToQuestion2( b, v, b1 );
			
			if not IsBool( a1 ) then
				f  := a1[3];
				z0 := a1[2];
				a1 := a1[1];
				h0 := SolveQuestion2( a1, z0 );

				if not IsBool(h0) then
					h0 := ConjugacyAutomorphismOfF2( F, h0 );
					return d*(b1*h0*f)^-1;
				fi;
			fi;
		od;
	od;

	return false;

end );