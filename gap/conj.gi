ExponentOfBraidWord := function( b )
	local e, i;

	if IsEmpty( b ) then
		return 0;
	fi;

	if b[1] = 4 then
		e := 6;
	else
		if b[1] > 0 then
			e := 1;
		else
			e := -1;
		fi;
	fi;

	for i in [2..Length(b)] do
		if b[i] > 0 then
			e := e + 1;
		else
			e := e - 1;
		fi;
	od;

	return e;

end;

ConjugacySA2 := function( a, b )
	local	e, b1, b2, j, conj, v, d, a1, f, z0, h0;

	b1 := WordOfAutomorphismOfF2( a );
	b2 := WordOfAutomorphismOfF2( b );
	b1 := WordOfSpecialAutomorphismOfF2ToBraidWord( b1 );
	b2 := WordOfSpecialAutomorphismOfF2ToBraidWord( b2 );
	e  := ExponentOfBraidWord( b1 )-ExponentOfBraidWord( b2 );

	if e mod 12 <> 0 then
		return false;
	else
		e := e/6;
		if e > 0 then 
			b2 := Concatenation( b2, ListWithIdenticalEntries( e, 4 ) );
		else
			b2 := Concatenation( b2, ListWithIdenticalEntries( AbsInt(e), -4 ) );
		fi;
	fi;
	
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

	cent := CentralizerAutomorphismOfF2InSA( a1^2 );
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

WordEqualToInverseBySigma := function( F, i, w )
	local sigma, r, s, sym;

	if i = 2 then
		sigma := AutomorphismOfF2( F, [ "s" ] );
	elif i = 3 then
		sigma := AutomorphismOfF2( F, [ -1, 2, 3, -1, 2, 3 ] );
	elif i = 4 then
		sigma := AutomorphismOfF2( F, [ "s", -1, 2, 3 ] );
	else
		Error( "i has to be 2,3 or 4.");
	fi;

	if ImageByAutomorphismOfF2( sigma, w ) <> w^-1 then
		return fail;
	fi;

	r := LetterRepAssocWord( w );
	s := NumberSyllables( w );

	if IsEvenInt( Length(r) ) then
		r := r{[1..Length(r)/2]};
		return rec( w1 := AssocWordByLetterRep( FamilyObj(w), r ), sym := One(F) );
	
	else
		if r = 1 then
			return rec( w1:= One(F), sym := w );
		fi;

		sym := r[(Length(r)+1)/2];
		r   := r{[1..(Length(r)-1)/2]};
		if sym < 0 then
			Add( r, sym );
			sym := AbsInt( sym );
		fi;
		return rec( w1 := AssocWordByLetterRep( FamilyObj(w), r ), sym := AssocWordByLetterRep( FamilyObj(w), [sym] ) );
	fi;
end;

WordConjugateToInverseBySigma := function( F, i, w )
	local sigma, u, u1, w0, w1, w2;

	if i = 2 then
		sigma := AutomorphismOfF2( F, [ "s" ] );
	elif i = 3 then
		sigma := AutomorphismOfF2( F, [ -1, 2, 3, -1, 2, 3 ] );
	elif i = 4 then
		sigma := AutomorphismOfF2( F, [ "s", -1, 2, 3 ] );
	else
		Error( "i has to be 2,3 or 4.");
	fi;

	u := RepresentativeAction( F, w^-1, ImageByAutomorphismOfF2( sigma, w ) );

	if IsBool( u ) then
		return fail;
	fi;

	u1 := WordEqualToInverseBySigma( F, i, u );
	w0 := w^(u1.w1)*u1.sym;
	w1 := WordEqualToInverseBySigma( F, i, w0 );
	w2 := w1.w1*w1.sym*ImageByAutomorphismOfF2( sigma, w1.w1^-1 )*u1.sym^-1;
	
	return rec( w2 := w2, w1 := w1.w1, sym := [ w1.sym, u1.sym ], u1 := u1.w1, u := u );
end;

SolveQuestion2Order2 := function( F, i, z )

	local sigma, w, r, p, w2, v, u, z1, z2, w1, s;

	if i = 2 then
		sigma := AutomorphismOfF2( F, [ "s" ] );
	elif i = 3 then
		sigma := AutomorphismOfF2( F, [ -1, 2, 3, -1, 2, 3 ] );
	elif i = 4 then
		sigma := AutomorphismOfF2( F, [ "s", -1, 2, 3 ] );
	else
		Error( "i has to be 2,3 or 4.");
	fi;

	w    := RootFreeGroup( z );
	r    := w[2];
	w    := w[1];

	w2 := WordConjugateToInverseBySigma( F, i, w );
	v  := w2.u1^-1;
	u  := w2.u^-1;
	z1 := w2.sym[2];
	z2 := w2.sym[1];
	w1 := w2.w1;
	w2 := w2.w2;

	if IsEvenInt( r ) then
		return rec( h := w^(-r/2), u := u );
	else
		p  := (r-1)/2;
		
		if z2 <> z1 then
			return false;
		fi;

		s := v^-1*z1*ImageByAutomorphismOfF2( sigma, v )*u;
		if s = s^0 then
			s := 0;
		else
			s := RootFreeGroup( s );
			s := s[2];
		fi;
		
		if IsEvenInt( s ) then
			return rec( h := v^-1*w2^(-s/2)*w1^-1*w2^(s/2-p)*v, u := u );
		else
			return rec( h := v^-1*w2^(-1*(s+1)/2)*w1*w2^((s-1)/2-p)*v, u := u );
		fi;
	fi;
end;

AutOrder2Conjugate := function( a )

	local F, s1, s2, s3, S, A, C, c, i;

	F  := a!.freeGroup;
	s1 := AutomorphismOfF2( F, [ "s" ] );
	s2 := AutomorphismOfF2( F, [ -1, 2, 3, -1, 2, 3 ] );
	s3 := AutomorphismOfF2( F, [ "s", -1, 2, 3 ] );

	S  := List( [s1,s2,s3], MatrixRepresentationOfAutomorphismOfF2 );
	A  := MatrixRepresentationOfAutomorphismOfF2( a );

	for i in [1,2,3] do
		C := ConjugacyGL2Z( A, S[i] );

		if not IsBool( C ) then
			c := AutomorphismOfF2ByMatrix( F, C );
			return rec( c := c, i := i+1 );
		fi;
	od;

end;

SolveQuestion2 := function( a, z )

	local F, c, i, h, A, t, d, r, C, p, b, w, a0, A0;

	if z = z^0 then
		return z;
	fi;

	if ImageByAutomorphismOfF2( a, z ) <> z^-1 then
		return false;
	fi;

	F := a!.freeGroup;

	if Order( a ) = 2 then
		c := AutOrder2Conjugate( a );
		i := c.i;
		c := c.c;
		h := SolveQuestion2Order2( F, i, ImageByAutomorphismOfF2( c, z ) );

		if IsBool(h) then
			return false;
		else
			return ImageByAutomorphismOfF2( c^-1, h.h );
		fi;
	fi;

	A := MatrixRepresentationOfAutomorphismOfF2( a );
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
		h := SolveQuestion2Order2( F, 3, ImageByAutomorphismOfF2( z ) );

		return ImageByAutomorphismOfF2( b^-1, h.h );
	fi;
end;	

InstallGlobalFunction( AreConjugateAutomorphismsOfF2, function( a, b )
	local	b1, v, d, F, s3, B, C, i, j, a1, f, z0, h0;

	if IsSpecialAutomorphismOfF2( a ) and IsSpecialAutomorphismOfF2( b ) then
		return ConjugacySA2( a, b );
	fi;

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
					h0 := ConjugacyAutomorphismOfF2( a!.freeGroup, h0 );
					return d*(b1*h0*f)^-1;
				fi;
			fi;
		od;
	od;

	return false;

end );