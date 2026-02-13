WhiteheadAutomorphismsOfF2 := function( F )
    local gens, type1, type2, typec, sigma1, sigma2;

    gens    := [];
    gens[1] := AutomorphismOfF2( F, [1] );
    gens[2] := AutomorphismOfF2( F, [2] );
    gens[3] := AutomorphismOfF2( F, [3] );
    gens[4] := AutomorphismOfF2( F, [-3, -2, 1, 2, 3] );
    type2   := ShallowCopy( gens );
    Add( type2, gens[1]^-1 );
    Add( type2, gens[2]^-1 );
    Add( type2, gens[3]^-1 );
    Add( type2, gens[4]^-1 );

    type1   := [];
    sigma2  := AutomorphismOfF2( F, ["s"] );
    sigma1  := gens[1]^-1*gens[2]*gens[3];
    Add( type1, sigma1 );
    Add( type1, sigma1^2 );
    Add( type1, sigma1^3 );
    Add( type1, sigma2 );
    Add( type1, sigma1*sigma2 );
    Add( type1, sigma1^2*sigma2 );
    Add( type1, sigma1^3*sigma2 );

	typec   := [];
	Add( typec, gens[1]*gens[3] );
	Add( typec, gens[2]*gens[4] );
	Add( typec, (gens[1]*gens[3])^-1 );
	Add( typec, (gens[2]*gens[4])^-1 );

    return rec( type1 := type1, type2 := type2, typec := typec, all := Concatenation( type1, type2, typec ) );

end;

WhiteheadLengthReduction := function( W, w )

    local i, img;

    for i in [1..Length(W)] do
        img := ImageByAutomorphismOfF2( W[i], w );
        if Length( img ) < Length(w) then
            return rec( aut := W[i], img := img );
        fi;
    od;

    return false;

end;

WhiteheadReductionAlgorithm := function( F, w )

    local W, u, aut;

    W   := WhiteheadAutomorphismsOfF2( F ).type2;
    u   := ShallowCopy( w );
    aut := AutomorphismOfF2( F, [] );
    u   := WhiteheadLengthReduction( W, u );

    while not IsBool(u) do
        aut := aut * u.aut;
        u   := u.img;
        u   := WhiteheadLengthReduction( W, u );
    od;

    return rec( aut := aut, min := ImageByAutomorphismOfF2( aut, w ) );

end;

InstallMethod( AreAutomorphicEquivalent,
	"For a free group and two elements",
	[ IsFreeGroup, IsAssocWordWithInverse, IsAssocWordWithInverse ],
	function( F, u, v )
		local red1, red2, W, ae, conj, w;

		if IsCyclicalyReducedWord( u ) and IsCyclicalyReducedWord(v) then
			red1 := WhiteheadReductionAlgorithm( F, u );
			red2 := WhiteheadReductionAlgorithm( F, v );

			if Length( red1.min ) = Length( red2.min ) then
				W := WhiteheadAutomorphismsOfF2( F ).type1;

				for w in W do
					if ImageByAutomorphismOfF2( w, red1.min ) = red2.min then
						return red1.aut*w*red2.aut^-1;
					fi;
				od;

				return false;
			else 
				return false;
			fi;

		else
			red1 := CyclicallyReducedWord( u );
			red2 := CyclicallyReducedWord( v );

			ae   := AreAutomorphicEquivalent( F, red1, red2 );

			if IsBool(ae) then
				return false;
			fi;

			conj := RepresentativeAction( F, ImageByAutomorphismOfF2( ae, u ), v );
			if IsBool( conj ) then
				return false;
			fi;

			conj := ConjugacyAutomorphismOfF2( F, conj );
			return ae*conj;
		fi;
end );

InstallMethod( AutomorphismOfF2ByImages,
    "for a free group and two elements",
    [ IsFreeGroup, IsAssocWordWithInverse, IsAssocWordWithInverse ],
	function( F, u, v )

		local	reduce, gens, ae, W, r;

		reduce := function( W, aut, w )
			local	a, img;

			img := ImageByAutomorphismOfF2( aut, w );
			for a in W{[2..Length(W)]} do 
				if Length( ImageByAutomorphismOfF2( a, img ) ) < Length( img ) then
					return a;
				fi;
			od;
			return false;
		end;

		gens := GeneratorsOfGroup( F );
		ae := AreAutomorphicEquivalent( F, gens[1], u );

		if IsBool( ae ) then
			return false;
		else
			ae := Inverse( ae );
			W := WhiteheadAutomorphismsOfF2( F ).all{[7,9,11,13,15,17,19]};
			r := AutomorphismOfF2( F, [] );
			
			while not IsBool( r ) do
				ae := ae*r;
				if ImageByAutomorphismOfF2( ae, v ) = gens[2] then
					return Inverse( ae );
				elif ImageByAutomorphismOfF2( ae, v ) = gens[2]^-1 then
					return Inverse( ae*W[1] );
				fi;
				r := reduce( W, ae, v );
			od;
		fi;

		return false;
end );

InstallMethod( AutomorphismOfF2ByMatrix,
    "for a free group and a matrix",
    [ IsFreeGroup, IsList ],
    function( F, M )

		local w, aut, t, s, v;

		w   := WordGL2ZinSU( M );
		aut := AutomorphismOfF2( F, [] );
		t   := AutomorphismOfF2( F, [-1, 2, 3, 1] );
		s   := AutomorphismOfF2( F, [-1, 2, 3] );

		if w.e = 1 then
			aut := AutomorphismOfF2( F, ["s"] );
		fi;

		aut := aut*s^(w.s);

		for v in w.w do
			if v = 1 then 
        		aut := aut*s;
			else
				aut := aut*t;
			fi;
    	od;
		
		return aut;
end );