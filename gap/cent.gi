## Solves if a and a^b are conjugate by an inner automorphism ##
SolveInnerConjugacyAutF2 := function( a, b )
    local Id, v, a1, f, z0, h0;
    
    Id := AutomorphismOfF2( a!.freeGroup, [ ] );

    v  := a^-1*(a^b);
    a1 := ReduceToQuestion2( a, v, Id );

    if not IsBool( a1 ) then
		f  := a1[3];
		z0 := a1[2];
						
		h0 := SolveQuestion2( a, z0 );

		if not IsBool(h0) then
			h0 := ConjugacyAutomorphismOfF2( a!.freeGroup, h0 );
			return b*(h0*f)^-1;
		fi;
	fi;

    return false;
end;



InstallGlobalFunction( "CentralizerAutomorphismOfF2", function( aut )
    local C, sigma, c, s2, px, py, A, CA, b, D, e, Fix, tmp;

    if IsSpecialAutomorphismOfF2( aut ) then
        C     := CentralizerAutomorphismOfF2InSA( aut );
        sigma := AutomorphismOfF2( aut!.freeGroup, ["s"] );
        c     := AreConjugateAutomorphismsOfF2( aut, aut^sigma);
        if not IsBool(c) then
            Add( C, sigma*c );
        fi;
    
    elif Order( aut ) = infinity then
        s2 := AutomorphismOfF2( aut!.freeGroup, [ 1, 2, 3, 1, 2, 3 ] );

        A  := MatrixRepresentationOfAutomorphismOfF2( aut );
        CA := CentralizerGL2Z( A );
        b  := AutomorphismOfF2ByMatrix( aut!.freeGroup, CA.gen );
        D  := DivisorsInt( CA.exponent );
        
        C  := [];
        tmp:= [];
        c := SolveInnerConjugacyAutF2( aut, s2 );

        if not IsBool(c) then
            Add( C, c );
        fi;
        
        for e in D do
            if IsEvenInt( e ) and IsEmpty(C) then
                c := SolveInnerConjugacyAutF2( aut, s2*b^e );
                
                if not IsBool(c) then
                    Add( C, c );
                fi;
            fi;


            if IsEmpty( tmp ) then
                c := SolveInnerConjugacyAutF2( aut, b^e );
                
                if not IsBool(c) then
                    Add( tmp, c );
                    Add( C, c );
                fi;
            fi; 
        od;

        Fix := FixedSubgroupAutomorphismOfF2( aut );
        if not IsEmpty( Fix ) then
            Add( C, ConjugacyAutomorphismOfF2( aut!.freeGroup, Fix[1] ) );
        fi;

        return C;
    else
        s2    := AutomorphismOfF2( aut!.freeGroup, [ 1, 2, 3, 1, 2, 3 ] );
        sigma := AutomorphismOfF2( aut!.freeGroup, ["s"] );
        c     := AreConjugateAutomorphismsOfF2( sigma, aut );
        if not IsBool(c) then
            C := [ aut, s2^c ];
        fi;

        sigma := AutomorphismOfF2( aut!.freeGroup, [ "s", 1, 2, 3 ] );
        c     := AreConjugateAutomorphismsOfF2( sigma, aut );
        if not IsBool(c) then
            px := AutomorphismOfF2( aut!.freeGroup, [ "d", 2, 3, 2, 2, 3, 2 ] );
            C := [ aut, s2^c, px^c ];
        fi;

        sigma := AutomorphismOfF2( aut!.freeGroup, [ "s", 1, 2, 3, -1, 3 ] );
        py    := AutomorphismOfF2( aut!.freeGroup, [-1,3] );
        c     := AreConjugateAutomorphismsOfF2( sigma, aut );
        if not IsBool(c) then
            C := [ aut, (s2*py)^c ];
        fi;
    fi;

    return C;
end );