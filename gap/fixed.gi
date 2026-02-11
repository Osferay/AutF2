FixedSubgroupSA2 := function( aut )
    local cent, gens, word, fix, i, j, w, t;

    if not IsSpecialAutomorphismOfF2( aut ) then
        Error( "input has to be a special automorphism." );
    fi;

    cent := CentralizerAutomorphismOfF2( aut );
    gens := List( cent, MatrixRepresentationOfAutomorphismOfF2 );
    word := TrivialWordsSL2Z( gens );
    fix  := [];
    
    for i in [1..Length(word)] do
        w := word[i];
        t := aut^0;
        for j in [1..Length(w)] do
		    t := t*cent[ AbsInt(w[j]) ]^SignInt(w[j]);
        od;
        Add( fix, t );
	od;
    
    fix := List( fix, ConjugacyElementConjugacyAutomorphismOfF2 );
    fix := NielsenReducedSet( fix );

    return fix;
end;

InstallGlobalFunction( FixedSubgroupAutomorphismOfF2, function( aut )
    local fix2, s2, s4, s4y, d;

    if IsSpecialAutomorphismOfF2( aut ) then
        return FixedSubgroupSA2( aut );
    fi;

    if Order( aut ) = infinity then
        fix2 := FixedSubgroupSA2( aut^2 );

        if Length( fix2 ) = 2 then
            if ImageAutomorphismOfF2( aut, fix2[1] ) = fix2[1] then
                return [ fix2[1] ];
            elif ImageAutomorphismOfF2( aut, fix2[2] ) = fix2[2] then
                return [ fix2[2] ];
            else
                return [];
            fi;

        elif Length( fix2 ) = 1 then
            if ImageAutomorphismOfF2( aut, fix2[1] ) = fix2[1] then
                return [ fix2[1] ];
            else
                return [];
            fi;

        else
            return [];
        fi;
    else
        s2  := AutomorphismOfF2( aut!.freeGroup, ["s"] );
        s4  := AutomorphismOfF2( aut!.freeGroup, ["s",-1,2,3] );
        s4y := AutomorphismOfF2( aut!.freeGroup, ["s",-1,2,3,1,3] );

        if ConjugacyAutomorphismOfF2( aut, s2 ) or ConjugacyAutomorphismOfF2( aut, s4y ) then
            return [];
        else
            d := ConjugacyAutomorphismOfF2( aut, s4 );
            return [ ImageAutomorphismOfF2( d^-1, aut!.freeGroup.1 ) ];
        fi;
    fi;
end );       