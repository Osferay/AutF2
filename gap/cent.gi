OrbitStabilizerCentralizer := function( aut, C )
    local dict, sdict, orbit, stab, o, y, j, todo, c, i, tmp;

    i     := AutomorphismOfF2( aut!.freeGroup, [] );
    dict  := NewDictionary( [], true );
    sdict := NewDictionary( [], true );
    AddDictionary( dict, aut!.lcf, 1 );
    orbit := [ [aut,i] ];
    stab  := [ ];
    todo  := [ [aut,i] ];

    while not IsEmpty(todo) do
        o := todo[1];
        Remove( todo, 1 );
        for c in C do
            y := o[1]^c;
            j := LookupDictionary( dict, y!.lcf );

            if IsBool(j) then
                AddDictionary( dict, y!.lcf, Length(orbit)+1 );
                Add( orbit, [ y, o[2]*c ] );
                Add( todo, [ y, o[2]*c ] );
            else
                tmp := (orbit[j][2])*(o[2]*c);
                if IsBool( LookupDictionary( sdict, tmp!.lcf ) ) then
                    AddDictionary( sdict, tmp!.lcf, Length(stab)+1 );
                    Add( stab, tmp );
                fi;
            fi;
        od;
    od;

    return [orbit,stab];
end;

InstallGlobalFunction( "CentralizerAutomorphismOfF2", function( aut )
    local C, C2, sigma, c, s2, px;

    if IsSpecialAutomorphismOfF2( aut ) then
        C     := CentralizerAutomorphismOfF2InSA( aut );
        sigma := AutomorphismOfF2( aut!.freeGroup, ["s"] );
        c     := AreConjugateAutomorphismsOfF2( aut, aut^sigma);
        if not IsBool(c) then
            Add( C, sigma*c );
        fi;
    
    elif Order( aut ) = infinity then
        C2 := CentralizerAutomorphismOfF2InSA( aut^2 );
    else
        s2    := AutomorphismOfF2( aut!.freeGroup, [ -1, 2, 3, -1, 2, 3 ] );
        sigma := AutomorphismOfF2( aut!.freeGroup, ["s"] );
        c     := AreConjugateAutomorphismsOfF2( aut, sigma );
        if not IsBool(c) then
            C := [ sigma, s2 ];
            C := List( C, x -> x^c );
        fi;

        sigma := AutomorphismOfF2( aut!.freeGroup, [ "s", -1, 2, 3 ] );
        c     := AreConjugateAutomorphismsOfF2( aut, sigma );
        if not IsBool(c) then
            px := AutomorphismOfF2( aut!.freeGroup, [ "d", 2, 3, 2, 2, 3, 2 ] );
            C := [ sigma, s2, px ];
            C := List( C, x -> x^c );
        fi;

        sigma := AutomorphismOfF2( aut!.freeGroup, [ "s", -1, 2, 3, 1, 3 ] );
        c     := AreConjugateAutomorphismsOfF2( aut, sigma );
        if not IsBool(c) then
            C := [ sigma, s2 ];
            C := List( C, x -> x^c );
        fi;
    fi;

    return C;
end );