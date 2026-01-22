InstallMethod( IsCyclicalyReducedWord,
	"for a associative word with inverse",
	[ IsAssocWordWithInverse ],
	function( u )
    local   rep;

    if IsOne( u ) then 
        return true;
    fi;

    rep := LetterRepAssocWord(u);

    return rep[1] <> -1*Last(rep);
end);

InstallMethod( IsSubword,
    "for a pair of words in a free group",
    [ IsAssocWordWithInverse, IsAssocWordWithInverse],
    function( u, v )
        local   l1, l2;

        l1 := LetterRepAssocWord( u );
        l2 := LetterRepAssocWord( v );

        return PositionSublist( l1, l2 );
    
end );

InstallMethod( IsPrefix,
    "for a pair of words in a free group",
    [ IsAssocWordWithInverse, IsAssocWordWithInverse],
    function( u, v )
        local   l1, l2;

        l1 := LetterRepAssocWord( u );
        l2 := LetterRepAssocWord( v );

        return PositionSublist( l1, l2 ) = 1;
    
end );

InstallMethod( IsSuffix,
    "for a pair of words in a free group",
    [ IsAssocWordWithInverse, IsAssocWordWithInverse],
    function( u, v )
        local   l1, l2;

        l1 := LetterRepAssocWord( u );
        l2 := LetterRepAssocWord( v );

        return PositionSublist( Reversed(l1), Reversed(l2) ) = 1;
    
end );

LexicographicOrderNSet := function( v, w )
    local l1, l2, r1, r2, i;

    if v = w then
        return false;
    fi;

    l1 := Length( v );
    l2 := Length( w );

    if l1 = l2 then 
        r1 := LetterRepAssocWord( v );
        r2 := LetterRepAssocWord( w );
        for i in [1..l1] do
            if r1[i] <> r2[i] then
                if r1[i] = -1*r2[i] then
                    if SignInt( r1[i] ) = 1 then
                        return true;
                    else
                        return false;
                    fi;
                elif AbsInt(r1[i]) < AbsInt(r2[i]) then
                    return true;
                else
                    return false;
                fi;
            fi;
        od;

    elif l1 < l2 then
        return true;
    
    else
        return false;
    fi;
end;

MinimumLexicographicOrderNSet := function( v, w )
    if LexicographicOrderNSet( v, w ) then
        return v;
    else
        return w;
    fi;
end;

MinimumLexicographicOrderNSetBacktrack := function( v, w )
    if LexicographicOrderNSet( v, w ) then
        return [ v, 1 ];
    else
        return [ w, 2 ];
    fi;
end;

NielsenReducedSetStep2 := function( V )

    local   one, W, i, j;

    one := V[1]^0;
    W   := ShallowCopy( V );

    for i in [1..Length(W)] do
        for j in [1..Length(W)] do
            if i <> j then
                if W[i]*W[j] = one or W[i]^-1*W[j] = one or W[i]*W[j]^-1 = one or W[i]^-1*W[j]^-1 = one then
                    Remove( W, i );
                    return W;
                elif LexicographicOrderNSet( W[i]*W[j], W[i] ) then
                    W[i] := MinimumLexicographicOrderNSet( W[i]*W[j], (W[i]*W[j])^-1 );
                    return W;
                elif LexicographicOrderNSet( W[i]^-1*W[j], W[i] ) then
                    W[i] := MinimumLexicographicOrderNSet( W[i]^-1*W[j], W[j]^-1*W[i] );
                    return W;
                elif LexicographicOrderNSet( W[i]*W[j]^-1, W[i] ) then
                    W[i] := MinimumLexicographicOrderNSet( W[i]*W[j]^-1, W[j]*W[i]^-1 );
                    return W;
                elif LexicographicOrderNSet( W[i]^-1*W[j]^-1, W[i] ) then
                    W[i] := MinimumLexicographicOrderNSet( W[i]^-1*W[j]^-1, W[j]*W[i] );
                    return W;
                fi;
            fi;
        od;
    od;

    return false;
end;

NielsenReducedSetStep3 := function( V )

    local   W, even, i, v, p, q, n, flag;

    W    := ShallowCopy( V );
    even := Filtered( W, x -> IsEvenInt(Length(x)) );

    for i in [1..Length( even )] do

        v := even[i];
        q := LetterRepAssocWord( v );
        q := q{[Length(v)/2+1..Length(v)]};
        q := AssocWordByLetterRep( FamilyObj(v), q );
        q := q^-1;

        p := PositionProperty( even, x -> x <> v and IsPrefix( x, q ) );
        
        if not IsBool( p ) then
            n    := Position( W, even[p] );
            W[n] := MinimumLexicographicOrderNSet( v*even[p], (v*even[p])^-1 );
            
            flag := true;
            while flag do
                flag := NielsenReducedSetStep2( W );
        
                if not IsBool( flag ) then
                    W := flag;
                    flag := true;
                fi;
            od;

            return W;
        fi; 

        p := PositionProperty( even, x -> x <> v and IsSuffix( x, q^-1 ) );

        if not IsBool( p ) then
            n    := Position( W, even[p] );
            W[n] := MinimumLexicographicOrderNSet( v*even[p]^-1, even[p]*v^-1 );
            
            flag := true;
            while flag do
                flag := NielsenReducedSetStep2( W );
        
                if not IsBool( flag ) then
                    W := flag;
                    flag := true;
                fi;
            od;

            return W;
        fi; 

    od;

    return false;
end;

InstallGlobalFunction( NielsenReducedSet, function( V )
    local   W, i, flag;

    W := ShallowCopy( V );

    #Duplicate free
    W   := Set( W );
    #Ordered
    W   := SortedList( W, LexicographicOrderNSet );
    #Without one
    if W[1] = W[1]^0 then
        Remove( W, 1 );
    fi;

    for i in [1..Length(W)] do
        if LexicographicOrderNSet( W[i]^-1, W[i]) then
            W[i] := W[i]^-1;
        fi;
    od;

    flag := true;
    while flag do
        flag := NielsenReducedSetStep2( W );
        
        if not IsBool( flag ) then
            W    := flag;
            flag := true;
        fi;
    od;

    flag := true;
    while flag do
        flag := NielsenReducedSetStep3( W );
        
        if not IsBool( flag ) then
            W    := flag;
            flag := true;
        fi;
    od;

    return W;
end );

NielsenReducedSetStep2Backtrack := function( V, words )

    local   one, new, W, i, j;

    one := V[1]^0;
    W   := ShallowCopy( V );
    new := ShallowCopy( words );

    for i in [1..Length(W)] do
        for j in [1..Length(W)] do
            if i <> j then
                if W[i]*W[j] = one or W[i]^-1*W[j] = one or W[i]*W[j]^-1 = one or W[i]^-1*W[j]^-1 = one then
                    Remove( W, i );
                    Remove( new, i );
                    return [ W, new ];
                elif LexicographicOrderNSet( W[i]*W[j], W[i] ) then
                    W[i] := MinimumLexicographicOrderNSetBacktrack( W[i]*W[j], (W[i]*W[j])^-1 )[1];
                    if MinimumLexicographicOrderNSetBacktrack( W[i]*W[j], (W[i]*W[j])^-1 )[2] = 1 then
                        new[i] := Concatenation( new[i], new[j] );
                    else
                        new[i] := -1*Reversed( Concatenation( new[i], new[j] ) );
                    fi;

                    return [ W, new ];
                elif LexicographicOrderNSet( W[i]^-1*W[j], W[i] ) then
                    W[i] := MinimumLexicographicOrderNSetBacktrack( W[i]^-1*W[j], W[j]^-1*W[i] )[1];
                    if MinimumLexicographicOrderNSetBacktrack( W[i]^-1*W[j], W[j]^-1*W[i] )[2] = 1 then
                        new[i] := Concatenation( -1*Reversed( new[i] ), new[j] );
                    else
                        new[i] := Concatenation( -1*Reversed( new[j] ), new[i] );
                    fi;

                    return [ W, new ];
                elif LexicographicOrderNSet( W[i]*W[j]^-1, W[i] ) then
                    W[i] := MinimumLexicographicOrderNSetBacktrack( W[i]*W[j]^-1, W[j]*W[i]^-1 )[1];
                    if MinimumLexicographicOrderNSetBacktrack( W[i]*W[j]^-1, W[j]*W[i]^-1 )[2] = 1 then
                        new[i] := Concatenation( new[i], -1*Reversed( new[j] ) );
                    else
                        new[i] := Concatenation( new[j], -1*Reversed( new[i] ) );
                    fi;

                    return [ W, new ];
                elif LexicographicOrderNSet( W[i]^-1*W[j]^-1, W[i] ) then
                    W[i] := MinimumLexicographicOrderNSetBacktrack( W[i]^-1*W[j]^-1, W[j]*W[i] );
                    if MinimumLexicographicOrderNSetBacktrack( W[i]^-1*W[j]^-1, W[j]*W[i] )[2] = 1 then
                        new[i] := -1*Reversed( Concatenation( new[j], new[i] ) );
                    else
                        new[i] := Concatenation( new[j], new[i] );
                    fi;

                    return [ W, new ];
                fi;
            fi;
        od;
    od;

    return false;
end;

NielsenReducedSetStep3Backtrack := function( V, words )

    local   W, new, even, ewor, i, v, y, p, q, n, flag;

    W    := ShallowCopy( V );
    new  := ShallowCopy( words );
    even := Filtered( W, x -> IsEvenInt(Length(x)) );
    ewor := Filtered( [1..Length(W) ], i -> IsEvenInt(Length(W[i])) );
    ewor := words{ ewor };

    for i in [1..Length( even )] do

        v := even[i];
        y := ewor[i];
        q := LetterRepAssocWord( v );
        q := q{[Length(v)/2+1..Length(v)]};
        q := AssocWordByLetterRep( FamilyObj(v), q );
        q := q^-1;

        p := PositionProperty( even, x -> x <> v and IsPrefix( x, q ) );
        
        if not IsBool( p ) then
            n    := Position( W, even[p] );
            W[n] := MinimumLexicographicOrderNSetBacktrack( v*even[p], (v*even[p])^-1 )[1];

            if MinimumLexicographicOrderNSetBacktrack( v*even[p], (v*even[p])^-1 )[2] = 1 then
                new[n] := Concatenation( y, new[n] );
            else
                new[n] := -1*Reversed( Concatenation( y, new[n] ) );
            fi;
            
            flag := true;
            while flag do
                flag := NielsenReducedSetStep2Backtrack( W, new );
        
                if not IsBool( flag ) then
                    W    := flag[1];
                    new  := flag[2];
                    flag := true;
                fi;
            od;

            return [ W, new ];
        fi; 

        p := PositionProperty( even, x -> x <> v and IsSuffix( x, q^-1 ) );

        if not IsBool( p ) then
            n    := Position( W, even[p] );
            W[n] := MinimumLexicographicOrderNSet( v*even[p]^-1, even[p]*v^-1 );
            
            if MinimumLexicographicOrderNSetBacktrack( v*even[p]^-1, even[p]*v^-1 )[2] = 1 then
                new[n] := Concatenation( y, -1*Reversed( new[n] ) );
            else
                new[n] := Concatenation( new[n], -1*Reversed( y ) );
            fi;

            flag := true;
            while flag do
                flag := NielsenReducedSetStep2Backtrack( W, new );
        
                if not IsBool( flag ) then
                    W    := flag[1];
                    new  := flag[2];
                    flag := true;
                fi;
            od;

            return [ W, new ];
        fi; 

    od;

    return false;
end;

NielsenReducedSetBacktrack := function( V, words )
    local   W, i, v, w, todo1, p, even, r, flag, todo2, new, n;

    W   := ShallowCopy( V );
    new := ShallowCopy( words );

    #Duplicate free
    i := 1;
    while not IsDuplicateFree( W ) do

        p := Position( W, W[i], i);
        while not IsBool( p ) do
            Remove( W, p );
            Remove( new, p );
            p := Position( W, W[i], i);
        od;

        i := i+1;
        
    od;
    #Ordered
    SortParallel( W, new, LexicographicOrderNSet );
    
    #Without one
    if W[1] = W[1]^0 then
        Remove( W, 1 );
        Remove( new, 1 );
    fi;
    
    for i in [1..Length(W)] do
        if LexicographicOrderNSet( W[i]^-1, W[i]) then
            W[i] := W[i]^-1;
            new[i] := -1*Reversed( new[i] );
        fi;
    od;

    flag := true;
    while flag do
        flag := NielsenReducedSetStep2Backtrack( W, new );
        
        if not IsBool( flag ) then
            W    := flag[1];
            new  := flag[2];
            flag := true;
        fi;
    od;

    flag := true;
    while flag do
        flag := NielsenReducedSetStep3Backtrack( W, new );
        
        if not IsBool( flag ) then
            W    := flag[1];
            new  := flag[2];
            flag := true;
        fi;
    od;

    return [ W, new ];

end;

InstallGlobalFunction( IsNielsenReducedSet, function( V )
    local one, i, j, k;

    if ForAny( V, x -> x^-1 in V ) then
        return false;
    fi;

    one := V[1]^0;

    for i in [1..Length(V)] do
        for j in [1..Length(V)] do
            if V[i]*V[j] <> one and Length( V[i]*V[j] ) < Length( V[i] ) then
                return false;
            fi;
        od;
    od;

    for i in [1..Length(V)] do
        for j in [1..Length(V)] do
            for k in [1..Length(V)] do
                if V[i]*V[j] <> one and V[j]*V[k] <> one and Length( V[i]*V[j]*V[k] ) < Length( V[i] ) - Length( V[j] ) + Length( V[k] ) then 
                    return false;
                fi;
            od;
        od;
    od;

    return true;
end );

CosetRepresentativeReducedNielsenSet := function( V, w )
    local v, u, W, flag, i, GP, half, p;

    v    := w^0;
    u    := ShallowCopy( w );
    W    := ShallowCopy( V );
    W    := Concatenation( W, List( V, x -> x^-1 ) );
    flag := true;
    half := Filtered( V, x -> IsEvenInt( Length( V ) ) );

    while flag do
        flag := false;
        for i in [1..Length(W)] do
            if Length( W[i]*u ) < Length( u ) then
                v := v*W[i]^-1;
                u := W[i]*u;
                flag := true;
            fi;
        od;

        for i in [1..Length(V)] do
            if Length( V[i]*u ) = Length( u ) then
                v := v*W[i]^-1;
                u := W[i]*u;
                flag := true;
            fi; 
        od;
    od;

    return [v,u];
end;

CosetRepresentativeReducedNielsenSetBacktrack := function( V, w )
    local v, u, W, flag, i, GP, new, word, p;

    v    := w^0;
    u    := ShallowCopy( w );
    W    := ShallowCopy( V );
    W    := Concatenation( W, List( V, x -> x^-1 ) );

    word := Concatenation( [1..Length(V)], -1*[1..Length(V)] );
    new  := [];

    flag := true;
    while flag do
        flag := false;
        for i in [1..Length(W)] do
            if Length( W[i]*u ) < Length( u ) then
                v := v*W[i]^-1;
                u := W[i]*u;
                flag := true;

                Add( new, -1*word[i] );
            fi;
        od;

        for i in [1..Length(V)] do
            if Length( V[i]*u ) = Length( u ) then
                v := v*W[i]^-1;
                u := W[i]*u;
                flag := true;

                Add( new, -i );
            fi; 
        od;
    od;
    
    return [v, u, new];
end;

InstallGlobalFunction( InFreeSubgroup, function( V, w )
    local W,r;

    W := NielsenReducedSet( V );
    r := CosetRepresentativeReducedNielsenSet( W, w );

    if r[2] = w^0 then
        return true;
    else
        return false;
    fi;
end );

InstallGlobalFunction( IsSubgroupOfFreeSubgroup, function( V, W )

    return ForAll( W, w -> InFreeSubgroup( V, w ) );

end );

InstallGlobalFunction( IsNormalSubgroupOfFreeSubgroup, function( V, W )

    return ForAll( W, w -> ForAll( V, v -> InFreeSubgroup( V, v^-1*w*v ) ) ) and IsSubgroupOfFreeSubgroup( V, W );

end );

InstallGlobalFunction( AreEquivalent, function( V, W ) 

    return IsSubgroupOfFreeSubgroup( V, W ) and IsSubgroupOfFreeSubgroup( W, V );

end );