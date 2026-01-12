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

InstallGlobalFunction( NielsenReducedSet, function( V )
    local   W, i, v, w, one, todo1, p, even, r, flag, todo2;

    W := ShallowCopy( V );

    one := W[1]^0;

    #Duplicate free
    W   := Set( W );
    #Ordered
    W   := SortedList( W, LexicographicOrderNSet );
    #Without one
    W   := Filtered( W, x -> x <> one );

    for i in [1..Length(W)] do
        if LexicographicOrderNSet( W[i]^-1, W[i]) then
            W[i] := W[i]^-1;
        fi;
    od;

    todo1 := ShallowCopy(W);
    while not IsEmpty( todo1 ) do
        v := todo1[1];
        Remove( todo1, 1 );
        p := Position( W, v );

        for i in [1..Length(W)] do
            if v <> W[i] then
                if v*W[i] = one or v^-1*W[i] = one or v*W[i]^-1 = one or v^-1*W[i]^-1 = one then
                    Remove( W, p );
                elif LexicographicOrderNSet( v*W[i], v ) then
                    W[p] := MinimumLexicographicOrderNSet( v*W[i], (v*W[i])^-1 );
                    Add( todo1, W[p] );
                elif LexicographicOrderNSet( v^-1*W[i], v ) then
                    W[p] := MinimumLexicographicOrderNSet( v^-1*W[i], W[i]^-1*v );
                    Add( todo1, W[p] );
                elif LexicographicOrderNSet( v*W[i]^-1, v ) then
                    W[p] := MinimumLexicographicOrderNSet( v*W[i]^-1, W[i]*v^-1 );
                    Add( todo1, W[p] );
                elif LexicographicOrderNSet( v^-1*W[i]^-1, v ) then
                    W[p] := MinimumLexicographicOrderNSet( v^-1*W[i]^-1, W[i]*v );
                    Add( todo1, W[p] );
                fi;
            fi;
        od;
    od;

    even := Filtered( W, x -> IsEvenInt(Length(x)) );
    flag := false;

    for i in [1..Length( even )] do

        v := even[i];
        r := LetterRepAssocWord( v );
        r := r{[1..Length(v)/2]};
        w := AssocWordByLetterRep( FamilyObj(v), r );
        w := w^-1;

        todo1 := Filtered( even, x -> x <> v and IsPrefix( x, w ) );
        todo2 := Filtered( even, x -> x <> v and IsSuffix( x, w^-1 ) );
        
        for i in [1..Length(todo1)] do
            p := Position( W, todo1[i] );
            W[p] := MinimumLexicographicOrderNSet( v*todo1[i], (v*todo1[i])^-1 );
            flag := true;
        od; 

        for i in [1..Length(todo2)] do
            p := Position( W, todo2[i] );
            W[p] := MinimumLexicographicOrderNSet( v*todo2[1]^-1, todo2[1]*v^-1 );
            flag := true;
        od;

        if flag then
            return NielsenReducedSet( W );
        else
            return W;
        fi;
    od;

end );

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