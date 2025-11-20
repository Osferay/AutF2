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