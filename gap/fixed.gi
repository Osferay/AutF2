FixedSubgroupSA2 := function( alpha )
    local cent, gens, word, fix, i, j, w, t;

    if not IsSpecialAutomorphismOfF2( alpha ) then
        Error( "input has to be a special automorphism." );
    fi;

    cent := CentralizerAutomorphismOfF2( alpha );
    gens := List( cent, MatrixRepresentationOfAutomorphismOfF2 );
    word := TrivialWordsSL2Z( gens );
    fix  := [];
    
    for i in [1..Length(word)] do
        w := word[i];
        t := alpha^0;
        for j in [1..Length(w)] do
		    t := t*cent[ AbsInt(w[j]) ]^SignInt(w[j]);
        od;
        Add( fix, t );
	od;
    
    fix := List( fix, ConjugacyElementConjugacyAutomorphismOfF2 );
    fix := NielsenReducedSet( fix );

    return fix;
end;