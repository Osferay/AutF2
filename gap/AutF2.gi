#
# AutF2: Computations for the automorphisms group of F2
#
# Implementations
#
InstallMethod(AutomorphismOfF2Family,
    "for a family",
    [ IsFamily ],
    function( Fam )
        local F;

        F:= NewFamily( "AutomorphismOfF2Family", IsAutomorphismOfF2 );

        if HasCharacteristic( Fam ) then
            SetCharacteristic( F, Characteristic( Fam ) );
        fi;

        return F;
    end 
);

GeneratorsOfGroupOfAutomorphismsOfF2 := function( F )
	local	autos, gens,
			AutF2IdentityFunction,
			ReplaceGen1ByInverse,
			ReplaceGen1ByGen1Gen2,
			ReplaceGen2ByGen2Gen1,
			ReplaceGen1ByInverseGen2Gen1,
			ReplaceGen2ByInverseGen1Gen2,
			ReplaceGen1ByGen1InverseGen2,
			ReplaceGen2ByGen2InverseGen1,
			ReplaceGen1ByGen2Gen1,
			ReplaceGen2ByGen1Gen2;

	autos := [];
	gens  := GeneratorsOfGroup( F );

	AutF2IdentityFunction := function( w )
		return w;
	end;
	Add( autos, AutF2IdentityFunction );
	#Word 1
	ReplaceGen1ByInverse := function( w )
		return EliminatedWord( w, gens[1], gens[1]^-1 );
	end;
	Add( autos, ReplaceGen1ByInverse );
	#Word 2
	ReplaceGen1ByGen1Gen2 := function( w )
		return EliminatedWord( w, gens[1], gens[1]*gens[2] );
	end;
	Add( autos, ReplaceGen1ByGen1Gen2 );
	#Word 3
	ReplaceGen2ByGen2Gen1 := function( w )
		return EliminatedWord( w, gens[2], gens[2]*gens[1] );
	end;
	Add( autos, ReplaceGen2ByGen2Gen1 );
	#Word 4
	ReplaceGen1ByInverseGen2Gen1 := function( w )
		return EliminatedWord( w, gens[1], gens[2]^-1*gens[1] );
	end;
	Add( autos, ReplaceGen1ByInverseGen2Gen1 );
	#Word 5
	ReplaceGen2ByInverseGen1Gen2 := function( w )
		return EliminatedWord( w, gens[2], gens[1]^-1*gens[2] );
	end;
	Add( autos, ReplaceGen2ByInverseGen1Gen2 );
	#Word 6
	ReplaceGen1ByGen1InverseGen2 := function( w )
		return EliminatedWord( w, gens[1], gens[1]*gens[2]^-1 );
	end;
	Add( autos, ReplaceGen1ByGen1InverseGen2 );
	#Word 7
	ReplaceGen2ByGen2InverseGen1 := function( w )
		return EliminatedWord( w, gens[2], gens[2]*gens[1]^-1 );
	end;
	Add( autos, ReplaceGen2ByGen2InverseGen1 );
	#Word 8
	ReplaceGen1ByGen2Gen1 := function( w )
		return EliminatedWord( w, gens[1], gens[2]*gens[1] );
	end;
	Add( autos, ReplaceGen1ByGen2Gen1 );
	#Word 9
	ReplaceGen2ByGen1Gen2 := function( w )
		return EliminatedWord( w, gens[2], gens[1]*gens[2] );
	end;
	Add( autos, ReplaceGen2ByGen1Gen2 );

	return autos;
end;

InstallMethod( AutomorphismOfF2, 
    "for a list of functions", 
    [ IsFreeGroup, IsList ], 
    function( F, word ) 
        local aut, v, w, i;

		

		v := GeneratorsOfGroup( F )[1];
		w := GeneratorsOfGroup( F )[2];
		for i in [1..Length(list[1])] do
			v := list[1][i]( v );
		od;
		for i in [1..Length(list[2])] do
			w := list[2][i]( w );
		od;
        
        aut := rec( functions := list, freeGroup := F, images := [v,w], word := word );

        aut := Objectify( NewType( AutomorphismOfF2Family( FamilyObj( aut ) ), IsAutomorphismOfF2 and RepAutomorphismOfF2 ),
                       aut ) ;
        return( aut );
    end 
);

InstallMethod( PrintObj, "for an automorphism of F2", [ IsAutomorphismOfF2 ],
    function( aut )
		local gens;

		gens := GeneratorsOfGroup( aut!.freeGroup );

		Print( gens[1], " -> ", aut!.images[1], "\n" );
		Print( gens[2], " -> ", aut!.images[2] );
    end 
);

InstallMethod( IsIdentityAutomorphismOfF2, 
	"for an automorphism of F2",
	[IsAutomorphismOfF2],
	function( aut )
		local gens, v, w;

		gens := GeneratorsOfGroup( aut!.freeGroup );
		v    := aut!.images[1];
		w    := aut!.images[2];
		return ( gens[1] = v and gens[2] = w );

end);

InstallMethod( MatrixRepresentationOfAutomorphismOfF2,
	"for an automorphism of F2",
	[IsAutomorphismOfF2],
	function( aut )
		local M, gens, v, w;

		M := [[0,0],[0,0]];
		gens := GeneratorsOfGroup( aut!.freeGroup );
		v    := aut!.images[1];
		w    := aut!.images[2];

		M[1][1] := ExponentSumWord( v, gens[1] );
		M[1][2] := ExponentSumWord( v, gens[2] );
		M[2][1] := ExponentSumWord( w, gens[1] );
		M[2][2] := ExponentSumWord( w, gens[2] );

		return M;
end );

InstallMethod( IsSpecialAutomorphismOfF2,
	"for an automorphism of F2",
	[IsAutomorphismOfF2],
	function( aut )
		local M;
		M := MatrixRepresentationOfAutomorphismOfF2(aut);
		return DeterminantIntMat(M) = 1;
end);

InstallOtherMethod( \*,
	"for automorphisms of F2",
	[ IsAutomorphismOfF2, IsAutomorphismOfF2],
	function( aut1, aut2 )
		local l1, l2, w1, w2, aut;

		if aut1!.freeGroup <> aut2!.freeGroup then
			Error( "The free groups of the automorphisms are different." );
		fi;

		l1 := Concatenation( aut1!.functions[1], aut2!.functions[1] );
		l2 := Concatenation( aut1!.functions[2], aut2!.functions[2] );
		w1 := Concatenation( aut1!.word[1], aut2!.word[1] );
		w2 := Concatenation( aut1!.word[2], aut2!.word[2] );
		aut := AutomorphismOfF2( aut1!.freeGroup, [l1,l2] );
		aut!.word := [w1, w2];

		return aut;	
end);

InstallOtherMethod( Inverse,
	"for automorphisms of F2",
	[ IsAutomorphismOfF2 ],
	function( aut )
		local 	autos,
				w1, w2, i, inv;

		if IsEmpty(aut!.word) then
			Error( "this automorphism does not belong an automorphism group." );
		fi;

		w1 := []; w2 := [];

		for i in [1..Length( aut!.word[1])/2 ] do
			if aut!.word[1][2*i-1] = 0 then
				w1[2*i-1] 	:= [0,0];
			elif aut!.word[1][2*i-1] = 1 then 
				w1[2*i-1] 	:= [1,1];
			elif aut!.word[1][2*i-1] = 2 and aut!.word[1][2*i] = 1 then 
				w1[2*i-1] 	:= [2,-1];
			elif aut!.word[1][2*i-1] = 3 and aut!.word[1][2*i] = 1 then 
				w1[2*i-1] 	:= [3,-1];
			elif aut!.word[1][2*i-1] = 4 and aut!.word[1][2*i] = 1 then 
				w1[2*i-1] 	:= [4,-1];
			elif aut!.word[1][2*i-1] = 5 and aut!.word[1][2*i] = 1 then 
				w1[2*i-1] 	:= [5,-1];
			elif aut!.word[1][2*i-1] = 6 and aut!.word[1][2*i] = -1 then 
				w1[2*i-1] 	:= [2,1];
			elif aut!.word[1][2*i-1] = 7 and aut!.word[1][2*i] = -1 then  
				w1[2*i-1] 	:= [3,1];
			elif aut!.word[1][2*i-1] = 8 and aut!.word[1][2*i] = -1 then  
				w1[2*i-1] 	:= [4,1];
			elif aut!.word[1][2*i-1] = 9 and aut!.word[1][2*i] = -1 then  
				w1[2*i-1] 	:= [5,1];
			fi;
		od;
		for i in [1..Length( aut!.word[2])/2 ] do
			if aut!.word[2][2*i-1] = 0 then
				w2[2*i-1] 	:= [0,0];
			elif aut!.word[2][2*i-1] = 1 then 
				w2[2*i-1] 	:= [1,1];
			elif aut!.word[2][2*i-1] = 2 and aut!.word[1][2*i] = 1 then  
				w2[2*i-1] 	:= [2,-1];
			elif aut!.word[2][2*i-1] = 3 and aut!.word[1][2*i] = 1 then  
				w2[2*i-1] 	:= [3,-1];
			elif aut!.word[2][2*i-1] = 4 and aut!.word[1][2*i] = 1 then  
				w2[2*i-1] 	:= [4,-1];
			elif aut!.word[2][2*i-1] = 5 and aut!.word[1][2*i] = 1 then  
				w2[2*i-1] 	:= [5,-1];
			elif aut!.word[2][2*i-1] = 2 and aut!.word[1][2*i] = -1 then  
				w2[2*i-1] 	:= [2,1];
			elif aut!.word[2][2*i-1] = 3 and aut!.word[1][2*i] = -1 then  
				w2[2*i-1] 	:= [3,1];
			elif aut!.word[2][2*i-1] = 4 and aut!.word[1][2*i] = -1 then  
				w2[2*i-1] 	:= [4,1];
			elif aut!.word[2][2*i-1] = 5 and aut!.word[1][2*i] = -1 then  
				w2[2*i-1] 	:= [5,1];
			fi;
		od;

		w1 := Flat( Reversed( w1 ) ); w2 := Flat( Reversed( w2 ) );

		inv := [ w1, w2 ];

		return inv;
end );

InstallOtherMethod( \^,
	"for automorphisms of F2",
	[ IsAutomorphismOfF2, IsInt],
	function( aut, e)
		local i, pow, mul;
		
		if e < 0 then
			pow := Inverse( aut );
			mul := Inverse( aut );
		else
			pow := aut;
			mul := aut;
		fi;

		for i in [2..e] do
			pow := pow * mul;
		od;
		return pow;
end);

InstallOtherMethod( Order,
	"for automorphisms of F2",
	[ IsAutomorphismOfF2],
	function( aut )
		local i, o;
		
		o := [1,2,3,4];
		for i in o do
			if IsIdentityAutomorphismOfF2(aut^i) then
				return i;
			fi;
		od;
		return infinity;
end);

GroupOfAutomorphismsOfF2 := function( F )

	local 	autos,
			e, t, f1, f2, f3, f4;

	autos := AllAutomorphismsOfF2( F );
	e := AutomorphismOfF2( F, [[autos[1]], [autos[1]]] );
	e!.word := [ [0,0], [0,0] ];
	t := AutomorphismOfF2( F, [[autos[2]], [autos[1]]] );
	t!.word := [ [1,1], [0,0] ];
	f1 := AutomorphismOfF2( F, [[autos[3]], [autos[1]]] );
	f1!.word := [ [2,1], [0,0] ];
	f2 := AutomorphismOfF2( F, [[autos[1]], [autos[4]]] );
	f2!.word := [ [0,0], [3,1] ];
	f3 := AutomorphismOfF2( F, [[autos[5]], [autos[1]]] );
	f3!.word := [ [4,1], [0,0] ];
	f4 := AutomorphismOfF2( F, [[autos[1]], [autos[6]]] );
	f4!.word := [ [0,0], [5,1] ];

	return [e,t,f1,f2,f3,f4];
end;