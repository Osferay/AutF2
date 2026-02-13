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
	local	gens, i,
			sigma, delta,
			phi1, phi2, phi3;

	gens  := GeneratorsOfGroup( F );

	sigma := function(w)
		local rep, r, new;
		rep := LetterRepAssocWord( w );
		new := [];
		for i in [1..Length( rep )] do
			r := rep[i];
			if AbsInt( r ) = 1 then Add( new, SignInt( r )*2 );
			else Add( new, SignInt( r )*1 );
			fi;
		od;
		return AssocWordByLetterRep( FamilyObj(w), new ); 
	end;
	phi1 := function( w )
		return EliminatedWord( w, gens[1], gens[1]*gens[2] );
	end;
	phi2 := function( w )
		return EliminatedWord( w, gens[2], gens[2]*gens[1] );
	end;
	phi3 := function( w )
		return EliminatedWord( w, gens[1], gens[2]^-1*gens[1] );
	end;
	delta := function( w )
		local new;

		new := ShallowCopy( w );
		new := EliminatedWord( new, gens[1], gens[1]^-1 );
		new := EliminatedWord( new, gens[2], gens[2]^-1 );
		new := EliminatedWord( new, gens[1], gens[2]*gens[1]*gens[2]^-1 );

		return new;
	end;

	return [ sigma, phi1, phi2, phi3, delta ];
end;

InverseOfGeneratorsOfGroupOfAutomorphismsOfF2 := function( F )
	local	gens,
			iphi1, iphi2, iphi3;

	gens  := GeneratorsOfGroup( F );

	iphi1 := function( w )
		return EliminatedWord( w, gens[1], gens[1]*gens[2]^-1 );
	end;
	iphi2 := function( w )
		return EliminatedWord( w, gens[2], gens[2]*gens[1]^-1 );
	end;
	iphi3 := function( w )
		return EliminatedWord( w, gens[1], gens[2]*gens[1] );
	end;

	return [ iphi1, iphi2, iphi3 ];
end;

ReduceWordOfAutomorphismOfF2 := function( word )
	local	red, i;

	red := ShallowCopy( word );

	for i in [1..Length( red )-1] do
		if IsInt( red[i] ) and IsInt( red[ i+1 ] ) then
			if red[i] = -red[i+1] then
				Remove( red, i+1 ); Remove( red, i );
				return ReduceWordOfAutomorphismOfF2( red );
			fi;
		elif red[i] = "s" and red[i+1] = "s" then
			Remove( red, i+1 ); Remove( red, i );
			return ReduceWordOfAutomorphismOfF2( red );
		elif red[i] = "d" and red[i+1] = "d" then
			Remove( red, i+1 ); Remove( red, i );
			return ReduceWordOfAutomorphismOfF2( red );
		fi;
	od;

	return red;
end;

MoveSigmaToLeft := function( word )
	local	pos, mov, tmp;

	mov := ReduceWordOfAutomorphismOfF2( word );
	pos := Positions( mov, "s" );

	if pos <> [1] and not IsEmpty( pos ) then 
		pos := Last( pos );
		if mov[ pos-1 ] = "d" then
			tmp := mov{[1..pos-2]};
			Append( tmp, [ "s", -1, 2, -1, -1, 2, 3] );
			Append( tmp, mov{ [pos+1..Length(mov) ] } );
			mov := ShallowCopy( tmp );
		elif AbsInt( mov[ pos-1 ] ) = 1 then
			mov[ pos ] := SignInt( mov[ pos-1 ] )*2;
			mov[ pos-1 ] := "s";
		elif AbsInt( mov[ pos-1 ] ) = 2 then
			mov[ pos ] := SignInt( mov[ pos-1 ] )*1;
			mov[ pos-1 ] := "s";
		elif AbsInt( mov[ pos-1 ] ) = 3 then
			tmp := mov{[1..pos-2]};
			Append( tmp, [ "s", -3, -2, SignInt( mov[ pos-1 ] )*1, 2, 3] );
			Append( tmp, mov{ [pos+1..Length(mov) ] } );
			mov := ShallowCopy( tmp );
		fi;
		
		mov := ReduceWordOfAutomorphismOfF2( mov );
		return MoveSigmaToLeft( mov );
	fi;

	return mov;
end;

WordOfSpecialAutomorphismOfF2ToBraidWord := function( word )
	local braid, i;

	braid := ShallowCopy( word );

	for i in [1..Length( braid )] do

		if braid[i] = "d" then
			braid[i] := 4;
		elif AbsInt( braid[i] ) = 1 then
			braid[i] := -1*braid[i];
		fi;
	
	od;
	
	return braid;
end;

BraidWordToWordOfSpecialAutomorphismOfF2 := function( braid )
	local word, i;

	word := ShallowCopy( braid );
	if word[1] = 0 then
		return [];
	fi;

	if word[1] = 4 then
		word[1] := "d";
	elif AbsInt( word[1] ) = 1 then
		word[1] := -1*word[1];
	fi;
	
	for i in [2..Length( word )] do

		if AbsInt( word[i] ) = 1 then
			word[i] := -1*word[i];
		fi;
	
	od;
	
	return word;
end;

LeftCanonicalFormAutomorphismOfF2 := function( word )
	local	lcf, i;

	if IsEmpty( word ) then
		return [];
	fi;

	if word[1] = "s" then
		lcf := word{[2..Length(word)]};
	else
		lcf := ShallowCopy( word );
	fi;
	lcf := WordOfSpecialAutomorphismOfF2ToBraidWord( lcf );
	
	AutF2WriteJSON( lcf );
	AutF2CallCpp( "lcf" );
	lcf := AutF2ReadJSON();
	
	lcf := BraidWordToWordOfSpecialAutomorphismOfF2( lcf );
	
	if word[1] = "s" then
		Add( lcf, "s", 1 );
	fi;

	return lcf;
end;

InstallMethod( AutomorphismOfF2, 
    "for a free group and a word of automorphisms", 
    [ IsFreeGroup, IsList ], 
    function( F, word ) 
        local aut, lcf, v, w, i, funs, gens, invs, identityfunction;

		lcf := MoveSigmaToLeft( word );
		lcf := LeftCanonicalFormAutomorphismOfF2( lcf );
        
        aut := rec( freeGroup := F, lcf := lcf );

        aut := Objectify( NewType( AutomorphismOfF2Family( FamilyObj( aut ) ), IsAutomorphismOfF2 and RepAutomorphismOfF2 ),
                       aut ) ;
        return( aut );
    end 
);

InstallMethod( WordOfAutomorphismOfF2, "for an automorphism of F2", [ IsAutomorphismOfF2 ],
	function( aut )
		return aut!.lcf;
end );

FunctionsAutomorphismOfF2 := function( aut )
	local F, autw, gens, invs, funs, i, identityfunction, v, w;

	F    := aut!.freeGroup;
	autw := aut!.lcf;		
	gens := GeneratorsOfGroupOfAutomorphismsOfF2( F );
	invs := InverseOfGeneratorsOfGroupOfAutomorphismsOfF2( F );
	funs := [];
	for i in [1..Length(autw)] do
		if autw[i] = 1 then
			Add( funs, gens[2] );
		elif autw[i] = 2 then
			Add( funs, gens[3] );
		elif autw[i] = 3 then
			Add( funs, gens[4] );
		elif autw[i] = -1 then
			Add( funs, invs[1] );
		elif autw[i] = -2 then
			Add( funs, invs[2] );
		elif autw[i] = -3 then
			Add( funs, invs[3] );
		elif autw[i] = "s" then
			Add( funs, gens[1] );
		elif autw[i] = "d" then
			Add( funs, gens[5] );
		else
			Error( "Word is not a word of automorphisms." );
		fi;
	od;

	if IsEmpty(funs) then
		identityfunction := function(w)
			return w;
		end;
		Add( funs, identityfunction );
	fi;

	return funs;
end;

InstallMethod( ImagesAutomorphismOfF2, "for an automorphism of F2", [ IsAutomorphismOfF2 ],
	function( aut )
		local F, v, w, funs, i;

		F    := aut!.freeGroup;
		funs := FunctionsAutomorphismOfF2( aut );

		v := GeneratorsOfGroup( F )[1];
		w := GeneratorsOfGroup( F )[2];
		for i in [1..Length(funs)] do
			v := funs[i]( v );
			w := funs[i]( w );
		od;

		return [ v, w ];
end );

InstallMethod( PrintObj, "for an automorphism of F2", [ IsAutomorphismOfF2 ],
    function( aut )
		Print( "Automorphism of F2 with word ", aut!.lcf );
    end 
);

InstallMethod( IsIdentityAutomorphismOfF2, 
	"for an automorphism of F2",
	[ IsAutomorphismOfF2 ],
	function( aut )
		return ( IsEmpty( aut!.lcf) );

end);

InstallMethod( MatrixRepresentationOfAutomorphismOfF2,
	"for an automorphism of F2",
	[ IsAutomorphismOfF2 ],
	function( aut )
		local M, lcf, i;

		M 	:= [[1,0],[0,1]];
		lcf := aut!.lcf;

		for i in [1..Length( lcf )] do
			if lcf[i] = "s"	then
				M := M*[[0,1],[1,0]];
			elif lcf[i] = "d" then
				M := M*[[-1,0],[0,-1]];
			elif AbsInt(lcf[i]) = 1 then
				M := M*[[ 1, SignInt( lcf[i] )*1 ],[ 0, 1 ]];
			elif AbsInt(lcf[i]) = 2 then
				M := M*[[ 1, 0 ],[ SignInt( lcf[i] )*1, 1 ]];
			elif AbsInt(lcf[i]) = 3 then
				M := M*[[ 1, -1*SignInt( lcf[i])*1 ],[ 0, 1 ]];
			fi;
		od;

		return M;
end );

InstallMethod( IsSpecialAutomorphismOfF2,
	"for an automorphism of F2",
	[ IsAutomorphismOfF2 ],
	function( aut )
		#local M;
		#M := MatrixRepresentationOfAutomorphismOfF2(aut);
		#return DeterminantIntMat(M) = 1;
		return WordOfAutomorphismOfF2( aut )[1] <> "s";
end);

InstallOtherMethod( \*,
	"for automorphisms of F2",
	[ IsAutomorphismOfF2, IsAutomorphismOfF2 ],
	function( aut1, aut2 )
		local w, aut;

		if aut1!.freeGroup <> aut2!.freeGroup then
			Error( "The free groups of the automorphisms are different." );
		fi;

		w := Concatenation( aut1!.lcf, aut2!.lcf );
		aut := AutomorphismOfF2( aut1!.freeGroup, w );

		return aut;	
end);

InstallOtherMethod( Inverse,
	"for automorphisms of F2",
	[ IsAutomorphismOfF2 ],
	function( aut )
		local 	word,
				inv,
				i;

		word := aut!.lcf;
		inv := [];

		for i in [1..Length( word )] do
			if IsInt(word[i]) then
				Add( inv, -word[i] );
			elif word[i] = "s" then
				Add( inv, "s" );
			else
				Add( inv, "d" );
			fi;
		od;
		
		return AutomorphismOfF2( aut!.freeGroup, Reversed(inv) );
end );

InstallOtherMethod( \^,
	"for automorphisms of F2",
	[ IsAutomorphismOfF2, IsInt ],
	function( aut, e)
		local i, pow, mul;
		
		if e = 0 then
			return AutomorphismOfF2( aut!.freeGroup, [] );
		fi;
		
		if e < 0 then
			pow := Inverse( aut );
			mul := Inverse( aut );
		else
			pow := aut;
			mul := aut;
		fi;

		for i in [2.. AbsInt(e)] do
			pow := pow * mul;
		od;
		return pow;
end);

ProductAutomorphismsOfF2ByWord := function( gens, word )
	local new, w, i, j, aut, inv, p;

	new := [];
	w   := List( gens, WordOfAutomorphismOfF2 );
	for i in [1..Length( word )] do
		p := ShallowCopy( w[ AbsInt(word[i]) ] );
		if word[i] > 0 then
			new := Concatenation( new, p );
		else
			for j in [1..Length( p )] do
				if IsInt( p[j] ) then
					p[j] := -1*p[j];
				fi;
			od;
			new := Concatenation( new, Reversed( p ) );
		fi;
	od;
	Error();
	aut := AutomorphismOfF2( gens[1]!.freeGroup, new );

	return aut;
end;

InstallOtherMethod( \^,
	"for two automorphisms of F2",
	[ IsAutomorphismOfF2, IsAutomorphismOfF2 ],
	function( a, b )
		local new;

		if a!.freeGroup <> b!.freeGroup then
			Error( "The free groups of the automorphisms are different." );
		fi;

		new := Inverse( b );
		new := ShallowCopy( WordOfAutomorphismOfF2( new ) );
		Append( new, WordOfAutomorphismOfF2( a ) );
		Append( new, WordOfAutomorphismOfF2( b ) );

		return AutomorphismOfF2( a!.freeGroup, new );
end);

InstallOtherMethod( Order,
	"for automorphisms of F2",
	[ IsAutomorphismOfF2 ],
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

InstallMethod( ImageByAutomorphismOfF2,
	"for an automorphism of F2 and an element of F2",
	[ IsAutomorphismOfF2, IsAssocWordWithInverse ],
	function( aut, w )
	local	funs, v, i;
	
	if not w in aut!.freeGroup then
		Error( "The given word has to be in the same free group as the automorphism." );
	fi;

	funs := FunctionsAutomorphismOfF2( aut );
	
	v := ShallowCopy( w );
	for i in [1..Length(funs)] do
		v := funs[i]( v );
	od;

	return v;
end);

InstallMethod( IsConjugacyAutomorphismOfF2,
	"for automorphisms of F2",
	[ IsAutomorphismOfF2 ],
	function(aut)
		local	F, gens, imgs, c1, c2;

		F    := aut!.freeGroup;
		gens := GeneratorsOfGroup( aut!.freeGroup );
		imgs := ImagesAutomorphismOfF2( aut );
		c1   := RepresentativeAction( F, gens[1], imgs[1] );
		c2   := RepresentativeAction( F, gens[2], imgs[2] );
		
		if not IsBool( c1 ) and not IsBool( c2 ) then
				if IsSuffix( c1, c2 ) then
					SetConjugacyElementConjugacyAutomorphismOfF2( aut, c1 );
					return true;
				elif IsSuffix( c2, c1 ) then
					SetConjugacyElementConjugacyAutomorphismOfF2( aut, c2 );
					return true;
				else
					return false;
				fi;
		fi;			

		return false;
	
end );

InstallMethod( ConjugacyElementConjugacyAutomorphismOfF2,
	"for a conjugacy automorphism in F2",
	[ IsAutomorphismOfF2 ],
	function( aut )
		if HasConjugacyElementConjugacyAutomorphismOfF2(aut) then
			return ConjugacyElementConjugacyAutomorphismOfF2( aut );
		elif IsConjugacyAutomorphismOfF2( aut ) then
			return ConjugacyElementConjugacyAutomorphismOfF2( aut );
		else
			return fail;
		fi;
end );

InstallMethod( ConjugacyAutomorphismOfF2,
	"for a conjugating element in F2",
	[ IsFreeGroup, IsAssocWordWithInverse ],
	function( F, word )
		local rep, r, aut, i;

		if Rank(F) <> 2 then
			Error( "This function is only for the free group of rank 2.");
		elif not word in F then
			Error( "The word has to be in the given free group." );
		fi;

		rep := LetterRepAssocWord( word );
		aut := [];
		for i in [1..Length(rep)] do
			r := rep[i];
			if AbsInt( r ) = 1 then
				Add( aut, SignInt( r )*2 );
				Add( aut, -3 );
				Add( aut, -2 );
				Add( aut, SignInt( r )*1 );
				Add( aut, 2 );
				Add( aut, 3 );
			elif AbsInt( r ) = 2 then
				Add( aut, SignInt( r )*1 );
				Add( aut, SignInt( r )*3 );
			fi;
		od;

		aut := AutomorphismOfF2( F, aut );
		SetIsConjugacyAutomorphismOfF2( aut, true );
		SetConjugacyElementConjugacyAutomorphismOfF2( aut, word );

		return aut;
end);

InstallOtherMethod( \=,
	"for two automorphisms of F2",
	[ IsAutomorphismOfF2, IsAutomorphismOfF2 ],
	function( a, b )

	return a!.lcf = b!.lcf;
end );

CentralizerAutomorphismOfF2 := function( a )
	local b, cent, i;

	if not IsSpecialAutomorphismOfF2( a ) then
		Error( "not yet." );
	fi;
	
	b := WordOfAutomorphismOfF2( a );
	b := WordOfSpecialAutomorphismOfF2ToBraidWord( b );

	AutF2WriteJSON( b );
	AutF2CallCpp( "cent" );
	cent := AutF2ReadJSON( );

	for i in [1..Length(cent)] do
		cent[i] := BraidWordToWordOfSpecialAutomorphismOfF2( cent[i] );
		cent[i] := AutomorphismOfF2( a!.freeGroup, cent[i] );
	od; 

	return cent;
end;