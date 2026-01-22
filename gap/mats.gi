NearestInteger := function( n )

    if IsInt( n ) then
        return n;
    elif n>0 then
        return Int( n+1/2 );
    elif n<0 then
        return Int( n-1/2 );
    fi;
end;

FloorRat := function( q )

    if IsInt(q) then
        return q;
    fi;

    if not IsRat( q ) then
        Error( "FloorRat is only for rational numbers." );
    fi;

    if q > 0 then
        return Int(q);
    else
        return Int(q)-1;
    fi;
end;

CeilRat := function( q )

    if IsInt(q) then
        return q;
    fi;

    if not IsRat( q ) then
        Error( "CeilRat is only for rational numbers." );
    fi;

    if q > 0 then
        return Int(q)+1;
    else
        return Int(q);
    fi;
end;

ReduceElitpicBQF := function( QF, P )
    local A, B, C, n, Q;

    A := QF[1];
    B := QF[2];
    C := QF[3];
    
    if A>C or ( A=C and -A <= B and B < 0 ) then
        Q := P*[[0,-1],[1,0]];
        return ReduceElitpicBQF( [ C, -B, A ], Q );
        
    elif A < AbsInt(B) then
        n := Int( B/(2*A) );
        Q := P*[[1,-n],[0,1]];
        return ReduceElitpicBQF( [ A, -2*A*n+B, A*n^2-B*n+C ], Q );
        
    else
        return rec( QF := QF, P := P );
    fi;
    
end;

ConjugacyClassElipticMatrix := function( A )

    local I, a, d, P, R, QF;

    I := [[1,0],[0,1]];

    if A = I or A = -I then
        return rec( rep := A, P := I );

    elif Order(A) = 2 then
        a := Eigenvectors( Rationals, A )[1];
        a := Product( List( a, DenominatorRat ) )*a;
        d := GcdInt( a[1], a[2] );
        if d <> 1 then
            a := a / d;
        fi;
        d := Gcdex( a[1], a[2] );
        P := [ [ a[1], a[2]], [ -d.coeff2, d.coeff1] ];
        P := P^-1;
        R := P^-1*A*P;
        if R[2][1] = 0 then

        elif IsEvenInt( R[2][1] ) then
            P := P*[[1,0], [ Int(R[2][1]/2), 1 ]];
            R := P^-1*A*P;
        else
            P := P*[[1,0], [ Int((R[2][1]-1)/2), 1 ]];
            R := P^-1*A*P;
        fi;
        
        return rec( rep := R, P := P );

    else
        if A[2][1] < 0 then 
            P := [[1,0],[0,-1]];
            R := P*A*P;
        else
            P := I;
            R := A;
        fi;
        
        QF := [ R[2][1], R[2][2]-R[1][1], -R[1][2] ];
        QF := ReduceElitpicBQF( QF, P );
        P  := P*QF.P;

        return  rec( rep:= QF.QF, P := P );
    fi;
                
end;

ConjugacyClassParabolicMatrix := function( A )

    local a, d, P, R;

    a := Eigenvectors( Rationals, A )[1];
    a := Product( List( a, DenominatorRat ) )*a;
    d := GcdInt( a[1], a[2] );
    if d <> 1 then
        a := a / d;
    fi;
    d := Gcdex( a[1], a[2] );
    P := [ [ a[1], a[2]], [ -d.coeff2, d.coeff1] ];
    P := P^-1;
    R := P^-1*A*P;

    if R[2][1] < 0 then
        P := P*[[1,0],[0,-1]];
        R := P^-1*A*P;
    fi;

    return rec( rep := R, P := P );

end;

ConjugacyClassHyperbolicMatrix := function( A )

    local R, P, Q, n;

    P := [[1,0],[0,1]];
    R := ShallowCopy( A );

    #Step 1
    if R[1][1] < R[2][2] then
        Q := [[0,1],[1,0]];
        P := P*Q;
        R := Q^-1*R*Q;
    fi;
    
    #Step 2
    if AbsInt( R[1][2] ) <= AbsInt(R[1][1]) and R[1][2]<0 then
        Q := [[1,0],[0,-1]];
        P := P*Q;
        R := Q^-1*R*Q;
    
    elif AbsInt( R[2][1] ) <= AbsInt(R[1][1]) and R[2][1]<0 then
        Q := [[1,0],[0,-1]];
        P := P*Q;
        R := Q^-1*R*Q;
    fi;
    
    #Step 3
    if 0 < R[1][2] and R[1][2] <= R[1][1] and R[2][2] < 0 then
        n := FloorRat( R[2][2]/R[1][2] );
        Q := [[1,0],[n,1]];
        P := P*Q;
        R := Q^-1*R*Q;

    elif 0 < R[2][1] and R[2][1] <= R[1][1] and R[2][2] < 0 then
        n := CeilRat( -R[2][2]/R[2][1] );
        Q := [[1,n],[0,1]];
        P := P*Q;
        R := Q^-1*R*Q;

    fi;

    #Step 4 and 5
    if R[1][2] < R[2][2] then
        n := FloorRat( (R[2][1]-R[2][2])/R[2][1] );
        Q := [[1,n],[0,1]];
        P := P*Q;
        R := Q^-1*R*Q;
    elif R[1][2] < R[2][2] then
        n := CeilRat( (R[2][2]-R[1][2])/R[1][2] );
        Q := [[1,0],[n,1]];
        P := P*Q;
        R := Q^-1*R*Q;
    fi;

    return rec( rep := R, P := P );
end;

ContinuedFraction := function( a, b )

    local x, y, r, Q;

    if a<b then 
        Error( "a should be greater than b" );
    fi;

    x := a; 
    y := b; 
    r := 1; 
    Q := [];

    while r <> 0 do
        Add( Q, QuoInt( x, y ) );
        r := RemInt( x, y );
        x := y;
        y := r;
    od;

    return Q;

end;

MinusContinuedFraction := function( a, b )

    local x, y, r, Q;

    x := a; 
    y := b; 
    r := 1; 
    Q := [];

    while r <> 0 do
        r := CeilRat( x/y );
        Add( Q, r );
        r := r*y-x;
        x := y;
        y := r;
    od;

    return Q;

end;

IsCyclicalyPermuted := function( c1, c2 )
    local n, p, perm, i;

    if Length( c1 ) <> Length( c2 ) then
        return false;
    fi;

    n := Length( c1 );
    p := [1..n];
    for i in [1..n] do
        perm := PermList( p );
        if Permuted( c2, perm ) = c1 then
            return p;
        else
            Add( p, Last(p), 1 );
            Remove( p, Length(p) );
        fi;
    od;

    return false;
end;

ConjugacyStandardHyperbolicMatrix := function( A, B )

    local c1, c2, p, pos, P, i;

    c1 := ContinuedFraction( A[1][1], A[2][1] );
    c2 := ContinuedFraction( B[1][1], B[2][1] );
    p  := IsCyclicalyPermuted( c1, c2 );

    if IsBool( p ) then
        return p;
    else
        pos := Position( p, 1 );
        P   := [[1,0],[0,1]];
        p   := p{[1..pos-1]};
        for i in c1{p} do
            P := P*[[i,1],[1,0]];
        od;
        return P^-1;
    fi;


end;

InstallGlobalFunction( ConjugacyGL2Z, function( A, B )

    local   t1, t2, d1, d2, o1, o2, R1, R2, C;

    t1 := Trace( A );
    t2 := Trace( B );
    d1 := DeterminantIntMat( A );
    d2 := DeterminantIntMat( B );
    o1 := Order(A);
    o2 := Order(B);

    if t1 <> t2 or d1 <> d2 or o1 <> o2 then
        return false;
    fi;

    if IsInt( o1 ) then
        #Elliptic case
        R1 := ConjugacyClassElipticMatrix( A );
        R2 := ConjugacyClassElipticMatrix( B );

        if R1.rep <> R2.rep then
            return false;
        else
            return (R1.P)*(R2.P)^-1;
        fi;
    elif AbsInt( t1 ) = 2 and d1 = 1 then
        #Parabolic case
        R1 := ConjugacyClassParabolicMatrix( A );
        R2 := ConjugacyClassParabolicMatrix( B );

        if R1.rep <> R2.rep then
            return false;
        else
            return (R1.P)*(R2.P)^-1;
        fi;
        
    else
        #Hyperbolic case
        R1 := ConjugacyClassHyperbolicMatrix( A );
        R2 := ConjugacyClassHyperbolicMatrix( B );

        C := ConjugacyStandardHyperbolicMatrix( R1.rep, R2.rep );

        if IsBool( C ) then 
            return C;
        else
            return (R1.P)*C*(R2.P)^-1;
        fi;
    fi;

end);

CentralizerElipticMatrix := function( A )
    local I;

    I := [[1,0],[0,1]];

    if A = I or A = -I then
        return [ [[0,-1],[1,1]], [[0,-1],[1,0]], [[0,1],[1,0]] ];
    else
        return rec( gen := A, exponent := 1 );
    fi;

end;

CentralizerParabolicMatrix := function( A )
    local   R, e, P, T;

    R := ConjugacyClassParabolicMatrix( A );
    e := R.rep[2][1];
    P := R.P;
    T := R.rep;
    T[2][1] := 1;

    return rec( gen := P*T*P^-1, exponent := e );
end;

DiscriminantMatrix22 := function( A )
    return Trace(A)^2-4*DeterminantIntMat(A);
end;

fundamentalUnit := function( delta )
        local a, a0, b, d, b0, q1, q2, A, u, v, tmp;

        a0 := 1;
        d  := RootInt( delta );
        if (delta - d) mod 2 = 0 then
            b0 := d;
        else
            b0 := d-1;
        fi;

        a  := a0; b  := b0;
        q1 := 0;  q2 := 1;

        repeat
            A   := Int( (b+d)/(2*a) );
            b   := 2*a*A-b;
            a   := (delta-b^2)/(4*a);
            tmp := ShallowCopy(q1);
            q1  := q1*A+q2;
            q2  := tmp;
        until a = a0 and b = b0;

        u := ( 2*a0*q2 + q1*b0 )/a0;
        v := q1/a0;

        return [u,v];
end;

CentralizerHyperbolicMatrix := function( A )
    local   t, delta, u, v, ad, b, c, C, e, U, V, Uu, Vv;

    t     := Trace(A);
    delta := DiscriminantMatrix22( A );
    u     := fundamentalUnit( delta );
    v     := u[2];
    u     := u[1];

    ad    := A[1][1]-A[2][2];
    b     := A[1][2];
    c     := A[2][1];

    C     := [ [ 1/2*(u+ad*v), b*v ], [ c*v, 1/2*(u-ad*v) ] ];

    e     := 1;
    U     := [ u ];
    V     := [ v ];
    while Last(U) <> t do
        Uu := u*Last(U)+v*Last(V)*delta;
        Vv := u*Last(V)+v*Last(U);
        Add( U, Uu );
        Add( V, Vv );
    od;

    return rec( gen := C, exponent := e );

end;

InstallGlobalFunction( CentralizerGL2Z, function( A )
    local o, t, d;

    o := Order( A );
    t := Trace( A );
    d := DeterminantIntMat( A );

    if IsInt( o ) then
        #Elliptic case
        return CentralizerElipticMatrix( A );

    elif AbsInt( t ) = 2 and d = 1 then
        #Parabolic case
        return CentralizerParabolicMatrix( A );
        
    else
        #Hyperbolic case
        return CentralizerHyperbolicMatrix( A );
    fi;

end );

MembershipCommutatorSL2Z := function( M )

    local   r,l,d, N, I, A, B, w, e;

    if M[1][1] mod 4 <> 1 or M[2][2] mod 4 <> 1 or M[1][2] mod 2 <> 0 or M[2][1] mod 2 <> 0 then
        return false;
    fi;

    l := function(M)
        return M[1][1]^2+M[2][1]^2;
    end;

    r := function(M)
        return M[2][2]^2+M[1][2]^2;
    end;

    d := function(M)
        return M[1][1]*M[1][2]+M[2][1]*M[2][2];
    end;

    N := M;
    I := M^0;
    A := [[1,2],[0,1]];
    B := [[1,0],[2,1]];
    w := [];

    while N <> I do
    
        if l(N) < AbsInt( d(N) ) then
            e := NearestInteger( -d(N)/(2*l(N)) );
            Append( w, ListWithIdenticalEntries( AbsInt(e), SignInt(e)*1 ) );
            N := N*A^e;
            
        elif r(N) < AbsInt( d(N) ) then
            e := NearestInteger( -d(N)/(2*r(N)) );
            Append( w, ListWithIdenticalEntries( AbsInt(e), SignInt(e)*2 ) );
            N := N*B^e;
        fi;
    od;
    
    return -1*Reversed(w);

end;

MatrixCommutatorSL2ZbyWordAB := function( w )

    local A, B, M, i;

    A := [[1,2],[0,1]];
    B := [[1,0],[2,1]];
    M := [[1,0],[0,1]];

    for i in [1..Length(w)] do
        if AbsInt(w[i]) = 1 then
            M := M*A^(SignInt(w[i]));
        else
            M := M*B^(SignInt(w[i]));
        fi;
    od;

    return M;
end;

TransversalRepresentativeCommutatorSL2Z := function( T, M )

    local t, N;

    for t in T do
        N := t^-1*M;
        if N[1][1] mod 4 = 1 and N[2][2] mod 4 = 1 and N[1][2] mod 2 = 0 and N[2][1] mod 2 = 0 then
            return t;
        fi;
    od;
end;

GeneratorsOfIntersectionCommutatorSL2Z := function( T, gens, F )

    local I, S, w, t, dict, orbit, stab, exps, s, y, j, e, A, b, U;

    I  := [[1,0],[0,1]];
    S  := ShallowCopy( gens );
    w  := ListWithIdenticalEntries( Length(gens), 1 );
    
    for t in T do
        #For each element in the transversal we compute the orbit-stabilizer
        dict  := NewDictionary( I, true );
        AddDictionary( dict, t, 1 );
        orbit := [ t ];
        stab  := [];
        exps  := [];

        for s in S do
            #Compute a new point
            y := TransversalRepresentativeCommutatorSL2Z( T, t*s );
            j := LookupDictionary( dict, y );
        
            e := 1;
            while IsBool( j ) do
                AddDictionary( dict, y, Length( orbit )+1 );
                Add( orbit, y ); 
                
                #If this is a new point, compute the block
                y := TransversalRepresentativeCommutatorSL2Z( T, s*y );
                e := e + 1;
                j := LookupDictionary( dict, y );
            od;
            
            Add( stab, s^e );
            Add( exps, e );
            
        od;
        S := stab;
        w := List( [1..Length(exps) ], i -> w[i]*exps[i] );
    od;
    
    S := List( S, MembershipCommutatorSL2Z );
    S := List( S, s -> AssocWordByLetterRep( FamilyObj( One(F) ), s ) );
    # We compute U = < <gens> \cap G' >
    b := [];
    for s in [1..Length(gens)] do
        Add( b, [ s ] );
    od;
    U := NielsenReducedSetBacktrack( S, b ); 

    A := [];
    for s in [1..Length( U[2] ) ] do
        S := [];
        for t in [1..Length( U[2][s] ) ] do
            S := Concatenation( S, ListWithIdenticalEntries( w[ U[2][s][t] ], U[2][s][t] ) );
        od;
        Add( A, S );
    od;
    
    return [ U[1], A ];

end;

SchreierVectorTransversalCommutatorSL2Z := function( T, gens, M )

    local I, t, dict, S, i, v, orbit, todo, j, o, y, h, new;

    I     := [[1,0],[0,1]];
    t     := TransversalRepresentativeCommutatorSL2Z( T, M );
    new   := [];

    dict  := NewDictionary( I, true );
    AddDictionary( dict, t, 1 );
    S     := NewDictionary( I, true );
    for i in [ 1.. 12 ] do
        AddDictionary( S, T[i], i );
    od;    
    v := ListWithIdenticalEntries( 12, 0 );

    orbit := [ t ];
    todo  := [ t ];
    j     := LookupDictionary( S, t );
    v[j]  := -1;

    #Compute the orbit of t
    while not IsEmpty( todo ) do

        o := todo[1];
        Remove( todo, 1 );

        for i in [1..Length(gens)] do

            y := TransversalRepresentativeCommutatorSL2Z( T, gens[i]*o );
            j := LookupDictionary( dict, y );
            
            if IsBool( j ) then
                AddDictionary( dict, y, Length( orbit )+1 );
                Add( orbit, y );
                Add( todo, y );

                j    := LookupDictionary( S, y );
                v[j] := i;
            fi;
        od;
    od;

    h := false;

    #Use the Schreier vector to obtain h
    while IsBool( h ) do 
        i := 1;

        t := T[i];
        j := v[i];
        if j <> 0 then
            
            h := I;
            y := t;

            while j <> -1 do
                h := gens[j]*h;
                Add( new, j );
                y := TransversalRepresentativeCommutatorSL2Z( T, gens[j]*y );
                i := LookupDictionary( S, y );
                j := v[i];
            od;
        fi;

        i := i+1;
    od;
    
    return [ t, h, new ];
end;

CosetRepresentativeSubgroupSL2Z := function( gens, M )

    local   S, U, F, T, I, wU, t, wh, h, u, wv;

    I  := [[1,0],[0,1]];
    S  := [[0,-1],[1,0]];
    U  := [[0,-1],[1,1]];
    F  := FreeGroup( 2 );

    #This is a transversal of the commutator for SL2Z
    T  := [ I, S, S^2, S^3, U, U^-1, S*U, S*U^-1, S^2*U, S^2*U^-1, S^3*U, S^3*U^-1 ];

    U  := GeneratorsOfIntersectionCommutatorSL2Z( T, gens, F );
    wU := U[2];
    U  := U[1];

    t  := SchreierVectorTransversalCommutatorSL2Z( T, gens, M );
    wh := t[3];
    h  := t[2];
    t  := t[1];

    u  := MembershipCommutatorSL2Z( h^-1*M*t^-1 );
    u  := AssocWordByLetterRep( FamilyObj( One(F) ), u );
    u  := CosetRepresentativeReducedNielsenSetBacktrack( U, u );
    wv := u[3];
    u  := u[2];

    u  := LetterRepAssocWord( u );
    u  := MatrixCommutatorSL2ZbyWordAB( u );
    
    return [ u*t, wU, wv, wh ];

end;

MembershipSubgroupSL2Z := function( gens, M )

    local r;

    r := CosetRepresentativeSubgroupSL2Z( gens, M );

    if r[1] = r[1]^0 then
        return true;
    else
        return false;
    fi;

end;

WordGL2ZinST := function( M )

    local A, S, T, C, N, c, w, d;

    A := M;
    S := [[0,-1],[1,0]];
    T := [[1,1],[0,1]];

    if DeterminantIntMat(A) = -1 then
        A := [[0,1],[1,0]]*A;
        d := -1;
    else 
        d := 1;
    fi;

    if A[2][1] = 0 and A[1][1] = 1 then
        return rec( det := d, wS := [0], wT := [ A[1][2] ] );
    elif A[2][1] = 0 and A[1][1] = -1 then
        return rec( det := d, wS := [2], wT := [ -A[1][2] ] );
    fi;

    C := MinusContinuedFraction( A[1][1], A[2][1] );

    N := [[1,0],[0,1]];
    for c in C do
        N := N*T^c*S;
    od;
    
    S := ListWithIdenticalEntries( Length(C), 1 );
    N := N^-1*A;
    Add( C, AbsInt(N[1][2]) );

    if N[1][1] = -1 then
        Add( S, 2 );
    else
        Add( S, 0 );
    fi; 
    
    return rec( det := d, wT := C, wS := S );

end;

MatrixGL2ZbyWordST := function( w )

    local M, i, T, S;

    M := [[1,0],[0,1]];
    S := [[0,-1],[1,0]];
    T := [[1,1],[0,1]];

    if w.det = -1 then
        M := M*[[0,1],[1,0]];
    fi;

    for i in [1..Length(w.wT)] do
        M := M*T^(w.wT[i])*S^(w.wS[i]);
    od;

    return M;

end;
