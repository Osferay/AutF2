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

WordSL2ZinST := function( M )

    local A, S, T, C, N, s, c;

    if DeterminantIntMat(M) <> 1 then
        Error( "M has to be in SL2Z" );
    fi;

    S := [[0,-1],[1,0]];
    T := [[1,1],[0,1]];

    if M = [[1,0],[0,1]] then
        return rec( s := 0, T := [] );
    fi;

    A := ShallowCopy( M );

    if A[2][1] = 0 and A[1][1] = 1 then
        return rec( s := 0, T := [ A[1][2] ] );
    elif A[2][1] = 0 and A[1][1] = -1 then
        return rec( s := 2, T := [ -1*( A[1][2] ) ] );
    fi;

    C := MinusContinuedFraction( A[1][1], A[2][1] );

    N := [[1,0],[0,1]];
    for c in C do
        N := N*T^c*S;
    od;
    N := N^-1*A;

    if N[1][1] = -1 then
        s := 2;
        Add( C, -1* N[1][2] );
    else
        s := 0;
        Add( C, N[1][2] );
    fi; 
    
    return rec( s := s, T := C );

end;

MatrixSL2ZbyWordST := function( w )

    local M, i, T, S;

    M := [[1,0],[0,1]];
    S := [[0,-1],[1,0]];
    T := [[1,1],[0,1]];

    if w.s = 2 then
        M := M*S^2;
    fi;

    for i in [1..Length(w.T)] do
        M := M*T^(w.T[i])*S;
    od;

    return M*S^-1;

end;

CanonicalWordInSU := function( w, s )
    local p, q;

    p := Position( w, 1 );
    q := Position( w, 2 );
    if p = 1 and IsBool(q) then
        return[ [], s + Length(w) ];
    elif p = 1 and not IsBool(q) then
        while w[p] = 1 do
            p := p + 1;
        od;
        return CanonicalWordInSU( w{[p..Length(w)]}, s+p-1 );
    fi;

    p := PositionSublist( w, [1,1,1,1] );
    if not IsBool( p ) then
        return CanonicalWordInSU( Concatenation( w{[1..p-1]}, w{[p+4..Length(w)]}), s );
    fi;

    p := PositionSublist( w, [2,2,2] );
    if not IsBool( p ) then
        return CanonicalWordInSU( Concatenation( w{[1..p-1]}, w{[p+3..Length(w)]}), s+2 );
    fi;

    p := PositionSublist( w, [1,1] );
    if not IsBool( p ) then
        return CanonicalWordInSU( Concatenation( w{[1..p-1]}, w{[p+2..Length(w)]}), s+2 );
    fi;

    return [ w, s mod 4 ];

end;

WordSL2ZinSU := function( M )

    local w, v, i, t;

    if DeterminantIntMat(M) <> 1 then
        Error( "M has to be in SL2Z" );
    fi;

    w := WordSL2ZinST( M );
    v := [];

    for t in w.T do
        if not IsEvenInt( t ) then
            v := Concatenation( v, [1,1] );
        fi;

        if t > 0 then
            for i in [1..t] do
                v := Concatenation( v, [1,2] );
            od;
        else
            for i in [1..-t] do
                v := Concatenation( v, [2,2,1] );
            od;
        fi;
        Add( v, 1 );
    od;

    if not IsEmpty( v ) then
        Remove( v, Length(v) );
    fi;

    v := CanonicalWordInSU( v, w.s );
    
    return rec( w := v[1], s := v[2] );
end;

MatrixSL2ZbyWordSU := function( w )

    local M, v, U, S;

    M := [[1,0],[0,1]];
    S := [[0,-1],[1,0]];
    U := [[0,-1],[1,1]];

    M := M*S^(w.s);

    for v in w.w do
        if v = 1 then
            M := M*S;
        else
            M := M*U;
        fi;
    od;

    return M;

end;

WordGL2ZinSU := function( M )

    local A, e, w;

    if M = [[1,0],[0,1]] then
        return rec( e := 0, s := 0, U := [] );
    fi;
    
    A := M;

    if DeterminantIntMat(A) = -1 then
        A := [[0,1],[1,0]]*A;
        e := 1;
    else 
        e := 0;
    fi;

    w := WordSL2ZinSU( A );

    return rec( e := e, s := w.s, w := w.w );

end;

MatrixGL2ZbyWordSU := function( w )

    local M, v, N;

    if w.e = 1 then
        M := [[0,1],[1,0]];
    else
        M := [[1,0],[0,1]];
    fi;
    
    v := rec( s := w.s, w := w.w );
    N := MatrixSL2ZbyWordSU( v );

    return M*N;

end;

LengthWordMatrixSL2ZinSU := function( M )
    local w;

    if DeterminantIntMat(M) <> 1 then
        Error( "M has to be in SL2Z" );
    fi;

    w := WordSL2ZinSU( M );

    return Length( w.w ) + w.s;
end;

HomomorphismSL2Zto12Z := function( M )
    local w, eU, eS, G, u, s;

    if DeterminantIntMat(M) <> 1 then
        Error( "M has to be in SL2Z" );
    fi;

    w := WordSL2ZinSU( M );

    eU := Length( Positions( w.w, 2 ) );
    eS := Length( Positions( w.w, 1 ) ) + w.s;
    
    return (2*eU + 3*eS) mod 12;

end;

MembershipCommutatorSL2Z := function( M )

    local   w, x, y, N, I, F, new, z;

    if DeterminantIntMat(M) <> 1 then
        Error( "M has to be in SL2Z" );
    fi;

    if HomomorphismSL2Zto12Z( M ) <> 0 then
        return false;
    fi;

    #[S,U] and [S^-1,U^-1]
    x   := [ [ 1, 1 ], [ 1, 2 ] ];
    y   := [ [ 2, 1 ], [ 1, 1 ] ];
    F   := [ [ x,1 ], [x^-1,-1], [y,2], [y^-1,-2] ];

    N   := ShallowCopy(M);
    new := [];
    I   := [[1,0],[0,1]];
    
    while N <> I do
        z := List( F, x -> LengthWordMatrixSL2ZinSU( N*x[1] ) );
        z := Position( z, Minimum(z) );
        Add( new, F[z][2] );
        N := N*F[z][1];
    od;

    return -1*Reversed(new);

end;

MatrixCommutatorSL2ZbyWordxy := function( w )

    local x, y, M, i;

    x := [ [ 1, 1 ], [ 1, 2 ] ];
    y := [ [ 2, 1 ], [ 1, 1 ] ];
    M := [[1,0],[0,1]];

    for i in [1..Length(w)] do
        if AbsInt(w[i]) = 1 then
            M := M*x^(SignInt(w[i]));
        else
            M := M*y^(SignInt(w[i]));
        fi;
    od;

    return M;
end;

TransversalRepresentativeCommutatorSL2Z := function( T, M )

    local t;

    if DeterminantIntMat(M) <> 1 then
        Error( "M has to be in SL2Z" );
    fi;

    t := HomomorphismSL2Zto12Z( M );

    return T[t+1];
end;

OrbitStabilizerMembership := function( T, t, gens )
    local I, S, w, dict, orbit, stab, o, i, y, j, todo, new, sdic, tmp;

    I     := [[1,0],[0,1]];
    dict  := NewDictionary( I, true );
    sdic  := NewDictionary( I, true );
    AddDictionary( dict, t, 1 );
    orbit := [ [ t, I, [] ] ];
    todo  := [ [ t, I, [] ] ];
    stab  := [];
    new   := [];

    while not IsEmpty(todo) do
        o := todo[1];
        Remove( todo, 1 );
        for i in [1..Length(gens[1])] do
            #Compute a new point
            y := TransversalRepresentativeCommutatorSL2Z( T, gens[1][i]*o[1] );
            j := LookupDictionary( dict, y );

            if IsBool(j) then
                AddDictionary( dict, y, Length(orbit)+1 );
                Add( orbit, [y, gens[1][i]*o[2], Concatenation( gens[2][i], o[3] ) ] );
                Add( todo, [y, gens[1][i]*o[2], Concatenation( gens[2][i], o[3] ) ] );
            else
                tmp := (orbit[j][2])^-1*(gens[1][i]*o[2]);
                if IsBool( LookupDictionary( sdic, tmp ) ) then
                    AddDictionary( sdic, tmp, Length(stab)+1 );
                    Add( stab, tmp );
                    Add( new, Concatenation( -1*Reversed( orbit[j][3] ), Concatenation( gens[2][i], o[3] ) ) );
                fi;
            fi;
        od;
    od;

    new := List( new, ReduceLetterRep );
    return [stab,new];
end;

GeneratorsOfIntersectionCommutatorSL2Z := function( T, gens, F )

    local S, t, b, i, j, U, A, w;

    S := [ gens, [] ];

    for i in [1..Length(gens)] do
        Add(S[2], [i]);
    od;
    
    for t in T do
        S := OrbitStabilizerMembership( T, t, S );
    od;
    
    w := S[2];
    S := List( S[1], MembershipCommutatorSL2Z );
    S := List( S, s -> AssocWordByLetterRep( FamilyObj( One(F) ), s ) );
    # We compute U = < <gens> \cap G' >

    b := [];
    for i in [1..Length(S)] do
        Add( b, [ i ] );
    od;
    
    U := NielsenReducedSetBacktrack( S, b ); 
    
    if IsEmpty(U[1]) then
        return U;
    fi;

    A := [];
    for i in [1..Length( U[2] ) ] do
        S := [];
        for j in [1..Length( U[2][i] ) ] do
            b := U[2][i][j];
            if b > 0 then
                S := Concatenation( S, w[b] );
            else
                S := Concatenation( S, -1*Reversed( w[-1*b] ) );
            fi;
        od;
        Add( A, ReduceLetterRep(S) );
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
    i := 1;
    while IsBool( h ) do 

        t := T[i];
        j := v[i];
        if j <> 0 then
            h := I;
            y := t;

            while j <> -1 do
                h := gens[j]*h;
                Add( new, j, 1 );
                y := TransversalRepresentativeCommutatorSL2Z( T, gens[j]^-1*y );
                i := LookupDictionary( S, y );
                j := v[i];
            od;
        fi;
        
        i := i+1;
    od;
    
    return [ t, h^-1, -1*Reversed(new) ];
end;

CosetRepresentativeSubgroupSL2Z := function( gens, M )

    local   S, U, F, T, I, wU, t, wh, h, u, wv;

    if ForAny( gens, x -> DeterminantIntMat(x) <> 1) or DeterminantIntMat(M) <> 1 then
        Error( "gens and M have to be in SL2Z" );
    fi;

    I  := [[1,0],[0,1]];
    S  := [[0,-1],[1,0]];
    U  := [[0,-1],[1,1]];
    F  := FreeGroup( 2 );

    #This is a transversal of the commutator for SL2Z
    T  := [ I, S^3*U^2, U, S, U^2, S*U, S^2, S*U^2, S^2*U, S^3, S^2*U^2, S^3*U ];

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
    u  := MatrixCommutatorSL2ZbyWordxy( u );
    
    return [ u*t, wU, wv, wh ];

end;

MembershipSubgroupSL2Z := function( gens, M )

    local r, w, i;

    r := CosetRepresentativeSubgroupSL2Z( gens, M );

    if r[1] = r[1]^0 then
        w := r[4];
        for i in [1..Length( r[3] ) ] do
            if r[3][i] > 0 then
                w := Concatenation( w, r[2][ r[3][i] ] );
            else 
                w := Concatenation( w, -1*Reversed( r[2][ -1*r[3][i] ] ) );
            fi;
        od;
        return ReduceLetterRep( w );
    else
        return false;
    fi;

end;

ReduceGeneratorSetSubgroupSL2Z := function( gens )
    local new, i, tmp;

    if ForAny( gens, x -> DeterminantIntMat(x) <> 1) then
        Error( "gens have to be in SL2Z" );
    fi;

    if Length(gens) = 1 then
        return gens;
    fi;

    new := ShallowCopy( gens );

    for i in [1..Length(new)] do
        tmp := Concatenation( new{[1..(i-1)]}, new{[(i+1).. Length(new)]} );
        if not IsBool( MembershipSubgroupSL2Z( tmp, new[i] ) ) then
            return ReduceGeneratorSetSubgroupSL2Z( tmp );
        fi;
    od;

    return new;
end;

ReduceParallelGeneratorSetSubgroupSL2Z := function( gens, w )
    local new, i, tmp, v;

    if Length(gens) = 1 then
        return [ gens, w ];
    fi;
    
    new := ShallowCopy( gens );
    v   := ShallowCopy( w );

    for i in [1..Length(new)] do
        tmp := Concatenation( new{[1..(i-1)]}, new{[(i+1).. Length(new)]} );
        if not IsBool( MembershipSubgroupSL2Z( tmp, new[i] ) ) then
            Remove( v, i );
            return ReduceParallelGeneratorSetSubgroupSL2Z( tmp, v );
        fi;
    od;

    return [ new, v ];
end;

MatrixbyWordInGens := function( gens, w )
    local v, M;

    M := [[1,0],[0,1]];
    for v in w do
        M := M*gens[ AbsInt(v) ]^SignInt(v);
    od;

    return M;
end;