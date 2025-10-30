DiscriminantMatrix22 := function( M )
    return Trace(M)^2-4*DeterminantIntMat(M);
end;

InstallGlobalFunction( IsReducedMatrix22, function( M )
    local D;

    D := DiscriminantMatrix22(M);
    if D < 0 then
        if AbsInt(M[2][2]-M[1][1]) = M[2][1] or M[2][1] = - M[1][2] then
            return (M[2][2]-M[1][1])>=0;
        else
            return AbsInt(M[2][2]-M[1][1]) <= M[2][1] and M[2][1] <= - M[1][2];
        fi;
    elif D>0 and not IsSquareInt( D ) then  
        return M[2][1]>0 and AbsoluteValue( Sqrt( Float(D) ) - 2*M[2][1] ) < Float( M[2][2]-M[1][1] ) and Float( M[2][2]-M[1][1] ) < Sqrt( Float( D ) );
    else
        if M[1][2] <> 0 then
            return false;
        elif M[1][1] = M[2][2] then
            return M[2][1] > 0;
        else
            return M[2][1] < M[1][1] - M[2][2];
        fi;
    fi;
end ); 

InstallGlobalFunction( ReduceMatrix22, function( M )
    local   ReduceMatrixPositiveDefinite,
            ReduceMatrixIndeterminate,
            ReduceMatrixSquare,
            D;
    
    ReduceMatrixPositiveDefinite := function( M )
        local   P, red, n;

        if IsReducedMatrix22( M ) then
            return rec( reduced := M, conj := [[1,0],[0,1]] );

        elif M[2][1] < 0 then
            P := [[1,0],[0,-1]];
            red := ReduceMatrixPositiveDefinite( P^-1*M*P );
            return rec( reduced := red.reduced, conj := red.conj*P );

        elif M[2][1]> -M[1][2] or (M[2][1]=M[1][2] and (-M[2][1] <= M[2][2]-M[1][1] and M[2][2]-M[1][1] < 0)) then
            P := [[0,-1],[1,0]];
            red := ReduceMatrixPositiveDefinite( P^-1*M*P );
            return rec( reduced := red.reduced, conj := red.conj*P );

        elif not IsReducedMatrix22(M) then
            n := (M[2][2]-M[1][1]-1)/(2*M[2][1]);
            n := Int( Ceil( Float(n) ) );
            P := [[1,-n],[0,1]];
            red := ReduceMatrixPositiveDefinite( P^-1*M*P );
            return rec( reduced := red.reduced, conj := red.conj*P );
        fi;
    end;

    ReduceMatrixIndeterminate := function( M )
        local cyc, red, n, D, P, R, conj, pos;

        D := DiscriminantMatrix22( M );
        red := function( b, a )
            local r;
            r := b mod (2*a);

            if Float( AbsInt(a) ) > Sqrt( Float( D ) ) then
                while r > AbsInt( a ) do
                    r := r-2*AbsInt(a);
                od;
                return r;
            else
                while Float(r) > Sqrt( Float( D ) ) do
                    r := r - 2*AbsInt(a);
                od;
                while Sqrt( Float( D ) ) - 2*AbsInt( a ) > Float( r ) do
                    r := r + 2*AbsInt(a);
                od;
                return r;
            fi;
        end;

        cyc := [];
        R   := M;
        conj := [ [[1,0],[0,1]] ];

        repeat 
            Add( cyc, R );
            n   := ( red( R[1][1]-R[2][2], R[1][2] ) + R[2][2] - R[1][1] ) / (2*R[1][2]);
            if R[1][2] < 0 then
                P := [ [ 0, -1 ], [ 1, -n ] ];
            else 
                P := [ [ 0, 1 ], [ 1, n ] ];
            fi;
            R := P^-1*R*P;
            Add( conj, Last(conj)*P );
        until R in cyc;

        pos := PositionsProperty( cyc, IsReducedMatrix22 );
        return rec( cycle := cyc{pos}, conj := conj{pos} );
    end;

    ReduceMatrixSquare := function( M )
        local a, d, P, R, n;

        if IsReducedMatrix22( M ) then
            return M;
        fi;

        a := Eigenvectors( Rationals, M )[1];
        d := GcdInt( a[1], a[2] );
        if d <> 1 then
            a := a / d;
        fi;
        d := Gcdex( a[1], a[2] );
        P := [ [ a[1], a[2]], [ -d.coeff2, d.coeff1] ];
        P := P^-1;
        R := P^-1*M*P;
        if R[2][1] < 0 then
            P := P*[[ 1, 0 ], [ 0, -1 ]];
            R := P^-1*M*P;
        fi;

        if R[1][1] = R[2][2] then
            return rec( reduced := R, conj := P );
        else
            if R[1][1] < R[1][1] - R[2][2] then
                if R[1][1] - R[2][2] > 0 then
                    n := Int(Floor( R[2][1]/( R[1][1] - R[2][2] ) ) ) + 1;
                else
                    n := Int(Ceil( R[2][1]/( R[1][1] - R[2][2] ) ) ) - 1;
                fi;
                P := P*[ [1, 0], [-n, 1] ];
                return rec( reduced := P^-1*M*P, conj := P );
            else
                return rec( reduced := R, conj := P ); 
            fi;
        fi;

    end;

    D := DiscriminantMatrix22( M );

    if D<0 then
        return ReduceMatrixPositiveDefinite( M );
    elif D>0 and not IsSquareInt( D ) then
        return ReduceMatrixIndeterminate( M );
    else
        return ReduceMatrixSquare( M );
    fi;

end );

InstallGlobalFunction( ConjugacyMatrix22, function( A, B )
    local   D1, D2,
            red1, red2,
            pos;

    D1 := DiscriminantMatrix22( A );
    D2 := DiscriminantMatrix22( B );

    red1 := ReduceMatrix22( A );
    red2 := ReduceMatrix22( B );

    if D1<0 and D2<0 then
    
        if red1.reduced = red2.reduced then
            return red2.conj*red1.conj^-1;
        else
            return false;
        fi;

    elif D1>0 and not IsSquareInt( D1 ) and D2>0 and not IsSquareInt( D2 ) then

        if red1.cycle[1] in red2.cycle then
            pos := Position( red2.cycle, red1.cycle[1] );
            return red2.conj[pos]*red1.conj[1]^-1;
        else
            return false;
        fi;

    elif IsSquareInt( D1 ) and IsSquareInt( D2 ) then

        if red1.reduced = red2.reduced then
            return red2.conj*red1.conj^-1;
        else
            return false;
        fi;
    else
        return false;
    fi;


end );

InstallGlobalFunction( CentralizerMatrix22, function( A )
    local   fundamentalUnit,
            delta, I2, U, a, b, c, d;

    fundamentalUnit := function( delta )
        local a0, d, b0, q1, q2, A, u, v;

        a0 := 1;
        d  := Int( Floor( Sqrt( Float( delta ) ) ) );
        if (delta - d) mod 2 = 0 then
            b0 := d;
        else
            b0 := d-1;
        fi;

        a  := a0; b  := b0;
        q1 := 1;  q2 := 2;

        repeat
            A  := Int( Floor( Float( (b+d)/(2*a) ) ) );
            b  := 2*a*A-b;
            a  := (delta-b^2)/(4*a);
            q1 := q1*A+q2;
            q2 := q1;
        until a = a0 and b = b0;

        u := ( 2*a0*q2 + q1*b0 )/a0;
        v := q1/a0;

        return [u,v];
    end;

    I2    := [ [1,0], [0,1] ];
    b     := A[1][1]-A[2][2];
    a     := A[1][2];
    c     := A[2][1];

    if a=0 and b=0 and c=0 then
        return [ [ [0, 1], [1, 0] ], [ [-1, 0], [0, 1] ], [ [1, 1], [0, 1] ] ];
    fi;

    d := Gcd( a, b, c );
    a := a/d; b := b/d; c := c/d;
    delta := b^2+4*a*c;
    
    if delta = -3 then
        return [ [ [ (1-b)/2, -a ], [ -c, (b+1)/2 ] ] ];
    elif delta = 4 or delta = -4 then
        return [ -I2, [ [ -b/2, -a ], [ -c, b/2 ] ] ];
    elif delta = 0 then
        return [ -I2, [ [ 1+b/2, a ], [ c, 1-b/2 ] ] ];
    elif delta = 1 then
        return [ -I2, [ [ -b, -2*a ], [ -2*c, b ] ] ];
    elif IsSquareInt( delta ) or delta < 0 then
        return [ -I2 ];
    else
        U := fundamentalUnit( delta );
        return [ -I2, [ [ (U[1]+b*U[2])/2, a*U[2] ], [ c*U[2], (U[1]-b*U[2])/2 ] ] ];
    fi;
end );