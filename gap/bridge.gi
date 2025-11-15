LoadPackage( "json" );

AutF2CallCpp := function( func )
    local filename, cpp;

    filename := DirectoriesPackageLibrary( "autf2", "src" );

    if func = "lcf" then
        cpp := Filename( filename[1], "lcf.o" );
    elif func = "conj" then
        cpp := Filename( filename[1], "conj.o" );
    elif func = "cent" then
        cpp := Filename( filename[1], "cent.o" );
    else
        Error( "No function implemented in c++." );
    fi;

    filename := Filename( filename[1], "todo.json" );
    filename := Concatenation( cpp, " \"", filename, "\"" );
    
    Exec( filename );

end;

AutF2WriteJSON := function( value )
    local filename, stream;

    filename := DirectoriesPackageLibrary( "autf2", "src" );
    filename := Filename( filename[1], "todo.json" );

    stream   := OutputTextFile( filename, false );
    GapToJsonStream( stream, value );
    CloseStream( stream );

    return true;
end;

AutF2ReadJSON := function( )
    local filename, stream, value;

    filename := DirectoriesPackageLibrary( "autf2", "src" );
    filename := Filename( filename[1], "todo.json" );

    stream   := InputTextFile( filename );
    value    := JsonStreamToGap( stream );
    CloseStream( stream );

    return value;
end;