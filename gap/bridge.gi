LoadPackage( "json" );

AutF2CallCpp := function( func )
    local filename, cpp;

    filename := DirectoriesPackageLibrary( "autf2", "src" );
    cpp := Filename( filename[1], "braidgap" );
    filename := Filename( filename[1], "todo.json" );     
    filename := Concatenation( cpp, " \"", filename, "\"", " \"", func, "\"" );       
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