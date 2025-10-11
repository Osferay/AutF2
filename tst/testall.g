#
# AutF2: Computations for the automorphisms group of F2
#
# This file runs package tests. It is also referenced in the package
# metadata in PackageInfo.g.
#
LoadPackage( "AutF2" );

TestDirectory(DirectoriesPackageLibrary( "AutF2", "tst" ),
  rec(exitGAP := true));

FORCE_QUIT_GAP(1); # if we ever get here, there was an error
