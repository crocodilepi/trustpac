L = 1.0;
Lc1 = 0.02;
Lc1 = 0.05;

Point( 1 ) = {0, 0, 0, Lc1}; 
Point( 2 ) = {L, 0, 0, Lc1};
Point( 3 ) = {L, L, 0, Lc1};
Point( 4 ) = {0, L, 0, Lc1};

Line( 11 )  = {1 , 2};
Line( 12 )  = {2 , 3};
Line( 13 )  = {3 , 4};
Line( 14 )  = {4 , 1};

Line Loop( 20 ) = {11, 12, 13, 14};

Plane Surface( 30 ) = { 20 };

Physical Line( "bas" ) = {11} ;
Physical Line( "droit" ) = {12} ;
Physical Line( "haut" ) = {13} ;
Physical Line( "gauche" ) = {14} ;

Physical Surface("dom") = { 30 } ;
