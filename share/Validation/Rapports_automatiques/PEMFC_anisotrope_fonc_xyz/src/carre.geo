//+
Dext = 100;
//+
Dint = 4;
//+
lc1 = Dint/4.;
//+
lc2 = 2*lc1;
//+
Point(1) = {0, 0, 0, lc2};
//+
Point(2) = {Dext, 0, 0, lc2};
//+
Point(3) = {Dext, Dext, 0, lc2};
//+
Point(4) = {0, Dext, 0, lc2};
//+
Point(5) = {0.5*(Dext-Dint), 0.5*(Dext-Dint), 0, lc1};
//+
Point(6) = {0.5*(Dext+Dint), 0.5*(Dext-Dint), 0, lc1};
//+
Point(7) = {0.5*(Dext+Dint), 0.5*(Dext+Dint), 0, lc1};
//+
Point(8) = {0.5*(Dext-Dint), 0.5*(Dext+Dint), 0, lc1};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 5};
//+
Line Loop(9) = {5, 6, 7, 8};
//+
Line Loop(10) = {1, 2, 3, 4};
//+
Plane Surface(11) = {9, 10};
//+
Physical Line("EXT") = {1, 2, 3, 4};
//+
Physical Line("INT") = {5, 6, 7, 8};
//+
Physical Surface("dome") = {11};
