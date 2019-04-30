//+
Point(1) = {0, 0, 0, 5.0};
//+
Point(2) = {0, 150, 0, 5.0};
//+
Point(3) = {100, 0, 0, 5.0};
//+
Point(4) = {100, 150, 0, 5.0};
//+
Line(1) = {2, 4};
//+
Line(2) = {4, 3};
//+
Line(3) = {3, 1};
//+
Line(4) = {1, 2};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Physical Surface("soil", 1) = {1};
//+
Physical Curve("top", 101) = {1};
//+
Physical Curve("bottom", 102) = {3};
//+
Physical Curve("sides", 103) = {4, 2};
//+
Transfinite Curve {1} = 50 Using Progression 1;
