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
Line(3) = {1, 1};
//+
Line(4) = {3, 1};
//+
Line(5) = {1, 2};
//+
Line Loop(1) = {1, 2, 4, 5};
//+
Plane Surface(1) = {1};
//+
Physical Surface("soil", 1) = {1};
//+
Physical Line("top", 101) = {1};
//+
Physical Line("bottom", 102) = {4};
//+
Physical Line("sides", 103) = {5, 2};
//+
Transfinite Line {1} = 50 Using Progression 1;
