mm=1e-3;
lc=7*mm;

Point(1) = {0, 0, 0, lc};
Point(2) = {600*mm, 0, 0, lc};
Point(3) = {600*mm, 240*mm, 0, lc};
Point(4) = {240*mm, 240*mm, 0, lc};
Point(5) = {240*mm, 600*mm, 0, lc};
Point(6) = {0, 600*mm, 0, lc};
Point(10) = {(600-10)*mm, 240*mm, 0, lc};
Point(11) = {(600-10)*mm, (240-30)*mm, 0, lc};
Point(12) = {600*mm, (240-30)*mm, 0, lc};

Line(1) = {6, 1};
Line(2) = {1, 2};
Line(3) = {2, 12};
Line(4) = {12, 3};
Line(5) = {10, 3};
Line(6) = {10, 4};
Line(7) = {4, 5};
Line(8) = {5, 6};
Line(9) = {12, 11};
Line(10) = {11, 10};

Line Loop(1) = {1, 2, 3, 9, 10, 6, 7, 8};
Plane Surface(1) = {1};

Line Loop(2) = -{9, 10, 5, -4};
Plane Surface(2) = {2};

Physical Surface(1000) = {1};
Physical Surface(1004) = {2};
Physical Line(1001) = {8};
Physical Line(1003) = {4};
Physical Point(1002) = {3};

Solver.AutoMesh = 0;
//Mesh.Algorithm = 1; 

