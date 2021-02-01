
/*
  All geometrical parameters of the model are used in both 
  Gmsh (CAO) and GetDP (shape optimisation).
  They are therefore defined in the file "ccore_common.pro".
*/

Include "ccore_common.pro";

// The mesh refinement factor is rather a global user parameter,
// therefore defined as a Onelab variable.
R = DefineNumber[1, Name "Model/Parameters/Mesh refinement factor",
		 Help "R=5 far raw mesh, R=1 normal mesh"];

// Center C-core in air box
CoreX=(L-A)/2; // X-position of C-core's bottom left corner
CoreY=(L-B)/2; // Y-position of C-core's bottom left corner

// external boundary
lc1=L/10.*R;

Point(1) = {0,0,0,lc1};
Point(2) = {L,0,0,lc1};
Point(3) = {L,L,0,lc1};
Point(4) = {0,L,0,lc1};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Curve Loop(1) = {1 ... 4};

// magnetic C-core

lc2=A/10.*R;
lc3=D/2.*R;

Point(5) = {CoreX,CoreY,0,lc2};
Point(6) = {CoreX+A,CoreY,0,lc2};
Point(7) = {CoreX+A,CoreY+(B-D)/2.,0,lc3};
Point(8) = {CoreX+A-E,CoreY+(B-D)/2.,0,lc3};
Point(9) = {CoreX+A-E,CoreY+E,0,lc2};
Point(10) = {CoreX+E,CoreY+E,0,lc2};
Point(11) = {CoreX+E,CoreY+B-E,0,lc2};
Point(12) = {CoreX+A-E,CoreY+B-E,0,lc2};
Point(13) = {CoreX+A-E,CoreY+(B+D)/2.,0,lc3};
Point(14) = {CoreX+A,CoreY+(B+D)/2.,0,lc3};
Point(15) = {CoreX+A,CoreY+B,0,lc2};
Point(16) = {CoreX,CoreY+B,0,lc2};

Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,9};
Line(9) = {9,10};
Line(10) = {10,11};
Line(11) = {11,12};
Line(12) = {12,13};
Line(13) = {13,14};
Line(14) = {14,15};
Line(15) = {15,16};
Line(16) = {16,5};

Curve Loop(2) = {5 ... 16};

// inductors

lc4=lc2; //F/2.*R;

Point(17) = {CoreX+E+F,CoreY+E+F,0,lc4};
Point(18) = {CoreX+E+F+G,CoreY+E+F,0,lc4};
Point(19) = {CoreX+E+F+G,CoreY+B-E-F,0,lc4};
Point(20) = {CoreX+E+F,CoreY+B-E-F,0,lc4};
Line(17) = {17,18};
Line(18) = {18,19};
Line(19) = {19,20};
Line(20) = {20,17};

Curve Loop(3) = {17 ... 20};

Point(21) = {CoreX-F-G,CoreY+E+F,0,lc4};
Point(22) = {CoreX-F,CoreY+E+F,0,lc4};
Point(23) = {CoreX-F,CoreY+B-E-F,0,lc4};
Point(24) = {CoreX-F-G,CoreY+B-E-F,0,lc4};
Line(21) = {21,22};
Line(22) = {22,23};
Line(23) = {23,24};
Line(24) = {24,21};

Curve Loop(4) = {21 ... 24};


Plane Surface(1) = {1,2,3,4};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};


Physical Surface("AIR", 1) = 1;
Physical Surface("CORE", 2) = 2;
Physical Surface("COILP", 3) = 3;
Physical Surface("COILN", 4) = 4;

Physical Line("DIR", 11) = {1 ... 4};

Include "onelab_optimize.geo";
