DefineConstant[
  mm = 1e-3
  deg = Pi/180
  R1 = 6*mm
  L2 = 18*mm
  L3 = 15*mm
  L4 = {L2/L3*Sqrt[L3^2-(12.5*mm)^2], Name "Optimization/parameters/L4", Closed 1}
  angleMeasure = 50*deg
  lc = 12*mm
  lcd = lc/15

  // mesh perturb
  PerturbMesh = 0
  VelocityTag = -1
  Perturbation = 1e-6
  modelpath = CurrentDir
];

Point(1) = {0, 0, 0, lcd};

For k In {0:10}
  theta = k * 90/10 * deg;
  DefineConstant[Rs~{k} = {R1, Name Sprintf["Optimization/parameters/spline radius %g",k+1], Closed 1}];
  Point(2020+k) = {Rs~{k}*Cos[theta], Rs~{k}*Sin[theta], 0, lcd};  
EndFor

BSpline(4020) = {2020,2021, 2022, 2023};
BSpline(4021) = {2023, 2024, 2025, 2026};
BSpline(4022) = {2026, 2027, 2028, 2029, 2030};

//Point(26) = {L2/L3*Sqrt[L3^2-((12.5-2)*mm)^2], (12.5-2)*mm, 0, lcd};

For k In {0:9}
  theta = k * 50/10 * deg;
  DefineConstant[Rso~{k} = {L2*L3/Sqrt[(L3*Cos[theta])^2+(L2*Sin[theta])^2], 
    Name Sprintf["Optimization/parameters/outer spline radius %g",k+1], Closed 1}];
  Point(6020+k) = {Rso~{k}*Cos[theta], Rso~{k}*Sin[theta], 0, lcd};  
EndFor

Point(24) = {20*mm-L4, 12.5*mm, 0, lcd};

BSpline(4023) = {6020, 6021, 6022, 6023};
BSpline(4024) = {6023, 6024, 6025, 6026};
BSpline(4025) = {6026, 6027, 6028};
BSpline(4026) = {6028, 6029, 24};

Point(6) = {20*mm, 0, 0, lcd};
Point(8) = {0, L3, 0, lc};
Point(9) = {25*mm, 0, 0, lc};
Point(10) = {163*mm, 0, 0, lc};
Point(11) = {25*mm, 50*mm, 0, lc};
Point(12) = {113*mm, 50*mm, 0, lc};
Point(13) = {113*mm, (50+80)*mm, 0, lc};
Point(14) = {0, (50+80)*mm, 0, lc};
Point(15) = {0, 180*mm, 0, lc};
Point(16) = {163*mm, 180*mm, 0, lc};
Point(17) = {25*mm, (50+7.5)*mm, 0, lc};
Point(18) = {113*mm, (50+7.5)*mm, 0, lc};
Point(19) = {25*mm, (50+7.5+39)*mm, 0, lc};
Point(20) = {113*mm, (50+7.5+39)*mm, 0, lc};

Point(1021) = {(9.5+2.25)*Cos[angleMeasure]*mm, (9.5+2.25)*Sin[angleMeasure]*mm, 0, lcd};
Point(1022) = {(9.5+2.25)*mm, 0, 0, lcd};

Point(23) = {20*mm, 12.5*mm, 0, lcd};

//DefineConstant[edge_x_1 = {20*mm-L4, Name "Optimization/parameters/edge radius 1"}];
//DefineConstant[edge_y_1 = {(12.5-2)*mm, Name "Optimization/parameters/edge radius 2"}];
//Point(25) = {edge_x_1, edge_y_1, 0, lcd};

Circle(1016) = {1022, 1, 1021};
Line(2) = {9, 10};
Line(3) = {10, 16};
Line(4) = {16, 15};
Line(5) = {15, 14};
Line(6) = {14, 13};
Line(8) = {12, 11};
Line(9) = {11, 9};
Line(10) = {17, 18};
Line(11) = {12, 18};
Line(12) = {18, 20};
Line(13) = {20, 19};
Line(14) = {19, 17};
Line(15) = {20, 13};
Line(18) = {1, 2020};
Line(24) = {2020, 6020};
Line(22) = {1, 2030};
Line(26) = {6, 23};
Line(27) = {24, 23};

Line(30) = {6020, 6};
Line(31) = {6, 9};
Line(32) = {2030, 14};

Curve Loop(1) = {-22, 18, 4020, 4021, 4022};
Plane Surface(1) = {1};
Curve Loop(2) = {4023,4024,4025,4026, 27, -26, -30};
Plane Surface(2) = {2};
Curve Loop(3) = {9, 2, 3, 4, 5, 6, -15, -12, -11, 8};
Plane Surface(3) = {3};
Curve Loop(4) = {13, 14, 10, 12};
Plane Surface(4) = {4};
//Curve Loop(6) = {24, 4023,4024,4025,4026, 27, 26, 31, -9, -8, 11, -10, -14, -13, 15, -6, -32, -4020, -4021, -4022};
//Plane Surface(6) = {6};
Curve Loop(6) = {24, 4023, 4024, 4025, 4026, 27, -26, 31, -9, -8, 11, -10, -14, -13, 15, -6, -32, -4022, -4021, -4020};
Plane Surface(6) = {6};

Physical Surface("pole", 1) = {3};
Physical Surface("inductor", 2) = {4};
Physical Surface("die molds", 3) = {1,2};
Physical Surface("air", 5) = {6};
Physical Line("dirichlet", 6) = {3,4,18,20,24,30,31,2};
Physical Line("neuman", 7) = {22,32};

If(PerturbMesh == 1)
  Printf("Computing velocity field ...");
  modelName = StrPrefix(StrRelative(General.FileName));  

  SyncModel;
  
  // Merge the original mesh
  Merge StrCat(modelpath, modelName, ".msh");

  // create a view with the original node coordinates as a vector dataset
  Plugin(NewView).NumComp = 3;
  Plugin(NewView).Run;
  Plugin(ModifyComponents).View = 0;
  Plugin(ModifyComponents).Expression0 = "x";
  Plugin(ModifyComponents).Expression1 = "y";
  Plugin(ModifyComponents).Expression2 = "z";
  Plugin(ModifyComponents).Run;

  // relocate the vertices of the original mesh on the perturbed geometry
  RelocateMesh Point "*";
  RelocateMesh Line "*";
  RelocateMesh Surface "*";

  // Create a view with the perturbed node coordinates as vector dataset
  Plugin(NewView).NumComp = 3;
  Plugin(NewView).Run;
  Plugin(ModifyComponents).View = 1;
  Plugin(ModifyComponents).Expression0 = "x";
  Plugin(ModifyComponents).Expression1 = "y";
  Plugin(ModifyComponents).Expression2 = "z";
  Plugin(ModifyComponents).Run;

  // compute the velocity field by subtracting the two vector datasets
  Plugin(ModifyComponents).View = 0;
  Plugin(ModifyComponents).OtherView = 1;
  Plugin(ModifyComponents).Expression0 = Sprintf("(w0 - v0)/(%g)", Perturbation);
  Plugin(ModifyComponents).Expression1 = Sprintf("(w1 - v1)/(%g)", Perturbation);
  Plugin(ModifyComponents).Expression2 = Sprintf("(w2 - v2)/(%g)", Perturbation);
  Plugin(ModifyComponents).Run;
  View.Name = "velocity";
  Delete View[1];
  SendToServer View[0] Sprintf("Optimization/Results/velocity_%g",VelocityTag); // Hidden
EndIf
