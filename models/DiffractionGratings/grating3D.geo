// Copyright (C) 2020 Guillaume Demésy
//
// This file is part of the model grating3D.pro.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with This program. If not, see <https://www.gnu.org/licenses/>.

Include "grating3D_data.geo";

SetFactory("OpenCASCADE");
e = 5*nm;
E = 10*(period_x+period_x/2);

Macro SetPBCs
  BlochXm() = Surface In BoundingBox{.5*(-period_x-dys)-e,-dyc/2-e, PML_bot_hh-e,(-period_x+dys)/2+e, dyc/2+e, PML_top_hh+PML_top+e};
  BlochXp() = Surface In BoundingBox{.5*( period_x-dys)-e,-dyc/2-e, PML_bot_hh-e,( period_x+dys)/2+e, dyc/2+e, PML_top_hh+PML_top+e};
  BlochYm() = Surface In BoundingBox{.5*(-period_x-dys)-e,-dyc/2-e, PML_bot_hh-e,( period_x+dys)/2+e,-dyc/2+e, PML_top_hh+PML_top+e};
  BlochYp() = Surface In BoundingBox{.5*(-period_x-dys)-e, dyc/2-e, PML_bot_hh-e,( period_x+dys)/2+e, dyc/2+e, PML_top_hh+PML_top+e};
  // For k In {1:#BlochXm()-1}
  //   Printf("BlochXm surf %g",BlochXm(k));
  // EndFor
  // For k In {1:#BlochYm()-1}
  //   Printf("BlochYm surf %g",BlochYm(k));
  // EndFor
  // For k In {1:#BlochXp()-1}
  //   Printf("BlochXp surf %g",BlochXp(k));
  // EndFor
  // For k In {1:#BlochYp()-1}
  //   Printf("BlochYp surf %g",BlochYp(k));
  // EndFor
  If (tag_geom==6) // bi-sin : BoundingBox does not catch BSpline Surfaces?
    BlochXm()+={17,12};
    BlochXp()+={19,14};
    BlochYm()+={15,20};
    BlochYp()+={13,18};
  EndIf
  Periodic Surface{BlochXp()} = {BlochXm()} Translate{period_x,        0,        0};
  Periodic Surface{BlochYp()} = {BlochYm()} Translate{    dys, dyc,        0};
Return

If (tag_geom==6)
  // Bi-sinusoidal profile by Spline surface filling
  // not the most general way to implement a freesurface height=f(x,y) and relies on multiple Coherence calls
  // several hacks are used, one of them is to do all this at the very begining
  nsin = 20;
  For i In {0:nsin}
      xx1~{i} = -period_x/2+i*period_x/2/(nsin-1);
      xx2~{i} =             i*period_x/2/(nsin-1);
  EndFor
  pp=newp;//pp=pp-1;
  For i In {0:nsin-1}
    yy=-period_y/2;
    Point(newp) = {xx1~{i},yy,rz/4*(Cos(2*Pi*xx1~{i}/period_x)+Cos(2*Pi*yy/period_y))+hh_L_3+thick_L_3/2};
    Point(newp) = {yy,xx1~{i},rz/4*(Cos(2*Pi*yy/period_x)+Cos(2*Pi*xx1~{i}/period_y))+hh_L_3+thick_L_3/2};
    yy=0;
    Point(newp) = {xx1~{i},yy,rz/4*(Cos(2*Pi*xx1~{i}/period_x)+Cos(2*Pi*yy/period_y))+hh_L_3+thick_L_3/2};
    Point(newp) = {yy,xx1~{i},rz/4*(Cos(2*Pi*yy/period_x)+Cos(2*Pi*xx1~{i}/period_y))+hh_L_3+thick_L_3/2};
    yy=period_y/2;
    Point(newp) = {xx1~{i},yy,rz/4*(Cos(2*Pi*xx1~{i}/period_x)+Cos(2*Pi*yy/period_y))+hh_L_3+thick_L_3/2};
    Point(newp) = {yy,xx1~{i},rz/4*(Cos(2*Pi*yy/period_x)+Cos(2*Pi*xx1~{i}/period_y))+hh_L_3+thick_L_3/2};
    yy=-period_y/2;
    Point(newp) = {xx2~{i},yy,rz/4*(Cos(2*Pi*xx2~{i}/period_x)+Cos(2*Pi*yy/period_y))+hh_L_3+thick_L_3/2};
    Point(newp) = {yy,xx2~{i},rz/4*(Cos(2*Pi*yy/period_x)+Cos(2*Pi*xx2~{i}/period_y))+hh_L_3+thick_L_3/2};
    yy=0;
    Point(newp) = {xx2~{i},yy,rz/4*(Cos(2*Pi*xx2~{i}/period_x)+Cos(2*Pi*yy/period_y))+hh_L_3+thick_L_3/2};
    Point(newp) = {yy,xx2~{i},rz/4*(Cos(2*Pi*yy/period_x)+Cos(2*Pi*xx2~{i}/period_y))+hh_L_3+thick_L_3/2};
    yy=period_y/2;
    Point(newp) = {xx2~{i},yy,rz/4*(Cos(2*Pi*xx2~{i}/period_x)+Cos(2*Pi*yy/period_y))+hh_L_3+thick_L_3/2};
    Point(newp) = {yy,xx2~{i},rz/4*(Cos(2*Pi*yy/period_x)+Cos(2*Pi*xx2~{i}/period_y))+hh_L_3+thick_L_3/2};
  EndFor
  For k In {1:11}
    BSpline(newl) = {pp+k:pp+12*nsin:12};
  EndFor
  BSpline(newl) = {pp  :pp+12*(nsin-1):12};
  Coherence;
  Box(1) = {-period_x/2,-period_y/2,     hh_L_3,period_x,period_y, thick_L_3};
  Delete{Volume{1}; Surface{1,2,3,4} ; Line{13,15,17,19};}
  Curve Loop(7) = {21, 20, -23, -16};
  Curve Loop(8) = {1, 2, -3, -12};
  Surface(7) = {8};
  Curve Loop(10) = {6, 5, -8, -3};
  Surface(8) = {10};
  Curve Loop(12) = {9, 10, -11, -8};
  Surface(9) = {12};
  Curve Loop(14) = {2, 9, -4, -7};
  Surface(10) = {14};
  Line(25) = {239, 229};
  Line(26) = {229, 238};
  Line(27) = {243, 235};
  Line(28) = {235, 242};
  Line(29) = {237, 244};
  Line(30) = {233, 240};
  Line(31) = {237, 245};
  Line(32) = {241, 233};
  Curve Loop(16) = {25, 12, 6, -27, -21};
  Surface(11) = {16};
  Curve Loop(18) = {20, -31, -11, -5, -27};
  Surface(12) = {18};
  Curve Loop(20) = {29, -18, -28, 5, 11};
  Surface(13) = {20};
  Curve Loop(22) = {10, 29, -24, -30, 4};
  Surface(14) = {22};
  Curve Loop(24) = {30, -14, -26, 1, 7};
  Surface(15) = {24};
  Curve Loop(26) = {25, 1, 7, -32, -16};
  Surface(16) = {26};
  Curve Loop(28) = {26, 22, -28, -6, -12};
  Surface(17) = {28};
  Curve Loop(30) = {31, -23, 32, 4, 10};
  Surface(18) = {30};
  Surface Loop(3) = {5, 16, 18, 12, 11, 9, 8, 7, 10};
  Volume(1) = {3};
  Surface Loop(2) = {6, 15, 14, 10, 7, 8, 9, 13, 17};
  Volume(2) = {2};
EndIf

For k In {7:0:-1}
  If (tag_geom==6)
    If (k!=3)
      p=newp; l=newl; ll=newll; s=news;
      Point(p)   = {0.5*(-period_x-dys), -dyc/2, hh_L~{k}};
      Point(p+1) = {0.5*( period_x-dys), -dyc/2, hh_L~{k}};
      Point(p+2) = {0.5*( period_x+dys),  dyc/2, hh_L~{k}};
      Point(p+3) = {0.5*(-period_x+dys),  dyc/2, hh_L~{k}};
      Line(l)={p,p+1};Line(l+1)={p+1,p+2};Line(l+2)={p+2,p+3};Line(l+3)={p+3,p};
      Curve Loop(ll) = {l,l+1,l+2,l+3};
      Surface(s) = {ll};
      Extrude {0, 0, thick_L~{k}} {Surface{s};}
    Else
      v=newv;
    EndIf
  Else
    p=newp; l=newl; ll=newll; s=news;
    Point(p)   = {0.5*(-period_x-dys), -dyc/2, hh_L~{k}};
    Point(p+1) = {0.5*( period_x-dys), -dyc/2, hh_L~{k}};
    Point(p+2) = {0.5*( period_x+dys),  dyc/2, hh_L~{k}};
    Point(p+3) = {0.5*(-period_x+dys),  dyc/2, hh_L~{k}};
    Line(l)={p,p+1};Line(l+1)={p+1,p+2};Line(l+2)={p+2,p+3};Line(l+3)={p+3,p};
    Curve Loop(ll) = {l,l+1,l+2,l+3};
    Surface(s) = {ll};
    Extrude {0, 0, thick_L~{k}} {Surface{s};}
  EndIf
EndFor
// Box(1) = {-period_x/2,-period_y/2, PML_bot_hh,period_x,period_y, PML_bot  };
// Box(2) = {-period_x/2,-period_y/2,     hh_L_6,period_x,period_y, thick_L_6};
// Box(3) = {-period_x/2,-period_y/2,     hh_L_5,period_x,period_y, thick_L_5};
// Box(4) = {-period_x/2,-period_y/2,     hh_L_4,period_x,period_y, thick_L_4};
// If (tag_geom!=6)
//   Box(5) = {-period_x/2,-period_y/2,     hh_L_3,period_x,period_y, thick_L_3};
// EndIf
// Box(6) = {-period_x/2,-period_y/2,     hh_L_2,period_x,period_y, thick_L_2};
// Box(7) = {-period_x/2,-period_y/2,     hh_L_1,period_x,period_y, thick_L_1};
// Box(8) = {-period_x/2,-period_y/2, PML_top_hh,period_x,period_y, PML_top  };

If (tag_geom==1)
  p=newp;
  Point(p) = {0,0,hh_L_3+rz};
  Line(97) = {29, 65};
  Line(98) = {65, 32};
  Line(99) = {65, 30};
  Line(100) = {65, 31};
  Curve Loop(92) = {97, 99, -43};
  Plane Surface(91) = {92};
  Curve Loop(93) = {99, 45, -100};
  Plane Surface(92) = {93};
  Curve Loop(94) = {100, 47, -98};
  Plane Surface(93) = {94};
  Curve Loop(95) = {97, 98, 48};
  Plane Surface(94) = {95};
  Surface Loop(9) = {91, 94, 93, 92, 42};
  Volume(9) = {9};
EndIf

If (tag_geom==2)
  Cylinder(9) = {0,0,hh_L_3,0,0,thick_L_3,rx};
EndIf

If (tag_geom==3)
  Torus(9) = {0,0,hh_L_3+ry-1*nm,rx,ry};
  Dilate { { 0,0,hh_L_3 }, { 1, 1, rz/ry } } { Volume{9}; }
  BooleanDifference{ Volume{9}; Delete; }{ Volume{4}; }
EndIf

If (tag_geom==4)
  Sphere(9) = {0,0,hh_L_3,rx};
  Dilate { { 0,0,hh_L_3 }, { 1, ry/rx, rz/rx } } { Volume{9}; }
  BooleanDifference{ Volume{9}; Delete; }{ Volume{4}; }
EndIf

If (tag_geom==5)
  Box(9) = {-rx/2,-ry/2, hh_L_2-rz, rx, ry, rz};
  Rotate {{0, 0, 1}, {0,0,0}, Pi/4} {Volume{9};}
EndIf

If (tag_geom==7)
  Box(9) = {-period_x/2,-ry/2, hh_L_3, period_x, ry, rz};
EndIf


Coherence;
Call SetPBCs;
If (tag_geom==6)
  Physical Volume("PMLBOT"        ,1) = {3};
  Physical Volume("LAYER_L6_SUBS" ,2) = {4};
  Physical Volume("LAYER_L5"      ,3) = {5};
  Physical Volume("LAYER_L4"      ,4) = {6};
  Physical Volume("LAYER_L3"      ,5) = {2};
  Physical Volume("LAYER_L2"      ,6) = {7};
  Physical Volume("LAYER_L1_SUPER",7) = {8};
  Physical Volume("PMLTOP"        ,8) = {9};
  Physical Volume("SCAT"          ,9) = {1};
Else
  Physical Volume("PMLBOT"        ,1) = {1};
  Physical Volume("LAYER_L6_SUBS" ,2) = {2};
  Physical Volume("LAYER_L5"      ,3) = {3};
  Physical Volume("LAYER_L4"      ,4) = {4};
  Physical Volume("LAYER_L3"      ,5) = {10};
  Physical Volume("LAYER_L2"      ,6) = {6};
  Physical Volume("LAYER_L1_SUPER",7) = {7};
  Physical Volume("PMLTOP"        ,8) = {8};
  Physical Volume("SCAT"          ,9) = {9};
EndIf
Physical Surface("BXM" ,101 ) = BlochXm();
Physical Surface("BXP" ,102 ) = BlochXp();
Physical Surface("BYM" ,201 ) = BlochYm();
Physical Surface("BYP" ,202 ) = BlochYp();
Physical Surface("STOP",301 ) = Surface In BoundingBox{-period_x/2-E,-period_y/2-E, PML_top_hh-e,period_x/2+E, period_y/2+E, PML_top_hh+e};
Physical Surface("SBOT",302 ) = Surface In BoundingBox{-period_x/2-E,-period_y/2-E, hh_L_6-e,period_x/2+E    , period_y/2+E, hh_L_6+e};
Physical Surface("SPMLTOP",401 ) = Surface In BoundingBox{-period_x/2-E,-period_y/2-E, PML_top_hh+PML_top-e,period_x/2+E, period_y/2+E, PML_top_hh+PML_top+e};
Physical Surface("SPMLBOT",402 ) = Surface In BoundingBox{-period_x/2-E,-period_y/2-E, PML_bot_hh-e,period_x/2+E, period_y/2+E, PML_bot_hh+e};
Physical Surface("SVIS",500 )    = Surface In BoundingBox{-period_x/2-E,-period_y/2-E, hh_L_3-e,period_x/2+E    , period_y/2+E, hh_L_3+e};

pts_PMLBOT()   = PointsOf{ Physical Volume{1}; };
pts_LAYER_L6() = PointsOf{ Physical Volume{2}; };
pts_LAYER_L5() = PointsOf{ Physical Volume{3}; };
pts_LAYER_L4() = PointsOf{ Physical Volume{4}; };
pts_LAYER_L3() = PointsOf{ Physical Volume{5}; };
pts_LAYER_L2() = PointsOf{ Physical Volume{6}; };
pts_LAYER_L1() = PointsOf{ Physical Volume{7}; };
pts_PMLTOP()   = PointsOf{ Physical Volume{8}; };
pts_ROD()      = PointsOf{ Physical Volume{9}; };

// if test_case==conv, we want to remesh the scatterer as well
// the following enforces to update lc_scat when paramaille changes
If (tag_geom==4)
  lc_scat = lambda_m/(paramaille*5);
EndIf

lc_PML = lambda_m/(paramaille*.5);
list_lc(0) = lc_PML;
For k In {1:6}
  n_L~{k}  = Sqrt(Abs(eps_re_L~{k}));
  lc_L~{k} = lambda_m/(paramaille*n_L~{k})/refine_mesh_L~{k};
  list_lc(7-k) = lc_L~{k};
EndFor
list_lc(7) = lc_PML;
list_lc(8) = lc_scat;

// This helps meshing: The default behavior of the PointsOf technique
// overides points belonging to several domains (by order of call of Characteristic Length)
If (tag_geom==7)
  meshing_sequence() = {1,8,2,3,5,6,7,4,9};
Else
    meshing_sequence() = {1,8,2,3,4,6,7,5,9};
EndIf

// Start with coarsest
Characteristic Length{:} = lc_PML;

For k In {0:8}
  which_dom = meshing_sequence(k);
  Characteristic Length{PointsOf{Physical Volume{which_dom};}} = list_lc(which_dom-1);
EndFor


If (tag_geom==3) // Split torus weird otherwise
  Mesh.Algorithm = 6;
EndIf
// Mesh.SurfaceEdges = 0;
Mesh.VolumeEdges = 0;
Mesh.ElementOrder = og;
