mm = 1.e-3;
deg = Pi/180.;
DefineConstant[
  NumMagnets = {2, Min 1, Max 20, Step 1, Name "Parameters/0Number of magnets"}
  Flag_InfiniteBox = {1, Choices{0,1}, Name "Infinite box/Add infinite box"}
  Flag_FullMenu = {0, Choices{0,1}, Name "Parameters/Show all parameters"}
];

For i In {1:NumMagnets}
  DefineConstant[
    X~{i} = {0, Min -100*mm, Max 100*mm, Step mm, Visible Flag_FullMenu,
      Name Sprintf("Parameters/Magnet %g/0X position [m]", i) },
    Y~{i} = { (i-1)*60*mm, Min -100*mm, Max 100*mm, Step mm,
      Name Sprintf("Parameters/Magnet %g/0Y position [m]", i) },
    Z~{i} = {0, Min -100*mm, Max 100*mm, Step mm, Visible Flag_FullMenu,
      Name Sprintf("Parameters/Magnet %g/0Z position [m]", i) },

    M~{i} = {(i==1)?0:1, Choices{0="Cylinder",1="Cube"},
      Name Sprintf("Parameters/Magnet %g/00Shape", i)},

    R~{i} = {20*mm, Min mm, Max 100*mm, Step mm,
      Name Sprintf("Parameters/Magnet %g/1Radius [m]", i),
      Visible (M~{i} == 0) },
    L~{i} = {50*mm, Min mm, Max 100*mm, Step mm,
      Name Sprintf("Parameters/Magnet %g/1Length [m]", i),
      Visible (M~{i} == 0) },

    Lx~{i} = {50*mm, Min mm, Max 100*mm, Step mm,
      Name Sprintf("Parameters/Magnet %g/1X length [m]", i),
      Visible (M~{i} == 1) },
    Ly~{i} = {50*mm, Min mm, Max 100*mm, Step mm, Visible Flag_FullMenu,
      Name Sprintf("Parameters/Magnet %g/1XY aspect ratio", i),
      Visible (M~{i} == 1) },
    Lz~{i} = {50*mm, Min mm, Max 100*mm, Step mm, Visible Flag_FullMenu,
      Name Sprintf("Parameters/Magnet %g/1XZ aspect ration", i),
      Visible (M~{i} == 1) },

    Rx~{i} = {0, Min -Pi, Max Pi, Step Pi/180, Visible Flag_FullMenu,
      Name Sprintf("Parameters/Magnet %g/2X rotation [deg]", i) },
    Ry~{i} = {0, Min -Pi, Max Pi, Step Pi/180, Visible Flag_FullMenu,
      Name Sprintf("Parameters/Magnet %g/2Y rotation [deg]", i) },
    Rz~{i} = {0, Min -Pi, Max Pi, Step Pi/180,
      Name Sprintf("Parameters/Magnet %g/2Z rotation [deg]", i) },

    MUR~{i} = {(i==1)?1.:1000.,
      Name Sprintf("Parameters/Magnet %g/3Mu relative []", i)},
    BR~{i} = {(i==1)?1.0:0.0,
      Name Sprintf("Parameters/Magnet %g/3Br [T]", i)}
  ];
EndFor

//The geometrical parameters of the Infinite box. 
DefineConstant[
  xInt = {1, Name "Infinite box/xInt", Visible 0}
  yInt = {1, Name "Infinite box/yInt", Visible 0}
  zInt = {1, Name "Infinite box/zInt", Visible 0}
  xExt = {xInt*2, Name "Infinite box/xExt", Visible 0}
  yExt = {yInt*2, Name "Infinite box/yExt", Visible 0}
  zExt = {zInt*2, Name "Infinite box/zExt", Visible 0}
  xCnt = {0, Name "Infinite box/xCenter", Visible 0}
  yCnt = {0, Name "Infinite box/yCenter", Visible 0}
  zCnt = {0, Name "Infinite box/zCenter", Visible 0}
];
