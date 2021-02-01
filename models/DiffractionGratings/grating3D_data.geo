
DefineConstant[
  test_case = {"halfellipsoid",
    Name "1Geometry",
    Choices {"hole", "pyramid", "torus", "2Dlamellar", "solarcell", "conv", "skew"},
    GmshOption "Reset", Autocheck 0
  }
];

Include StrCat["grating3D_data_",test_case,".geo"];
