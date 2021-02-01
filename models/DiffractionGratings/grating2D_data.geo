

DefineConstant[
  test_case = {"LamellarGrating",
    Name "0Test case",
    Choices {"AnisotropicGrating", "LamellarGrating", "PhotonicCrystalSlab", "ResonantGrating", "plasmonics"},
    GmshOption "Reset", Autocheck 0
  }
];

nb_plot_periods = 3;

Include StrCat["grating2D_data_",test_case,".geo"];
