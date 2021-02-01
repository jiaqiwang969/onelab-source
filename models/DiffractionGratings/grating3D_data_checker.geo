nm  = 10000;
pp1 = "1Incident Plane Wave";
pp2 = "2Layers Thicknesses";
pp3 = "3Scatterer Properties";
pp4 = "4Layer Materials";
pp5 = "5Computational Parameters";
pp6 = "6Output";
DefineConstant[
    lambda0       = {500   , Name StrCat[pp1,"/1lambda0 [nm]"]},
    thetadeg      = {0     , Name StrCat[pp1,"/2theta0 [deg]"]},
    phideg        = {0     , Name StrCat[pp1,"/3phi0 [deg]"]},
    psideg        = {45    , Name StrCat[pp1,"/4psi0 [deg]"]},
    period_x      = {1.25*Sqrt[2]*lambda0, Name StrCat[pp2,"/1X period [nm]"]},
    period_y      = {1.25*Sqrt[2]*lambda0, Name StrCat[pp2,"/2Y period [nm]"]},
    thick_L_1     = {50    , Name StrCat[pp2,"/3thickness layer 1 [nm] (superstrate)"]},
    thick_L_2     = {50    , Name StrCat[pp2,"/4thickness layer 2 [nm]"]},
    thick_L_3     = {1.2*lambda0, Name StrCat[pp2,"/5thickness layer 3 [nm]"]},
    thick_L_4     = {50    , Name StrCat[pp2,"/6thickness layer 4 [nm]"]},
    thick_L_5     = {50    , Name StrCat[pp2,"/7thickness layer 5 [nm]"]},
    thick_L_6     = {50    , Name StrCat[pp2,"/8thickness layer 6 [nm] (substrate)"]},
    xsideg        = {0    , Name StrCat[pp2,"/9skew angle [deg]"],Visible 1},

    tag_geom      = {  5      , Name StrCat[pp3,"/0Shape"], Choices {1="Pyramid",2="Cylindrical Hole",3="Torus",4="HalfEllipspoid",5="Checkerboard",6="bi-sinusoidal",7="2D lamellar"}},
    rx            = {1.25*lambda0, Name StrCat[pp3,"/1rx"]},
    ry            = {1.25*lambda0, Name StrCat[pp3,"/2ry"]},
    rz            = {lambda0     , Name StrCat[pp3,"/3rz"]},
    flag_mat_scat = { 0    , Name StrCat[pp3,"/4Scatterer permittivity model"], Choices {0="Custom (Value Below)",1="SiO2",2="Ag (palik)",3="Al (palik)",4="Au (johnson)",5="Nb2O5",6="ZnSe",7="MgF2",8="TiO2",9="PMMA",10="Si",11="ITO",12="Cu (palik)"} },
    eps_re_Scat   = { 2.25 , Name StrCat[pp3,"/7eps_re_Scat"]},
    eps_im_Scat   = { 0    , Name StrCat[pp3,"/8eps_im_Scat"]},

    flag_mat_1    = { 0    , Name StrCat[pp4,"/1Layer 1"], Choices {0="Custom (Value Below)",1="SiO2",2="Ag (palik)",3="Al (palik)",4="Au (johnson)",5="Nb2O5",6="ZnSe",7="MgF2",8="TiO2",9="PMMA",10="Si",11="ITO",12="Cu (palik)"} },
    flag_mat_2    = { 0    , Name StrCat[pp4,"/2Layer 2"], Choices {0="Custom (Value Below)",1="SiO2",2="Ag (palik)",3="Al (palik)",4="Au (johnson)",5="Nb2O5",6="ZnSe",7="MgF2",8="TiO2",9="PMMA",10="Si",11="ITO",12="Cu (palik)"} },
    flag_mat_3    = { 0    , Name StrCat[pp4,"/3Layer 3"], Choices {0="Custom (Value Below)",1="SiO2",2="Ag (palik)",3="Al (palik)",4="Au (johnson)",5="Nb2O5",6="ZnSe",7="MgF2",8="TiO2",9="PMMA",10="Si",11="ITO",12="Cu (palik)"} },
    flag_mat_4    = { 0    , Name StrCat[pp4,"/4Layer 4"], Choices {0="Custom (Value Below)",1="SiO2",2="Ag (palik)",3="Al (palik)",4="Au (johnson)",5="Nb2O5",6="ZnSe",7="MgF2",8="TiO2",9="PMMA",10="Si",11="ITO",12="Cu (palik)"} },
    flag_mat_5    = { 0    , Name StrCat[pp4,"/5Layer 5"], Choices {0="Custom (Value Below)",1="SiO2",2="Ag (palik)",3="Al (palik)",4="Au (johnson)",5="Nb2O5",6="ZnSe",7="MgF2",8="TiO2",9="PMMA",10="Si",11="ITO",12="Cu (palik)"} },
    flag_mat_6    = { 0    , Name StrCat[pp4,"/6Layer 6"], Choices {0="Custom (Value Below)",1="SiO2",2="Ag (palik)",3="Al (palik)",4="Au (johnson)",5="Nb2O5",6="ZnSe",7="MgF2",8="TiO2",9="PMMA",10="Si",11="ITO",12="Cu (palik)"} },
    eps_re_L_1    = {2.25  , Name StrCat[pp4,"/layer 1: real part of relative permittivity"]},
    eps_im_L_1    = {0     , Name StrCat[pp4,"/layer 1: imag part of relative permittivity"]},
    eps_re_L_2    = {2.25  , Name StrCat[pp4,"/layer 2: real part of relative permittivity"]},
    eps_im_L_2    = {0     , Name StrCat[pp4,"/layer 2: imag part of relative permittivity"]},
    eps_re_L_3    = {1     , Name StrCat[pp4,"/layer 3: real part of relative permittivity"]},
    eps_im_L_3    = {0     , Name StrCat[pp4,"/layer 3: imag part of relative permittivity"]},
    eps_re_L_4    = {1     , Name StrCat[pp4,"/layer 4: real part of relative permittivity"]},
    eps_im_L_4    = {0     , Name StrCat[pp4,"/layer 4: imag part of relative permittivity"]},
    eps_re_L_5    = {1     , Name StrCat[pp4,"/layer 5: real part of relative permittivity"]},
    eps_im_L_5    = {0     , Name StrCat[pp4,"/layer 5: imag part of relative permittivity"]},
    eps_re_L_6    = {1     , Name StrCat[pp4,"/layer 6: real part of relative permittivity"]},
    eps_im_L_6    = {0     , Name StrCat[pp4,"/layer 6: imag part of relative permittivity"]},
    
    og            = {0          , Name StrCat[pp5,"/0geometrical order [-]"]  , Choices {0="1",1="2"} },
    oi            = {1          , Name StrCat[pp5,"/0interpolation order [-]"], Choices {0="1",1="2"} },
    paramaille    = {6     , Name StrCat[pp5,"/1Number of mesh elements per wavelength [-]"]},
    lc_scat       = {lambda0/(1.5*paramaille)    , Name StrCat[pp5,"/2Scatterer absolute mesh size [nm]"]},
    PML_top       = {lambda0, Name StrCat[pp5,"/4PML top thickness [nm]"]},
    PML_bot       = {lambda0, Name StrCat[pp5,"/5PML bot thickness [nm]"]},
    Nmax          = {2      , Name StrCat[pp5,"/6Number of non specular order to output [-]"]},
    refine_mesh_L_1= {1.5    , Name StrCat[pp5,"/7refine layers/1refine mesh layer 1 [-]"]},
    refine_mesh_L_2= {1.5    , Name StrCat[pp5,"/7refine layers/2refine mesh layer 2 [-]"]},
    refine_mesh_L_3= {1      , Name StrCat[pp5,"/7refine layers/3refine mesh layer 3 [-]"]},
    refine_mesh_L_4= {1      , Name StrCat[pp5,"/7refine layers/4refine mesh layer 4 [-]"]},
    refine_mesh_L_5= {1      , Name StrCat[pp5,"/7refine layers/5refine mesh layer 5 [-]"]},
    refine_mesh_L_6= {1      , Name StrCat[pp5,"/7refine layers/6refine mesh layer 6 [-]"]},
    FlagLinkFacets = {0      , Name StrCat[pp5,"/8FlagLinkFacets? [-]"], Choices {0,1}, Visible 0},
    
    InterpSampling     = { 30   , Name StrCat[pp6,"/0Interpolation grid step [nm]"]},
    Flag_interp_cubic  = { 1    , Name StrCat[pp6,"/1Interpolate on cubic grid?"], Choices {0,1} },
    FlagOutEtotCuts    = { 1    , Name StrCat[pp6,"/2Output Total Electric Field cuts?"] , Choices {0,1} },
    FlagOutHtotCuts    = { 0    , Name StrCat[pp6,"/3Output Total Magnetic Field cuts?"] , Choices {0,1} },
    FlagOutEscaCuts    = { 1    , Name StrCat[pp6,"/4Output Scattered Electric Field cuts?"] , Choices {0,1} },
    FlagOutPoyCut      = { 1    , Name StrCat[pp6,"/5Output Poynting cuts?"] , Choices {0,1} },
    FlagOutEtotFull    = { 0    , Name StrCat[pp6,"/6Total Electric Field Full Output?"] , Choices {0,1} },
    FlagOutEscaFull    = { 0    , Name StrCat[pp6,"/7Scattered Electric Field Full Output?"] , Choices {0,1} },
    FlagOutPoyFull     = { 0    , Name StrCat[pp6,"/8Poynting Full Output?"] , Choices {0,1} }
];

lambda0       = nm * lambda0;
period_x      = nm * period_x;
period_y      = nm * period_y;
thick_L_1     = nm * thick_L_1;
thick_L_2     = nm * thick_L_2;
thick_L_3     = nm * thick_L_3;
thick_L_4     = nm * thick_L_4;
thick_L_5     = nm * thick_L_5;
thick_L_6     = nm * thick_L_6;
rx            = nm * rx;
ry            = nm * ry;
rz            = nm * rz;
lc_scat       = nm * lc_scat;
PML_top       = nm * PML_top;
PML_bot       = nm * PML_bot;
InterpSampling= nm * InterpSampling;

lambda_m = lambda0;
og+=1;
oi+=1;

hh_L_6 = -thick_L_6;
For k In {1:5}
    hh_L~{6-k} = hh_L~{7-k}+thick_L~{7-k};
EndFor
PML_bot_hh = hh_L_6-PML_bot;
PML_top_hh = hh_L_1+thick_L_1;
hh_L_7     = PML_bot_hh;
hh_L_0     = PML_top_hh;
thick_L_7  = PML_bot;
thick_L_0  = PML_top;

theta0 = thetadeg*Pi/180;
phi0   = phideg*Pi/180;
psi0   = psideg*Pi/180;
xsi    = xsideg*Pi/180;

dyc = period_y*Cos[xsi];
dys = period_y*Sin[xsi];

DomainZsizeSca  = PML_top_hh+PML_bot-(hh_L_6-PML_bot);
DomainZsizeTot  = PML_top_hh-hh_L_6;
npts_interpX    = period_x/InterpSampling;
npts_interpY    = period_y/InterpSampling;
npts_checkpoyX  = 50;
npts_checkpoyY  = 50;
npts_interpZSca = DomainZsizeSca/InterpSampling;
npts_interpZTot = DomainZsizeTot/InterpSampling;
