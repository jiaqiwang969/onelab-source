
Include "transfo_common.pro";

Solver.AutoShowLastStep = 0; // don't automatically show the last time step

// CORE

p_Leg_1_L_0=newp; Point(newp) = {-width_Core/2., 0, 0, c_Core};
p_Leg_1_R_0=newp; Point(newp) = {-width_Core/2.+width_Core_Leg, 0, 0, c_Core};

p_Leg_1_L_1=newp; Point(newp) = {-width_Core/2., height_Core/2., 0, c_Core};
p_Leg_1_R_1=newp; Point(newp) = {-width_Core/2.+width_Core_Leg, height_Core/2.-width_Core_Leg, 0, c_Core};


p_Leg_2_L_0=newp; Point(newp) = {width_Core/2.-width_Core_Leg, 0, 0, c_Core};
p_Leg_2_R_0=newp; Point(newp) = {width_Core/2., 0, 0, c_Core};

p_Leg_2_L_1=newp; Point(newp) = {width_Core/2.-width_Core_Leg, height_Core/2.-width_Core_Leg, 0, c_Core};
p_Leg_2_R_1=newp; Point(newp) = {width_Core/2., height_Core/2., 0, c_Core};


l_Core_In[]={};
l_Core_In[]+=newl; Line(newl) = {p_Leg_1_R_0, p_Leg_1_R_1};
l_Core_In[]+=newl; Line(newl) = {p_Leg_1_R_1, p_Leg_2_L_1};
l_Core_In[]+=newl; Line(newl) = {p_Leg_2_L_1, p_Leg_2_L_0};

l_Core_Out[]={};
l_Core_Out[]+=newl; Line(newl) = {p_Leg_2_R_0, p_Leg_2_R_1};
l_Core_Out[]+=newl; Line(newl) = {p_Leg_2_R_1, p_Leg_1_L_1};
l_Core_Out[]+=newl; Line(newl) = {p_Leg_1_L_1, p_Leg_1_L_0};

l_Core_Leg_1_Y0[]={};
l_Core_Leg_1_Y0[]+=newl; Line(newl) = {p_Leg_1_L_0, p_Leg_1_R_0};

l_Core_Leg_2_Y0[]={};
l_Core_Leg_2_Y0[]+=newl; Line(newl) = {p_Leg_2_L_0, p_Leg_2_R_0};


ll_Core=newll; Curve Loop(newll) = {l_Core_In[], l_Core_Leg_2_Y0[], l_Core_Out[], l_Core_Leg_1_Y0[]};
s_Core=news; Plane Surface(news) = {ll_Core};

Physical Surface("CORE", CORE) = {s_Core};

// COIL_1_PLUS

x_[]=Point{p_Leg_1_R_0};
p_Coil_1_p_int_0=newp; Point(newp) = {x_[0]+gap_Core_Coil_1, 0, 0, c_Coil_1};
p_Coil_1_p_ext_0=newp; Point(newp) = {x_[0]+gap_Core_Coil_1+width_Coil_1, 0, 0, c_Coil_1};
p_Coil_1_p_int_1=newp; Point(newp) = {x_[0]+gap_Core_Coil_1, height_Coil_1/2, 0, c_Coil_1};
p_Coil_1_p_ext_1=newp; Point(newp) = {x_[0]+gap_Core_Coil_1+width_Coil_1, height_Coil_1/2, 0, c_Coil_1};

l_Coil_1_p[]={};
l_Coil_1_p[]+=newl; Line(newl) = {p_Coil_1_p_ext_0, p_Coil_1_p_ext_1};
l_Coil_1_p[]+=newl; Line(newl) = {p_Coil_1_p_ext_1, p_Coil_1_p_int_1};
l_Coil_1_p[]+=newl; Line(newl) = {p_Coil_1_p_int_1, p_Coil_1_p_int_0};

l_Coil_1_p_Y0[]={};
l_Coil_1_p_Y0[]+=newl; Line(newl) = {p_Coil_1_p_int_0, p_Coil_1_p_ext_0};

ll_Coil_1_p=newll; Curve Loop(newll) = {l_Coil_1_p[], l_Coil_1_p_Y0[]};
s_Coil_1_p=news; Plane Surface(news) = {ll_Coil_1_p};

Physical Surface("COIL_1_PLUS", COIL_1_PLUS) = {s_Coil_1_p};


// COIL_1_MINUS

x_[]=Point{p_Leg_1_L_0};
p_Coil_1_m_int_0=newp; Point(newp) = {x_[0]-gap_Core_Coil_1, 0, 0, c_Coil_1};
p_Coil_1_m_ext_0=newp; Point(newp) = {x_[0]-(gap_Core_Coil_1+width_Coil_1), 0, 0, c_Coil_1};
p_Coil_1_m_int_1=newp; Point(newp) = {x_[0]-gap_Core_Coil_1, height_Coil_1/2, 0, c_Coil_1};
p_Coil_1_m_ext_1=newp; Point(newp) = {x_[0]-(gap_Core_Coil_1+width_Coil_1), height_Coil_1/2, 0, c_Coil_1};

l_Coil_1_m[]={};
l_Coil_1_m[]+=newl; Line(newl) = {p_Coil_1_m_ext_0, p_Coil_1_m_ext_1};
l_Coil_1_m[]+=newl; Line(newl) = {p_Coil_1_m_ext_1, p_Coil_1_m_int_1};
l_Coil_1_m[]+=newl; Line(newl) = {p_Coil_1_m_int_1, p_Coil_1_m_int_0};

l_Coil_1_m_Y0[]={};
l_Coil_1_m_Y0[]+=newl; Line(newl) = {p_Coil_1_m_int_0, p_Coil_1_m_ext_0};

ll_Coil_1_m=newll; Curve Loop(newll) = {l_Coil_1_m[], l_Coil_1_m_Y0[]};
s_Coil_1_m=news; Plane Surface(news) = {ll_Coil_1_m};

Physical Surface("COIL_1_MINUS", COIL_1_MINUS) = {s_Coil_1_m};


// COIL_2_PLUS

x_[]=Point{p_Leg_2_R_0};
p_Coil_2_p_int_0=newp; Point(newp) = {x_[0]+gap_Core_Coil_2, 0, 0, c_Coil_2};
p_Coil_2_p_ext_0=newp; Point(newp) = {x_[0]+gap_Core_Coil_2+width_Coil_2, 0, 0, c_Coil_2};
p_Coil_2_p_int_1=newp; Point(newp) = {x_[0]+gap_Core_Coil_2, height_Coil_2/2, 0, c_Coil_2};
p_Coil_2_p_ext_1=newp; Point(newp) = {x_[0]+gap_Core_Coil_2+width_Coil_2, height_Coil_2/2, 0, c_Coil_2};

l_Coil_2_p[]={};
l_Coil_2_p[]+=newl; Line(newl) = {p_Coil_2_p_ext_0, p_Coil_2_p_ext_1};
l_Coil_2_p[]+=newl; Line(newl) = {p_Coil_2_p_ext_1, p_Coil_2_p_int_1};
l_Coil_2_p[]+=newl; Line(newl) = {p_Coil_2_p_int_1, p_Coil_2_p_int_0};

l_Coil_2_p_Y0[]={};
l_Coil_2_p_Y0[]+=newl; Line(newl) = {p_Coil_2_p_int_0, p_Coil_2_p_ext_0};

ll_Coil_2_p=newll; Curve Loop(newll) = {l_Coil_2_p[], l_Coil_2_p_Y0[]};
s_Coil_2_p=news; Plane Surface(news) = {ll_Coil_2_p};

Physical Surface("COIL_2_PLUS", COIL_2_PLUS) = {s_Coil_2_p};


// COIL_2_MINUS

x_[]=Point{p_Leg_2_L_0};
p_Coil_2_m_int_0=newp; Point(newp) = {x_[0]-gap_Core_Coil_2, 0, 0, c_Coil_2};
p_Coil_2_m_ext_0=newp; Point(newp) = {x_[0]-(gap_Core_Coil_2+width_Coil_2), 0, 0, c_Coil_2};
p_Coil_2_m_int_1=newp; Point(newp) = {x_[0]-gap_Core_Coil_2, height_Coil_2/2, 0, c_Coil_2};
p_Coil_2_m_ext_1=newp; Point(newp) = {x_[0]-(gap_Core_Coil_2+width_Coil_2), height_Coil_2/2, 0, c_Coil_2};

l_Coil_2_m[]={};
l_Coil_2_m[]+=newl; Line(newl) = {p_Coil_2_m_ext_0, p_Coil_2_m_ext_1};
l_Coil_2_m[]+=newl; Line(newl) = {p_Coil_2_m_ext_1, p_Coil_2_m_int_1};
l_Coil_2_m[]+=newl; Line(newl) = {p_Coil_2_m_int_1, p_Coil_2_m_int_0};

l_Coil_2_m_Y0[]={};
l_Coil_2_m_Y0[]+=newl; Line(newl) = {p_Coil_2_m_int_0, p_Coil_2_m_ext_0};

ll_Coil_2_m=newll; Curve Loop(newll) = {l_Coil_2_m[], l_Coil_2_m_Y0[]};
s_Coil_2_m=news; Plane Surface(news) = {ll_Coil_2_m};

Physical Surface("COIL_2_MINUS", COIL_2_MINUS) = {s_Coil_2_m};


// AIR_WINDOW

l_Air_Window_Y0[]={};
l_Air_Window_Y0[]+=newl; Line(newl) = {p_Leg_1_R_0, p_Coil_1_p_int_0};
l_Air_Window_Y0[]+=newl; Line(newl) = {p_Coil_1_p_ext_0, p_Coil_2_m_ext_0};
l_Air_Window_Y0[]+=newl; Line(newl) = {p_Coil_2_m_int_0, p_Leg_2_L_0};

ll_Air_Window=newll; Curve Loop(newll) = {-l_Core_In[], -l_Coil_1_p[], l_Coil_2_m[], l_Air_Window_Y0[]};
s_Air_Window=news; Plane Surface(news) = {ll_Air_Window};

Physical Surface("AIR_WINDOW", AIR_WINDOW) = {s_Air_Window};


// AIR_EXT

x_[]=Point{p_Leg_2_R_1};
p_Air_Ext_1_R_0=newp; Point(newp) = {x_[0]+gap_Core_Box_X, 0, 0, c_Box};
p_Air_Ext_1_R_1=newp; Point(newp) = {x_[0]+gap_Core_Box_X, x_[1]+gap_Core_Box_Y, 0, c_Box};

x_[]=Point{p_Leg_1_L_1};
p_Air_Ext_1_L_0=newp; Point(newp) = {x_[0]-gap_Core_Box_X, 0, 0, c_Box};
p_Air_Ext_1_L_1=newp; Point(newp) = {x_[0]-gap_Core_Box_X, x_[1]+gap_Core_Box_Y, 0, c_Box};

l_Air_Ext[]={};
l_Air_Ext[]+=newl; Line(newl) = {p_Air_Ext_1_R_0, p_Air_Ext_1_R_1};
l_Air_Ext[]+=newl; Line(newl) = {p_Air_Ext_1_R_1, p_Air_Ext_1_L_1};
l_Air_Ext[]+=newl; Line(newl) = {p_Air_Ext_1_L_1, p_Air_Ext_1_L_0};


l_Air_Ext_Y0[]={};
l_Air_Ext_Y0[]+=newl; Line(newl) = {p_Leg_2_R_0, p_Coil_2_p_int_0};
l_Air_Ext_Y0[]+=newl; Line(newl) = {p_Coil_2_p_ext_0, p_Air_Ext_1_R_0};

l_Air_Ext_Y0[]+=newl; Line(newl) = {p_Air_Ext_1_L_0, p_Coil_1_m_ext_0};
l_Air_Ext_Y0[]+=newl; Line(newl) = {p_Coil_1_m_int_0, p_Leg_1_L_0};


ll_Air_Ext=newll; Curve Loop(newll) = {-l_Core_Out[], -l_Coil_2_p[], l_Coil_1_m[], l_Air_Ext[], l_Air_Ext_Y0[]};
s_Air_Ext=news; Plane Surface(news) = {ll_Air_Ext};

Physical Surface("AIR_EXT", AIR_EXT) = {s_Air_Ext};
Physical Curve("SUR_AIR_EXT", SUR_AIR_EXT) = {l_Air_Ext[]};
