// Copyright (C) 2020 Guillaume Demésy
//
// This file is part of the model grating2D.pro.
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

myDir = "res2D_conical/";
DefineConstant[
  lambda0   = {lambda_min , Min lambda_min, Max lambda_max, Step (lambda_max-lambda_min)/(nb_lambdas-1), Name StrCat[pp2, "1wavelength [nm]"] , Loop 1, Highlight Str[colorpp2],Graph "200000200020", ServerAction "Reset GetDP/T0, GetDP/R0, GetDP/Lambda_step, GetDP/Omegaal absorption"},
  phi_deg = {  45 , Name StrCat[pp2, "3incident plane wave angle (phi) [deg]"] , Highlight Str[colorpp2], Closed close_menu},
  psi_deg = {  45 , Name StrCat[pp2, "4incident plane wave angle (psi) [deg]"] , Highlight Str[colorpp2], Closed close_menu}
];
lambda0    = lambda0    * nm;
lambda_min = lambda_min * nm;
lambda_max = lambda_max * nm;
A          = 1.;
nb_orders  = 2;

Group {
  // Boundaries
  SurfBlochLeft   = Region[SURF_BLOCH_X_LEFT];
  SurfBlochRight  = Region[SURF_BLOCH_X_RIGHT];
  SurfCutSubs1    = Region[SURF_INTEG_SUB1];
  SurfCutSuper1   = Region[SURF_INTEG_SUP1];
  SurfCutSubs2    = Region[SURF_INTEG_SUB2];
  SurfCutSuper2   = Region[SURF_INTEG_SUP2];

  // Domains
  pmlbot          = Region[PMLBOT];
  sub             = Region[{SUB,SURF_INTEG_SUB1,SURF_INTEG_SUB2}];
  layer_dep       = Region[LAYERDEP];
  rod_out         = Region[RODOUT];
  For i In {0:N_rods-1:1}
    rod~{i} = Region[{ROD~{i}}];
    rods   += Region[{ROD~{i}}];
  EndFor
  layer_cov       = Region[LAYERCOV];
  sup             = Region[{SUP,SURF_INTEG_SUP1,SURF_INTEG_SUP2}];
  pmltop          = Region[PMLTOP];
  Omega_source    = Region[{layer_dep,rod_out,rods,layer_cov}];
  Omega_nosource  = Region[{pmltop,pmlbot,sup,sub}];
  Omega           = Region[{Omega_source,Omega_nosource}];
  Omega_top       = Region[{layer_dep,rod_out,rods,layer_cov,sup}];
  Omega_bot       = Region[{sub}];
  Omega_pml       = Region[{pmltop,pmlbot}];

  Plot_domain     = Region[{Omega,-pmltop,-pmlbot}];
  Plot_bnd        = Region[SURF_PLOT];
  Print_point     = Region[PRINT_POINT];

  Omega_DefTr     = Region[{SUB,SUP}];  
  Gama = Region[{SurfCutSubs1,SurfCutSuper1}];
  Tr   = ElementsOf[Omega_DefTr, ConnectedTo Gama];
}

Function{
  I[] = Complex[0.0,1.0];
  EZ[]= Vector[0,0,1];
  bndCol[Plot_bnd] = Complex[1,1];


  epsr_re_interp_mat_1[]  = InterpolationLinear[$1]{ListAlt[lambdamat_1 ,epsr_mat_re_1 ]};
  epsr_im_interp_mat_1[]  = InterpolationLinear[$1]{ListAlt[lambdamat_1 ,epsr_mat_im_1 ]};
  epsr_re_interp_mat_2[]  = InterpolationLinear[$1]{ListAlt[lambdamat_2 ,epsr_mat_re_2 ]};
  epsr_im_interp_mat_2[]  = InterpolationLinear[$1]{ListAlt[lambdamat_2 ,epsr_mat_im_2 ]};
  epsr_re_interp_mat_3[]  = InterpolationLinear[$1]{ListAlt[lambdamat_3 ,epsr_mat_re_3 ]};
  epsr_im_interp_mat_3[]  = InterpolationLinear[$1]{ListAlt[lambdamat_3 ,epsr_mat_im_3 ]};
  epsr_re_interp_mat_4[]  = InterpolationLinear[$1]{ListAlt[lambdamat_4 ,epsr_mat_re_4 ]};
  epsr_im_interp_mat_4[]  = InterpolationLinear[$1]{ListAlt[lambdamat_4 ,epsr_mat_im_4 ]};
  epsr_re_interp_mat_5[]  = InterpolationLinear[$1]{ListAlt[lambdamat_5 ,epsr_mat_re_5 ]};
  epsr_im_interp_mat_5[]  = InterpolationLinear[$1]{ListAlt[lambdamat_5 ,epsr_mat_im_5 ]};
  epsr_re_interp_mat_6[]  = InterpolationLinear[$1]{ListAlt[lambdamat_6 ,epsr_mat_re_6 ]};
  epsr_im_interp_mat_6[]  = InterpolationLinear[$1]{ListAlt[lambdamat_6 ,epsr_mat_im_6 ]};
  epsr_re_interp_mat_7[]  = InterpolationLinear[$1]{ListAlt[lambdamat_7 ,epsr_mat_re_7 ]};
  epsr_im_interp_mat_7[]  = InterpolationLinear[$1]{ListAlt[lambdamat_7 ,epsr_mat_im_7 ]};
  epsr_re_interp_mat_8[]  = InterpolationLinear[$1]{ListAlt[lambdamat_8 ,epsr_mat_re_8 ]};
  epsr_im_interp_mat_8[]  = InterpolationLinear[$1]{ListAlt[lambdamat_8 ,epsr_mat_im_8 ]};
  epsr_re_interp_mat_9[]  = InterpolationLinear[$1]{ListAlt[lambdamat_9 ,epsr_mat_re_9 ]};
  epsr_im_interp_mat_9[]  = InterpolationLinear[$1]{ListAlt[lambdamat_9 ,epsr_mat_im_9 ]};
  epsr_re_interp_mat_10[] = InterpolationLinear[$1]{ListAlt[lambdamat_10,epsr_mat_re_10]};
  epsr_im_interp_mat_10[] = InterpolationLinear[$1]{ListAlt[lambdamat_10,epsr_mat_im_10]};
  epsr_re_interp_mat_11[] = InterpolationLinear[$1]{ListAlt[lambdamat_11,epsr_mat_re_11]};
  epsr_im_interp_mat_11[] = InterpolationLinear[$1]{ListAlt[lambdamat_11,epsr_mat_im_11]};
  epsr_re_interp_mat_12[] = InterpolationLinear[$1]{ListAlt[lambdamat_12,epsr_mat_re_12]};
  epsr_im_interp_mat_12[] = InterpolationLinear[$1]{ListAlt[lambdamat_12,epsr_mat_im_12]};
  epsr_re_interp_mat_13[] = InterpolationLinear[$1]{ListAlt[lambdamat_13,epsr_mat_re_13]};
  epsr_im_interp_mat_13[] = InterpolationLinear[$1]{ListAlt[lambdamat_13,epsr_mat_im_13]};
  epsr_re_interp_mat_14[] = InterpolationLinear[$1]{ListAlt[lambdamat_14,epsr_mat_re_14]};
  epsr_im_interp_mat_14[] = InterpolationLinear[$1]{ListAlt[lambdamat_14,epsr_mat_im_14]};
  epsr_re_interp_mat_15[] = InterpolationLinear[$1]{ListAlt[lambdamat_15,epsr_mat_re_15]};
  epsr_im_interp_mat_15[] = InterpolationLinear[$1]{ListAlt[lambdamat_15,epsr_mat_im_15]};
  epsr_re_interp_mat_16[] = InterpolationLinear[$1]{ListAlt[lambdamat_16,epsr_mat_re_16]};
  epsr_im_interp_mat_16[] = InterpolationLinear[$1]{ListAlt[lambdamat_16,epsr_mat_im_16]};

  For i In {1:5}
    For j In {1:nb_available_materials}
      If(flag_mat~{i}==j-1)
        epsr_re_dom~{i}[] = epsr_re_interp_mat~{j}[lambda0/nm*1e-9];
        epsr_im_dom~{i}[] = epsr_im_interp_mat~{j}[lambda0/nm*1e-9];
      EndIf
    EndFor
  EndFor
  For k In {0:nb_available_lossless_materials-1}
    If(flag_mat_6==lossless_material_list(k)-1)
      epsr_re_dom_6[] = epsr_re_interp_mat~{lossless_material_list(k)}[lambda0/nm*1e-9];
      epsr_im_dom_6[] = epsr_im_interp_mat~{lossless_material_list(k)}[lambda0/nm*1e-9];
    EndIf
  EndFor

  epsr2_re[]       = epsr_re_dom_1[];
  epsr2_im[]       = epsr_im_dom_1[];
  epsr_layer_dep_re[] = epsr_re_dom_2[];
  epsr_layer_dep_im[] = epsr_im_dom_2[];
  epsr_rod_out_re[]   = epsr_re_dom_3[];
  epsr_rod_out_im[]   = epsr_im_dom_3[];
  epsr_rods_re[]      = epsr_re_dom_4[];
  epsr_rods_im[]      = epsr_im_dom_4[];
  epsr_layer_cov_re[] = epsr_re_dom_5[];
  epsr_layer_cov_im[] = epsr_im_dom_5[];
  epsr1_re[]       = epsr_re_dom_6[];
  epsr1_im[]       = epsr_im_dom_6[];
  
  om0      = 2.*Pi*cel/lambda0;
  k0       = 2.*Pi/lambda0;
  epsr1[]  = Complex[epsr1_re[],epsr1_im[]];
  epsr2[]  = Complex[epsr2_re[],epsr2_im[]];
  n1[]     = Sqrt[epsr1[]];
  n2[]     = Sqrt[epsr2[]];
  k1norm[] = k0*n1[];
  k2norm[] = k0*n2[];

  phi0    =  phi_deg   * deg2rad;
  theta0  =  theta_deg * deg2rad;
  psi0    =  psi_deg   * deg2rad;
  // Ae      =  1/Sqrt[ep0/mu0];
  Ae      = 1;

  alpha[] = -k0*n1[]*Sin[theta0]*Sin[phi0];
  beta1[] = -k0*n1[]*Cos[theta0];
  gamma[] = -k0*n1[]*Sin[theta0]*Cos[phi0];
  beta2[] = -Sqrt[k0^2*epsr2[]-alpha[]^2-gamma[]^2];
  k1[]  = Vector[alpha[], beta1[], gamma[]];
  k2[]  = Vector[alpha[], beta2[], gamma[]];
  k1r[] = Vector[alpha[],-beta1[], gamma[]];

  rs[] = (beta1[]-beta2[])/(beta1[]+beta2[]);
  ts[] =    2.*beta1[] /(beta1[]+beta2[]);
  rp[] = (beta1[]*epsr2[]-beta2[]*epsr1[])/(beta1[]*epsr2[]+beta2[]*epsr1[]);
  tp[] =            (2.*beta1[]*epsr2[])/(beta1[]*epsr2[]+beta2[]*epsr1[]);

  spol[]    = Vector[-Cos[phi0],0,Sin[phi0]];
  AmplEis[] =      spol[];
  AmplErs[] = rs[]*spol[];
  AmplEts[] = ts[]*spol[];
  AmplHis[] = Sqrt[ep0*epsr1[]/mu0]     *spol[];
  AmplHrs[] = Sqrt[ep0*epsr1[]/mu0]*rp[]*spol[];
  AmplHts[] = Sqrt[ep0*epsr1[]/mu0]*tp[]*spol[];

  Eis[]     = AmplEis[] * Exp[I[]*k1[] *XYZ[]];
  Ers[]     = AmplErs[] * Exp[I[]*k1r[]*XYZ[]];
  Ets[]     = AmplEts[] * Exp[I[]*k2[] *XYZ[]];
  His[]     = AmplHis[] * Exp[I[]*k1[] *XYZ[]];
  Hrs[]     = AmplHrs[] * Exp[I[]*k1r[]*XYZ[]];
  Hts[]     = AmplHts[] * Exp[I[]*k2[] *XYZ[]];
  Eip[]     = -1/(om0*ep0*epsr1[]) * k1[]  /\ His[];
  Erp[]     = -1/(om0*ep0*epsr1[]) * k1r[] /\ Hrs[];
  Etp[]     = -1/(om0*ep0*epsr2[]) * k2[]  /\ Hts[];

  BlochX_phase_shift[] = Exp[I[]*alpha[]*d];

  For i In {0:2*nb_orders}
    alpha~{i}[]      = -CompX[k1[]]  + 2*Pi/d*(i-nb_orders);
    expialphax~{i}[] = Exp[I[]*alpha~{i}[]*X[]];
  EndFor
  For i In {0:2*nb_orders}
    beta_super~{i}[] = Sqrt[k0^2*epsr1[] - alpha~{i}[]^2 - gamma[]^2];
    beta_subs~{i}[]  = Sqrt[k0^2*epsr2[] - alpha~{i}[]^2 - gamma[]^2];
  EndFor

  a_pml        = 1.;
  b_pml        = 1.;
  sx           = 1.;
  sy[]         = Complex[a_pml,b_pml]; 
  sz           = 1.;
  PML_Tensor[] = TensorDiag[sz*sy[]/sx,sx*sz/sy[],sx*sy[]/sz];

  epsilonr[sup]           = epsr1[] * TensorDiag[1,1,1];
  epsilonr[sub]           = epsr2[] * TensorDiag[1,1,1];
  epsilonr[layer_dep]     = Complex[epsr_layer_dep_re[],epsr_layer_dep_im[]] * TensorDiag[1,1,1];
  epsilonr[rod_out]       = Complex[epsr_rod_out_re[],epsr_rod_out_im[]    ] * TensorDiag[1,1,1];
  If (Flag_aniso==0)
    epsilonr[rods]        = Complex[epsr_rods_re[],epsr_rods_im[]] * TensorDiag[1,1,1];
  Else
    epsilonr[rods]        = Tensor[Complex[epsr_custom_anisoXX_re, epsr_custom_anisoXX_im],
           Complex[epsr_custom_anisoXY_re, epsr_custom_anisoXY_im],
           0,
           Complex[epsr_custom_anisoXY_re,-epsr_custom_anisoXX_im],
           Complex[epsr_custom_anisoYY_re,epsr_custom_anisoYY_im],
           0,
           0,
           0,
           Complex[epsr_custom_anisoZZ_re,epsr_custom_anisoZZ_im]
           ];
  EndIf
  epsilonr[layer_cov]       = Complex[epsr_layer_cov_re[],epsr_layer_cov_im[]] * TensorDiag[1,1,1];
  epsilonr[pmltop]          = epsr1_re[]*PML_Tensor[];
  epsilonr[pmlbot]          = epsr2_re[]*PML_Tensor[]; 
  
  epsilonr_annex[sub]       = epsr2[] * TensorDiag[1,1,1];
  epsilonr_annex[sup]       = epsr1[] * TensorDiag[1,1,1];
  epsilonr_annex[layer_dep] = epsr1[] * TensorDiag[1,1,1];
  epsilonr_annex[rod_out]   = epsr1[] * TensorDiag[1,1,1];
  epsilonr_annex[rods]      = epsr1[] * TensorDiag[1,1,1];
  epsilonr_annex[layer_cov] = epsr1[] * TensorDiag[1,1,1];
  epsilonr_annex[pmltop]    = epsr1_re[]*PML_Tensor[];
  epsilonr_annex[pmlbot]    = epsr2_re[]*PML_Tensor[];

  mur[pmltop]     = PML_Tensor[];
  mur[pmlbot]     = PML_Tensor[];
  mur[sub]        = TensorDiag[1,1,1];
  mur[sup]        = TensorDiag[1,1,1];
  mur[layer_dep]  = TensorDiag[1,1,1];
  mur[rods]       = TensorDiag[1,1,1];
  mur[rod_out]    = TensorDiag[1,1,1];
  mur[layer_cov]  = TensorDiag[1,1,1];

  Ei[] = Ae*(Cos[psi0]*Eip[]-Sin[psi0]*Eis[]);
  Er[] = Ae*(Cos[psi0]*Erp[]-Sin[psi0]*Ers[]);
  Et[] = Ae*(Cos[psi0]*Etp[]-Sin[psi0]*Ets[]);
  Hi[] =  1/(om0*mu0*mur[]) * k1[]  /\ Ei[];
  Hr[] =  1/(om0*mu0*mur[]) * k1r[] /\ Er[];
  Ht[] =  1/(om0*mu0*mur[]) * k2[]  /\ Et[];

  
  E1[pmltop]     = Ei[]+Er[];
  E1[Omega_top]  = Ei[]+Er[];
  E1[Omega_bot]  = Et[];  
  E1[pmlbot]     = Et[];  
  E1d[pmltop]    = Er[];
  E1d[Omega_top] = Er[];
  E1d[Omega_bot] = Et[];
  E1d[pmlbot]    = Et[];
  E1t[]  = Vector[CompX[E1[]], CompY[E1[]]  ,    0        ];
  E1l[]  = Vector[     0     ,      0       , CompZ[E1[]] ];
  E1dt[] = Vector[CompX[E1d[]], CompY[E1d[]],    0        ];
  E1dl[] = Vector[     0      ,      0      , CompZ[E1d[]]];


  H1[pmltop]     = Hi[]+Hr[];
  H1[Omega_top]  = Hi[]+Hr[];
  H1[Omega_bot]  = Ht[];
  H1[pmlbot]     = Ht[];
  H1d[Omega_top] = Hr[];
  H1d[Omega_bot] = Ht[];
  Pinc[]         = 0.5*Ae^2*Sqrt[epsr1[]*ep0/mu0] * Cos[theta0];
}

Constraint {
  {Name Bloch;
    Case {
      { Region SurfBlochRight; Type LinkCplx ; RegionRef SurfBlochLeft;
               Coefficient BlochX_phase_shift[]; Function Vector[$X-d,$Y,$Z] ;
      }
    }
  }
}

Jacobian {
  { Name JVol ;
    Case { 
      { Region All ; Jacobian Vol ; }
    }
  }
  { Name JSur ;
    Case {
      { Region All ; Jacobian Sur ; }
    }
  }
}

Integration {
  If (Flag_o2_geom==0)
    { Name Int_1 ;
      Case { 
        { Type Gauss ;
          Case { 
            { GeoElement Point    ; NumberOfPoints  1 ; }
            { GeoElement Line     ; NumberOfPoints  4 ; }
            { GeoElement Triangle ; NumberOfPoints  12 ; }
          }
        }
      }
    }
  EndIf
  If (Flag_o2_geom==1)
    { Name Int_1 ;
      Case { 
        { Type Gauss ;
          Case { 
            { GeoElement Point     ; NumberOfPoints  1 ; }
            { GeoElement Line2     ; NumberOfPoints  4 ; }
            { GeoElement Triangle2 ; NumberOfPoints  12 ; }
          }
        }
      }
    }
  EndIf
}

FunctionSpace {
  { Name Hcurl; Type Form1;
    BasisFunction {
      { Name se ; NameOfCoef he ; Function BF_Edge     ; Support Omega; Entity EdgesOf[All]; }
      { Name se2; NameOfCoef he2; Function BF_Edge_2E  ; Support Omega; Entity EdgesOf[All]; }
      If (Flag_o2_interp==1)
        { Name se3; NameOfCoef he3; Function BF_Edge_3F_b; Support Omega; Entity FacetsOf[All]; }
        { Name se4; NameOfCoef he4; Function BF_Edge_3F_c; Support Omega; Entity FacetsOf[All]; }
        { Name se5; NameOfCoef he5; Function BF_Edge_4E  ; Support Omega; Entity EdgesOf[All]; }
      EndIf
    }
    Constraint {
      { NameOfCoef he ; EntityType EdgesOf ; NameOfConstraint Bloch; }
      { NameOfCoef he2; EntityType EdgesOf ; NameOfConstraint Bloch; }
      If (Flag_o2_interp==1)
        { NameOfCoef he5; EntityType EdgesOf ; NameOfConstraint Bloch; }
      EndIf
    }
  }

  { Name Hgrad_perp; Type Form1P;
    BasisFunction {
      { Name sn ; NameOfCoef hn ; Function BF_PerpendicularEdge   ; Support Omega; Entity NodesOf[All]; }
      { Name sn2; NameOfCoef hn2; Function BF_PerpendicularEdge_2E; Support Omega; Entity EdgesOf[All]; }
      If (Flag_o2_interp==1)
        { Name sn3; NameOfCoef hn3; Function BF_PerpendicularEdge_2F; Support Omega; Entity FacetsOf[All]; }
        { Name sn4; NameOfCoef hn4; Function BF_PerpendicularEdge_3E; Support Omega; Entity EdgesOf[All]; }
        { Name sn5; NameOfCoef hn5; Function BF_PerpendicularEdge_3F; Support Omega; Entity FacetsOf[All]; }
      EndIf
      // BF_PerpendicularEdge_2F, 1., QUA|HEX|PRI, 0 },
      // BF_PerpendicularEdge_2V, 1., QUA|HEX, 0 },
      // BF_PerpendicularEdge_3E, 2., ALL, 0 },
      // BF_PerpendicularEdge_3F, 2., TRI|QUA|TET|HEX|PRI, 0 },
      // BF_PerpendicularEdge_3V, 2., HEX|PRI, 0 },
    }
    Constraint {
      { NameOfCoef hn ; EntityType NodesOf ; NameOfConstraint Bloch; }
      { NameOfCoef hn2; EntityType EdgesOf ; NameOfConstraint Bloch; }
      If (Flag_o2_interp==1)
        { NameOfCoef hn4; EntityType EdgesOf ; NameOfConstraint Bloch; }
      EndIf
    }
  }



  { Name L2_lambda; Type Form0;
    BasisFunction{
      { Name ln ; NameOfCoef lambda_n ; Function BF_Node   ; Support Region[Gama]; Entity NodesOf[All];}
      // { Name ln2; NameOfCoef lambda_n2; Function BF_Node_2E; Support Region[Gama]; Entity EdgesOf[All];}
    }
  }
}

Formulation {
  {Name helmoltz_conical; Type FemEquation;
    Quantity {
      { Name Et; Type Local; NameOfSpace Hcurl; }
      { Name El; Type Local; NameOfSpace Hgrad_perp; }
      { Name uy; Type Local; NameOfSpace L2_lambda;}
    }
    Equation {
      Galerkin {[ 1/mur[] * Dof{d Et}                          , {d Et}     ] ; In Omega; Integration Int_1; Jacobian JVol; }
      Galerkin {[ 1/mur[] * Dof{d El}                          , {d Et}     ] ; In Omega; Integration Int_1; Jacobian JVol; }
      Galerkin {[ 1/mur[] * Dof{d El}                          , {d El}     ] ; In Omega; Integration Int_1; Jacobian JVol; }
      Galerkin {[ 1/mur[] * Dof{d Et}                          , {d El}     ] ; In Omega; Integration Int_1; Jacobian JVol; }
      Galerkin {[-Complex[0,gamma[]]  /mur[] * Dof{d Et}       , EZ[]*^{Et} ] ; In Omega; Integration Int_1; Jacobian JVol; }
      Galerkin {[-Complex[0,gamma[]]  /mur[] * Dof{d El}       , EZ[]*^{Et} ] ; In Omega; Integration Int_1; Jacobian JVol; }
      Galerkin {[ Complex[0,gamma[]]  /mur[] * (EZ[]*^Dof{Et}) , {d El}     ] ; In Omega; Integration Int_1; Jacobian JVol; }
      Galerkin {[ Complex[0,gamma[]]  /mur[] * (EZ[]*^Dof{Et}) , {d Et}     ] ; In Omega; Integration Int_1; Jacobian JVol; }
      Galerkin {[           gamma[]^2 /mur[] * (EZ[]*^Dof{Et}) , EZ[]*^{Et} ] ; In Omega; Integration Int_1; Jacobian JVol; }
      
      Galerkin {[ -k0^2 * epsilonr[] * Dof{Et}                 , {Et}       ] ; In Omega; Integration Int_1; Jacobian JVol; }
      Galerkin {[ -k0^2 * epsilonr[] * Dof{Et}                 , {El}       ] ; In Omega; Integration Int_1; Jacobian JVol; }
      Galerkin {[ -k0^2 * epsilonr[] * Dof{El}                 , {Et}       ] ; In Omega; Integration Int_1; Jacobian JVol; }
      Galerkin {[ -k0^2 * epsilonr[] * Dof{El}                 , {El}       ] ; In Omega; Integration Int_1; Jacobian JVol; }

      Galerkin {[ -k0^2 * (epsilonr[]-epsilonr_annex[]) * E1t[] , {Et}      ] ; In Omega_source; Integration Int_1; Jacobian JVol; }
      Galerkin {[ -k0^2 * (epsilonr[]-epsilonr_annex[]) * E1t[] , {El}      ] ; In Omega_source; Integration Int_1; Jacobian JVol; }
      Galerkin {[ -k0^2 * (epsilonr[]-epsilonr_annex[]) * E1l[] , {Et}      ] ; In Omega_source; Integration Int_1; Jacobian JVol; }
      Galerkin {[ -k0^2 * (epsilonr[]-epsilonr_annex[]) * E1l[] , {El}      ] ; In Omega_source; Integration Int_1; Jacobian JVol; }

      Galerkin{ [                         Dof{uy}   , {uy}]; In Gama; Jacobian JSur; Integration Int_1;}
      Galerkin{ [ Trace[Dof{Et}, Tr]*Vector[0,-1,0] , {uy}]; In Gama; Jacobian JSur; Integration Int_1;}
    }
  }
}

Resolution {
  { Name helmoltz_scalar;
    System {
      { Name S; NameOfFormulation helmoltz_conical; Type ComplexValue; }
    }
    Operation {
      CreateDir[Str[myDir]];
      Printf["%f",lambda0/nm];
      Generate[S];
      Solve[S];
    }
  }
}

PostProcessing {
  { Name postpro_energy; NameOfFormulation helmoltz_conical;
    Quantity {
      { Name lambda_step   ; Value { Local { [ lambda0/nm ]; In Omega ; Jacobian JVol; } } }
      { Name E2d    ; Value { Local { [ {Et}+{El}      ]; In Omega; Jacobian JVol; } } }
      { Name Etot   ; Value { Local { [ {Et}+{El}+E1[] ]; In Omega; Jacobian JVol; } } }
      { Name boundary    ; Value { Local { [ bndCol[]    ]; In Plot_bnd ; Jacobian JVol ; } } }
      For i In {-nb_plot_periods:nb_plot_periods}
        { Name Etot~{i}     ; Value { Local { [ ({Et}+{El}+E1[])*Exp[I[]*i*CompX[k1[]]*d] ]; In Omega; Jacobian JVol; } } }
      EndFor
      { Name H2d    ; Value { Local { [ (-I[]/(mur[]*mu0*om0))*({Curl Et}+{Curl El}) ]; In Omega; Jacobian JVol; } } }
      { Name Htot   ; Value { Local { [ (-I[]/(mur[]*mu0*om0))*({Curl Et}+{Curl El})+H1[] ]; In Omega; Jacobian JVol; } } }
      { Name E2d_z  ; Value { Local { [ CompZ[{Et}+{El}] ]; In Omega; Jacobian JVol; } } }
      { Name E2d_t  ; Value { Local { [    {Et}          ]; In Omega; Jacobian JVol; } } }
      { Name Etot_z ; Value { Local { [ CompZ[{Et}+{El}+E1[]] ]; In Omega; Jacobian JVol; } } }
      { Name Htot_z ; Value { Local { [ CompZ[(-I[]/(mur[]*mu0*om0))*({Curl Et}+{Curl El})+H1[]] ]; In Omega; Jacobian JVol; } } }
      { Name H2d_z  ; Value { Local { [ CompZ[(-I[]/(mur[]*mu0*om0))*({Curl Et}+{Curl El})] ]; In Omega; Jacobian JVol; } } }

      { Name E1   ; Value { Local { [ E1[]  ]; In Omega; Jacobian JVol; } } }
      { Name E1t  ; Value { Local { [ E1t[] ]; In Omega; Jacobian JVol; } } }
      { Name E1l  ; Value { Local { [ E1l[] ]; In Omega; Jacobian JVol; } } }

      { Name Q_subs      ; Value{ Integral{ [  0.5 * ep0*om0*Fabs[epsr2_im[]      ]        * ( SquNorm[{Et}+{El}+E1[]] ) / (Pinc[]*d) ] ; In sub         ; Integration Int_1 ; Jacobian JVol ; } } }
      { Name Q_rod_out   ; Value{ Integral{ [  0.5 * ep0*om0*Fabs[epsr_rod_out_im[]]       * ( SquNorm[{Et}+{El}+E1[]] ) / (Pinc[]*d) ] ; In rod_out   ; Integration Int_1 ; Jacobian JVol ; } } }
      { Name Q_layer_dep ; Value{ Integral{ [  0.5 * ep0*om0*Fabs[epsr_layer_dep_im[]]     * ( SquNorm[{Et}+{El}+E1[]] ) / (Pinc[]*d) ] ; In layer_dep   ; Integration Int_1 ; Jacobian JVol ; } } }
      { Name Q_layer_cov ; Value{ Integral{ [  0.5 * ep0*om0*Fabs[epsr_layer_cov_im[]]     * ( SquNorm[{Et}+{El}+E1[]] ) / (Pinc[]*d) ] ; In layer_cov   ; Integration Int_1 ; Jacobian JVol ; } } }
      { Name Q_tot       ; Value{ Integral{ [  0.5 * ep0*om0*Fabs[Im[CompZZ[epsilonr[]]]]  * ( SquNorm[{Et}+{El}+E1[]] ) / (Pinc[]*d) ] ; In Plot_domain ; Integration Int_1 ; Jacobian JVol ; } } }

      For i In {0:2*nb_orders}
        { Name int_x_t~{i} ; Value{ Integral{ [ CompX[{Et}+E1t[] ] * expialphax~{i}[]/d ] ; In SurfCutSubs1  ; Integration Int_1 ; Jacobian JSur ; } } }
        { Name int_y_t~{i} ; Value{ Integral{ [ ({uy}+CompY[E1t[]])* expialphax~{i}[]/d ] ; In SurfCutSubs1  ; Integration Int_1 ; Jacobian JSur ; } } }
        { Name int_z_t~{i} ; Value{ Integral{ [ CompZ[{El}+E1l[] ] * expialphax~{i}[]/d ] ; In SurfCutSubs1  ; Integration Int_1 ; Jacobian JSur ; } } }

        { Name int_x_r~{i} ; Value{ Integral{ [ CompX[{Et}+E1dt[]] * expialphax~{i}[]/d ] ; In SurfCutSuper1 ; Integration Int_1 ; Jacobian JSur ; } } }
        { Name int_y_r~{i} ; Value{ Integral{ [({uy}+CompY[E1dt[]])* expialphax~{i}[]/d ] ; In SurfCutSuper1 ; Integration Int_1 ; Jacobian JSur ; } } }
        { Name int_z_r~{i} ; Value{ Integral{ [ CompZ[{El}+E1dl[]] * expialphax~{i}[]/d ] ; In SurfCutSuper1 ; Integration Int_1 ; Jacobian JSur ; } } }
      EndFor

      For i In {0:2*nb_orders}
        { Name eff_t~{i}   ; Value { Term{Type Global; [
                1/(Ae^2*-beta1[]) * ( beta_subs~{i}[] * SquNorm[$int_x_t~{i}]+
                                      beta_subs~{i}[] * SquNorm[$int_y_t~{i}]+
                                      beta_subs~{i}[] * SquNorm[$int_z_t~{i}] ) ] ; In SurfCutSubs1 ; } } }
        { Name eff_r~{i}   ; Value { Term{Type Global; [
                1/(Ae^2*-beta1[]) * ( beta_super~{i}[] * SquNorm[$int_x_r~{i}]+
                                      beta_super~{i}[] * SquNorm[$int_y_r~{i}]+
                                      beta_super~{i}[] * SquNorm[$int_z_r~{i}] ) ] ; In SurfCutSuper1 ; } } }                                  
        // // BUGGY
        // { Name eff_t~{i}   ; Value{ Term{Type Global; [
        //         1/(Ae^2*beta_subs~{i}[]*-beta1[]) * ((beta_subs~{i}[]^2+gamma[]^2    )*SquNorm[$int_z_t~{i}]+
        //                                              (beta_subs~{i}[]^2+alpha~{i}[]^2)*SquNorm[$int_x_t~{i}]+
        //                                             2*alpha~{i}[]*gamma[]*Re[$int_z_t~{i}*Conj[$int_x_t~{i}]] ) ] ; In SurfCutSubs1 ; } } }
        // { Name eff_r~{i}   ; Value { Term{Type Global; [
        //         1/(Ae^2*beta_super~{i}[]*-beta1[]) * ((beta_super~{i}[]^2+gamma[]^2    )*SquNorm[$int_z_r~{i}]+
        //                                               (beta_super~{i}[]^2+alpha~{i}[]^2)*SquNorm[$int_x_r~{i}]+
        //                                               2*alpha~{i}[]*gamma[]*Re[$int_z_r~{i}*Conj[$int_x_r~{i}]]) ] ; In SurfCutSuper1 ; } } }
      EndFor
    }
  }
}

PostOperation {
  { Name postop_energy; NameOfPostProcessing postpro_energy ;
    Operation {
      Print[ lambda_step, OnPoint{0,0,0}, Format Table, File > StrCat[myDir, "temp_lambda_step.txt"], SendToServer "GetDP/Lambda_step" ] ;
      Print [ Etot , OnElementsOf Omega, File "Etot.pos" ];
      For i In {0:2*nb_orders}
        Print[ int_x_t~{i}[SurfCutSubs1] , OnGlobal, StoreInVariable $int_x_t~{i}, File > StrCat[myDir, "int_x_t.txt"], Format Table];
        Print[ int_y_t~{i}[SurfCutSubs1] , OnGlobal, StoreInVariable $int_y_t~{i}, File > StrCat[myDir, "int_y_t.txt"], Format Table];
        Print[ int_z_t~{i}[SurfCutSubs1] , OnGlobal, StoreInVariable $int_z_t~{i}, File > StrCat[myDir, "int_z_t.txt"], Format Table];
        Print[ int_x_r~{i}[SurfCutSuper1], OnGlobal, StoreInVariable $int_x_r~{i}, File > StrCat[myDir, "int_x_r.txt"], Format Table];
        Print[ int_y_r~{i}[SurfCutSuper1], OnGlobal, StoreInVariable $int_y_r~{i}, File > StrCat[myDir, "int_y_r.txt"], Format Table];
        Print[ int_z_r~{i}[SurfCutSuper1], OnGlobal, StoreInVariable $int_z_r~{i}, File > StrCat[myDir, "int_z_r.txt"], Format Table];
      EndFor
      For i In {0:2*nb_orders}
        Print[ eff_t~{i}[SurfCutSubs1]  , OnRegion SurfCutSubs1 , File > StrCat[myDir, "eff_t.txt"], Format Table ];
        Print[ eff_r~{i}[SurfCutSuper1] , OnRegion SurfCutSuper1, File > StrCat[myDir, "eff_r.txt"], Format Table ];
      EndFor
      Print[ Q_tot[Plot_domain]     , OnGlobal, Format FrequencyTable, File > StrCat[myDir, "absorption-Q_tot.txt"]];
      Print[ Q_subs[sub]            , OnGlobal, Format FrequencyTable, File > StrCat[myDir, "absorption-Q_subs.txt"]];
      Print[ Q_rod_out[rod_out]     , OnGlobal, Format FrequencyTable, File > StrCat[myDir, "absorption-Q_rod_out.txt"]];
      Print[ Q_layer_dep[layer_dep] , OnGlobal, Format FrequencyTable, File > StrCat[myDir, "absorption-Q_layer_dep.txt"]];
      Print[ Q_layer_cov[layer_cov] , OnGlobal, Format FrequencyTable, File > StrCat[myDir, "absorption-Q_layer_cov.txt"]];
      If(multiplot)
        Echo[ Str["For i In {PostProcessing.NbViews-1:0:-1}",
                    "  If(!StrCmp(View[i].Name, 'boundary') || !StrCmp(View[i].Name, 'boundary_Combine'))",
                    "    Delete View[i];",
                    "  EndIf",
                    "EndFor"], File StrCat[myDir,"tmp1.geo"]] ;
        For i In {-nb_plot_periods:nb_plot_periods}
          Print[ Etot~{i}, OnElementsOf Plot_domain, File StrCat[myDir, Sprintf("u_tot_lambda%.2fnm_%g.pos", lambda0/nm,i+nb_plot_periods)], ChangeOfCoordinates {$X+i*d,$Y,$Z}, Name Sprintf("E_tot_%.2fnm.pos", lambda0/nm) ] ;
          Print[ boundary, OnElementsOf Plot_bnd, File StrCat[myDir,Sprintf("boundary_%g.pos",i+nb_plot_periods)], ChangeOfCoordinates {$X+i*d,$Y,$Z} ] ;
        EndFor
        Echo[ Str["Combine ElementsByViewName;",
                  "Hide {",
                  "Point{1,2,7,8,9,10,20,22};",
                  "Line{1,7,8,9,10,30,32,34,2,3,4,5,6,12,16,20,24,28};",
                  "Surface{36,48};}",
                  "Geometry.Color.Lines = {0,0,0};",
                  "l=PostProcessing.NbViews-1; View[l].ColorTable={Black}; ",
                  "View[l-1].Visible=1; View[l-1].ShowScale=0;",
                  "View[l].ShowScale=0; View[l].LineWidth=1.5; View[l].LineType=1;Geometry.LineWidth=0;"],
                  File StrCat[myDir,"tmp3.geo" ]] ;
      EndIf
    }
  }
}

DefineConstant[
  R_ = {"helmoltz_scalar", Name "GetDP/1ResolutionChoices", Visible 1},
  C_ = {"-solve -pos -petsc_prealloc 200 -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps", Name "GetDP/9ComputeCommand", Visible 1},
  P_ = {"postop_energy", Name "GetDP/2PostOperationChoices", Visible 1}
];

If(plotRTgraphs)
  DefineConstant[
    refl_  = {0, Name "GetDP/R0", ReadOnly 1, Graph "02000000", Visible 1},
    abs_   = {0, Name "GetDP/Omegaal absorption", ReadOnly 1, Graph "00000002", Visible 1},
    trans_ = {0, Name "GetDP/T0", ReadOnly 1, Graph "000000000002", Visible 1}
  ];
EndIf
