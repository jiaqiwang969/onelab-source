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

myDir = "res2D_scalar/";
DefineConstant[
  lambda0 = {lambda_min , Min lambda_min, Max lambda_max, Step (lambda_max-lambda_min)/(nb_lambdas-1), Name StrCat[pp2, "1wavelength [nm]"] , Loop 1, Highlight Str[colorpp2],Graph "200000200020", ServerAction "Reset GetDP/T0, GetDP/R0, GetDP/Lambda_step, GetDP/total absorption"}
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
}

Function{
  I[] = Complex[0.0,1.0];
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

  epsr2_re[]          = epsr_re_dom_1[];
  epsr2_im[]          = epsr_im_dom_1[];
  epsr_layer_dep_re[] = epsr_re_dom_2[];
  epsr_layer_dep_im[] = epsr_im_dom_2[];
  epsr_rod_out_re[]   = epsr_re_dom_3[];
  epsr_rod_out_im[]   = epsr_im_dom_3[];
  epsr_rods_re[]      = epsr_re_dom_4[];
  epsr_rods_im[]      = epsr_im_dom_4[];
  epsr_layer_cov_re[] = epsr_re_dom_5[];
  epsr_layer_cov_im[] = epsr_im_dom_5[];
  epsr1_re[]          = epsr_re_dom_6[];
  epsr1_im[]          = epsr_im_dom_6[];
  
  om0     =  2.*Pi*cel/lambda0;
  k0      =  2.*Pi/lambda0;
  epsr1[] =  Complex[epsr1_re[],epsr1_im[]];
  epsr2[] =  Complex[epsr2_re[],epsr2_im[]];
  n1[]    =  Sqrt[epsr1[]];
  n2[]    =  Sqrt[epsr2[]];
  k1norm[]=  k0*n1[];
  k2norm[]=  k0*n2[];
  alpha[] =  k0*n1[]*Sin[theta_deg*deg2rad];
  beta1[] = -k0*n1[]*Cos[theta_deg*deg2rad];
  beta2[] = -Sqrt[k0^2*epsr2[]-alpha[]^2];
  
  k1[] = Vector[alpha[], beta1[],0];
  k1r[]= Vector[alpha[],-beta1[],0];
  k2[] = Vector[alpha[], beta2[],0];

  If (flag_polar==1)
    r[]    = (CompY[k1[]]*epsr2[]-CompY[k2[]]*epsr1[])/(CompY[k1[]]*epsr2[]+CompY[k2[]]*epsr1[]);
    t[]    =                  (2.*CompY[k1[]]*epsr2[])/(CompY[k1[]]*epsr2[]+CompY[k2[]]*epsr1[]);
    Pinc[] = 0.5*A^2*Sqrt[mu0/(ep0*epsr1_re[])] * Cos[theta_deg*deg2rad];
  EndIf
  If (flag_polar==0)
    r[]    = (CompY[k1[]]-CompY[k2[]])/(CompY[k1[]]+CompY[k2[]]);
    t[]    =           2.*CompY[k1[]] /(CompY[k1[]]+CompY[k2[]]);
    Pinc[] =  0.5*A^2*Sqrt[ep0*epsr1_re[]/mu0] * Cos[theta_deg*deg2rad]; 
  EndIf

  BlochX_phase_shift[] = Exp[I[]*alpha[]*d];

  For i In {0:2*nb_orders}
    alpha~{i}[]      = CompX[k1[]]  + 2.*Pi/d*(i-nb_orders);
    expmialpha~{i}[] = Exp[-I[]*alpha~{i}[]*X[]];
    beta_super~{i}[] = -Sqrt[k0^2*epsr1[]-alpha~{i}[]^2];
    beta_subs~{i}[]  = -Sqrt[k0^2*epsr2[]-alpha~{i}[]^2];
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

  ui[pmltop]    = 0.;
  ui[pmlbot]    = 0.;
  ui[sup]       = A*Exp[I[]*k1[]*XYZ[]];
  ui[layer_cov] = A*Exp[I[]*k1[]*XYZ[]];
  ui[rod_out]   = A*Exp[I[]*k1[]*XYZ[]];
  ui[rods]      = A*Exp[I[]*k1[]*XYZ[]];
  ui[layer_dep] = A*Exp[I[]*k1[]*XYZ[]];
  ui[sub]       = 0.;
  
  ur[pmltop]    = 0.;
  ur[pmlbot]    = 0.;
  ur[sup]       = r[]*Exp[I[]*k1r[]*XYZ[]];
  ur[layer_dep] = r[]*Exp[I[]*k1r[]*XYZ[]];
  ur[rod_out]   = r[]*Exp[I[]*k1r[]*XYZ[]];
  ur[rods]      = r[]*Exp[I[]*k1r[]*XYZ[]];
  ur[layer_cov] = r[]*Exp[I[]*k1r[]*XYZ[]];
  ur[sub]       = 0;

  ut[pmltop]    = 0.;
  ut[pmlbot]    = 0.;
  ut[sup]       = 0.;
  ut[layer_dep] = 0.;
  ut[rod_out]   = 0.;
  ut[rods]      = 0.;
  ut[layer_cov] = 0.;
  ut[sub]       = t[]*Exp[I[]*k2[]*XYZ[]];
  
  u1[]          = ui[]+ur[]+ut[];
  u1d[]         = ur[]+ut[];
  
  If (flag_polar==1)
    E1i[]  = -1/(om0*ep0*epsilonr_annex[]) * k1[]  /\ Vector[0,0,ui[]];
    E1r[]  = -1/(om0*ep0*epsilonr_annex[]) * k1r[] /\ Vector[0,0,ur[]];
    E1t[]  = -1/(om0*ep0*epsilonr_annex[]) * k2[]  /\ Vector[0,0,ut[]];
    E1[]   = E1i[]+E1r[]+E1t[];
    Ex1[]  = CompX[E1[]];
    Ey1[]  = CompY[E1[]];

    detepst[]       = CompXX[epsilonr[]]*CompYY[epsilonr[]]-CompXY[epsilonr[]]*Conj[CompYX[epsilonr[]]];
    detepst_annex[] = CompXX[epsilonr_annex[]]*CompYY[epsilonr_annex[]]-CompXY[epsilonr_annex[]]*Conj[CompYX[epsilonr_annex[]]];
    xsi[]           = Transpose[epsilonr[]]/detepst[];
    xsi_annex[]     = Transpose[epsilonr_annex[]]/detepst_annex[];
    chi[]           = CompZZ[mur[]];
    source_xsi_r[]      = (xsi[]-xsi_annex[]) * -k1[]  * I[] * ui[];
    source_xsi_i[]      = (xsi[]-xsi_annex[]) * -k1r[] * I[] * ur[];
    source_chi_r[]      = 0;
    source_chi_i[]      = 0;

  Else
    detmut[]    = CompXX[mur[]]*CompYY[mur[]]-CompXY[mur[]]*Conj[CompYX[mur[]]];
    xsi[]       = Transpose[mur[]]/detmut[];
    chi[]       = CompZZ[epsilonr[]];
    chi_annex[] = CompZZ[epsilonr_annex[]];
    source_xsi_r[] = Vector[0.,0.,0.];
    source_xsi_i[] = Vector[0.,0.,0.];
    source_chi_r[] = k0^2 * (chi[]-chi_annex[]) * ur[];
    source_chi_i[] = k0^2 * (chi[]-chi_annex[]) * ui[];
  EndIf
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
  { Name Hgrad; Type Form0;
    BasisFunction {
      { Name sn;  NameOfCoef un;  Function BF_Node;    Support Region[Omega]; Entity NodesOf[Omega]; }
      { Name sn2; NameOfCoef un2; Function BF_Node_2E; Support Region[Omega]; Entity EdgesOf[Omega]; }
      If (Flag_o2_interp==1)
        { Name un3; NameOfCoef un3; Function BF_Node_2F; Support Region[Omega]; Entity FacetsOf[Omega]; }
        { Name un4; NameOfCoef un4; Function BF_Node_3E; Support Region[Omega]; Entity EdgesOf[Omega]; }
        { Name un5; NameOfCoef un5; Function BF_Node_3F; Support Region[Omega]; Entity FacetsOf[Omega]; }
      EndIf
      
    }
    Constraint {
      { NameOfCoef un;  EntityType NodesOf ; NameOfConstraint Bloch; }
      { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Bloch; }
      If (Flag_o2_interp==1)
        { NameOfCoef un3; EntityType EdgesOf ; NameOfConstraint Bloch; }
      EndIf
    }
  }
}

Formulation {
  {Name helmoltz_scalar; Type FemEquation;
    Quantity {
      { Name u2d; Type Local; NameOfSpace Hgrad;}}
    Equation {
      Galerkin { [-xsi[]*Dof{Grad u2d} , {Grad u2d}]; In Omega; Jacobian JVol; Integration Int_1; }
      Galerkin { [k0^2*chi[]*Dof{u2d}  ,      {u2d}]; In Omega; Jacobian JVol; Integration Int_1; }
      Galerkin { [source_xsi_i[]       , {Grad u2d}]; In Omega; Jacobian JVol; Integration Int_1; }
      Galerkin { [source_xsi_r[]       , {Grad u2d}]; In Omega; Jacobian JVol; Integration Int_1; }        
      Galerkin { [source_chi_r[]       ,      {u2d}]; In Omega; Jacobian JVol; Integration Int_1; }
      Galerkin { [source_chi_i[]       ,      {u2d}]; In Omega; Jacobian JVol; Integration Int_1; }
    }
  }
}

Resolution {
  { Name helmoltz_scalar;
    System {
      { Name S; NameOfFormulation helmoltz_scalar; Type ComplexValue; }
    }
    Operation {
      CreateDir[Str[myDir]];
      Printf["%f",lambda0/nm];
      Generate[S] ;
      Solve[S] ;
    }
  }
}

PostProcessing {
  { Name postpro_energy; NameOfFormulation helmoltz_scalar;
    Quantity {
      { Name u           ; Value { Local { [ {u2d}       ]; In Omega; Jacobian JVol; } } }
      { Name u_diff      ; Value { Local { [ {u2d}+u1d[] ]; In Omega; Jacobian JVol; } } }
      { Name u_tot       ; Value { Local { [ {u2d}+u1[]  ]; In Omega; Jacobian JVol; } } }
      { Name u1          ; Value { Local { [ u1[]        ]; In Omega; Jacobian JVol; } } }
      { Name lambda_step ; Value { Local { [ lambda0/nm  ]; In Omega; Jacobian JVol; } } }
      { Name boundary    ; Value { Local { [ bndCol[]    ]; In Plot_bnd ; Jacobian JVol ; } } }
      For i In {-nb_plot_periods:nb_plot_periods}
        { Name u_tot~{i}     ; Value { Local { [ ({u2d}+u1[])*Exp[I[]*i*CompX[k1[]]*d] ]; In Omega; Jacobian JVol; } } }
      EndFor
      For i In {0:2*nb_orders}
        { Name s_r~{i} ; Value { Integral{ [ expmialpha~{i}[] * ({u2d}+u1d[])/d ] ; In SurfCutSuper1 ; Jacobian JSur ; Integration Int_1 ; } } }
        { Name s_t~{i} ; Value { Integral{ [ expmialpha~{i}[] * ({u2d}+u1d[])/d ] ; In SurfCutSubs1  ; Jacobian JSur ; Integration Int_1 ; } } }
        { Name order_t_angle~{i} ; Value { Local{ [-Atan2[Re[alpha~{i}[]],Re[beta_subs~{i}[]]]/deg2rad ] ; In Omega; Jacobian JVol; } } }
        { Name order_r_angle~{i} ; Value { Local{ [ Atan2[Re[alpha~{i}[]],Re[beta_super~{i}[]]]/deg2rad ] ; In Omega; Jacobian JVol; } } }
      EndFor
      If (flag_polar==1)
        For i In {0:2*nb_orders}
          { Name eff_r~{i} ; Value { Term{ Type Global; [ SquNorm[#i]*beta_super~{i}[]/beta1[] ] ; In SurfCutSuper1 ; } } }
          { Name eff_t~{i} ; Value { Term{ Type Global; [ SquNorm[#(2*nb_orders+1+i)]*(beta_subs~{i}[]/beta1[])*(epsr1[]/epsr2[]) ] ; In SurfCutSubs1 ; } } }
        EndFor
        For i In {0:N_rods-1:1}
          { Name Q_rod~{i}   ; Value { Integral { [  0.5 * ep0*om0*Fabs[epsr_rods_im[]]          * ( SquNorm[CompY[{Grad u2d}]*I[]/(om0*ep0*CompXX[epsilonr[]])+Ex1[]/CompXX[epsilonr[]]*CompXX[epsilonr_annex[]]] + SquNorm[-CompX[{Grad u2d}]*I[]/(om0*ep0*CompXX[epsilonr[]])+Ey1[]/CompXX[epsilonr[]]*CompXX[epsilonr_annex[]]] ) / (Pinc[]*d) ] ; In rod~{i}   ; Integration Int_1 ; Jacobian JVol ; } } }
        EndFor
        { Name Q_subs        ; Value { Integral { [  0.5 * ep0*om0*Fabs[epsr2_im[]      ]        * ( SquNorm[CompY[{Grad u2d}]*I[]/(om0*ep0*CompXX[epsilonr[]])+Ex1[]/CompXX[epsilonr[]]*CompXX[epsilonr_annex[]]] + SquNorm[-CompX[{Grad u2d}]*I[]/(om0*ep0*CompYY[epsilonr[]])+Ey1[]/CompYY[epsilonr[]]*CompYY[epsilonr_annex[]] ] ) / (Pinc[]*d) ] ; In sub         ; Integration Int_1 ; Jacobian JVol ; } } }
        { Name Q_rod_out     ; Value { Integral { [  0.5 * ep0*om0*Fabs[epsr_rod_out_im[]]       * ( SquNorm[CompY[{Grad u2d}]*I[]/(om0*ep0*CompXX[epsilonr[]])+Ex1[]/CompXX[epsilonr[]]*CompXX[epsilonr_annex[]]] + SquNorm[-CompX[{Grad u2d}]*I[]/(om0*ep0*CompYY[epsilonr[]])+Ey1[]/CompYY[epsilonr[]]*CompYY[epsilonr_annex[]] ] ) / (Pinc[]*d) ] ; In rod_out   ; Integration Int_1 ; Jacobian JVol ; } } }
        { Name Q_layer_dep   ; Value { Integral { [  0.5 * ep0*om0*Fabs[epsr_layer_dep_im[]]     * ( SquNorm[CompY[{Grad u2d}]*I[]/(om0*ep0*CompXX[epsilonr[]])+Ex1[]/CompXX[epsilonr[]]*CompXX[epsilonr_annex[]]] + SquNorm[-CompX[{Grad u2d}]*I[]/(om0*ep0*CompYY[epsilonr[]])+Ey1[]/CompYY[epsilonr[]]*CompYY[epsilonr_annex[]] ] ) / (Pinc[]*d) ] ; In layer_dep   ; Integration Int_1 ; Jacobian JVol ; } } }
        { Name Q_layer_cov   ; Value { Integral { [  0.5 * ep0*om0*Fabs[epsr_layer_cov_im[]]     * ( SquNorm[CompY[{Grad u2d}]*I[]/(om0*ep0*CompXX[epsilonr[]])+Ex1[]/CompXX[epsilonr[]]*CompXX[epsilonr_annex[]]] + SquNorm[-CompX[{Grad u2d}]*I[]/(om0*ep0*CompYY[epsilonr[]])+Ey1[]/CompYY[epsilonr[]]*CompYY[epsilonr_annex[]] ] ) / (Pinc[]*d) ] ; In layer_cov   ; Integration Int_1 ; Jacobian JVol ; } } }
        { Name Q_tot         ; Value { Integral { [  0.5 * ep0*om0*Fabs[Im[CompZZ[epsilonr[]]]]  * ( SquNorm[CompY[{Grad u2d}]*I[]/(om0*ep0*CompXX[epsilonr[]])+Ex1[]/CompXX[epsilonr[]]*CompXX[epsilonr_annex[]]] + SquNorm[-CompX[{Grad u2d}]*I[]/(om0*ep0*CompYY[epsilonr[]])+Ey1[]/CompYY[epsilonr[]]*CompYY[epsilonr_annex[]] ] ) / (Pinc[]*d) ] ; In Plot_domain ; Integration Int_1 ; Jacobian JVol ; } } }
      EndIf
      If (flag_polar==0)        
        For i In {0:2*nb_orders}
          { Name eff_r~{i} ; Value { Term{ Type Global; [ SquNorm[#i]*beta_super~{i}[]/beta1[] ] ; In SurfCutSuper1 ; } } }
          { Name eff_t~{i} ; Value { Term{ Type Global; [ SquNorm[#(2*nb_orders+1+i)]*(beta_subs~{i}[]/beta1[])] ; In SurfCutSubs1 ; } } }
        EndFor
        For i In {0:N_rods-1:1}
          { Name Q_rod~{i}  ; Value { Integral { [0.5 * ep0*om0*Fabs[epsr_rods_im[]]          * (SquNorm[{u2d}+u1[]]) / (Pinc[]*d) ] ; In rod~{i}    ; Integration Int_1 ; Jacobian JVol ; } } }
        EndFor
        { Name Q_subs       ; Value { Integral { [  0.5 * ep0*om0*Fabs[epsr2_im[]]            *(SquNorm[{u2d}+u1[]]) / (Pinc[]*d) ] ; In sub         ; Integration Int_1 ; Jacobian JVol ; } } }
        { Name Q_rod_out    ; Value { Integral { [  0.5 * ep0*om0*Fabs[epsr_rod_out_im[]]     *(SquNorm[{u2d}+u1[]]) / (Pinc[]*d) ] ; In rod_out     ; Integration Int_1 ; Jacobian JVol ; } } }
        { Name Q_layer_dep  ; Value { Integral { [  0.5 * ep0*om0*Fabs[epsr_layer_dep_im[]]   *(SquNorm[{u2d}+u1[]]) / (Pinc[]*d) ] ; In layer_dep   ; Integration Int_1 ; Jacobian JVol ; } } }
        { Name Q_layer_cov  ; Value { Integral { [  0.5 * ep0*om0*Fabs[epsr_layer_cov_im[]]   *(SquNorm[{u2d}+u1[]]) / (Pinc[]*d) ] ; In layer_cov   ; Integration Int_1 ; Jacobian JVol ; } } }
        { Name Q_tot        ; Value { Integral { [  0.5 * ep0*om0*Fabs[Im[CompXX[epsilonr[]]]]*(SquNorm[{u2d}+u1[]]) / (Pinc[]*d) ] ; In Plot_domain ; Integration Int_1 ; Jacobian JVol ; } } }
      EndIf
    }
  }
}

PostOperation {
  { Name postop_energy; NameOfPostProcessing postpro_energy ;
    Operation {
      Print[ lambda_step, OnPoint{0,0,0}, Format Table, File > StrCat[myDir, "temp_lambda_step.txt"], SendToServer "GetDP/Lambda_step" ] ;
      For i In {0:2*nb_orders}
        Print[ s_r~{i}[SurfCutSuper1], OnGlobal, Store i                , Format Table , File > StrCat[myDir, Sprintf("temp_s_r_%g.txt", i-nb_orders)]];
        Print[ s_t~{i}[SurfCutSubs1] , OnGlobal, Store (2*nb_orders+1+i), Format Table , File > StrCat[myDir, Sprintf("temp_s_t_%g.txt", i-nb_orders)]];
      EndFor
      For i In {0:2*nb_orders}
        Print[ eff_r~{i}[SurfCutSuper1], OnRegion SurfCutSuper1,Store (4*nb_orders+1+i), Format FrequencyTable, File > StrCat[myDir, Sprintf("efficiency_r_%g.txt", i-nb_orders)]];
        Print[ eff_t~{i}[SurfCutSubs1] , OnRegion SurfCutSubs1 ,Store (6*nb_orders+1+i), Format FrequencyTable, File > StrCat[myDir, Sprintf("efficiency_t_%g.txt", i-nb_orders)]];
        Print[ order_r_angle~{i}     , OnPoint{0,0,0}, Format Table , File > StrCat[myDir, Sprintf("order_r_angle_%g.txt", i-nb_orders)]];
        Print[ order_t_angle~{i}     , OnPoint{0,0,0}, Format Table , File > StrCat[myDir, Sprintf("order_t_angle_%g.txt", i-nb_orders)]];          
      EndFor
      Print[ eff_r~{nb_orders}[SurfCutSuper1], OnRegion SurfCutSuper1, Format Table, SendToServer "GetDP/R0", File StrCat[myDir, "temp_R0.txt"]];
      Print[ eff_t~{nb_orders}[SurfCutSubs1] , OnRegion SurfCutSubs1 , Format Table, SendToServer "GetDP/T0", File StrCat[myDir, "temp_T0.txt"]];
      Print[ Q_tot[Plot_domain]             , OnGlobal  , Format FrequencyTable ,SendToServer "GetDP/total absorption", File > StrCat[myDir, "absorption-Q_tot.txt"]];
      For i In {0:N_rods-1:1}
        Print[ Q_rod~{i}[rod~{i}] , OnGlobal, Format FrequencyTable, File > StrCat[myDir, Sprintf("absorption-Q_rod_%g.txt", i+1) ]];           
      EndFor  
      Print[ Q_tot[Plot_domain]     , OnGlobal, Format FrequencyTable, File > StrCat[myDir, "absorption-Q_tot.txt"]];
      Print[ Q_subs[sub]            , OnGlobal, Format FrequencyTable, File > StrCat[myDir, "absorption-Q_subs.txt"]];
      Print[ Q_rod_out[rod_out]     , OnGlobal, Format FrequencyTable, File > StrCat[myDir, "absorption-Q_rod_out.txt"]];
      Print[ Q_layer_dep[layer_dep] , OnGlobal, Format FrequencyTable, File > StrCat[myDir, "absorption-Q_layer_dep.txt"]];
      Print[ Q_layer_cov[layer_cov] , OnGlobal, Format FrequencyTable, File > StrCat[myDir, "absorption-Q_layer_cov.txt"]];
      If (flag_polar==1)
        Print[ u_tot    , OnElementsOf Plot_domain, File StrCat[myDir, Sprintf("Hz_tot_lambda%.2fnm_1.pos", lambda0/nm)] , Name Sprintf("Hz_tot_%.2fnm.pos", lambda0/nm)];
      EndIf
      If (flag_polar==0)
        Print[ u_tot  , OnElementsOf Plot_domain, File StrCat[myDir, Sprintf("Ez_tot_lambda%.2fnm_1.pos", lambda0/nm)] , Name Sprintf("Ez_tot_%.2fnm.pos", lambda0/nm)];
        // Print[ u, OnElementsOf Plot_domain, File StrCat[myDir, Sprintf("u2d_%.2fnm_1.pos", lambda0/nm)] , Name Sprintf("u2d_%.2fnm.pos", lambda0/nm)];
      EndIf
      If(multiplot)
        Echo[ Str["For i In {PostProcessing.NbViews-1:0:-1}",
                    "  If(!StrCmp(View[i].Name, 'boundary') || !StrCmp(View[i].Name, 'boundary_Combine'))",
                    "    Delete View[i];",
                    "  EndIf",
                    "EndFor"], File StrCat[myDir,"tmp1.geo"]] ;
        For i In {-nb_plot_periods:nb_plot_periods}
          Print[ u_tot~{i}, OnElementsOf Plot_domain, File StrCat[myDir, Sprintf("u_tot_lambda%.2fnm_%g.pos", lambda0/nm,i+nb_plot_periods)], ChangeOfCoordinates {$X+i*d,$Y,$Z}, Name Sprintf("u_tot_%.2fnm.pos", lambda0/nm) ] ;
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
  C_ = {"-solve -pos -petsc_prealloc 100 -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps", Name "GetDP/9ComputeCommand", Visible 1},
  P_ = {"postop_energy", Name "GetDP/2PostOperationChoices", Visible 1}
];

If(plotRTgraphs)
  DefineConstant[
    refl_  = {0, Name "GetDP/R0", ReadOnly 1, Graph "02000000", Visible 1},
    abs_   = {0, Name "GetDP/total absorption", ReadOnly 1, Graph "00000002", Visible 1},
    trans_ = {0, Name "GetDP/T0", ReadOnly 1, Graph "000000000002", Visible 1}
  ];
EndIf
