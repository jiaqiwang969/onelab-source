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
Include "grating3D_materials.pro";

myDir = "res3D/";

Group {
	// SubDomains
	PMLbot   = Region[1];
	L_6_temp = Region[2];
	L_5      = Region[3];
	L_4      = Region[4];
	L_3      = Region[5];
	L_2      = Region[6];
	L_1_temp = Region[7];
	PMLtop   = Region[8];
	Scat     = Region[9];

	// Boundaries
	SurfBlochXm    = Region[101];
	SurfBlochXp    = Region[102];
	SurfBlochYm    = Region[201];
	SurfBlochYp    = Region[202];
  SurfIntTop     = Region[301];
  SurfIntBot     = Region[302];


	SurfDirichlet  = Region[{401,402}];
  SurfBloch      = Region[{SurfBlochXm,SurfBlochXp,SurfBlochYm,SurfBlochYp}];

  If (FlagLinkFacets==1)
    SurfExcludeFacets  = Region[{}];
  Else
    SurfExcludeFacets  = Region[{SurfBloch}];
  EndIf

  L_1 = Region[{L_1_temp,SurfIntTop}];
  L_6 = Region[{L_6_temp,SurfIntBot}];

	Omega          = Region[{PMLbot,L_6,L_5,L_4,L_3,L_2,L_1,PMLtop,Scat}];
	Omega_nosource = Region[{PMLbot,L_6,L_1,PMLtop}];
	Omega_source   = Region[{Scat,L_2,L_3,L_4,L_5}];
	Omega_super    = Region[{Scat,L_1,L_2,L_3,L_4,L_5,PMLtop}];
	Omega_subs     = Region[{L_6,PMLbot}];
	Omega_plot     = Region[{L_6,L_5,L_4,L_3,L_2,L_1,Scat}];

	// SurfNeumann    = Region[{SurfBlochXm,SurfBlochXp,SurfBlochYm,SurfBlochYp}];
}



Function{
  I[] = Complex[0,1];
  zhat[] = Vector[0,0,1];

  small_delta = 0.0*nm;
	mu0         = 4*Pi*1.e2*nm;
	ep0         = 8.854187817e-3*nm;
	cel         = 1/Sqrt[ep0 * mu0];
	om0         = 2*Pi*cel/lambda0;
	k0          = 2.*Pi/lambda0;
	Ae          = 1;
	Pinc        = 0.5*Ae^2*Sqrt[eps_re_L_1*ep0/mu0] * Cos[theta0];

  // permittivities
  For i In {1:6}
    If (flag_mat~{i}==0)
      epsr[L~{i}]  = Complex[eps_re_L~{i} , eps_im_L~{i}] * TensorDiag[1,1,1];
      If (i==1)
        epsr1[] = Complex[eps_re_L~{i} , eps_im_L~{i}];
      EndIf
      If (i==6)
        epsr2[] = Complex[eps_re_L~{i} , eps_im_L~{i}];
      EndIf
    Else
      For j In {1:nb_available_materials}
        If(flag_mat~{i}==j)
          epsr[L~{i}]  = Complex[epsr_re_interp_mat~{j}[lambda0/nm*1e-9] , epsr_im_interp_mat~{j}[lambda0/nm*1e-9]] * TensorDiag[1,1,1];
          If (i==1)
            epsr1[] = Complex[epsr_re_interp_mat~{j}[lambda0/nm*1e-9] , epsr_im_interp_mat~{j}[lambda0/nm*1e-9]];
          EndIf
          If (i==6)
            epsr2[] = Complex[epsr_re_interp_mat~{j}[lambda0/nm*1e-9] , epsr_im_interp_mat~{j}[lambda0/nm*1e-9]];
          EndIf
        EndIf
      EndFor
    EndIf
  EndFor

  If (flag_mat_scat==0)
    epsr[Scat] = Complex[eps_re_Scat , eps_im_Scat] * TensorDiag[1,1,1];
  Else
    For j In {flag_mat_scat:flag_mat_scat}
      epsr[Scat] = Complex[epsr_re_interp_mat~{j}[lambda0/nm*1e-9] , epsr_im_interp_mat~{j}[lambda0/nm*1e-9]] * TensorDiag[1,1,1];
    EndFor
  EndIf

	For i In {1:6}
		mur[L~{i}]   = TensorDiag[1,1,1];
	EndFor
	mur[Scat]  = TensorDiag[1,1,1];

  ////// PMLS
	a_pml           = 1.;
	b_pml           = 1.;
  // bermu
  n1[]    = Sqrt[epsr1[]];
  n2[]    = Sqrt[epsr2[]];
  k1norm[]= k0*n1[];
  k2norm[]= k0*n2[];
  Zmax = PML_top_hh;
  Zmin = hh_L_6;
  Damp_pml_top[] = 1/Abs[Zmax + PML_top - Z[]];
  Damp_pml_bot[] = 1/Abs[Zmin - PML_bot - Z[]];
  // DampingProfileX[] = 1/(Xmax + SizePMLX - Fabs[X[]]) - 1/(SizePMLX);

  sx            = 1.;
	sz[]          = Complex[a_pml,b_pml];
	sz_bermutop[] = Complex[1,Damp_pml_top[]/k1norm[]];
	sz_bermubot[] = Complex[1,Damp_pml_bot[]/k2norm[]];
	sy            = 1.;

	epsr[PMLtop]  = Re[epsr1[]]*TensorDiag[sz[]*sy/sx,sx*sz[]/sy,sx*sy/sz[]];
	mur[PMLtop]   =             TensorDiag[sz[]*sy/sx,sx*sz[]/sy,sx*sy/sz[]];
  epsr[PMLbot]  = Re[epsr2[]]*TensorDiag[sz[]*sy/sx,sx*sz[]/sy,sx*sy/sz[]];
  mur[PMLbot]   =             TensorDiag[sz[]*sy/sx,sx*sz[]/sy,sx*sy/sz[]];

  // epsr[PMLbot] = epsr2[];
  // mur[PMLbot]  = TensorDiag[1,1,1];

	epsr_annex[PMLbot]       = epsr[];
	epsr_annex[PMLtop]       = epsr[];
	epsr_annex[Omega_source] = epsr1[] * TensorDiag[1,1,1];
	epsr_annex[L_1]          = epsr[];
	epsr_annex[L_6]          = epsr[];

  //// Reference Field solution of annex problem (simple diopter)
  k1x[] = -k0*n1[]*Sin[theta0]*Cos[phi0];
  k1y[] = -k0*n1[]*Sin[theta0]*Sin[phi0];
  k1z[] = -k0*n1[]*Cos[theta0];
  k2x[] =  k1x[];
  k2y[] =  k1y[];
  k2z[] = -Sqrt[k0^2*epsr2[]-k1x[]^2-k1y[]^2];
  k1[]  = Vector[k1x[],k1y[], k1z[]];
  k2[]  = Vector[k2x[],k2y[], k2z[]];
  k1r[] = Vector[k1x[],k1y[],-k1z[]];

  rs[] = (k1z[]-k2z[])/(k1z[]+k2z[]);
  ts[] =    2.*k1z[] /(k1z[]+k2z[]);
  rp[] = (k1z[]*epsr2[]-k2z[]*epsr1[])/(k1z[]*epsr2[]+k2z[]*epsr1[]);
  tp[] =            (2.*k1z[]*epsr2[])/(k1z[]*epsr2[]+k2z[]*epsr1[]);

  spol[]    = Vector[Sin[phi0],-Cos[phi0],0];
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

  Ei[] = Ae*(Cos[psi0]*Eip[]-Sin[psi0]*Eis[]);
  Er[] = Ae*(Cos[psi0]*Erp[]-Sin[psi0]*Ers[]);
  Et[] = Ae*(Cos[psi0]*Etp[]-Sin[psi0]*Ets[]);
  Hi[] =  1/(om0*mu0*mur[]) * k1[]  /\ Ei[];
  Hr[] =  1/(om0*mu0*mur[]) * k1r[] /\ Er[];
  Ht[] =  1/(om0*mu0*mur[]) * k2[]  /\ Et[];

  E1[Omega_super]  = Ei[]+Er[];
  E1[Omega_subs]   = Et[];
  E1d[Omega_super] = Er[];
  E1d[Omega_subs]  = Et[];
  H1[Omega_super]  = Hi[]+Hr[];
  H1[Omega_subs]   = Ht[];
  H1d[Omega_super] = Hr[];
  H1d[Omega_subs]  = Ht[];

  source[] = (om0/cel)^2*(epsr[]-epsr_annex[])*E1[];

  // Bloch phase shifts
  skx1[] =  k1x[];
  // sky1[] = -k0*n1[]*Sin[theta0]*Sin[phi0+xsi];
  sky1[] = -k0*n1[]*Sin[theta0]*Sin[phi0+xsi];

  dephX[] = Exp[I[]*skx1[]*period_x];
	dephY[] = Exp[I[]*sky1[]*period_y];

	// Fourier coefficients variables
  Nb_ordre = 2*Nmax+1;
  For i In {0:Nb_ordre-1}
    For j In {0:Nb_ordre-1}
      alpha~{i}~{j}[] = -k1x[] + 2*Pi/period_x*(i-Nmax);
      beta~{i}~{j}[]  = -k1y[] + 2*Pi/period_y*(j-Nmax)/Cos[xsi] - 2*Pi/period_x*(i-Nmax)*Tan[xsi];
      expialphaxy~{i}~{j}[] = Exp[I[]*(alpha~{i}~{j}[]*X[]+beta~{i}~{j}[]*Y[])];
    EndFor
  EndFor
  For i In {0:Nb_ordre-1}
    For j In {0:Nb_ordre-1}
      gammar~{i}~{j}[] = Sqrt[k0^2*epsr1[] - alpha~{i}~{j}[]^2 - beta~{i}~{j}[]^2];
      gammat~{i}~{j}[] = Sqrt[k0^2*epsr2[] - alpha~{i}~{j}[]^2 - beta~{i}~{j}[]^2];
    EndFor
  EndFor

}

Constraint {
  { Name Dirichlet; Type Assign;
    Case {
      { Region SurfDirichlet; Value 0.; }
    }
  }
  { Name BlochX;
    Case {
    { Region SurfBlochXp; Type LinkCplx ; RegionRef SurfBlochXm;
      Coefficient dephX[]; Function Vector[$X-period_x,$Y,$Z] ; }
    }
  }
  { Name BlochY;
    Case {
      { Region SurfBlochYp; Type LinkCplx ; RegionRef SurfBlochYm;
      Coefficient dephY[]; Function Vector[$X-dys,$Y-dyc,$Z] ; }
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
  { Name JLin ;
    Case {
      { Region All ; Jacobian Lin ; }
    }
  }
}

Integration {
  { Name I1 ;
    Case {
      { Type Gauss ;
        Case {
          { GeoElement Point       ; NumberOfPoints   1 ; }
          { GeoElement Line        ; NumberOfPoints   4 ; }
          { GeoElement Triangle    ; NumberOfPoints  12 ; }
          { GeoElement Triangle2   ; NumberOfPoints  12 ; }
          { GeoElement Tetrahedron ; NumberOfPoints  16 ; }
          { GeoElement Tetrahedron2; NumberOfPoints  16 ; }
        }
      }
    }
  }
}

FunctionSpace {
  { Name Hcurl; Type Form1;
    BasisFunction {
      { Name sn; NameOfCoef un; Function BF_Edge; Support Region[{Omega}]; Entity EdgesOf[All]; }
      { Name sn2; NameOfCoef un2; Function BF_Edge_2E;Support Region[{Omega}]; Entity EdgesOf[All]; }
      If(oi==2)
        { Name sn3; NameOfCoef un3; Function BF_Edge_3F_b; Support Region[Omega]; Entity FacetsOf[Omega, Not SurfExcludeFacets]; }
        { Name sn4; NameOfCoef un4; Function BF_Edge_3F_c; Support Region[Omega]; Entity FacetsOf[Omega, Not SurfExcludeFacets]; }
        { Name sn5; NameOfCoef un5; Function BF_Edge_4E  ; Support Region[Omega]; Entity  EdgesOf[Omega, Not SurfExcludeFacets]; }
      EndIf
    }
    Constraint {
      { NameOfCoef un;  EntityType EdgesOf ; NameOfConstraint BlochX; }
      { NameOfCoef un;  EntityType EdgesOf ; NameOfConstraint BlochY; }
      { NameOfCoef un;  EntityType EdgesOf ; NameOfConstraint Dirichlet; }
      { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint BlochX; }
      { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint BlochY; }
      { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Dirichlet; }
      If (FlagLinkFacets==1)
        { NameOfCoef un3; EntityType FacetsOf ; NameOfConstraint BlochX; }
        { NameOfCoef un3; EntityType FacetsOf ; NameOfConstraint BlochY; }
        { NameOfCoef un4; EntityType FacetsOf ; NameOfConstraint BlochX; }
        { NameOfCoef un4; EntityType FacetsOf ; NameOfConstraint BlochY; }
        { NameOfCoef un5; EntityType EdgesOf  ; NameOfConstraint BlochX; }
        { NameOfCoef un5; EntityType EdgesOf  ; NameOfConstraint BlochY; }
      EndIf
      { NameOfCoef un3; EntityType FacetsOf ; NameOfConstraint Dirichlet; }
      { NameOfCoef un4; EntityType FacetsOf ; NameOfConstraint Dirichlet; }
      { NameOfCoef un5; EntityType EdgesOf  ; NameOfConstraint Dirichlet; }
    }
  }
}

Formulation {
  { Name helmholtz_vector; Type FemEquation;
    Quantity {
      { Name u; Type Local; NameOfSpace Hcurl;}
    }
		Equation{
      Galerkin { [-1/mur[]    * Dof{Curl u} , {Curl u}]; In Omega       ; Jacobian JVol; Integration I1; }
      Galerkin { [k0^2*epsr[] * Dof{u}      ,      {u}]; In Omega       ; Jacobian JVol; Integration I1; }
      Galerkin { [source[]                  ,      {u}]; In Omega_source; Jacobian JVol; Integration I1; }
    }
  }
}

Resolution {
  { Name helmholtz_vector;
    System {
      { Name M; NameOfFormulation helmholtz_vector; Type ComplexValue; }
    }
    Operation {
      CreateDir[Str[myDir]];
      Generate[M];
      Solve[M]; //SaveSolution[M];
    }
  }
}

PostProcessing {
  { Name postpro_helmholtz_vector; NameOfFormulation helmholtz_vector; NameOfSystem M;
    Quantity {
      { Name u      ; Value { Local { [ {u}       ]; In Omega; Jacobian JVol; } } }
      { Name Etot   ; Value { Local { [ {u}+E1[]  ]; In Omega; Jacobian JVol; } } }
      { Name Htot   ; Value { Local { [ H1[]-I[]/(mur[]*mu0*om0)*{Curl u}]; In Omega; Jacobian JVol; } } }
      { Name Edif   ; Value { Local { [ {u}+E1d[] ]; In Omega; Jacobian JVol; } } }
      { Name E1     ; Value { Local { [     E1[]  ]; In Omega; Jacobian JVol; } } }
      { Name Em     ; Value { Local { [     E1d[] ]; In Omega; Jacobian JVol; } } }
      { Name source ; Value { Local { [  source[] ]; In Omega; Jacobian JVol; } } }
      { Name epsr_xx; Value { Local { [  CompXX[epsr[]] ]; In Omega; Jacobian JVol; } } }
      { Name Damp_pml_top; Value { Local { [Damp_pml_top[]  ]; In Omega; Jacobian JVol; } } }
      { Name Poy_tot; Value { Local { [ 0.5*Re[Cross[{u}+E1[] , Conj[ H1[]-I[]/(mur[]*mu0*om0)*{Curl u}]]] ]; In Omega; Jacobian JVol; } } }
      { Name Poy_ref; Value { Local { [ 0.5*Re[Cross[{u}+E1d[], Conj[H1d[]-I[]/(mur[]*mu0*om0)*{Curl u}]]] ]; In Omega; Jacobian JVol; } } }
      { Name Poy_inc; Value { Local { [ 0.5*Re[Cross[    Ei[] , Conj[Hi[]]]                              ] ]; In Omega; Jacobian JVol; } } }
      { Name lambda_step   ; Value { Local { [ lambda0/nm ]; In Omega ; Jacobian JVol; } } }

      For k In {2:6}
        { Name Abs_L~{k} ; Value { Integral { [ ep0*om0 * 0.5*Im[CompXX[epsr[]]]*(SquNorm[{u}+E1[]]) / (Pinc*period_x*dyc) ] ; In L~{k} ; Integration I1 ; Jacobian JVol ; } } }
      EndFor
      { Name Abs_scat ; Value { Integral { [ ep0*om0 * 0.5*Im[CompXX[epsr[]]]*(SquNorm[{u}+E1[]]) / (Pinc*period_x*dyc) ] ; In Scat ; Integration I1 ; Jacobian JVol ; } } }
      { Name Abs_scat2; Value { Integral { [ ep0*om0 * 0.5*Im[CompXX[epsr[]]]*(SquNorm[{u}+E1[]]) ] ; In Scat ; Integration I1 ; Jacobian JVol ; } } }

      For i In {0:Nb_ordre-1}
        For j In {0:Nb_ordre-1}
          { Name int_x_t~{i}~{j} ; Value { Integral { [ CompX[{u}+E1[] ]*expialphaxy~{i}~{j}[]/(period_x*period_y) ] ; In SurfIntBot ; Integration I1 ; Jacobian JSur ; } } }
          { Name int_y_t~{i}~{j} ; Value { Integral { [ CompY[{u}+E1[] ]*expialphaxy~{i}~{j}[]/(period_x*period_y) ] ; In SurfIntBot ; Integration I1 ; Jacobian JSur ; } } }
          { Name int_x_r~{i}~{j} ; Value { Integral { [ CompX[{u}+E1d[]]*expialphaxy~{i}~{j}[]/(period_x*period_y) ] ; In SurfIntTop ; Integration I1 ; Jacobian JSur ; } } }
          { Name int_y_r~{i}~{j} ; Value { Integral { [ CompY[{u}+E1d[]]*expialphaxy~{i}~{j}[]/(period_x*period_y) ] ; In SurfIntTop ; Integration I1 ; Jacobian JSur ; } } }
          { Name eff_t~{i}~{j}   ; Value { Term{Type Global; [
                 1/(gammat~{i}~{j}[]*-k1z[]*Cos[xsi]^2) * ( (gammat~{i}~{j}[]^2+alpha~{i}~{j}[]^2)*SquNorm[$int_x_t~{i}~{j}]+
                                                            (gammat~{i}~{j}[]^2+ beta~{j}~{j}[]^2)*SquNorm[$int_y_t~{i}~{j}]+
                                                             2*alpha~{i}~{j}[]*beta~{i}~{j}[]*Re[$int_x_t~{i}~{j}*Conj[$int_y_t~{i}~{j}]] ) ] ; In SurfIntBot ; } } }
          { Name eff_r~{i}~{j}   ; Value { Term{Type Global; [
                 1/(gammar~{i}~{j}[]*-k1z[]*Cos[xsi]^2) * ((gammar~{i}~{j}[]^2+alpha~{i}~{j}[]^2)*SquNorm[$int_x_r~{i}~{j}]+
                                                           (gammar~{i}~{j}[]^2+ beta~{i}~{j}[]^2)*SquNorm[$int_y_r~{i}~{j}]+
                                                            2*alpha~{i}~{j}[]*beta~{i}~{j}[]*Re[$int_x_r~{i}~{j}*Conj[$int_y_r~{i}~{j}]]) ] ; In SurfIntTop ; } } }
          { Name numbering_ij~{i}~{j}   ; Value { Term{Type Global; [Vector[i-Nmax,j-Nmax,0]] ; In SurfIntBot ; } } }
        EndFor
      EndFor
    }
  }
}

PostOperation {
  { Name postop_helmholtz_vector; NameOfPostProcessing postpro_helmholtz_vector ;
    Operation {
      If (FlagOutEscaFull==1)
        If (Flag_interp_cubic==1)
          Print [ u , OnBox { {-period_x/2,-period_y/2,hh_L_6-PML_bot} {period_x/2,-period_y/2,hh_L_6-PML_bot} {-period_x/2,period_y/2,hh_L_6-PML_bot} {-period_x/2,-period_y/2,hh_L_1+thick_L_1+PML_top} } {npts_interpX,npts_interpY,npts_interpZSca} , File StrCat[myDir,"u_grid.pos"], Name "u_grid"];
        Else
          Print [ u , OnElementsOf Omega, File StrCat[myDir,"Etot.pos"]];
        EndIf
      EndIf
      If (FlagOutEtotFull==1)
        If (Flag_interp_cubic==1)
          Print [ Etot , OnBox { {-period_x/2,-period_y/2,hh_L_6} {period_x/2,-period_y/2,hh_L_6} {-period_x/2,period_y/2,hh_L_6} {-period_x/2,-period_y/2,hh_L_1+thick_L_1} } {npts_interpX,npts_interpY,npts_interpZTot} , File StrCat[myDir,"Etot_grid.pos"], Name "Etot_grid"];
        Else
          Print [ Etot , OnElementsOf Omega_plot, File StrCat[myDir,"Etot.pos"]];
        EndIf
      EndIf
      If (FlagOutPoyFull==1)
        If (Flag_interp_cubic==1)
          Print [ Poy_tot , OnBox { {-period_x/2,-period_y/2,hh_L_6} {period_x/2,-period_y/2,hh_L_6} {-period_x/2,period_y/2,hh_L_6} {-period_x/2,-period_y/2,hh_L_1+thick_L_1} } {npts_interpX,npts_interpY,npts_interpZTot} , File StrCat[myDir,"Poy_tot_grid.pos"], Name "Poy_tot_grid"];
        Else
          Print [ Poy_tot , OnElementsOf Omega_plot, File StrCat[myDir,"Poy_tot.pos"]];
        EndIf
      EndIf
      If (FlagOutEscaCuts==1)
        Print [ u , OnPlane { {-period_x/2,0,hh_L_6-PML_bot} {period_x/2,0,hh_L_6-PML_bot} {-period_x/2,0,hh_L_1+thick_L_1+PML_top} } {npts_interpX,npts_interpZSca} , File StrCat[myDir,"u_cut_Y=0.pos"], Name "u_cut_Y=0"];
        Print [ u , OnPlane { {0,-period_y/2,hh_L_6-PML_bot} {0,period_y/2,hh_L_6-PML_bot} {0,-period_y/2,hh_L_1+thick_L_1+PML_top} } {npts_interpY,npts_interpZSca} , File StrCat[myDir,"u_cut_X=0.pos"], Name "u_cut_X=0"];
      EndIf
      If (FlagOutEtotCuts==1)
        Print [ Etot , OnPlane { {-period_x/2,0,hh_L_6} {period_x/2,0,hh_L_6} {-period_x/2,0,hh_L_1+thick_L_1} } {npts_interpX,npts_interpZTot} , File StrCat[myDir,"Etot_cut_Y=0.pos"], Name "Etot_cut_Y=0"];
        Print [ Etot , OnPlane { {0,-period_y/2,hh_L_6} {0,period_y/2,hh_L_6} {0,-period_y/2,hh_L_1+thick_L_1} } {npts_interpY,npts_interpZTot} , File StrCat[myDir,"Etot_cut_X=0.pos"], Name "Etot_cut_X=0"];
      EndIf
      If (FlagOutHtotCuts==1)
        Print [ Htot , OnPlane { {-period_x/2,0,hh_L_6} {period_x/2,0,hh_L_6} {-period_x/2,0,hh_L_1+thick_L_1} } {npts_interpX,npts_interpZTot} , File StrCat[myDir,"Htot_cut_Y=0.pos"], Name "Htot_cut_Y=0"];
        Print [ Htot , OnPlane { {0,-period_y/2,hh_L_6} {0,period_y/2,hh_L_6} {0,-period_y/2,hh_L_1+thick_L_1} } {npts_interpY,npts_interpZTot} , File StrCat[myDir,"Htot_cut_X=0.pos"], Name "Htot_cut_X=0"];
      EndIf
      If (FlagOutPoyCut==1)
        Print [ Poy_tot , OnPlane { {-period_x/2,0,hh_L_6} {period_x/2,0,hh_L_6} {-period_x/2,0,hh_L_1+thick_L_1} } {npts_interpX,npts_interpZTot} , File StrCat[myDir,"Poy_tot_cut_Y=0.pos"], Name "Poy_tot_cut_Y=0"];
        Print [ Poy_tot , OnPlane { {0,-period_y/2,hh_L_6} {0,period_y/2,hh_L_6} {0,-period_y/2,hh_L_1+thick_L_1} } {npts_interpY,npts_interpZTot} , File StrCat[myDir,"Poy_tot_cut_X=0.pos"], Name "Poy_tot_cut_X=0"];
      EndIf
      // // // // For DEBUG
      // Print [ epsr_xx , OnElementsOf Omega, File StrCat[myDir,"epsr_xx.pos"]];
      // Print [ H1x , OnElementsOf Omega, File StrCat[myDir,"H1x.pos"]];
      // Print [ Etot , OnElementsOf Omega, File StrCat[myDir,"Etot.pos"]];
      // Print [ u , OnElementsOf Omega, File StrCat[myDir,"u.pos"]];
      // Print [ E1 , OnElementsOf Omega, File StrCat[myDir,"E1.pos"]];
      // Print [ Htotx , OnElementsOf Omega, File StrCat[myDir,"Htotx.pos"]];
      // Print [ uper   , OnElementsOf Omega, File StrCat[myDir,"uper.pos"]];
      // Print [ uperx  , OnElementsOf Omega, File StrCat[myDir,"uperx.pos"]];
      // Print [ upery  , OnElementsOf Omega, File StrCat[myDir,"upery.pos"]];
      // Print [ uperz  , OnElementsOf Omega, File StrCat[myDir,"uperz.pos"]];

      // Print [ E1perx , OnElementsOf Omega, File StrCat[myDir,"E1perx.pos"]];
      // Print [ E1pery , OnElementsOf Omega, File StrCat[myDir,"E1pery.pos"]];
      // Print [ Etotpery , OnElementsOf Omega, File StrCat[myDir,"Etotpery.pos"]];
      // Print [ Poy_tot , OnPlane { {0.5*(-period_x-dys), -dyc/2,(hh_L_6+hh_L_5)/2}
      //                             {0.5*( period_x-dys), -dyc/2,(hh_L_6+hh_L_5)/2}
      //                             {0.5*(-period_x+dys),  dyc/2,(hh_L_6+hh_L_5)/2} }
      //                             {npts_checkpoyX-1,npts_checkpoyY-1} , File StrCat[myDir,"Poy_tot.pos"]];
      // Print [ Poy_ref , OnPlane { {0.5*(-period_x-dys), -dyc/2, hh_L_1+thick_L_1/2}
      //                             {0.5*( period_x-dys), -dyc/2, hh_L_1+thick_L_1/2}
      //                             {0.5*(-period_x+dys),  dyc/2, hh_L_1+thick_L_1/2} }
      //                             {npts_checkpoyX-1,npts_checkpoyY-1} , File StrCat[myDir,"Poy_ref.pos"]];
      // Print [ Poy_inc , OnPlane { {0.5*(-period_x-dys), -dyc/2, hh_L_1+thick_L_1/2}
      //                             {0.5*( period_x-dys), -dyc/2, hh_L_1+thick_L_1/2}
      //                             {0.5*(-period_x+dys),  dyc/2, hh_L_1+thick_L_1/2} }
      //                             {npts_checkpoyX-1,npts_checkpoyY-1} , File StrCat[myDir,"Poy_inc.pos"]];
      //////////// END DEBUG  /////////

      Print [ Poy_tot , OnPlane { {0.5*(-period_x-dys), -dyc/2,(hh_L_6+hh_L_5)/2}
                                  {0.5*( period_x-dys), -dyc/2,(hh_L_6+hh_L_5)/2}
                                  {0.5*(-period_x+dys),  dyc/2,(hh_L_6+hh_L_5)/2} }
                                  {npts_checkpoyX-1,npts_checkpoyY-1} , File StrCat[myDir,"Poy_tot_gd.pos"], Format Table];
      Print [ Poy_ref , OnPlane { {0.5*(-period_x-dys), -dyc/2, hh_L_1+thick_L_1/2}
                                  {0.5*( period_x-dys), -dyc/2, hh_L_1+thick_L_1/2}
                                  {0.5*(-period_x+dys),  dyc/2, hh_L_1+thick_L_1/2} }
                                  {npts_checkpoyX-1,npts_checkpoyY-1} , File StrCat[myDir,"Poy_ref_gd.pos"], Format Table];
      Print [ Poy_inc , OnPlane { {0.5*(-period_x-dys), -dyc/2, hh_L_1+thick_L_1/2}
                                  {0.5*( period_x-dys), -dyc/2, hh_L_1+thick_L_1/2}
                                  {0.5*(-period_x+dys),  dyc/2, hh_L_1+thick_L_1/2} }
                                  {npts_checkpoyX-1,npts_checkpoyY-1} , File StrCat[myDir,"Poy_inc_gd.pos"], Format Table];

      For k In {2:6}
        Print[ Abs_L~{k}[L~{k}], OnGlobal, File > StrCat[myDir,Sprintf["temp-Q_L_%g.txt",k]], Format Table ];
      EndFor
      Print[ Abs_scat[Scat]  , OnGlobal, File > StrCat[myDir,"temp-Q_scat.txt"], Format Table ];

      For i In {0:Nb_ordre-1}
        For j In {0:Nb_ordre-1}
          Print[ int_x_t~{i}~{j}[SurfIntBot], OnGlobal, StoreInVariable $int_x_t~{i}~{j}, Format Table];
          Print[ int_y_t~{i}~{j}[SurfIntBot], OnGlobal, StoreInVariable $int_y_t~{i}~{j}, Format Table];
          Print[ int_x_r~{i}~{j}[SurfIntTop], OnGlobal, StoreInVariable $int_x_r~{i}~{j}, Format Table];
          Print[ int_y_r~{i}~{j}[SurfIntTop], OnGlobal, StoreInVariable $int_y_r~{i}~{j}, Format Table];
        EndFor
      EndFor

      For i In {0:Nb_ordre-1}
        For j In {0:Nb_ordre-1}
          Print[ eff_t~{i}~{j}[SurfIntBot], OnRegion SurfIntBot, File > StrCat[myDir, "eff_t.txt"], Format Table ];
          Print[ eff_r~{i}~{j}[SurfIntTop], OnRegion SurfIntTop, File > StrCat[myDir, "eff_r.txt"], Format Table ];
        EndFor
      EndFor
      For i In {0:Nb_ordre-1}
        For j In {0:Nb_ordre-1}
          Print[ numbering_ij~{i}~{j}[SurfIntBot], OnRegion SurfIntBot, File > StrCat[myDir,"numbering_ij.txt"], Format Table ];
        EndFor
      EndFor
      Print[ lambda_step, OnPoint{0,0,0}, Format Table, File > StrCat[myDir, "temp_lambda_step.txt"], SendToServer "GetDP/Lambda_step" ] ;
    }
  }
}

DefineConstant[
  R_ = {"helmholtz_vector", Name "GetDP/1ResolutionChoices", Visible 1},
  C_ = {"-solve -pos -petsc_prealloc 500 -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -ksp_error_if_not_converged", Name "GetDP/9ComputeCommand", Visible 1},
  P_ = {"postop_helmholtz_vector", Name "GetDP/2PostOperationChoices", Visible 1}
];

// Notes on possible mumps craches (NaN's when large nDOF)
// -ksp_view
// -ksp_error_if_not_converged
// -ksp_monitor_true_residual
// related to? https://lists.mcs.anl.gov/mailman/htdig/petsc-users/2018-October/036377.html

// N.B. pastix works with petsc 3.12
// -pc_factor_mat_solver_type pastix
