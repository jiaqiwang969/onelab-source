/* -------------------------------------------------------------------
   Tutorial: Optimization of ferromagnetic C-core

   Features:
   - Nonlinear magnetostatic model of a C-core
   - Imposed current density or imposed current
   - Air gap and core width as design parameters
   - Prescribed By in the air gap as objective function 
   - Direct approach of the shape sensitivity analysis
   - Sensivity analysis with Lie derivative

   To compute the solution in a terminal:
       getdp shape.pro -solve GetPerformancesAndGradients

   To compute the solution interactively from the Gmsh GUI:
       File > Open > shape.pro
       Run (button at the bottom of the left panel)
   ------------------------------------------------------------------- */

/*
  This tutorial deals with shape sensitivity analysis in Magnetostatics.
  The model consists of a ferromagnetic C-core powered
  by a coil with 'NbTurns' windings suplied with a DC current. 
  The ferromagnetic behavior of the C-core material can be 
  represented by either a linear of a nonlinear saturable law. 
  The power supply is set to maintain either a constant total current
  or a constant current density as the design parameters vary. 
  The optimisation problem aims at having a prescribed value
  for the By component of the magnetic flux density
  in the middle of the air gap.
  The two design variables selected for achieving this objective are
  - the air gap width 'D',
  - the core width 'E'.

  The 'shp.py' is organised so as to allow either 
  - to check the correctness of the sensitivity computation,
  - effectively solve the optimisation problem.

  The computation of sensitivities is a delicate step of the optimisation,
  which can however be fully validated by finite difference.
  This tutorial is thus dedicated to both the explanation and verification
  of the direct computation of sensitivities in Magnetostatics.

  Communication between the optimisation loop ('shp.py')
  and the Onelab model is done via the Onelab database. 
  Design variables (here 'D' and 'E') are defined as Onelab variables 
  whose name has the reserved prefix 'Optimization/Parameters'.
  Similarly, the unknown field 'a', velocities 'velocity_*'
  and sensibilities 'dwdtau_*' are also exchanged in memory
  between different functions or steps of the model
  as Onelab variables with reserved prefixes 'Optimization/Results/*'.

  Finally, an extra mechanism is implemented 
  to apply a temporary small perturbation of a design parameter 
  in order to compute the associated velocity field
  without interfering with the value stored in the Onelab database.
  This affectation is done with a '-setnumber' line argument 
  when calling Gmsh. 
  In consequence, the definition of the design variables
  as Onelab parameters must be enclosed in a 'DefineConstant[]' statement.
  Remember that affectations within a 'DefineConstant[]'
  are ignored if the variable already exists because,
  and in particular if it has been set by a '-setnumber' line argument.
 */

Include "ccore_common.pro";

DefineConstant[
  R_ = {"GetPerformances", Name "GetDP/1ResolutionChoices", Visible 0}
  C_ = {"-solve -v 5 -v2", Name "GetDP/9ComputeCommand", Visible 0}
  P_ = {"", Name "GetDP/2PostOperationChoices", Visible 0}
  VelocityTag = -1
  OptiIterNumber = 0
];

// Onelab parameters of the C-core model 

NL_tol_abs = 1e-8; 	// absolute tolerance on residual for noninear iterations
NL_tol_relax = 1.0; 	// relaxation on residual for noninear iterations
NL_iter_max = 50; 	// maximum number of noninear iterations
Flag_NL = 
  DefineNumber[ 1, Name "Model/Parameters/Non linear core material", Choices {0,1}];
murCore = 
  DefineNumber[300, Name "Model/Parameters/Mur core", Visible 1,//!Flag_NL,
			 Help "Magnetic relative permeability of Core"];
Flag_Jfixed = 
  DefineNumber[ 1, Name "Model/Parameters/Fixed current density", Choices {0,1}];
CurrentDensity = 1e6 *
  DefineNumber[1, Name "Model/Parameters/Current Density", Visible Flag_Jfixed,
	       Help "Current density in coil [A/mm2]"];
NbTurns = 300;
Current = 
  DefineNumber[CurrentDensity*CoilSection_ref/NbTurns, 
	       Name "Model/Parameters/Current [A]", Visible !Flag_Jfixed,
	       Help "Current in coil [A]"];

CoilSection = G*(B-2*E-2*F);
If(Flag_Jfixed)
  Mmf = CurrentDensity * CoilSection;
  Current = Mmf / NbTurns;
Else
  Mmf = NbTurns * Current;
  CurrentDensity = Mmf / CoilSection;
EndIf



Group {
  // Physical regions (in capital letters):
  AIR    = Region[ 1 ];   
  CORE   = Region[ 2 ];
  COILP  = Region[ 3 ];   
  COILN  = Region[ 4 ];   
  NOFLUX = Region[ 11 ];

  // Abstract regions
  Vol_Mag     = Region[ {AIR, CORE, COILP, COILN} ];
  Vol_S_Mag   = Region[ {COILP, COILN} ];
  Sur_Dir_Mag = Region[ {NOFLUX} ];
  Sur_Neu_Mag = Region[ {} ];

  Vol_NL_Mag = Region[ {} ];
  If(Flag_NL)
    Vol_NL_Mag = Region[ {CORE} ];
  EndIf
  Vol_L_Mag  = Region[ {Vol_Mag,-Vol_NL_Mag} ];
}

Function {
  mu0 = 4.e-7 * Pi ;
  nu0 = 1. / mu0 ;

  nu [ Region[{AIR, COILP, COILN}] ] = nu0;
 
  Mat_h_r = { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 250, 
	      300, 400, 500, 600,  700, 800, 900, 1000, 1250, 1500, 2000, 2500, 5000, 
	      7500,  10000, 15000, 20000, 59000, 174000, 514000, 1520000, 4470000, 
	      13200000, 38900000, 115000000, 339000000, 1000000000 } ;
  Mat_b_r = { 0.0, 0.194880829963, 0.377143018857, 0.537767739762, 0.672888260835, 
	      0.783043000477, 0.871342430831,0.941778611986, 0.998183303557, 1.04378111223, 
	      1.08110469369, 1.14963767549, 1.19607212343, 1.22964695907, 1.25515221835,
	      1.29162498935, 1.31678879432, 1.35015120537, 1.37220092877, 1.38859114656, 
	      1.4017440574, 1.41287024565, 1.42264180514, 1.43146158921, 1.45082466146, 
	      1.46784549989, 1.49819370601, 1.52578650709, 1.64314027719, 1.73458485332, 
	      1.8039068939,1.89568786291, 1.95213815187, 2.1390774927, 2.45827909293, 
	      3.32303272825, 5.85485500678, 13.2701832298, 35.2114648741, 99.8027446541, 
	      291.062951228, 854.036370229, 2515.3105707 } ;

  Mat_b2_r = Mat_b_r()^2;
  Mat_nu_r = Mat_h_r()/Mat_b_r();
  Mat_nu_r(0) = Mat_nu_r(1);
  Mat_nu_b2_r = ListAlt[Mat_b2_r(), Mat_nu_r()];
  nu_interp_r[] = InterpolationLinear[ SquNorm[$1] ]{Mat_nu_b2_r()};
  dnudb2_interp_r[] = dInterpolationLinear[SquNorm[$1]]{Mat_nu_b2_r()};
  h_interp_r[] = nu_interp_r[$1] * $1 ;

  If(!Flag_NL)
    nu [ CORE ]  = nu0/murCore;
    dhdb_NL[ CORE ] = 0;
    mu_analytic = mu0*murCore;
  Else
    nu [ Vol_NL_Mag ] = nu_interp_r[$1] ;
    dhdb_NL[ Vol_NL_Mag ] = 2*dnudb2_interp_r[$1#2]*SquDyadicProduct[#2];
    mu_analytic = 1/Mat_nu_r(0);
  EndIf

 // shape function of J : js[] = nI * Jshape[]
  Jshape[ COILP ] = Vector[0,0,+1]/CoilSection;
  Jshape[ COILN ] = Vector[0,0,-1]/CoilSection;
  js[] = Mmf * Jshape[];

  // Analytical solution
  RelGap = D/(mu0*E);
  RelCore = (2*(A+B-2*E)-D)/(mu_analytic*E);
  Flux = NbTurns*Mmf/(RelGap+RelCore); // analytic approximation
  SetNumber("Model/Results/Flux analytic", Flux);
  //Flux = 0;
  If( SweepOnParameter == 1 ) // D
    dFluxdTau = Flux/E*Flux/(NbTurns*Mmf)*(1/mu_analytic-1/mu0);
  ElseIf( SweepOnParameter == 2 ) // E
    dFluxdTau = Flux/E*(1.+4*Flux/(NbTurns*Mmf*mu_analytic));
    If( Flag_Jfixed )
      dFluxdTau -= 2*Flux*G/CoilSection;
    EndIf
  ElseIf( SweepOnParameter == 3 ) // F
    dFluxdTau = 0;
  Else
    dFluxdTau = 0;
  EndIf


  VV[] = Vector[$1,$2,0*$3] ;
  dV[] = TensorV[$1,$2,0*$3] ;
  Lie2form[] = TTrace[$2]*$1 - Transpose[$2]*$1; // $1=2-form (vector field), $2=dV (tensor)

  If(Flag_Jfixed)
    LieOf_js[] = Lie2form[js[], $1]; // $1=dV
  Else
    LieOf_js[] = 0;
  EndIf

  // Lie derivative of H(B) where B is a 2-form and H a 1-form ($1:{d a}, $2:dV)
  LieOf_HB[] = nu[$1] * (Transpose[$2] * $1 - TTrace[$2] * $1 + $2 * $1) ;
  LieOf_HB_NL[] = dhdb_NL[$1] * (Transpose[$2] - TTrace[$2] * TensorDiag[1,1,1]) * $1;
}



Group {
  Dom_Hcurl_a_Mag_2D = Region[ {Vol_Mag, Sur_Neu_Mag} ];
}

Function{
  l_a = ListFromServer["Optimization/Results/a"];
  aFromServer[] = ValueFromIndex[]{l_a()};
  For i In {1:3}
    l_v~{i} = ListFromServer[Sprintf["Optimization/Results/velocity_%g_%g",VelocityTag,i]];
    velocity~{i}[] = ValueFromIndex[]{l_v~{i}()};
  EndFor
}

Constraint {
  { Name Dirichlet_a_Mag;
    Case {
      { Region NOFLUX; Type Assign; Value 0; }
    }
  }
  { Name aFromServer;
    Case {
      { Region Vol_Mag; Type Assign; Value aFromServer[]; }
      { Region NOFLUX; Type Assign; Value 0; }
    }
  }
  For i In {1:3}
    { Name velocity~{i} ;
      Case {
        { Region Vol_Mag ; Value velocity~{i}[]; }
      }
    }
  EndFor
}



FunctionSpace {
  { Name Hcurl_a_2D; Type Form1P;
    BasisFunction {
      { Name se; NameOfCoef ae; Function BF_PerpendicularEdge;
        Support Dom_Hcurl_a_Mag_2D ; Entity NodesOf[ All ]; }
    }
    Constraint {
      { NameOfCoef ae; EntityType NodesOf;
        NameOfConstraint Dirichlet_a_Mag; }
    }
  }

  { Name Hcurl_a_2D_fullyfixed; Type Form1P;
    BasisFunction{
      { Name se1; NameOfCoef ae1; Function BF_PerpendicularEdge;
        Support Vol_Mag; Entity NodesOf[All]; }
    }
    Constraint{
      { NameOfCoef ae1; EntityType NodesOf; NameOfConstraint aFromServer; }
    }
  }
  For i In {1:3}
    { Name H_v~{i} ; Type Form0;
      BasisFunction {
        { Name sn ; NameOfCoef un ; Function BF_Node ;
          Support Vol_Mag; Entity NodesOf[ All ] ; }
      }
      Constraint {
        { NameOfCoef un ; EntityType NodesOf ; NameOfConstraint velocity~{i}; }
      }
    }
  EndFor
}

Jacobian {
  { Name Vol ;
    Case { 
      { Region All ; Jacobian Vol ; }
    }
  }
  { Name Sur;
    Case {
      { Region All; Jacobian Sur; }
    }
  }
}

Integration {
  { Name Int ;
    Case { { Type Gauss ;
	Case {
          { GeoElement Point; NumberOfPoints  1; }
          { GeoElement Line; NumberOfPoints  5; }
          { GeoElement Triangle; NumberOfPoints  4; }
          { GeoElement Quadrangle; NumberOfPoints  4; } 
	}
      }
    }
  }
}


// -------------------------------------------------------------------------
// This resolution solves the direct problem to compute the sensitivity
// of the induction flux with respect to a given design variable.


Formulation {
  { Name MagSta_a; Type FemEquation;
    Quantity {
      { Name a; Type Local; NameOfSpace Hcurl_a_2D; }
    }
    Equation {
      Integral { [ nu[{d a}] * Dof{d a} , {d a} ];
	In Vol_Mag ; Jacobian Vol; Integration Int; }
      Integral { JacNL [ dhdb_NL[{d a}] * Dof{d a} , {d a} ];
	In Vol_NL_Mag; Jacobian Vol; Integration Int; }
      Integral { [ -js[] , {a} ];
	In Vol_S_Mag; Jacobian Vol; Integration Int; }
    }
  }
}

Resolution {
  { Name GetPerformances;
    System {
      { Name SYS; NameOfFormulation MagSta_a;}
    }
    Operation {
      If( OptiIterNumber == 1 )
        CreateDir[ResDir];
        DeleteFile[StrCat[ResDir,"w.txt"]];
        DeleteFile[StrCat[ResDir,"Lie_w.txt"]];
      EndIf
      InitSolution[SYS];
      IterativeLoop[NL_iter_max, NL_tol_abs, NL_tol_relax]{
        GenerateJac[SYS];
        SolveJac[SYS];
      }
      PostOperation[Get_ObjectiveConstraints];
    }
  }
}

PostProcessing {
  { Name MagSta_a_2D; NameOfFormulation MagSta_a;
    Quantity {
      { Name az;
        Value {
          Term { [ CompZ[{a}] ]; In Dom_Hcurl_a_Mag_2D; Jacobian Vol; }
        }
      }
      { Name b;
        Value {
          Term { [ {d a} ]; In Dom_Hcurl_a_Mag_2D; Jacobian Vol; }
        }
      }
      { Name by;
        Value {
          Term { [ CompY[{d a}] ]; In Dom_Hcurl_a_Mag_2D; Jacobian Vol; }
        }
      }
      { Name js;
        Value {
          Term { [ js[] ]; In Vol_S_Mag; Jacobian Vol; }
        }
      }
      { Name Flux ; Value { Integral { 
	    [ {a} * NbTurns * Jshape[] ] ;
	    In Vol_S_Mag  ; Jacobian Vol ; Integration Int ; } }
      }
      { Name fobj_b_at_point; 
        Value {
          Term { [ SquNorm[ CompY[{d a}] - By_target ] ]; 
	    In Dom_Hcurl_a_Mag_2D; Jacobian Vol; }
        }
      }
    }
  }
}

PostOperation {
  { Name Get_ObjectiveConstraints; NameOfPostProcessing MagSta_a_2D;
    Operation{

      CreateDir["res"];
      //Print[js, OnElementsOf Vol_S_Mag, File "res/js.pos"];
      //Print[b, OnElementsOf Vol_Mag, File "res/b.pos"];

      Print[az, OnElementsOf Vol_Mag, File "res/az.pos"];
      Echo[ Str["k= PostProcessing.NbViews-1;",
		"View[k].IntervalsType = 3;",
		"View[0].NbIso = 30;"], 
	    File "tmp.geo", LastTimeStepOnly] ;

      // In memory transfer of the a field for assembling sensitivity equations
      Print[az, OnElementsOf Vol_Mag, Format NodeTable, File "",
	    SendToServer "Optimization/Results/a", Hidden 1];

     Print[ Flux[ Vol_S_Mag ], OnGlobal, 
	     Format Table, File > "", 
	     SendToServer "Model/Results/Flux computed"];

     Print[ by, OnPoint { Px_target, Py_target, 0}, 
	    Format Table, File "", 
	    SendToServer "Model/Results/b probe"];
     
     Print[ fobj_b_at_point, OnPoint { Px_target, Py_target, 0}, 
	    Format Table, File > "res/w.txt", 
	    SendToServer "Optimization/Results/w"];
    }
  }
}

// -------------------------------------------------------------------------
// This resolution solves the direct problem to compute the sensitivity
// of the induction flux with respect to a given design variable.

Formulation {
  { Name DirectFormulationOf_MagSta_a_2D; Type FemEquation ;
    Quantity {
      { Name Lva; Type Local; NameOfSpace Hcurl_a_2D; }
      { Name a; Type Local; NameOfSpace Hcurl_a_2D_fullyfixed; }
      For i In {1:3}
        { Name v~{i} ; Type Local ; NameOfSpace H_v~{i}; }
      EndFor
    }
    Equation {
      Integral { [ nu[{d a}] * Dof{d Lva}, {d Lva} ];
        In Vol_Mag; Jacobian Vol; Integration Int; }
      Integral { [ dhdb_NL[{d a}] * Dof{d Lva}, {d Lva} ];
        In Vol_NL_Mag; Jacobian Vol; Integration Int; }
      Integral { [LieOf_HB[{d a},dV[{d v_1},{d v_2},{d v_3}]], {d Lva}];
        In Vol_Mag; Jacobian Vol; Integration Int; }
      Integral { [LieOf_HB_NL[{d a},dV[{d v_1},{d v_2},{d v_3}]], {d Lva}];
        In Vol_NL_Mag; Jacobian Vol; Integration Int; }

      Integral { [ -LieOf_js[dV[{d v_1},{d v_2},{d v_3}]], {Lva} ] ;
        In Vol_S_Mag ; Jacobian Vol; Integration Int; }

      Integral { [ 0*Dof{a}, {a} ] ;
        In Vol_Mag; Jacobian Vol ; Integration Int ; }
      For i In {1:3}
        Integral { [ 0*Dof{v~{i}}, {v~{i}} ] ;
          In Vol_Mag; Jacobian Vol ; Integration Int ; }
      EndFor
    }
  }
}

Resolution {
  { Name GetGradient_wrt_dv;
    System {
      { Name DIR; NameOfFormulation DirectFormulationOf_MagSta_a_2D; }
    }
    Operation {
      InitSolution[DIR];
      Generate[DIR];
      Solve[DIR];
      PostOperation[Get_GradOf_w];
    }
  }
}

PostProcessing{
  { Name Direct_MagSta; NameOfFormulation DirectFormulationOf_MagSta_a_2D;
    PostQuantity {
      { Name az; 
	Value{ 
	  Term{ [ CompZ[{a}] ]; In Vol_Mag; Jacobian Vol; } 
	}
      }
      { Name VV; 
	Value{ 
	  Term{ [ Vector[{v_1},{v_2},{v_3}] ]; In Vol_Mag ; Jacobian Vol;}
	}
      }
      { Name Lie_az; 
	Value { 
	  Term { [ CompZ[ {Lva}] ] ; In Vol_Mag ; Jacobian Vol; }
	}
      }
      { Name Lie_b_at_point ; Value {
          Term {  [ 2 * (CompY[{d a}] - By_target) * CompY[{d Lva}] ] ; 
	    In Vol_Mag; Jacobian Vol;  }
	}
      }
    }
  }
}

PostOperation{
  { Name Get_GradOf_w; NameOfPostProcessing Direct_MagSta;
    Operation {
      CreateDir["res"];

      Print[ VV, OnElementsOf Vol_Mag,File "res/velocity.pos" ];
      Echo[ Str["View[PostProcessing.NbViews-1].GlyphLocation = 2;"], 
	    File "tmp.geo", LastTimeStepOnly] ;

      //Print[az, OnElementsOf Vol_Mag, File "res/azstar.pos"];
      // Echo[ Str["k= PostProcessing.NbViews-1;",
      // 		"View[k].IntervalsType = 3;",
      // 		"View[0].NbIso = 30;"], 
      // 	    File "tmp.geo", LastTimeStepOnly] ;

      Print[ Lie_az, OnElementsOf Vol_Mag, File "res/Lie_az.pos" ];
      Echo[ Str["k= PostProcessing.NbViews-1;",
		"View[k].IntervalsType = 3;",
		"View[0].NbIso = 30;"], 
	    File "tmp.geo", LastTimeStepOnly] ;
    
      Print[ Lie_b_at_point, OnPoint{ Px_target, Py_target, 0},
	     Format Table, File > "res/Lie_w.txt",
	     SendToServer Sprintf["Optimization/Results/dwdtau_%g",VelocityTag]];

    }
  }
}
