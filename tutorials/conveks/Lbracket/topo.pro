Struct OPT::YOUNG_LAW [ Enum, simp, modifiedSimp, ramp, none];

DefineConstant [
  Opt_ResDir = StrCat[CurrentDir,"res_opt/"]
  Opt_ResDir_Onelab = "Optimization/Results/"
  Flag_PrintLevel = {1, Name "General/Verbosity", Visible 1}
  
  // Optimization parameters
  Optimizer = {0,Choices {0="Conveks - MMA",1="Conveks - GCMMA"},Name "Optimization/00Optimizer"}
  Opt_maxIter = {1600, Name "Optimization/01Max iterations"}
  Opt_maxChange = {0.01, Name "Optimization/02Max change"}
  Opt_isInnerIter = {0, Name "Optimization/03Inner Iteration", Visible 0, Choices {0, 1}}

  // Density field parameters 
  densityFieldInit = {0.5, Name "Optimization/3Density/0Inital value"}
  Opt_filter_exists = {1, Name "Optimization/3Density/1filter", Choices{0,1}}
  Opt_filter_radius = {0.008, Name "Optimization/3Density/2filter radius", Visible Opt_filter_exists}

  // Interpolation parameters
  Flag_opt_matlaw = {OPT::YOUNG_LAW.modifiedSimp,
    Choices {
      OPT::YOUNG_LAW.simp="simp",
      OPT::YOUNG_LAW.modifiedSimp="modified-simp",
      OPT::YOUNG_LAW.ramp="ramp",
      OPT::YOUNG_LAW.none="none"},
    Name "Optimization/5Young Interpolation/0Material law"}
  Flag_opt_matlaw_penal = {3.0,
    Name "Optimization/5Young Interpolation/1penalization",
    Visible (Flag_opt_matlaw == OPT::YOUNG_LAW.simp
       || Flag_opt_matlaw == OPT::YOUNG_LAW.ramp
       || Flag_opt_matlaw == OPT::YOUNG_LAW.modifiedSimp) }
  Flag_opt_stress_penal = {2.5,
    Name "Optimization/6Von-Mises parameters/1penalization of stress",
    Visible Flag_opt_matlaw != OPT::YOUNG_LAW.none}
  PARAM_pnorm = {8,Name "Optimization/6Von-Mises parameters/2Pnorm penalization"}

  // Model parameters
  Flag_EPC = 1
  iP = 1
  Flag_getGradient = 1
];

fobj = DefineNumber[0, Name "Optimization/Results/objective", 
	Graph "0300000000000000000000000000000000000"];
maxConstr = DefineNumber[0, Name "Optimization/Results/max(|Constraints|)", 
	Graph "0000000003000000000000000000000000000"];

Group {
  Sur_Clamp~{iP} = Region[{1001}];
  Vol_Force~{iP} = Region[{}];
  Sur_Force~{iP} = Region[{1003}];
  Vol_Elast~{iP} = Region[{1000,1004}];
  Domain~{iP} = Region[ { Vol_Elast~{iP},
			  Vol_Force~{iP},
                          Sur_Force~{iP},
			  Sur_Clamp~{iP}} ];

  // Define the group for topology optimization
  DomainDensity = Region[1000];
}

Function {
  YoungFull = 1.0; // N/m2
  PoissonFull = 0.33;
  rhoElastFull = 7874.0; // kg/m3

  Poisson[] = PoissonFull;
  rhoElast[Region[{Vol_Elast~{iP},-DomainDensity}]] = rhoElastFull;
  Young[Region[{Vol_Elast~{iP},-DomainDensity}]] = YoungFull;

  // Handling of the modified Young law - allowing to set up a continuous
  // problem with smoothly varying material properties between the "material"
  // and "void" regions.
  If(Flag_opt_matlaw == OPT::YOUNG_LAW.simp)
    Young[DomainDensity] = YoungFull * $1 ^ Flag_opt_matlaw_penal;
    dYoung[DomainDensity] = Flag_opt_matlaw_penal * YoungFull * $1 ^ (Flag_opt_matlaw_penal - 1);
  ElseIf(Flag_opt_matlaw == OPT::YOUNG_LAW.modifiedSimp)
    YoungMin = 1e-6;
    Young[DomainDensity] = YoungMin + (YoungFull-YoungMin) * $1 ^ Flag_opt_matlaw_penal;
    dYoung[DomainDensity] = Flag_opt_matlaw_penal * (YoungFull-YoungMin) * $1 ^ (Flag_opt_matlaw_penal - 1);
  ElseIf(Flag_opt_matlaw == OPT::YOUNG_LAW.ramp)
    YoungMin = 1.0e-6;
    Young[DomainDensity] = YoungMin + $1*(YoungFull-YoungMin)/(1.0+Flag_opt_matlaw_penal*(1.0-$1));
    dYoung[DomainDensity] = (YoungFull-YoungMin)*(1.0+Flag_opt_matlaw_penal)/(1.0+Flag_opt_matlaw_penal*(1.0-$1))^2.0;
  EndIf

  // Force for elastic model
  pressure_x[Region[1003]] = 0;
  pressure_y[Region[1003]] = -2.5;

  // Plane stress
  a[] = Young[$1]/(1.-Poisson[]^2); //$1:{xe}
  c[] = Young[$1]*Poisson[]/(1.-Poisson[]^2);
  a1[] = 1/(1.-Poisson[]^2); //$1:{xe}
  c1[] = Poisson[]/(1.-Poisson[]^2);
  b1[] = 1./2./(1.+Poisson[]);

  C_xx[] = Tensor[ (Young[$1])#9000*a1[],0  ,0  ,    0  ,#9000*b1[],0 , 0 ,0,#9000*b1[]];
  C_xy[] = Tensor[ 0  ,#9000*c1[],0  ,    #9000*b1[],0  ,0  ,    0  ,0  ,0   ];

  C_yx[] = Tensor[ 0  ,#9000*b1[],0  ,    #9000*c1[],0  ,0  ,    0  ,0  ,0   ];
  C_yy[] = Tensor[ #9000*b1[],0  ,0  ,    0  ,#9000*a1[],0  ,    0  ,0  ,#9000*b1[]];

  // macroscopic stress (voigt notation),
  // with $1:{d ux}, $2:{d uy}, $3:{xe}
  C_2d[] = a[$1] * TensorSym[1,Poisson[],0,  1,0,  (1-Poisson[])/2];
  epsElast[] = Vector[CompX[$1],CompY[$2],CompY[$1]+CompX[$2]];
  sigM[] = C_2d[1] * epsElast[$1,$2];

  // Von-Mises stress (scalar)
  // with $1:{d ux}, $2:{d uy}, $3:{xe}
  VMmat[] = TensorSym[1., -0.5, 0., 1., 0., 3.];
  sigVM[] = Sqrt[sigM[$1,$2]#1015 * (VMmat[] * #1015)];

  sigVMqp[#{Vol_Elast~{iP},-DomainDensity}] = sigVM[$1,$2];
  sigVMqp[DomainDensity] = ($3 ^ (Flag_opt_matlaw_penal-Flag_opt_stress_penal)) * sigVM[$1,$2];

  // $1:{xe}, $2:{d ux}, $3:{d uy}
  InternalEnergy[] = 0.5 * (C_xx[$1] * $2 + C_xy[$1] * $3) * $2
   + 0.5 * (C_yx[$1] * $2 + C_yy[$1] * $3) * $3;

  da1[] = 1 / (1 - Poisson[]^2);
  dc1[] = Poisson[] / (1 - Poisson[]^2);
  db1[] = 0.5 / (1 + Poisson[]);

  dC_xx[] = Tensor[ (dYoung[$1]#6001)* da1[],0,0,0, #6001*db1[],0, 0,0,#6001*db1[] ];
  dC_xy[] = Tensor[ 0,#6001*dc1[],0 ,#6001*db1[],0,0 ,0,0,0 ];

  dC_yx[] = Tensor[ 0, #6001*db1[],0, #6001*dc1[],0,0, 0,0,0 ];
  dC_yy[] = Tensor[ #6001*db1[],0,0, 0, #6001*da1[], 0,    0,0,#6001*db1[] ];

  // $1:{xe}, $2:{d ux}, $3:{d uy}, $4:{d lamx}, $5:{d lamy}
  GradPartialElast[] = (dC_xx[$1] * $2 + dC_xy[$1] * $3) * $4
    + (dC_yx[$1] * $2 + dC_yy[$1] * $3) * $5;

  uElast[] = Vector[$1,$2,0];
}

Jacobian {
  { Name Vol; Case { { Region All ; Jacobian Vol; } } }
  { Name Sur; Case { { Region All ; Jacobian Sur; } } }
}

Integration {
  { Name I1 ; Case {
      { Type Gauss ;
        Case {
          { GeoElement Triangle   ; NumberOfPoints  6 ; }
	  { GeoElement Quadrangle ; NumberOfPoints  4 ; }
	  { GeoElement Line       ; NumberOfPoints  13 ; }
        }
      }
    }
  }
}

// -------------------------------------------------------------------------

// Density field handling (through a fully-fixed function space defined on
// DomainDensity - a subset of ConductingDomain) and filtering. This is done on
// densities, defined ayt the barycenter of elements.
//
// Beware of the special 1-point Gauss integration rule, which is required to
// evaluate the densities precisely where the "Depth 0" post-operations generate
// the data!

Integration {
  { Name GaussOnePoint ; Case {
      { Type Gauss ;
        Case {
          { GeoElement Point       ; NumberOfPoints  1; }
          { GeoElement Line        ; NumberOfPoints  1; }
          { GeoElement Triangle    ; NumberOfPoints  1; }
          { GeoElement Quadrangle  ; NumberOfPoints  1; }
          { GeoElement Prism       ; NumberOfPoints  1; }
          { GeoElement Tetrahedron ; NumberOfPoints  1; }
          { GeoElement Hexahedron  ; NumberOfPoints  1; }
          { GeoElement Pyramid     ; NumberOfPoints  1; }
        }
      }
    }
  }
}

Function {
  l_xe() = ListFromServer["Optimization/Results/designVariable"];
  If(#l_xe() > 1)
    densityField[] = ValueFromIndex[]{l_xe()};
  Else
    densityField[] = densityFieldInit;
  EndIf
}

Constraint {
  { Name ConstraintDensity;
    Case {
      { Region DomainDensity; Value (densityField[]*ElementVol[]); }
    }
  }
}

FunctionSpace {
  { Name H_DensityField; Type Form3; // filtered - L2
    BasisFunction {
      { Name sxn ; NameOfCoef uxn ; Function BF_Volume ;
        Support DomainDensity; Entity VolumesOf[ All ] ; }
    }
    Constraint {
      { NameOfCoef uxn; EntityType VolumesOf; NameOfConstraint ConstraintDensity; }
    }
  }
}

Formulation {
  { Name DensityField ; Type FemEquation ;
    Quantity {
      { Name xe ; Type Local ; NameOfSpace H_DensityField;}
    }
    Equation {
      Galerkin { [ 0*Dof{xe}, {xe} ] ;
        In DomainDensity; Jacobian Vol; Integration GaussOnePoint; }
    }
  }
}

Resolution {
  { Name DensityField;
    System {
      { Name D; NameOfFormulation DensityField; }
    }
    Operation{
      Generate[D];
      // Solve[D]; // no need: no unknowns :-)
      PostOperation[DensityField];
    }
  }
}

PostProcessing {
  { Name DensityField; NameOfFormulation DensityField;
    PostQuantity {
      { Name xe; Value { Term { [{xe}]; In DomainDensity; Jacobian Vol;} } }
    }
  }
}

PostOperation {
  { Name DensityField; NameOfPostProcessing DensityField;
    Operation {
      Print[ xe, OnElementsOf DomainDensity, Depth 0,
        Format ElementTable, File "", Hidden 1,
        SendToServer StrCat[Opt_ResDir_Onelab,"designVariable"]];
      If(Flag_PrintLevel > 0)
        CreateDir[Opt_ResDir];
        Print[ xe, OnElementsOf DomainDensity,
          File StrCat[Opt_ResDir,"xeCurrentMesh.pos"]];
      EndIf
    }
  }
}

Function {
  l_filt() = ListFromServer["Optimization/Results/filterInput"];
  If(#l_filt() > 1)
    filtSource[] = ValueFromIndex[]{l_filt()};
  Else
    filtSource[] = ***;
  EndIf
}

Constraint {
  { Name psi0;
    Case {
      { Region DomainDensity; Value filtSource[] * ElementVol[]; }
    }
  }
}

FunctionSpace {
  { Name H_psi0; Type Form3;
    BasisFunction {
      { Name sn ; NameOfCoef un ; Function BF_Volume ;
        Support DomainDensity; Entity VolumesOf[ All ] ; }
    }
    Constraint {
      { NameOfCoef un; EntityType VolumesOf; NameOfConstraint psi0; }
    }
  }
  { Name H_psif ; Type Form0;
    BasisFunction {
      { Name sn ; NameOfCoef un ; Function BF_Node ;
        Support DomainDensity; Entity NodesOf[ All ] ; }
    }
  }
}

Formulation {
  { Name FilterDensity ; Type FemEquation ;
    Quantity {
      { Name psi0 ; Type Local ; NameOfSpace H_psi0;}
      { Name psif ; Type Local ; NameOfSpace H_psif;}
    }
    Equation {
      Galerkin { [ Dof{d psif} * Opt_filter_radius^2.0, {d psif} ] ;
        In DomainDensity; Jacobian Vol ; Integration GaussOnePoint; }
      Galerkin { [ Dof{psif}, {psif} ] ;
        In DomainDensity; Jacobian Vol; Integration GaussOnePoint; }

      Galerkin { [ -{psi0}, {psif} ] ;
        In DomainDensity; Jacobian Vol; Integration GaussOnePoint; }

      Galerkin { [ 0*Dof{psi0}, {psi0} ] ;
        In DomainDensity; Jacobian Vol; Integration GaussOnePoint ; }
    }
  }
}

Resolution {
  { Name FilterDensity;
    System {
      { Name F; NameOfFormulation FilterDensity; }
    }
    Operation{
      InitSolution[F];
      Generate[F];
      Solve[F];
      PostOperation[FilterDensity];
    }
  }
}

PostProcessing {
  { Name FilterDensity; NameOfFormulation FilterDensity;
    PostQuantity {
      { Name filtSource; Value { Term { [ {psi0} ]; In DomainDensity;Jacobian Vol;} } }
      { Name filtOut; Value { Term { [ {psif} ]; In DomainDensity; Jacobian Vol; } } }
    }
  }
}

PostOperation {
  { Name FilterDensity; NameOfPostProcessing FilterDensity;
    Operation {
      Print[ filtOut, OnElementsOf DomainDensity, Depth 0,
        Format ElementTable, File "", Hidden 1,
        SendToServer StrCat[Opt_ResDir_Onelab,"filterOutput"]];
      If(Flag_PrintLevel > 1)
        CreateDir[Opt_ResDir];
        Print[ filtOut, OnElementsOf DomainDensity, Depth 0, Format Table,
          File StrCat[Opt_ResDir,"filtOut.pos"] ];
      EndIf
      If(Flag_PrintLevel > 5)
        Print[ filtSource, OnElementsOf DomainDensity,
          File StrCat[Opt_ResDir,"filtSource.pos"]];
        Print[ filtOut, OnElementsOf DomainDensity,
          File StrCat[Opt_ResDir,"filtOut_v.pos"]];
      EndIf
    }
  }
}

// -------------------------------------------------------------------------

// Handling of the direct linear elastic problem, using the densities to
// interpolate the hook material law as well as the material density.

Constraint {
  { Name Displacement_x;
    Case {
      { Region Sur_Clamp~{iP} ; Type Assign ; Value 0; }
    }
  }
  { Name Displacement_y;
    Case {
       { Region Sur_Clamp~{iP} ; Type Assign ; Value 0; }
    }
  }
  { Name Displacement_z;
    Case {
       { Region Sur_Clamp~{iP} ; Type Assign ; Value 0; }
    }
  }
}

FunctionSpace {
  { Name H_ux_Mec ; Type Form0 ;
    BasisFunction {
      { Name sxn ; NameOfCoef uxn ; Function BF_Node ;
        Support Domain~{iP} ; Entity NodesOf[ All ] ; }
    }
    Constraint {
      { NameOfCoef uxn ;
        EntityType NodesOf ; NameOfConstraint Displacement_x ; }
    }
  }
  { Name H_uy_Mec ; Type Form0 ;
    BasisFunction {
      { Name syn ; NameOfCoef uyn ; Function BF_Node ;
        Support Domain~{iP} ; Entity NodesOf[ All ] ; }
    }
    Constraint {
      { NameOfCoef uyn ;
        EntityType NodesOf ; NameOfConstraint Displacement_y ; }
    }
  }
}

Formulation {
  { Name Elast_u ; Type FemEquation ;
    Quantity {
      { Name ux  ; Type Local ; NameOfSpace H_ux_Mec ; }
      { Name uy  ; Type Local ; NameOfSpace H_uy_Mec ; }
      { Name xe ; Type Local ; NameOfSpace H_DensityField;}
    }
    Equation {
      Galerkin { [ C_xx[{xe}] * Dof{d ux}, {d ux} ] ;
        In Vol_Elast~{iP} ; Jacobian Vol ; Integration I1 ; }
      Galerkin { [ C_xy[{xe}] * Dof{d uy}, {d ux} ] ;
        In Vol_Elast~{iP} ; Jacobian Vol ; Integration I1 ; }
      Galerkin { [ C_yx[{xe}] * Dof{d ux}, {d uy} ] ;
        In Vol_Elast~{iP} ; Jacobian Vol ; Integration I1 ; }
      Galerkin { [ C_yy[{xe}] * Dof{d uy}, {d uy} ] ;
        In Vol_Elast~{iP} ; Jacobian Vol ; Integration I1 ; }

      Galerkin { [ -pressure_x[] , {ux} ];
        In Sur_Force~{iP} ; Jacobian Sur ; Integration I1 ; }
      Galerkin { [ -pressure_y[] , {uy} ];
        In Sur_Force~{iP} ; Jacobian Sur ; Integration I1 ; }

      // DO NOT REMOVE!!!
      Galerkin { [ 0*Dof{xe}, {xe} ] ;
        In DomainDensity; Jacobian Vol ; Integration GaussOnePoint ; }
    }
  }
}

// -------------------------------------------------------------------------

// Handling of the adjoint linear elastic problems, depending on the chosen
// objective function.

FunctionSpace {
  { Name H_lamx_Mec ; Type Form0 ;
    BasisFunction {
      { Name sxn ; NameOfCoef uxn ; Function BF_Node ;
        Support Domain~{iP} ; Entity NodesOf[ All ] ; }
    }
    Constraint {
      { NameOfCoef uxn ;
        EntityType NodesOf ; NameOfConstraint Displacement_x ; }
    }
  }
  { Name H_lamy_Mec ; Type Form0 ;
    BasisFunction {
      { Name syn ; NameOfCoef uyn ; Function BF_Node ;
        Support Domain~{iP} ; Entity NodesOf[ All ] ; }
    }
    Constraint {
      { NameOfCoef uyn ;
        EntityType NodesOf ; NameOfConstraint Displacement_y ; }
    }
  }
}

Macro AdjointFormulationOf_Elast_u
  /*
    Generate the adjoint formulation of EM_Formulations
    for a given performance function, without the Frechet derivative terms.
    These are added by the user latter on.

    Input: _AdjointName_ -> the name of the adjoint formulation.
  */

  Formulation {
    { Name S2N[_AdjointName_] ; Type FemEquation ;
      Quantity {
        { Name ux  ; Type Local ; NameOfSpace H_ux_Mec ; }
        { Name uy  ; Type Local ; NameOfSpace H_uy_Mec ; }
        { Name lamx  ; Type Local ; NameOfSpace H_lamx_Mec ; }
        { Name lamy  ; Type Local ; NameOfSpace H_lamy_Mec ; }
        { Name xe ; Type Local ; NameOfSpace H_DensityField;}
      }
      Equation {
        Galerkin { [ Transpose[C_xx[{xe}]] * Dof{d lamx}, {d lamx} ] ;
          In Vol_Elast~{iP} ; Jacobian Vol ; Integration I1 ; }
        Galerkin { [ Transpose[C_xy[{xe}]] * Dof{d lamx}, {d lamy} ] ;
          In Vol_Elast~{iP} ; Jacobian Vol ; Integration I1 ; }
        Galerkin { [ Transpose[C_yx[{xe}]] * Dof{d lamy}, {d lamx} ] ;
          In Vol_Elast~{iP} ; Jacobian Vol ; Integration I1 ; }
        Galerkin { [ Transpose[C_yy[{xe}]] * Dof{d lamy}, {d lamy} ] ;
          In Vol_Elast~{iP} ; Jacobian Vol ; Integration I1 ; }

        // DO NOT REMOVE!!!
        Galerkin { [ 0*Dof{xe}, {xe} ] ;
          In DomainDensity; Jacobian Vol ; Integration GaussOnePoint ; }
      }
    }
  }
Return

Macro CreateAdjointFormulation
  // This Macro creates the adjoint formulation for a given perfromance function
  // and a given physical problem.
  //
  // Input:
  // * OPT_AdjointWithoutFrechet: name of the Macro which generates the adjoint
  //   problem without the Frechet derivative terms.
  // * OPT_FrechetDerivative: name of Macro which writes the Frechet
  //   derivative Galerkin terms of the performance function
  // * OPT_PartialDerivative: name of Macro which writes the partial
  //   derivative post-op terms of the performance function
  // Output:
  // * The adjoint formulation for the given perfromance function with a name,
  //   AdjointFormulationOf_OPT_PhysFormulationName_OPT_FrechetDerivative
  // * The corresponding post-processing and post-operation
  _AdjointName_ = StrCat[OPT_AdjointWithoutFrechet,"_",OPT_FrechetDerivative];

  Printf[StrCat["Add the Frechet term ",OPT_FrechetDerivative,
    " to the adjoint formulation ",OPT_AdjointWithoutFrechet]];

  // Load the template of the adjoint formulation for a given physical system
  // (specified by the string OPT_PhysFormulationName).
  Call S2N[OPT_AdjointWithoutFrechet];

  Formulation {
    { // Adjoint formulation of OPT_FUNC_NAME
      Append;
      Name S2N[_AdjointName_]; Type FemEquation;

      // Add the Fréchet derivative terms of the perfromanec function
      Equation {
        Call S2N[OPT_FrechetDerivative];
      }
    }
  }

  PostProcessing {
    { // Adjoint postpro of OPT_FUNC_NAME
      Append;
      Name S2N[_AdjointName_]; NameOfFormulation S2N[_AdjointName_];
      PostQuantity {
        // Add the Partial derivative of the performance function
        If(Exists(OPT_PartialDerivative))
          Call S2N[OPT_PartialDerivative];
        EndIf
      }
    }
  }
Return

// -------------------------------------------------------------------------

// In the following, we create the post-operations necessary
// for the computation of objective and constraints.
// We also make use of the adjoint method to compute their
// sensitivity wit respect to the density field.

PostProcessing{
  { Name ObjectiveConstraints; NameOfFormulation Elast_u;
    Quantity {
      { Name density ; Value { Term { [ {xe} ]; In DomainDensity; Jacobian Vol; } } }
      { Name Felast_pressure; Value {
        Term{ [ Vector[pressure_x[{xe}],pressure_y[{xe}],0] ];In Sur_Force~{iP}; }}}
      { Name uElast ; Value { Term { [ uElast[ {ux},{uy}] ];
        In Vol_Elast~{iP}; Jacobian Vol; } } }
      { Name Young; Value { Term { [ Young[{xe}] ];
        In Vol_Elast~{iP} ; Jacobian Vol; }}}

      { Name VolumeDomainFull; Value { Integral { [ 1 ];
        In DomainDensity; Jacobian Vol; Integration GaussOnePoint;}}}
      { Name VolumeMaterial; Value { Integral { [ {xe} ];
        In DomainDensity; Jacobian Vol; Integration GaussOnePoint;}}}

      // internal energy
      { Name Compliance; Value {
        Integral { [ InternalEnergy[{xe}, {d ux}, {d uy}] ];
          In DomainDensity ; Jacobian Vol; Integration I1; }}}

      // macroscopic Von-Mises stresses
      { Name vonMises; Value { Term { [ sigVM[{d ux},{d uy}] ];
        In Vol_Elast~{iP} ; Jacobian Vol; }}}

      // qp-relaxed Von-Mises stresses
      { Name vonMises_qp; Value { Term { [ sigVMqp[{d ux},{d uy},{xe}] ];
        In Vol_Elast~{iP} ; Jacobian Vol; }}}

      // P-norm of qp-relaxed Von-Mises stresses
      { Name vonMises_qp_sum; Value { Integral {
        [ (sigVMqp[{d ux},{d uy},{xe}]^PARAM_pnorm)/ElementVol[] ];
        In DomainDensity;Jacobian Vol;Integration I1; }}}

      { Name vonMises_qp_pnorm; Value { Term { Type Global;
        [ $vonMises_qp_sum ^ (1/PARAM_pnorm) ]; In Vol_Elast~{iP}; }}}
    }
  }
}

PostOperation {
  { Name Get_DensityField; NameOfPostProcessing ObjectiveConstraints;
    Operation {
      CreateDir[Opt_ResDir];
      Print[density, OnElementsOf DomainDensity,
        File StrCat[Opt_ResDir,"densityField.pos"]];
    }
  }
}

PostOperation {
  { Name Get_ObjectiveConstraints; NameOfPostProcessing ObjectiveConstraints;
    Operation {
      // Volume
      Print[ VolumeDomainFull[DomainDensity], OnGlobal,
        Format Table, File "", Color "LightYellow",
        SendToServer StrCat[Opt_ResDir_Onelab,"VolumeDomainFull"] ];

      Print[ VolumeMaterial[DomainDensity], OnGlobal,
        Format Table, File "", Color "LightYellow",
        SendToServer StrCat[Opt_ResDir_Onelab,"VolumeDomain"] ];

      // p-norm of Von-Mises
      Print[ vonMises_qp_sum[DomainDensity], OnGlobal,
        Format Table, File "", StoreInVariable $vonMises_qp_sum ];

      Print[ vonMises_qp_pnorm, OnRegion DomainDensity,
        Format Table, StoreInVariable $vonMises_qp_pnorm, File "", Color "LightYellow",
        SendToServer StrCat[Opt_ResDir_Onelab,"pnorm of Von-Mises"] ] ;

      If(Flag_PrintLevel > 2 && !Opt_isInnerIter)
        CreateDir[Opt_ResDir];
        Print[ vonMises_qp, OnElementsOf Vol_Elast~{iP},
          File StrCat[Opt_ResDir,"vm_qp.pos"]];
      EndIf
      If(Flag_PrintLevel > 4 && !Opt_isInnerIter)
        Print[ vonMises, OnElementsOf Vol_Elast~{iP},
          File StrCat[Opt_ResDir,"vm.pos"]];
      EndIf
      If(Flag_PrintLevel > 5 && !Opt_isInnerIter)
        Print[ Young, OnElementsOf Vol_Elast~{iP},
          File StrCat[Opt_ResDir, "young.pos"] ];
        //Print[ Felast_pressure, OnElementsOf Sur_Force~{iP},
        //  File StrCat[Opt_ResDir, "force_pressure.pos"] ];
        Print[ uElast, OnElementsOf Vol_Elast~{iP},
          File StrCat[Opt_ResDir, "u.pos"] ];
      EndIf
    }
  }
}

// -------------------------------------------------------------------------

// Sensitivity of volume

PostProcessing{
  { Name GradOf_Volume; NameOfFormulation Elast_u;
    Quantity {
      { Name SensVolume; Value { Integral{ [ 1. ];
        In DomainDensity; Jacobian Vol; Integration GaussOnePoint; }}}
    }
  }
}

PostOperation{
  { Name Get_GradOf_Volume; NameOfPostProcessing GradOf_Volume;
    Operation {
      // evaluate the sensitivity of the volume fraction at the barycenter
      // of the elements and store it with a reference "SensVolumeFraction"
      Print[ SensVolume, OnElementsOf DomainDensity,
        Format ElementTable, Depth 0, Name "", File "", Hidden 1,
        SendToServer StrCat[Opt_ResDir_Onelab,"GradOfVolume"] ];
      If (Flag_PrintLevel > 5)
        CreateDir[Opt_ResDir];
        Print[ SensVolume, OnElementsOf DomainDensity, Depth 1,
          File StrCat[Opt_ResDir, StrCat["SensVolume.pos"]]] ;
      EndIf
    }
  }
}

// -------------------------------------------------------------------------

// Sensitivity of the p-norm of the Von-Mises

Function {
  pnormVMCoef_qp[] = ( $vonMises_qp_pnorm^(1-PARAM_pnorm) )
    * ( sigVMqp[$1,$2,$3]^(PARAM_pnorm-1) ) / ElementVol[];

  Frechet_VM_qp[] = pnormVMCoef_qp[$1,$2,$3]
    * ($3 ^ (Flag_opt_matlaw_penal-Flag_opt_stress_penal))
    * (Transpose[C_2d[1]] * (VMmat[] * sigM[$1,$2]))
    / sigVM[$1,$2];

  Partial_VM_qp[] = pnormVMCoef_qp[$1,$2,$3]
    * (Flag_opt_matlaw_penal-Flag_opt_stress_penal)
    * sigVMqp[$1,$2,$3] / $3 ;

  Frechet_VM_wrt_dux[] = Vector[CompX[Frechet_VM_qp[$1,$2,$3]#105], CompZ[#105], 0];
  Frechet_VM_wrt_duy[] = Vector[CompZ[Frechet_VM_qp[$1,$2,$3]#105], CompY[#105], 0];
}

Macro OPT_FRECHET_PnormVonMises
  // Galerkin Term of the Fréchet derivative of torque at a given u*
  Galerkin { [ - Frechet_VM_wrt_dux[{d ux},{d uy},{xe}], {d lamx} ];
    In DomainDensity ; Jacobian Vol ; Integration I1 ; }
  Galerkin { [ - Frechet_VM_wrt_duy[{d ux},{d uy},{xe}], {d lamy} ];
    In DomainDensity ; Jacobian Vol ; Integration I1 ; }
Return

// Finally, we create the adjoint formulation of CurrentMismatchNorm
OPT_AdjointWithoutFrechet = "AdjointFormulationOf_Elast_u";
OPT_FrechetDerivative = "OPT_FRECHET_PnormVonMises";
Call CreateAdjointFormulation;
nameAdjoint_vonMises = Str(_AdjointName_);

PostProcessing {
  { Name GradOf_VonMises; NameOfFormulation S2N[nameAdjoint_vonMises];
    PostQuantity {
      { Name Grad_vm;
        Value {
          Integral {[ Partial_VM_qp[{d ux},{d uy},{xe}] ];
            In DomainDensity; Jacobian Vol; Integration GaussOnePoint; }
          Integral {[ -GradPartialElast[{xe},{d ux},{d uy},{d lamx},{d lamy}] ];
            In DomainDensity; Jacobian Vol; Integration GaussOnePoint; }
        }
      }

      { Name lamElast ; Value { Term { [ uElast[ {lamx},{lamy}] ];
        In Vol_Elast~{iP}; Jacobian Vol; } } }
    }
  }
}

PostOperation{
  { Name Get_GradOf_VonMises; NameOfPostProcessing GradOf_VonMises;
    Operation {
      Print[ Grad_vm, OnElementsOf DomainDensity,
        Name "", Depth 0, Format ElementTable, File "", Hidden 1,
        SendToServer StrCat[Opt_ResDir_Onelab,"GradOfVonMises"] ];

      If(Flag_PrintLevel > 5)
        CreateDir[Opt_ResDir];
        Print[ Grad_vm, OnElementsOf DomainDensity,
          File StrCat[Opt_ResDir,"GradOfVonMises.pos"] ];
        Print[ lamElast, OnElementsOf Vol_Elast~{iP},
          File StrCat[Opt_ResDir,"lamElast.pos"] ];
      EndIf
    }
  }
}

// -------------------------------------------------------------------------

// Let us define the optimization problem.

Resolution {
  { Name ComputePerformancesAndGradient;
    System {
      { Name SYS; NameOfFormulation Elast_u;}
      { Name LAM; NameOfFormulation S2N[nameAdjoint_vonMises];}
    }
    Operation {
      // Initialization of the systems
      InitSolution[SYS];
      InitSolution[LAM];

      // Show the density field
      If(Flag_PrintLevel > 0 && !Opt_isInnerIter)
        PostOperation[Get_DensityField];
      EndIf

      // Provide the density field to {xe} of the physical system
      UpdateConstraint[SYS];

      // Solve the physical problem
      Generate[SYS];
      Solve[SYS];

      // Compute the objective function and the constraints
      PostOperation[Get_ObjectiveConstraints];

      // Compute the sensitivity of objective/constraints
      If(Flag_getGradient)
        PostOperation[Get_GradOf_Volume];
        Generate[LAM];
        SolveAgainWithOther[LAM,SYS];
        PostOperation[Get_GradOf_VonMises];
      EndIf
    }
  }
}
