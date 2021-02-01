/* -------------------------------------------------------------------
   Tutorial 4 : floating potential of a microstrip electrode

   Features:
   - Global quantities and their special shape functions
   - Computation of the energy dual, i.e. of the armature charge of the electrode
   - More on ONELAB parameters (flags, model options, check boxes, menus, ...)

   To compute the solution interactively from the Gmsh GUI:
       File > Open > floating.pro
       Run (button at the bottom of the left panel)

   ------------------------------------------------------------------- */

/* GetDP is pretty good at the management of global (non-local) basis
   functions. Finite element expansions typically associate basis functions to
   individual nodes or edges in the mesh. But consider the situation where a
   scalar field is set to be uniform over a region of the problem (say a
   floating potential electrode in an electrostatic problem). By factorizing the
   identical nodal value "v_electrode", a global (non-local) basis function
   "BF_electrode" is obtained as factor which is the sum of the shape functions
   of all the nodes in the electrode region. This basis function "BF_electrode"
   - is a continuous function, scalar in this case;
   - is equal to 1 at the nodes of the electrode, and to 0 at all other nodes;
   - decreases from 1 to 0 over the one-element-thick layer of elements sharing
     at least one node with the electrode region.

   One such global basis function can be associated with each electrode in the
   system, so that the finite element expansion of the electric scalar potential
   reads:

   v = Sum_k sn_k vn_k + Sum_electrode v_electrode BF_electrode

   with the the sum_k running over all nodes except those of the electrode
   regions.

   This tutorial shows how GetDP takes advantage of global quantities and
   the associated global basis functions
   - to reduce the number of unknowns;
   - to efficiently compute the electrode charges "Q_electrode", which are
     precisely the energy duals of the global "v_electrode" quantities;
   - to deal with floating potentials, which are the computed electrode
     potentials when the electrode charge is imposed;
   - to provide output quantities (charges, armature voltages, capacitances,
     ...) that can be immediately used in a external circuit. */

Group {
  /* Geometrical regions: */

  Air = Region[101];
  Diel1 = Region[111];

  Ground = Region[120];
  Microstrip = Region[121];
  SurfInf = Region[130];

  /* Abstract regions:

     Vol_Ele            : volume where -div(epsilon grad v) = 0 is solved
     Sur_Neu_Ele        : surface with imposed non homogeneous Neumann boundary conditions
                          (on n.d = -n . (epsilon grad v))
     Sur_Electrodes_Ele : electrode regions */

  Vol_Ele = Region[ {Air, Diel1} ];
  Sur_Neu_Ele = Region[ {} ];
  Sur_Electrodes_Ele = Region [ {Ground, Microstrip} ];
}

/* A number of ONELAB parameters are defined to provide model parameters or model
   options interactively. */

MicrostripTypeBC = DefineNumber[0, Name "1Microstrip excitation/Type",
                                Choices{0="Fixed voltage", 1="Fixed charge"}] ;
MicrostripValueBC = DefineNumber[1e-3, Name "1Microstrip excitation/Value"] ;
EpsilonRelDiel = DefineNumber[9.8, Name "2Dielectric/Relative permittivity"] ;
DisplayGlobalBF = DefineNumber[0, Name "3Options/Display global basis functions",
                               Choices {0,1} ] ;
OverwriteOutput = DefineNumber[1, Name "3Options/Overwrite output.txt file",
                               Choices {0,1} ] ;

Function {
  eps0 = 8.854187818e-12;  // permittivity of empty space
  epsilon[Air] = eps0;
  epsilon[Diel1] = EpsilonRelDiel * eps0;
}

Constraint {
  /* The Dirichlet boundary condition on the local electric potential is only
     used for the homogeneous condition on the top and right boundaries. The
     microstrip and the ground are now treated as electrodes, whose voltage is
     imposed with the "SetGlobalPotential" constraint below. */
  { Name Dirichlet_Ele; Type Assign;
    Case {
      { Region SurfInf; Value 0; }
    }
  }

  { Name SetGlobalPotential; Type Assign;
    Case {
      /* Impose the potential regionwise on the different parts of
	 "Sur_Electrodes_Ele". The voltage is imposed to the Microstrip
	 electrode when the "Fixed voltage" option is enabled, i.e. if
	 MicrostripTypeBC == 0. */
      { Region Ground; Value 0; }
      If(MicrostripTypeBC == 0)
        { Region Microstrip; Value MicrostripValueBC; }
      EndIf
    }
  }
  { Name SetArmatureCharge; Type Assign;
    Case {
      /* Impose the charge on the Microstrip electrode when the "Fixed charged"
      option is enabled, i.e. if MicrostripTypeBC == 1. */
      If(MicrostripTypeBC == 1)
        { Region Microstrip; Value MicrostripValueBC; }
      EndIf
    }
  }
}

Group{
  /* The domain of definition comprises all regions on which the field "v" is
     defined.*/
  Dom_Hgrad_v_Ele =  Region[ {Vol_Ele, Sur_Neu_Ele, Sur_Electrodes_Ele} ];
}

FunctionSpace {
  /* The magic in the treatment of global quantitities by GetDP lies in the fact
     that nearly all the work is done at the level of the FunctionSpace
     definition. The finite element expansion of "v" is

     v = Sum_k sn_k vn_k + Sum_electrode v_electrode BF_electrode

     with the sum_k running over all nodes except those of the electrode
     regions. This is exactly what one finds in the FunctionSpace definition
     below with "sf" standing for "BF_electrode" and "vf" for "v_electrode".

     The global quantities are attributed an explicit and meaningful name in the
     "GlobalQuantity" section; these names are used in the corresponding
     "GlobalTerm" in the Formulation. Such global terms are the equivalent of a
     "Integral" term, but where no integration needs to be performed. The
     "AssociatedWith" statement refers to the fact that the global potential of
     an electrode is the (electrostatic) energy dual of the electric charge
     carried by that electrode. Indeed, let us consider the electrostatic weak
     formulation derived in Tutorial 1: find v in Hgrad_v_Ele such that

     (epsilon grad v, grad v')_Vol_Ele - (n . (epsilon grad v), v')_Bnd_Vol_Ele = 0

     holds for all test functions v'. When the test-function v' is BF_electrode,
     the boundary term reduces to

     (n . (epsilon grad v), BF_electrode)_Sur_Electrodes_Ele.

     Since BF_electrode == 1 on the electrode, the boundary term is actually
     simply equal to the integral of (n . epsilon grad v) on the electrode,
     i.e. the flux of the displacement field, which is by definition the
     charge Q_electrode carried by the electrodes.

     By checking the "Display global basis functions" checkbox and running the
     model, you can take a look at how the two "BF_electrode" basis functions in
     this model look like.  Constraints can then be set on either component of
     the FunctionSpace.  Besides the usual Dirichlet boundary condition on the
     local field (left here for the sake of completeness but not
     used in this model), one may fix either the
     GlobalPotential or the ArmatureCharge of each indidual electrode (never both,
     of course). When the ArmatureCharge is fixed, the computed GlobalPotential
     for that electrode is the so-called floating potential. */

  { Name Hgrad_v_Ele; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef vn; Function BF_Node;
        Support Dom_Hgrad_v_Ele; Entity NodesOf[ All, Not Sur_Electrodes_Ele ]; }
      { Name sf; NameOfCoef vf; Function BF_GroupOfNodes;
        Support Dom_Hgrad_v_Ele; Entity GroupsOfNodesOf[ Sur_Electrodes_Ele ]; }
    }
    GlobalQuantity {
      { Name GlobalPotential; Type AliasOf       ; NameOfCoef vf; }
      { Name ArmatureCharge ; Type AssociatedWith; NameOfCoef vf; }
    }
    Constraint {
      { NameOfCoef vn; EntityType NodesOf;
        NameOfConstraint Dirichlet_Ele; }
      { NameOfCoef GlobalPotential; EntityType GroupsOfNodesOf;
	NameOfConstraint SetGlobalPotential; }
      { NameOfCoef ArmatureCharge; EntityType GroupsOfNodesOf;
	NameOfConstraint SetArmatureCharge; }
    }
    // Subspace definition only needed to display BF_electrode in PostProcessing
    SubSpace {
      { Name vf; NameOfBasisFunction sf; }
    }
  }
}

Jacobian {
  { Name Vol ;
    Case {
      { Region All ; Jacobian Vol ; }
    }
  }
}

Integration {
  { Name Int ;
    Case { {Type Gauss ;
            Case { { GeoElement Triangle    ; NumberOfPoints  4 ; }
                   { GeoElement Quadrangle  ; NumberOfPoints  4 ; } }
      }
    }
  }
}

Formulation {
  /* The formulation contains only minor changes compared to formulation from
     the first tutorial. The global quantities are declared as "Global" in the
     "Quantity" section, and a "GlobalTerm" is added that triggers the assembly
     of the additional equation per electrode (the "pre-integrated" boundary
     term) in the system to compute the charge Q_electrode, which
     satisfies (just consider the equation corresponding to the test function
     BF_electrode):

     Q_electrode = (-epsilon grad v, grad BF_electrode)_Vol_Ele */
  { Name Electrostatics_v; Type FemEquation;
    Quantity {
      { Name v; Type Local; NameOfSpace Hgrad_v_Ele; }
      { Name U; Type Global; NameOfSpace Hgrad_v_Ele [GlobalPotential]; }
      { Name Q; Type Global; NameOfSpace Hgrad_v_Ele [ArmatureCharge]; }
      // next line only needed to display the BF_electrode in PostProcessing
      { Name vf; Type Local; NameOfSpace Hgrad_v_Ele [vf]; }
    }
    Equation {
      Integral { [ epsilon[] * Dof{d v} , {d v} ];
        In Vol_Ele; Jacobian Vol; Integration Int; }
      GlobalTerm { [ -Dof{Q} , {U} ]; In Sur_Electrodes_Ele; }
    }
  }
}

Resolution {
  { Name EleSta_v;
    System {
      { Name Sys_Ele; NameOfFormulation Electrostatics_v; }
    }
    Operation {
      Generate[Sys_Ele]; Solve[Sys_Ele]; SaveSolution[Sys_Ele];
    }
  }
}

PostProcessing {
 { Name EleSta_v; NameOfFormulation Electrostatics_v;
    Quantity {
      { Name v; Value {
          Term { [ {v} ]; In Vol_Ele; Jacobian Vol; }
        }
      }
      { Name e; Value {
          Term { [ -{d v} ]; In Vol_Ele; Jacobian Vol; }
        }
      }
      { Name d; Value {
          Term { [ -epsilon[] * {d v} ]; In Vol_Ele; Jacobian Vol; }
        }
      }
      { Name Q; Value {
          Term { [ {Q} ]; In Sur_Electrodes_Ele; }
        }
      }
      { Name U; Value {
          Term { [ {U} ]; In Sur_Electrodes_Ele; }
        }
      }
      { Name C; Value {
          Term { [ {Q}/{U} ]; In Sur_Electrodes_Ele; }
        }
      }
      { Name energy; Value {
          Integral { Type Global; // global integral over Vol_Ele
            [ epsilon[] / 2. * SquNorm[{d v}] ];
            In Vol_Ele; Jacobian Vol; Integration Int;
          }
	}
      }
      // next lines only needed to display global BF in PostProcessing
      { Name BFGround; Value {
          Term { [ BF{vf} ]; In Dom_Hgrad_v_Ele; SubRegion Ground; Jacobian Vol; }
        }
      }
      { Name BFMicrostrip; Value {
          Term { [ BF{vf} ]; In Dom_Hgrad_v_Ele; SubRegion Microstrip; Jacobian Vol; }
        }
      }
    }
  }
}

/* Several output results are generated, which are both displayed in the
   graphical user interface, and stored in disk files. In particular, global
   quantities related results are stored in the "output.txt" file. A user option
   allows to chose to not overwrite the "output.txt" file when running a new
   simulation. */

PostOperation {
  { Name Map; NameOfPostProcessing EleSta_v;
     Operation {
       If( DisplayGlobalBF )
         Print[ BFGround, OnElementsOf Dom_Hgrad_v_Ele, File "BFGround.pos" ];
         Echo[ StrCat["l=PostProcessing.NbViews-1;",
		      "View[l].IntervalsType = 1;",
		      "View[l].NbIso = 40;",
		      "View[l].ShowElement = 1;"],
	       File "BFGround.opt", LastTimeStepOnly] ;

         Print[ BFMicrostrip, OnElementsOf Dom_Hgrad_v_Ele,
		File "BFMicrostrip.pos" ];
	 Echo[ StrCat["l=PostProcessing.NbViews-1;",
		      "View[l].IntervalsType = 1;",
		      "View[l].NbIso = 40;",
		      "View[l].ShowElement = 1;"],
	       File "BFMicrostrip.opt", LastTimeStepOnly] ;
       EndIf

       Print[ v, OnElementsOf Dom_Hgrad_v_Ele, File "v.pos" ];
       Echo[ StrCat["l=PostProcessing.NbViews-1;",
                    "View[l].IntervalsType = 3;",
                    "View[l].NbIso = 40;"],
             File "v.opt", LastTimeStepOnly] ;

       If(OverwriteOutput)
         DeleteFile[ "output.txt" ];
       EndIf

       Echo[ "Microstrip charge [C]:", Format Table, File > "output.txt"] ;
       Print[ Q, OnRegion Microstrip, File > "output.txt", Color "AliceBlue",
	      Format Table, SendToServer "Output/Microstrip/Charge [C]" ];
       Echo[ "Microstrip potential [V]:", Format Table, File > "output.txt"] ;
       Print[ U, OnRegion Microstrip, File > "output.txt", Color "AliceBlue",
	      Format Table, SendToServer "Output/Microstrip/Potential [V]" ];
       Echo[ "Ground charge [C]:", Format Table, File > "output.txt"] ;
       Print[ Q, OnRegion Ground, File > "output.txt", Color "AliceBlue",
	      Format Table, SendToServer "Output/Ground/Charge [C]" ];
       Echo[ "Microstrip capacitance [F]:", Format Table, File > "output.txt"] ;
       Print[ C, OnRegion Microstrip, File > "output.txt", Color "AliceBlue",
	      Format Table, SendToServer "Output/Global/Capacitance [F]" ];
       Echo[ "Electrostatic energy [J]:", Format Table, File > "output.txt"] ;
       Print[ energy, OnGlobal, File > "output.txt",
	      Color "AliceBlue",
	      Format Table, SendToServer "Output/Global/Energy [J]" ];
     }
  }
}
