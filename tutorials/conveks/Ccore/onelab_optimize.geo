/*
  This file is a generic block of commands
  for the automatic computation of velocity fields.
  This velocity is used for the computation of sensitivities [1].
  It should be included at the bottom 
  of all shape optimisation .geo files.

   [1]:Erin Kuci, François Henrotte, Pierre Duysinx, and Christophe Geuzaine.
   ``Design sensitivity analysis for shape optimization based on the Lie derivative'',
   Computer Methods in Applied Mechanics and Engineering 317 (2017), pp. 702 –722.
   issn: 0045-7825.

 */

DefineConstant[
  PerturbMesh = 0
  VelocityTag = -1
  Perturbation = 1e-10
  modelpath = CurrentDir
];

If(PerturbMesh == 1)
  Printf("Computing velocity field ...");
  modelName = StrPrefix(StrRelative(General.FileName));  

  SyncModel;
  
  // Merge the original mesh
  Merge StrCat(modelpath, modelName, ".msh");

  // create a view with the original node coordinates as a vector dataset
  Plugin(NewView).NumComp = 3;
  Plugin(NewView).Run;
  Plugin(ModifyComponents).View = 0;
  Plugin(ModifyComponents).Expression0 = "x";
  Plugin(ModifyComponents).Expression1 = "y";
  Plugin(ModifyComponents).Expression2 = "z";
  Plugin(ModifyComponents).Run;
  //Save View[0] StrCat(modelpath, Sprintf("view_%g.msh",0));

  // relocate the vertices of the original mesh on the perturbed geometry
  RelocateMesh Point "*";
  RelocateMesh Line "*";
  RelocateMesh Surface "*";
  Save StrCat(modelpath, modelName, Sprintf("Perturb_%g.msh",VelocityTag));

  // Create a view with the perturbed node coordinates as vector dataset
  Plugin(NewView).NumComp = 3;
  Plugin(NewView).Run;
  Plugin(ModifyComponents).View = 1;
  Plugin(ModifyComponents).Expression0 = "x";
  Plugin(ModifyComponents).Expression1 = "y";
  Plugin(ModifyComponents).Expression2 = "z";
  Plugin(ModifyComponents).Run;
  //Save View[1] StrCat(modelpath, Sprintf("view_%g.msh",1));

  // compute the velocity field by subtracting the two vector datasets
  Plugin(ModifyComponents).View = 0;
  Plugin(ModifyComponents).OtherView = 1;
  Plugin(ModifyComponents).Expression0 = Sprintf("(w0 - v0)/(%g)", Perturbation);
  Plugin(ModifyComponents).Expression1 = Sprintf("(w1 - v1)/(%g)", Perturbation);
  Plugin(ModifyComponents).Expression2 = Sprintf("(w2 - v2)/(%g)", Perturbation);
  Plugin(ModifyComponents).Run;
  View.Name = "velocity";
  //Save View[0] StrCat(modelpath, Sprintf("view_%g.msh",2));
  Delete View[1];
  Save View[0] StrCat(modelpath, Sprintf("velocity_%g.msh",VelocityTag));
  SendToServer View[0] Sprintf("Optimization/Results/velocity_%g",VelocityTag);
EndIf
