# Open this file with Gmsh (interactively with File->Open->shape.py, or on the
# command line with 'gmsh shape.py')
#

import numpy as np
import onelab
import conveks

c = onelab.client(__file__)

# get gmsh and getdp locations from Gmsh options
mygmsh = c.getString('General.ExecutableFileName')
mygetdp = ''
for s in range(9):
   n = c.getString('Solver.Name' + str(s))
   if(n == 'GetDP'):
      mygetdp = c.getString('Solver.Executable' + str(s))
      break
if(not len(mygetdp)):
   c.sendError('This appears to be the first time you are trying to run GetDP')
   c.sendError('Please run a GetDP model interactively once with Gmsh to ' +
               'initialize the solver location')
   exit(0)
c.sendInfo('Will use gmsh={0} and getdp={1}'.format(mygmsh, mygetdp))

# append correct path to model file names
modelName = 'shape'
file_geo = c.getPath(modelName+'.geo')
file_msh = c.getPath(modelName+'.msh')
file_pro = c.getPath(modelName+'.pro')

# build command strings with appropriate parameters and options
mygetdp += ' -p 0 ' + file_pro + ' -msh ' + file_msh + ' '
mygmsh += ' -p 0 -v 2 -parametric ' + file_geo + ' '

# load geometry in GUI to have something to look at
c.openProject(file_geo)

# dry getdp run (without -solve or -pos option) to get model parameters in the GUI
c.runSubClient('myGetDP', mygetdp)

# read optimization parameters as Onelab parameter (editable in the GUI)
Optimizer = c.getNumber('Optimization/00Optimizer')
maxIter = c.getNumber('Optimization/01Max iterations')
maxChange = c.getNumber('Optimization/02Max change')
maxKKTNorm = c.getNumber('Optimization/03Max KKT norm')

# end of check phase. We're done if we do not start the optimization
if c.action == 'check' :exit(0)

# Let us define now the design variables as a list 
# for which the entries are defiend as follows:
# ['onelab/path', value, lower-bound, upper-bound] 
allParameters = []
xe = 6.0e-3
for k in xrange(11):
    allParameters.append(['Optimization/parameters/spline radius '+str(k+1),xe, xe-4.0e-3, xe+4.0e-3, 'Rs_'+str(k)])

L2 = 18.0e-3;L3 = 15.0e-3
for k in xrange(10):
    theta = float(k) * 50.0/10.0 * np.pi/180.0
    xe = L2*L3/((L3*np.cos(theta))**2.0+(L2*np.sin(theta))**2.0)**0.5
    allParameters.append(['Optimization/parameters/outer spline radius '+str(k+1), xe, xe-4.0e-3, xe+1.5e-3, 'Rso_'+str(k)])

allParameters.append(['Optimization/parameters/L4', 0.00994987, 0.00994987-4.0e-3, xe+4.0e-3, 'L4'])

# we create a dictionnary based on the list of design variables
x = {}
for k in xrange(len(allParameters)):
    x[k] = allParameters[k]

# some utility functions:
# (1) to read the values from a table  
def readSimpleTable(path):
    with open(path) as FileObj:
        return np.array([float(lines.split()[-1]) for lines in FileObj])

# (2) to set the values of a node table
def setNodeTable(var, nodes, val):
    data = np.ravel(np.column_stack((nodes, val)))
    data = np.insert(data, 0, float(len(nodes)))
    c.setNumber(var, values=data, visible=0)

# number of design variables and number of constraints
numVariables = len(x)

# initial value of the variables in the design space,
# the lower bound for design variables,
# as well as the upper bound for design variables.
initialPoint = np.zeros(numVariables)
lowerBound = np.zeros(numVariables)
upperBound = np.zeros(numVariables)
for label, var in x.iteritems():
    initialPoint[label] = var[1]
    lowerBound[label] = var[2]
    upperBound[label] = var[3]

# We define the function which returns the objective and 
# constraints values (as a list) for a given point 'xFromOpti'. 
def get_objective(xFromOpti):
    # send the current point to GetDP model
    for label, var in x.iteritems():
        x[label][1] = xFromOpti[label]
        c.setNumber(var[0],value=x[label][1])

    # reload the geometry in GUI to see the geometrical changes
    c.openProject(file_geo)

    # mesh the geometry
    c.runSubClient('myGmsh', mygmsh + ' -2')

    # run Fem and evaluate the objective function
    c.runSubClient('myGetDP', mygetdp + '-solve GetPerformances')
    objective = np.sum(readSimpleTable(c.getPath('res_opt/w.txt')))
    constraints = np.array([np.sum(xFromOpti)/np.sum(upperBound)-1]);  
    return objective,constraints

# We also define the function which returns the gradient of the 
# objective and constraints (as a C-style array) for a given point. 
def get_grad(perturb=1e-6):
    # generate the velocity field of each design variable
    for id, var in x.iteritems():
        c.runNonBlockingSubClient('myGmsh', mygmsh \
          + ' -setnumber PerturbMesh 1'\
          + ' -setnumber VelocityTag '+str(id)\
          + ' -setnumber '+var[-1]+' '+str(var[1]+perturb)\
          + ' -run')
    c.waitOnSubClients()
    # send the x,y,z components of each velocity field
    for label, var in x.iteritems():
        d = c.getNumbers('Optimization/Results/velocity_'+str(label))
        if label==0:nodes = np.array(d[1::4])
        setNodeTable('Optimization/Results/velocity_'+str(label)+'_1',nodes,np.array(d[2::4]))
        setNodeTable('Optimization/Results/velocity_'+str(label)+'_2',nodes,np.array(d[3::4]))
        setNodeTable('Optimization/Results/velocity_'+str(label)+'_3',nodes,np.array(d[4::4]))
        c.setNumber(var[0],value=x[label][1])

    # compute the sensitivity with respect to design variables at `x'
    for dv, var in x.iteritems():
        c.runNonBlockingSubClient('myGetDP', mygetdp \
          + ' -setnumber VelocityTag '+str(dv)\
          + ' -solve GetGradient_wrt_dv')
    c.waitOnSubClients()
    grad_objective = np.asarray([np.sum(readSimpleTable(c.getPath('res_opt/Grad_w_wrt_dv_'+str(dv)+'.txt')))\
        for dv in xrange(numVariables)])
    grad_constraints = np.ones(numVariables)/np.sum(upperBound);  

    return grad_objective,grad_constraints

# Initialize conveks
conveks.initialize()

# Initialize the MMA optimizer
conveks.mma.initialize(initialPoint, lowerBound, upperBound)

# Set some options for MMA
conveks.mma.option.setNumber('General.Verbosity', 0*4)
if Optimizer != 1:
    conveks.mma.option.setNumber('General.Verbosity', 0) # def 4
    conveks.mma.option.setNumber('SubProblem.isContinuation', 1)
    conveks.mma.option.setNumber('SubProblem.cInitial'  , 1e2)
    conveks.mma.option.setNumber('SubProblem.cFinal'    , 1e5)
    conveks.mma.option.setNumber('SubProblem.cMultiply' , 10)
    conveks.mma.option.setNumber('SubProblem.cIteration', 50)

# Get iteration count (here it will be 1 - could be different in case of restart)
it = conveks.mma.getOuterIteration()
innerit = 0; inneritTot = 0; change = 1.0; kkt_norm = 1.0

# we kkep track of the objective values in all_objectives 
all_objectives = []

# get (copy of) current point
xFromMMA = conveks.mma.getCurrentPoint()

# compute objective function and constraints
# as well as their sensitivity with respect to design variables at `x'
objective,constraints = get_objective(xFromMMA) 
grad_objective,grad_constraints = get_grad()

# send data to onelab server
c.setNumber('Optimization/Results/objective', value=objective) 
c.addNumberChoice('Optimization/Results/objective', value=objective) 
c.setNumber('Optimization/Results/max(|Constraints|)', value=np.max(np.abs(constraints))) 
c.addNumberChoice('Optimization/Results/max(|Constraints|)', value=np.max(np.abs(constraints))) 

# show current iteration
print 'iter.  inner iter.   obj.   max(constr.)   L2-norm(kkt)     change   point'
print '%3i  %3i   %.4e   %.4e   %.4e   %.4e'%(it-1, innerit, objective, np.max(np.abs(constraints)), kkt_norm, change),
print xFromMMA[0:5]

while it <= maxIter and change > maxChange \
    and c.getString('shape/Action') != 'stop':

    # Solve MMA
    if Optimizer==0:conveks.mma.setMoveLimits(lowerBound, upperBound, 0.1)
    conveks.mma.updateCurrentPoint(objective,grad_objective,constraints,grad_constraints)

    # Let us check now if the MMA approximation is conservative 
    # at the new point; and adapt the latter if it is not the case.
    if Optimizer == 1:
        # Evaluate the objective at the new point and check if 
        # the approximation is conservative at the current point
        obj_new,c_new = get_objective(conveks.mma.getCurrentPoint())
        conserv = conveks.mma.isConservativeSubProblem(obj_new, c_new) 

        innerit = 0
        print '\t it.  new-point (xn)       f(xc)       f(xn)      max(c(xc))     max(c(xn))'
        print '\t%3i'%(innerit),
        print conveks.mma.getCurrentPoint()[0:5],
        print '  %.4e  %.4e  %.4e  %.4e'%(objective,obj_new,np.max(constraints), np.max(c_new))
    
        while innerit <= 15 and conserv==0:
            innerit += 1

            # adapt the subproblem to make it conservative and generate a new point
            conveks.mma.updateCurrentPoint(objective,grad_objective,
                objectiveAtUpdatedPoint=obj_new,constraintsAtUpdatedPoint=c_new)

            # evaluate the objective and the constraints at the new point
            # and check if the subproblem is conservative at this point
            obj_new,c_new = get_objective(conveks.mma.getCurrentPoint())
            conserv = conveks.mma.isConservativeSubProblem(obj_new, c_new) 

            # show inner iteration
            print '\t%3i'%(innerit),
            print conveks.mma.getCurrentPoint()[0:5],
            print '  %.4e  %.4e  %.4e  %.4e'%(objective,obj_new,np.max(constraints), np.max(c_new))

        inneritTot += innerit

    # get (copy of) current point
    xFromMMA = conveks.mma.getCurrentPoint()

    # compute objective function and constraints
    # as well as their sensitivity with respect to design variables at `x'
    objective,constraints = get_objective(xFromMMA) 
    grad_objective,grad_constraints = get_grad()
    all_objectives.append(objective)

    # evaluate the L2-norm of KKT 
    kkt_norm = conveks.mma.getKKTNorm(constraints,grad_objective,grad_constraints)

    # get the change between two successive points 
    change = conveks.mma.getDesignChange()

    # show the current state
    print '%3i  %3i   %.4e   %.4e   %.4e   %.4e'%(it, innerit, objective, np.max(np.abs(constraints)), kkt_norm, change),
    print xFromMMA[0:5]

    # get the outer iteration count here
    it = conveks.mma.getOuterIteration()

    # send data to the onelab server
    c.setNumber('Optimization/Results/01Current iteration', value=it, readOnly=1)
    c.setNumber('Optimization/Results/02Current change', value=change, readOnly=1)
    c.setNumber('Optimization/Results/03Current KKT norm', value=kkt_norm, readOnly=1)
    c.setNumber('Optimization/Results/objective', value=objective) 
    c.addNumberChoice('Optimization/Results/objective', value=objective) 
    c.setNumber('Optimization/Results/max(|Constraints|)', value=np.max(np.abs(constraints))) 
    c.addNumberChoice('Optimization/Results/max(|Constraints|)', value=np.max(np.abs(constraints))) 

# show the problem at the converged state 
print '*'*80
print 'Optimization has converged after %i iterations.'%(it-1)
print '  > point           :',xFromMMA[0:5]
print '  > objective       : %.4e'%(objective)
print '  > max(constraints): %.4e'%(np.max(np.abs(constraints)))
print '  > L2-norm of kkt  : %.4e'%(kkt_norm)
print '  > inner iterations: %i'%(inneritTot)
print '  > change          : %.4e'%(change)
print '*'*80

# This should be called at the end of MMA
conveks.mma.finalize()

# This should be called at the end
conveks.finalize()