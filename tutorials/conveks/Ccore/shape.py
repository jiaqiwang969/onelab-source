# Open this file with Gmsh
# - interactively with File->Open->shp.py,
# - in the command line with 'gmsh shape.py'


import numpy as np
import onelab
import conveks
from scipy.optimize import minimize

modelName = 'ccore'

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
file_geo = c.getPath(modelName+'.geo')
file_msh = c.getPath(modelName+'.msh')
file_pro = c.getPath(modelName+'.pro')

# build command strings with appropriate parameters and options
mygmsh += ' -p 0 -v 2 -parametric ' + file_geo + ' '
mygetdp += ' -p 0 ' + file_pro + ' -msh ' + file_msh + ' '

# load geometry in GUI to have something to look at
c.openProject(file_geo)

c.sendCommand('Solver.AutoMesh = -1;')

# dry getdp run (without -solve or -pos option) to get model parameters in the GUI
c.runSubClient('myGetDP', mygetdp)

# end of check phase. We're done if we do not start the optimization
if c.action == 'check' :
   exit(0)

## def readSimpleTable(path):
##     with open(path) as FileObj:
##         return np.array([float(lines.split()[-1]) for lines in FileObj])

def setNodeTable(var, nodes, val):
    data = np.ravel(np.column_stack((nodes, val)))
    data = np.insert(data, 0, float(len(nodes)))
    c.setNumber(var, values=data, visible=0)

def getVelocityField(xvar, perturb=1e-7):
    for id, var in xvar.iteritems():
        c.runNonBlockingSubClient('myGmsh', mygmsh \
          + ' -setnumber PerturbMesh 1'\
          + ' -setnumber Perturbation '+str(perturb)\
          + ' -setnumber VelocityTag '+str(id)\
          + ' -setnumber '+var[-1]+' '+str(var[1]+perturb)\
          + ' -run')
    c.waitOnSubClients()
    # send the x,y,z components of each velocity field
    for label, var in x.iteritems():
        d = c.getNumbers('Optimization/Results/velocity_'+str(label))
        if label==0:
            nodes = np.array(d[1::4])
        setNodeTable('Optimization/Results/velocity_'+str(label)+'_1',nodes,np.array(d[2::4]))
        setNodeTable('Optimization/Results/velocity_'+str(label)+'_2',nodes,np.array(d[3::4]))
        setNodeTable('Optimization/Results/velocity_'+str(label)+'_3',nodes,np.array(d[4::4]))
        c.setNumber(var[0],value=x[label][1])

def fobj(xFromOpti):
    # send the current point to GetDP model
    for label, var in x.iteritems():
        x[label][1] = xFromOpti[label]
        c.setNumber(var[0],value=x[label][1])

    # reload the geometry in GUI to see the geometrical changes
    c.openProject(file_geo)

    c.runSubClient('myGmsh', mygmsh + ' -2')
    c.runSubClient('myGetDP', mygetdp + '-solve GetPerformances' \
                   + ' - setnumber OuterIteration ' + str(0)) # it !!
    return c.getNumber('Optimization/Results/w')


def grad(xFromOpti):
    # send the current point to GetDP model
    for label, var in x.iteritems():
        x[label][1] = xFromOpti[label]
        c.setNumber(var[0],value=x[label][1])
        
    # as well as their sensitivity with respect to design variables at `x'
    getVelocityField(x)
    for dv, var in x.iteritems():
        c.runNonBlockingSubClient('myGetDP', mygetdp \
          + ' -setnumber VelocityTag '+str(dv)\
          + ' -solve GetGradient_wrt_dv')
    c.waitOnSubClients()
    return np.asarray([c.getNumber('Optimization/Results/dwdtau_'+str(dv)) for dv, var in x.iteritems()])
    
Nfeval = 1
def callbackF(Xi):
    global Nfeval
    print '{0:4d} {1: 3.6f} {2: 3.6f} {3: 3.6f}'.format(Nfeval, Xi[0], fobj(Xi), grad(Xi)[0])
    Nfeval += 1

def test(n):
    tab = np.zeros((n+1,4))
    for step in range(n+1):
        xstep = x[0][2] + (x[0][3]-x[0][2])/float(n)*step
        f,g,d = fobj([xstep]), grad([xstep])[0], 0
        tab[step] = [xstep, f, g, d]
        if step >= 2 :
            tab[step-1][3] = (tab[step][1]-tab[step-2][1])/(tab[step][0]-tab[step-2][0])
        if c.getString(c.name+'/Action') == 'stop':
            exit(1)
    print
    for step in range(n+1):
        ind1, ind2 = max(step-1,0), min(step+1,n)
        tab[step][3] = (tab[ind2][1]-tab[ind1][1])/(tab[ind2][0]-tab[ind1][0])
        print "%4d" % step,
        for item in tab[step] : print "%3.6f" % item,
        print
    
    return tab

# List of design parameters with format
# [ "Onelab name", value, Min, Max, "variable name" ]
designParameters={
    0 : ['Optimization/Parameters/D', 0.02, 0.001, 0.050, 'D'],
    1 : ['Optimization/Parameters/E', 0.01 , 0.01, 0.15 , 'E']
}
Optimizers={
    0 : 'scipy',
    1 : 'mma'
}

x = {}
index = int(c.getNumber('Optimization/Sweep on parameter'))-1
Optimizer = Optimizers[int(c.getNumber('Optimization/Optimizer'))]

# define optimization parameters as Onelab parameter (editable in the GUI)
maxIter = c.defineNumber('Optimization/00Max iterations', value=100, visible=(Optimizer=='mma'))
maxChange = c.defineNumber('Optimization/01Max change', value=1.0e-3,  visible=(Optimizer=='mma'))

if index >= 0: # parameter sweep and optimisation over one selected design parameter
    x[0] = designParameters[index][:] # make a copy
    test(int(c.getNumber('Optimization/Steps in range')))
else: # optimisation over full parameter space
    x = designParameters.copy()

# initial value of the variables in the design space,
# the lower bound for design variables,
# as well as the upper bound for design variables.
numVariables = 1 if index >= 0 else len(designParameters)
initialPoint = np.zeros(numVariables)
lowerBound = np.zeros(numVariables)
upperBound = np.zeros(numVariables)
if index >= 0: # single parameter
    x[0] = designParameters[index]
    initialPoint[0] = designParameters[index][1]
    lowerBound[0] = designParameters[index][2]
    upperBound[0] = designParameters[index][3]
else: # all the parameters
    for label, var in designParameters.iteritems():
        x[label] = var
        initialPoint[label] = var[1]
        lowerBound[label] = var[2]
        upperBound[label] = var[3]

if Optimizer=='scipy':
    bounds = [(lowerBound[i], upperBound[i]) for i in range(numVariables)]
    res = minimize(fobj, initialPoint.tolist(), jac=grad, bounds=bounds, \
        callback=callbackF, method = 'L-BFGS-B', tol=None, \
        options={'ftol': 1e-5, 'gtol': 1e-3, 'disp':True} )
    print 'Scipy result:', res

elif Optimizer=='mma':
    # Initialize the MMA optimizer
    conveks.mma.initialize(initialPoint, lowerBound, upperBound)

    # Set some options for MMA
    conveks.mma.option.setNumber('General.Verbosity', 0) 

    # Get iteration count (here it will be 1 - could be different in case of restart)
    it = conveks.mma.getOuterIteration()
    change = 1.0
    
    while it <= maxIter and c.getString(c.name+'/Action') != 'stop' and change > maxChange:
        c.setNumber('Optimization/01Current iteration', value=it, readOnly=1,
            attributes={'Highlight':'LightYellow'})

        xFromMMA = np.array(conveks.mma.getCurrentPoint())
        objective = fobj(xFromMMA)
        grad_objective = grad(xFromMMA)
        # since mma is tailored for constrained optimization,
        # we add a dummy constrained which is always fulfilled 
        constraints = np.array([np.sum(xFromMMA)/10.0-1.0])
        grad_constraints = np.ones(numVariables)/10.0
    
        print it, change, xFromMMA, objective, grad_objective
        c.sendInfo('Optimization: it. {}, obj. {}, constr. {}, ch. {}'.format(it,objective,constraints[0],change))

        # call MMA and update the current point
        #conveks.mma.setMoveLimits(lowerBound, upperBound, 1e-0) # parametre interessant
        conveks.mma.updateCurrentPoint(grad_objective,constraints,grad_constraints)
        change = conveks.mma.getDesignChange()
        it = conveks.mma.getOuterIteration()

    print 'Mma result: ', conveks.mma.getCurrentPoint()

    # This should be called at the end
    conveks.mma.finalize()
