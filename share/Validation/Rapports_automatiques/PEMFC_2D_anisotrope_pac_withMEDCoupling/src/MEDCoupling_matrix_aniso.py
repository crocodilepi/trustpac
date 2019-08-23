import MEDLoader as ml
import medcoupling as mc

'''
This script aims to calculate the matrix of anisotropic conductivity of GDL crushed by the bipolar plate. 
 matrix D = [D00 D01; D10 D11] with D = M.D_ortho.M^-1
  M is rotation matrix, M = [cos -sin; sin cos] with an angle theta = (gradT, uy), theta positive in the anti-clockwise
   gradT is a temporary vector assuming that it is the same with the through-plane axe of GDL, gradT is calculated in a previous step and stored in the .med file
   uy is the unit vector of Oy
  D_ortho is the matrix of anisotropie in the cas that the direction in-phane and through-plane of GDL are the same of Ox and Oy of the cartesien Oxy, 
   D_ortho = diag(Dxx, Dyy) is a diagonal matrix
    Dxx is the conductivity in the horizontal direction (Ox)
    Dyy is the conductivity in the  vertical direction (Oy)
    Dxx et Dyy depend on the thickness of GDL and the initial conductivity (normal GDL non crushed)

to check: M^-1 = transpose(M), det(M) = 1

'''

# INPUT
medfilename = "gradT_0000.med"
crushedGDLgroupname = "Haut"
ymin = 0.		# limit of geometry (got by TRUST)
ymax = 1.

# INPUT med file, read mesh of GDL
meshMEDFileRead = ml.MEDFileMesh.New(medfilename) # MEDFileUMesh
mesh2d = meshMEDFileRead.getMeshAtLevel(0)	# MEDCouplingUMesh
mesh2d.setName("mesh2D")
mesh2d.getCoords().setInfoOnComponents(["X [m]","Y [m]"])
#print "MESH2D ",mesh2d

# read mesh of the boundary of crushed GDL
mesh1d = meshMEDFileRead.getGroup(-1, crushedGDLgroupname)	# MEDFileUMesh
mesh1d.zipCoords() # important remove inutil points
mesh = mesh1d.deepCopy()
mesh.setName("mesh1D")
mesh.changeSpaceDimension(1)
mesh.checkConsistencyLight()
#print "MESH1D", mesh

# conducvitiy in-plane and through-plan
Cxx1 = 100  # in-plane, not crushed
Cxx2 = 200 # through-plane, crushed
Cyy1 = 1 # in-plane, not crushed
Cyy2 = 2 # through-plane, crushed
ep1 = 1. # max thickness = ymax - ymin
ep2 = 0.5 # min thickness

# creer un champ sur la maille surfacique, just for TEST
time = 4.22		# ms
valArray = mesh1d.getCoords()[:,1]-ymin	# epaisseur
#print "MIN MAX VALEUR D'EPAISSEUR EN METRE ", valArray.getMinMaxPerComponent()

# creer un champ sur le maillage du bord
f = ml.MEDCouplingFieldDouble(ml.ON_NODES, ml.ONE_TIME)	# impossible pour un champ no_time ???
f.setTimeUnit("ms") # Time unit is ms.
f.setTime(time,1,-1)	 # Time attached is 4.22 ms, iteration id is 2 and order id (or sub iteration id) is -1
f.setArray(valArray)
f.setMesh(mesh)
f.setName("champ_epaisseur_1D")
#print "CHAMP DEFINED ON THE PLATE MESH FOR INTERPOLATION", f
#ml.WriteField("ep1d_mod.med",f, True)

# JUST FOR TEST
bary = mesh.computeCellCenterOfMass()
valArray = f.getValueOnMulti(bary)
#print "GETTING VALUES (INTERPOLATED) IN BARYCENTRE OF CELLS OF MESH ", valArray

# creer un champ sur le maillage volumique
coo2Dx = mesh2d.computeCellCenterOfMass()[:, 0]
val2dArr = f.getValueOnMulti(coo2Dx)
f2 = ml.MEDCouplingFieldDouble(ml.ON_CELLS, ml.ONE_TIME)	# impossible pour un champ no_time ???
f2.setTimeUnit("ms") # Time unit is ms.
f2.setTime(time,1,-1)	 # Time attached is 4.22 ms, iteration id is 2 and order id (or sub iteration id) is -1
f2.setArray(val2dArr)
f2.setMesh(mesh2d)
f2.setName("champ_epaisseur_2D")
#print "CHAMP DEFINED IN THE MESH2D ", f2
#ml.WriteField("ep2d_mod.med",f2, True)

# get champ gradT at the last time step in the .med file
timeStepsIds=mc.GetCellFieldIterations(medfilename,"dom","GRADT_ELEM_dom")
lastTimeStepIds = timeStepsIds[len(timeStepsIds)-1]
lastTimeStep = mc.GetTimeAttachedOnFieldIteration(medfilename,"GRADT_ELEM_dom",lastTimeStepIds[0] , lastTimeStepIds[1])
#print "LAST TIME STEP ITERATOR AND ORDER ", lastTimeStepIds
#print "LAST TIME STEP IN SECONDS ",lastTimeStep
gradTField = mc.ReadField(medfilename, "GRADT_ELEM_dom", lastTimeStepIds[0] , lastTimeStepIds[1])
#print "GRADT FIELD ", gradTField

# vector unit of gradT
gradT = gradTField.getArray()		# 2 components
gradTmagn = gradT.magnitude()		# 1 component
gradT[:,0] /= gradTmagn
gradT[:,1] /= gradTmagn
#print "vector unit of gradT", gradT
#print "min max per component vector unit of gradT ", gradT.getMinMaxPerComponent()
#print "magnitude of vector unit of gradT ", gradT.magnitude()

# number of Cells of mesh3D
n = gradTField.getMesh().getNumberOfCells()
#print "number of Cells", n

# coefficient de conductivite dans le plan
conduc_xx = mc.DataArrayDouble(n)
conduc_xx[:] = Cxx1

# coefficient de conductivite hors le plan
conduc_yy = mc.DataArrayDouble(n)
conduc_yy[:] = Cyy1

for i in range(n):
    alpha = (val2dArr[i] - ep2)/(ep1 - ep2)
    conduc_xx[i] = alpha * Cxx1 + (1-alpha)*Cxx2
    conduc_yy[i] = alpha * Cyy1 + (1-alpha)*Cyy2
    
#print "min max conductivity in-plan xx ", conduc_xx.getMinMaxPerComponent()
#print "min max conductivity through-plan yy ", conduc_yy.getMinMaxPerComponent()

# matrix inverse of rotation M^-1
matM = mc.DataArrayDouble(n,4)	# 4 comp or 2x2 comp
matM[:, 1] = gradT[:, 0]
matM[:, 3] = gradT[:, 1]
matM[:, 0] = gradT[:,1]
matM[:, 2] = -gradT[:,0]

# matrix of rotation = transpose of matrix inverse of rotation
# matM = mc.DataArrayDouble(n,9)	# 9 comp or 3x3 comp
# matM[:,0] = invM[:,0]
# matM[:,1] = invM[:,3]
# matM[:,2] = invM[:,6]
# matM[:,3] = invM[:,1]
# matM[:,4] = invM[:,4]
# matM[:,5] = invM[:,7]
# matM[:,6] = invM[:,2]
# matM[:,7] = invM[:,5]
# matM[:,8] = invM[:,8]

# calcul du coefficient de conductivite
conduc = mc.DataArrayDouble(n,4)		# 4 components
conduc[:,:] = 0.
#print "CONDUCTIVITY INITIAL", conduc
conduc[:,0] += matM[:,0]*matM[:,0]*conduc_xx
conduc[:,0] += matM[:,1]*matM[:,1]*conduc_yy

conduc[:,1] += matM[:,0]*matM[:,2]*conduc_xx
conduc[:,1] += matM[:,1]*matM[:,2]*conduc_yy

conduc[:,2] = conduc[:,1]

conduc[:,3] += matM[:,2]*matM[:,2]*conduc_xx
conduc[:,3] += matM[:,3]*matM[:,3]*conduc_yy

#print "MIN MAX CONDUCTIVITY PER COMPONENT", conduc.getMinMaxPerComponent()

# export champ to a .med file
fconduc = ml.MEDCouplingFieldDouble(ml.ON_CELLS, ml.ONE_TIME)	# impossible pour un champ no_time ???
fconduc.setTimeUnit("s") # Time unit is second
fconduc.setTime(lastTimeStep,lastTimeStepIds[0],lastTimeStepIds[1])	 # Time attached is 4.22 ms, iteration id is 2 and order id (or sub iteration id) is -1
fconduc.setArray(conduc)
fconduc.setMesh(mesh2d)
fconduc.setName("CONDUCTIVITY_ELEM_dom")
#print "CHAMP CONDUCTIVITY DEFINED IN THE MESH2D ", fconduc
ml.WriteField("conduc2d_mod.med",fconduc, True)
