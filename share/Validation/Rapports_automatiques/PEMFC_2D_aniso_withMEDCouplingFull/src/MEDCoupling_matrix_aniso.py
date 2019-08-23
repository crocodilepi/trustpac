import MEDLoader as ml
import medcoupling as mc

# INPUT
medfilename = "Mesh_1.med"
crushedGDLgroupname = "Haut"
ymin = 0.
ymax = 1.

# INPUT med file, read mesh of GDL
meshMEDFileRead = ml.MEDFileMesh.New(medfilename) # MEDFileUMesh
mesh2d = meshMEDFileRead.getMeshAtLevel(0)	# MEDCouplingUMesh
mesh2d.setName("mesh2D")
mesh2d.getCoords().setInfoOnComponents(["X [m]","Y [m]"])
#print "MESH2D ",mesh2d

# read mesh of the boundary of crushed GDL
mesh1d = meshMEDFileRead.getGroup(-1, crushedGDLgroupname)	# MEDFileUMesh
mesh1d.zipCoords() # important to remove inutil points
mesh = mesh1d.deepCopy()
mesh.setName("mesh1D")
mesh.changeSpaceDimension(1) # magic function!
mesh.checkConsistencyLight() # check
#print "MESH1D", mesh

# thickness of GDL non crushed, assuming that GDL is crushed in vertical axe, e.g Oy
Cxx1 = 100  # in-plane, not crushed
Cxx2 = 200 # in-plane, crushed
Cyy1 = 1 # through-plane, not crushed
Cyy2 = 2 # through-plane, crushed
ep1 = 1. # max thickness = ymax - ymin , not crushed
ep2 = 0.5 # min thickness, crushed

# creer un champ sur la maille surfacique
time = 4.22 # random ms
epArr = mesh1d.getCoords()[:,1]-ymin # epaisseur calculated from the curve bord y-coordinate
#print "MIN MAX VALEUR D'EPAISSEUR EN METRE ", epArr.getMinMaxPerComponent()

# creer un champ sur le maillage du bord
f = ml.MEDCouplingFieldDouble(ml.ON_NODES, ml.ONE_TIME)	# impossible pour un champ no_time ???
f.setTimeUnit("ms") # Time unit is ms.
f.setTime(time,1,-1)	 # Time attached is 4.22 ms, iteration id is 2 and order id (or sub iteration id) is -1
f.setArray(epArr)
f.setMesh(mesh)
f.setName("champ_epaisseur_1D")
#print "CHAMP DEFINED ON THE PLATE MESH FOR INTERPOLATION", f
ml.WriteField("ep1d.med",f, True)

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
ml.WriteField("ep2d.med",f2, True)

# ortho field interpolated
n = mesh2d.getNumberOfCells()
gradTField = ml.MEDCouplingFieldDouble(ml.ON_CELLS, ml.ONE_TIME)	# impossible pour un champ no_time ???
gradTData = mc.DataArrayDouble(n, 2) # nb_cells x 2 comp
gradTData[:,0] = [0.]
gradTData[:,1] = [1.]
gradTField.setTimeUnit("ms") # Time unit is ms.
gradTField.setTime(time,1,-1) # Time attached is 4.22 ms, iteration id is 2 and order id (or sub iteration id) is -1
gradTField.setArray(gradTData)
gradTField.setMesh(mesh2d)
gradTField.setName("champ_normal_2D")
# create the normal vector on the boundary 'haut'
normaField = mesh1d.buildOrthogonalField()
normaField.setMesh(mesh)
normaData = normaField.getArray() # ATTENTION: ortho field not correct!!! vector [-0 -1] for elem #0 -> #79
for i in range(len(normaData)):
	if(normaData[i,1] < 0):
		normaData[i,0] = -normaData[i, 0]
		normaData[i,1] = -normaData[i, 1]

#print "normaField", normaField
#print "normaData", normaData 

bary = gradTField.getMesh().computeCellCenterOfMass() # [x, y]
norm = normaField.getValueOnMulti(bary[:,0])
epai = f.getValueOnMulti(bary[:,0])
uy = mc.DataArrayDouble(n, 2)
uy[:,0] = 0.
uy[:,1] = 1.
dist = norm - uy
eps = 1e-8
for i in range(n):
	# get epaisseur
	epi = epai[i]
	# get y
	yi = bary[i,1]
	# alpha
	alpha = (yi - ymin) / (epi - ymin)
	#if(dist[i,0] > eps):
	if(dist.magnitude()[i] > eps):
		gradTData[i] = alpha*norm[i] + (1. - alpha)*uy[i] # ATTENTION: gradTData not identity

gradTData /= gradTData.magnitude()

#print "gradTData" , gradTData
ml.WriteField("ortho.med",gradTField, True)

# number of Cells
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
matM[:, 1] = gradTData[:, 0]
matM[:, 3] = gradTData[:, 1]
matM[:, 0] = gradTData[:,1]
matM[:, 2] = -gradTData[:,0]

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
fconduc.setTime(time,1,-1) # Time attached is 4.22 ms, iteration id is 2 and order id (or sub iteration id) is -1
fconduc.setArray(conduc)
fconduc.setMesh(mesh2d)
fconduc.setName("CONDUCTIVITY_ELEM_dom")
#print "CHAMP CONDUCTIVITY DEFINED IN THE MESH2D ", fconduc
ml.WriteField("conduc2d.med",fconduc, True)
