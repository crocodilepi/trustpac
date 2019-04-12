# -*- coding: utf-8 -*-

###
### This file is generated automatically by SALOME v8.5.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
theStudy = salome.myStudy

# import salome_notebook
# notebook = salome_notebook.NoteBook(theStudy)
# sys.path.insert( 0, r'/home/vd256574/Formation_TRUST/poisson_anisotrope/salome')

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New(theStudy)

# all elements of the geometry
O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
O_1 = geompy.MakeVertex(0, 0, 0)
OX_1 = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY_1 = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ_1 = geompy.MakeVectorDXDYDZ(0, 0, 1)
Vertex_1 = geompy.MakeVertex(0, 0, 0)
Vertex_2 = geompy.MakeVertex(10, 0, 0)
Vertex_3 = geompy.MakeVertex(10, 0.5, 0)
Vertex_4 = geompy.MakeVertex(6, 0.5, 0)
Vertex_5 = geompy.MakeVertex(4, 1, 0)
Vertex_6 = geompy.MakeVertex(0, 1, 0)
Line_1 = geompy.MakeLineTwoPnt(Vertex_1, Vertex_2)
Line_2 = geompy.MakeLineTwoPnt(Vertex_2, Vertex_3)
Line_3 = geompy.MakeLineTwoPnt(Vertex_3, Vertex_4)
Line_4 = geompy.MakeLineTwoPnt(Vertex_5, Vertex_6)
Line_5 = geompy.MakeLineTwoPnt(Vertex_6, Vertex_1)
Curve_1 = geompy.MakeCurveParametric("10*t", "0.75-0.25*sin(5*(t-0.5)*pi)", "0", 0.4, 0.6, 10, GEOM.Polyline, True)
Face_1 = geompy.MakeFaceWires([Line_1, Line_2, Line_3, Line_4, Line_5, Curve_1], 1)

# groups of the geometry
Bas = geompy.CreateGroup(Face_1, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Bas, [3])
Droit = geompy.CreateGroup(Face_1, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Droit, [6])
Haut = geompy.CreateGroup(Face_1, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Haut, [8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30])
Gauche = geompy.CreateGroup(Face_1, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Gauche, [32])

# add to study all elements
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( O_1, 'O' )
geompy.addToStudy( OX_1, 'OX' )
geompy.addToStudy( OY_1, 'OY' )
geompy.addToStudy( OZ_1, 'OZ' )
geompy.addToStudy( Vertex_1, 'Vertex_1' )
geompy.addToStudy( Vertex_2, 'Vertex_2' )
geompy.addToStudy( Vertex_3, 'Vertex_3' )
geompy.addToStudy( Vertex_4, 'Vertex_4' )
geompy.addToStudy( Vertex_5, 'Vertex_5' )
geompy.addToStudy( Vertex_6, 'Vertex_6' )
geompy.addToStudy( Line_1, 'Line_1' )
geompy.addToStudy( Line_2, 'Line_2' )
geompy.addToStudy( Line_3, 'Line_3' )
geompy.addToStudy( Line_4, 'Line_4' )
geompy.addToStudy( Line_5, 'Line_5' )
geompy.addToStudy( Curve_1, 'Curve_1' )
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudyInFather( Face_1, Bas, 'Bas' )
geompy.addToStudyInFather( Face_1, Droit, 'Droit' )
geompy.addToStudyInFather( Face_1, Haut, 'Haut' )
geompy.addToStudyInFather( Face_1, Gauche, 'Gauche' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

# create a new mesher with MG_CADSurf
smesh = smeshBuilder.New(theStudy)
Mesh_1 = smesh.Mesh(Face_1)
MG_CADSurf = Mesh_1.Triangle(algo=smeshBuilder.MG_CADSurf)
MG_CADSurf_Parameters_1 = MG_CADSurf.Parameters()
MG_CADSurf_Parameters_1.SetPhySize( 0.05 )
MG_CADSurf_Parameters_1.SetMinSize( 0.05 )
MG_CADSurf_Parameters_1.SetMaxSize( 0.05 )
MG_CADSurf_Parameters_1.SetChordalError( 0 )

# compute the meshing
isDone = Mesh_1.Compute()
smesh.SetName(Mesh_1, 'Mesh_1')

# group of the meshing
Bas_1 = Mesh_1.GroupOnGeom(Bas,'Bas',SMESH.EDGE)
Droit_1 = Mesh_1.GroupOnGeom(Droit,'Droit',SMESH.EDGE)
Haut_1 = Mesh_1.GroupOnGeom(Haut,'Haut',SMESH.EDGE)
Gauche_1 = Mesh_1.GroupOnGeom(Gauche,'Gauche',SMESH.EDGE)

## Set names of Mesh objects
smesh.SetName(MG_CADSurf.GetAlgorithm(), 'MG-CADSurf')
smesh.SetName(MG_CADSurf_Parameters_1, 'MG-CADSurf Parameters_1')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Gauche_1, 'Gauche')
smesh.SetName(Haut_1, 'Haut')
smesh.SetName(Droit_1, 'Droit')
smesh.SetName(Bas_1, 'Bas')

# export to MED
try:
  Mesh_1.ExportMED( r'./Mesh_1.med', 0, SMESH.MED_MINOR_0, 1, None ,1)
  pass
except:
  print 'ExportToMEDX() failed. Invalid file name?'

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(True)
