""" Last updated: 2023/05/03

To be run in Abaqus standard/explicit

To run within Linux Terminal, Type: abaqus python SCRIPTNAME.py
        NOTE: To run, you must first initial abaqus/standard through:
                    module load abaqus intel
To run within ABAQUS/CAE, type (into Command Line): execfile('SCRIPTNAME.py')
        NOTE: To run through ABAQUS/CAE Kernal, you must add:
                    from abaqus import *
                    from abaqusConstants import *

Script for creation of 1/8 ellipsoid models within Consolini et al. 2023.
Appropriate modifications must be made for creating models to examine the effect
of initial cut length and curvature. 
1) To investigate initial cut length simply change:
    self.cutSize1 = (self.b1/2) - (self.b1/8.086) --> self.cutSize1 = (self.b1/2) - (self.b1/4.000)
    self.cutSize2 = (self.b1/2) + (self.b1/8.086) --> self.cutSize1 = (self.b1/2) + (self.b1/4.000)
2) To investigate curvature change:
    self.a1 = 9.315 --> self.a1 = 4.000
    self.b1 = 6.506 --> self.b1 = 2.794
    self.c1 = 6.506 --> self.c1 = 2.794
   and,
    self.cutSize1 = (self.b1/2) - (self.b1/7.423) --> self.cutSize1 = (self.b1/2) - (self.b1/8.000)
    self.cutSize2 = (self.b1/2) + (self.b1/7.423) --> self.cutSize1 = (self.b1/2) + (self.b1/8.000)
"""
from abaqus import *
from abaqusConstants import *
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
from mesh import ElemType
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
import numpy as np

#///// Create transverse model geometry and rigid constraint
def ICRMouseGeometryTransCut(modelName, duraMaterPart, rigidPart, a1, b1, c1, a2, b2, c2, cutSize1, cutSize2):
    # A) Creation of ellipsoid (a = y, b = z, c = x) (a > b = c)
    mdb.Model(name=modelName, modelType=STANDARD_EXPLICIT) 
    dS = mdb.models[modelName].ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    g = dS.geometry
    dS.setPrimaryObject(option=STANDALONE)
    # A.1) Define geometry
    dS.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
    dS.FixedConstraint(entity=g[2])
    dS.EllipseByCenterPerimeter(center=(0.0, 0.0), axisPoint1=(0.0, a1), 
        axisPoint2=(b1, 0.0))
    dS.EllipseByCenterPerimeter(center=(0.0, 0.0), axisPoint1=(0.0, a2), 
        axisPoint2=(b2, 0.0))
    dS.Line(point1=(0.0, 0.0), point2=(0.0, a1))
    dS.Line(point1=(0.0, 0.0), point2=(b1, 0.0))
    dS.autoTrimCurve(curve1=g[3], point1=(-b1/1.5, a1/2))
    dS.autoTrimCurve(curve1=g[5], point1=(-b1/1.5, a1/2))
    dS.autoTrimCurve(curve1=g[9], point1=(b1/1.5, -a1/2))
    dS.autoTrimCurve(curve1=g[10], point1=(b1/1.5, -a1/2))
    dS.autoTrimCurve(curve1=g[7], point1=(0.10, a1/2))
    dS.autoTrimCurve(curve1=g[8], point1=(1.0, -0.10))
    # A.2) Create part
    dP = mdb.models[modelName].Part(name=duraMaterPart, 
        dimensionality=THREE_D, type=DEFORMABLE_BODY)
    dP = mdb.models[modelName].parts[duraMaterPart]
    dP.BaseSolidRevolve(sketch=dS, angle=90.0, flipRevolveDirection=OFF)
    dS.unsetPrimaryObject() 

    # B) Create two parts split along the midline
    # B.1) Create Datum Plane
    # B.1.1) First partition along x-y midline
    dP = mdb.models[modelName].parts[duraMaterPart]
    f, e, v = dP.faces, dP.edges, dP.vertices
    t = dP.MakeSketchTransform(sketchPlane=f[3], sketchUpEdge=e[6], 
        sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 0.0))
    dPPS = mdb.models[modelName].ConstrainedSketch(name='__profile__', 
        sheetSize=19.01, gridSpacing=0.47, transform=t)
    dPPS.setPrimaryObject(option=SUPERIMPOSE)
    dP.projectReferencesOntoSketch(sketch=dPPS, filter=COPLANAR_EDGES)
    dPPS.Line(point1=(a1/2, 0.0), point2=(a1/2, -b1))
    dPPF = f.getSequenceFromMask(mask=('[#8 ]', ), )
    dP.PartitionFaceBySketch(sketchUpEdge=e[6], faces=dPPF, sketch=dPPS)
    dPPS.unsetPrimaryObject()
    # B.1.2) Second partition along y midline
    t = dP.MakeSketchTransform(sketchPlane=f[5], sketchUpEdge=e[10], 
        sketchPlaneSide=SIDE1, sketchOrientation=BOTTOM, origin=(0.0, 0.0, 
        0.0))
    dPPS = mdb.models[modelName].ConstrainedSketch(name='__profile__', 
        sheetSize=20.1, gridSpacing=0.5, transform=t)
    g = dPPS.geometry
    dPPS.setPrimaryObject(option=SUPERIMPOSE)
    dP.projectReferencesOntoSketch(sketch=dPPS, filter=COPLANAR_EDGES)
    dPPS.Line(point1=(0.0, -(a1/2)), point2=(b1, -(a1/2)))
    dPPF = f.getSequenceFromMask(mask=('[#20 ]', ), )
    dP.PartitionFaceBySketch(sketchUpEdge=e[10], faces=dPPF, 
        sketchOrientation=BOTTOM, sketch=dPPS)
    dPPS.unsetPrimaryObject()
    # B.1.3) Datum axis
    dP.DatumAxisByTwoPoint(point1=v[0], point2=v[4])
    # B.1.4) Partition through face along y midline
    dPPF = f.getSequenceFromMask(mask=('[#8 ]', ), )
    dP.PartitionFaceByShortestPath(point1=v[0], point2=v[4], faces=dPPF)
    # B.1.5) Datum plane
    d = dP.datums
    dP.DatumPlaneByLinePoint(line=d[4], point=dP.InterestingPoint(edge=e[0], 
        rule=MIDDLE))
    # B.2) Create two parts
    parts = ['1', '2']
    flip = [OFF, ON]
    cutEdge = [3, 1]
    for i, j, k in zip(parts, flip, cutEdge):
    # B.2.1) copy and cut duraMaterPart
        duraP = mdb.models[modelName].Part(name=duraMaterPart+i, 
            objectToCopy=mdb.models[modelName].parts[duraMaterPart])
        t = duraP.MakeSketchTransform(sketchPlane=d[6], sketchUpEdge=e[7], 
            sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0, a1/2, 0))
        duraPS = mdb.models[modelName].ConstrainedSketch(name='__profile__', 
            sheetSize=23.52, gridSpacing=0.58, transform=t)
        duraP.projectReferencesOntoSketch(sketch=duraPS, filter=COPLANAR_EDGES)
    # B.2.2) rectangular cut
        duraPS.rectangle(point1=(0, 0), point2=(-40, -40))
        duraP.CutExtrude(sketchPlane=d[6], sketchUpEdge=e[7], sketchPlaneSide=SIDE1, 
            sketchOrientation=RIGHT, sketch=duraPS, flipExtrudeDirection=j)
        duraPS.unsetPrimaryObject()
    # B.2.3) create cut partitions
        fDuraP, eDuraP = duraP.faces, duraP.edges
        t = duraP.MakeSketchTransform(sketchPlane=fDuraP[0], sketchUpEdge=eDuraP[k], 
            sketchPlaneSide=SIDE1, sketchOrientation=BOTTOM, origin=(0.0, a1/2, 
            0.0))
        duraPSPartition = mdb.models[modelName].ConstrainedSketch(name='__profile__', 
            sheetSize=7.99, gridSpacing=0.19, transform=t)
        g = duraPSPartition.geometry
        duraPSPartition.setPrimaryObject(option=SUPERIMPOSE)
        duraP.projectReferencesOntoSketch(sketch=duraPSPartition, filter=COPLANAR_EDGES)
    # B.2.3) Create rectangle partition along y-z surface for cut
        if i == '1':
            duraPSPartition.rectangle(point1=(0.0, cutSize1), point2=(b1, cutSize2))
            # duraPSPartition.rectangle(point1=(cutSize1, 0.0), point2=(cutSize2, b1))
        elif i == '2':
            duraPSPartition.rectangle(point1=(cutSize1, 0.0), point2=(cutSize2, b1))
            # duraPSPartition.rectangle(point1=(0.0, cutSize1), point2=(b1, cutSize2))
        duraPSPartitionFace = fDuraP.getSequenceFromMask(mask=('[#1 ]', ), )
        duraP.PartitionFaceBySketch(sketchUpEdge=eDuraP[k], faces=duraPSPartitionFace, 
        sketchOrientation=BOTTOM, sketch=duraPSPartition)
        duraPSPartition.unsetPrimaryObject()

    # C) Create rigid part to constrain geometry
    rS = mdb.models[modelName].ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    g = rS.geometry
    rS.setPrimaryObject(option=STANDALONE)
    # C.1) Define geometry
    rS.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
    rS.FixedConstraint(entity=g[2])
    rS.EllipseByCenterPerimeter(center=(0.0, 0.0), axisPoint1=(0.0, a2), 
        axisPoint2=(b2, 0.0))
    rS.EllipseByCenterPerimeter(center=(0.0, 0.0), axisPoint1=(0.0, a2-0.01), 
        axisPoint2=(b2-0.01, 0.0))
    rS.Line(point1=(0.0, 0.0), point2=(0.0, a2))
    rS.Line(point1=(0.0, 0.0), point2=(b2, 0.0))
    rS.autoTrimCurve(curve1=g[3], point1=(-a2/3.5, a2/2))
    rS.autoTrimCurve(curve1=g[5], point1=(-a2/3.5, a2/2))
    rS.autoTrimCurve(curve1=g[9], point1=(a2/3.5, -a2/2))
    rS.autoTrimCurve(curve1=g[10], point1=(a2/3.5, -a2/2))
    rS.autoTrimCurve(curve1=g[7], point1=(0.05, a2/1.25))
    rS.autoTrimCurve(curve1=g[8], point1=(b2/3.75, -0.10))
    # C.2) Create part
    rP = mdb.models[modelName].Part(name=rigidPart, dimensionality=THREE_D, 
        type=DISCRETE_RIGID_SURFACE)
    rP.BaseSolidRevolve(sketch=rS, angle=90.0, flipRevolveDirection=OFF)
    rS.unsetPrimaryObject()
    rC, rE = rP.cells, rP.edges
    rP.ReferencePoint(point=rP.InterestingPoint(edge=rE[2], rule=MIDDLE))
    rP.RemoveCells(cellList = rC[0:1])
    rP.checkGeometry()

    # D) Delete duraMaterPart
    del mdb.models[modelName].parts[duraMaterPart]

#///// Create longitudinal model geometry and rigid constraint
def ICRMouseGeometryLongCut(modelName, duraMaterPart, rigidPart, a1, b1, c1, a2, b2, c2, cutSize1, cutSize2):
    # A) Creation of ellipsoid (a = y, b = z, c = x) (a > b = c)
    mdb.Model(name=modelName, modelType=STANDARD_EXPLICIT) 
    dS = mdb.models[modelName].ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    g = dS.geometry
    dS.setPrimaryObject(option=STANDALONE)
    # A.1) Define geometry
    dS.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
    dS.FixedConstraint(entity=g[2])
    dS.EllipseByCenterPerimeter(center=(0.0, 0.0), axisPoint1=(0.0, a1), 
        axisPoint2=(b1, 0.0))
    dS.EllipseByCenterPerimeter(center=(0.0, 0.0), axisPoint1=(0.0, a2), 
        axisPoint2=(b2, 0.0))
    dS.Line(point1=(0.0, 0.0), point2=(0.0, a1))
    dS.Line(point1=(0.0, 0.0), point2=(b1, 0.0))
    dS.autoTrimCurve(curve1=g[3], point1=(-b1/1.5, a1/2))
    dS.autoTrimCurve(curve1=g[5], point1=(-b1/1.5, a1/2))
    dS.autoTrimCurve(curve1=g[9], point1=(b1/1.5, -a1/2))
    dS.autoTrimCurve(curve1=g[10], point1=(b1/1.5, -a1/2))
    dS.autoTrimCurve(curve1=g[7], point1=(0.10, a1/2))
    dS.autoTrimCurve(curve1=g[8], point1=(1.0, -0.10))
    # A.2) Create part
    dP = mdb.models[modelName].Part(name=duraMaterPart, 
        dimensionality=THREE_D, type=DEFORMABLE_BODY)
    dP = mdb.models[modelName].parts[duraMaterPart]
    dP.BaseSolidRevolve(sketch=dS, angle=90.0, flipRevolveDirection=OFF)
    dS.unsetPrimaryObject()
    
    # B) Create two parts split along the midline
    # B.1) Create Datum Plane
    # B.1.1) First partition along x-z midline
    f, e, v = dP.faces, dP.edges, dP.vertices
    t = dP.MakeSketchTransform(sketchPlane=f[2], sketchUpEdge=e[7], 
        sketchPlaneSide=SIDE1, sketchOrientation=BOTTOM, origin=(0.0, 0.0, 
        0.0))
    dPP1 = mdb.models[modelName].ConstrainedSketch(name='__profile__', 
        sheetSize=9.22, gridSpacing=0.23, transform=t)
    g = dPP1.geometry
    dPP1.setPrimaryObject(option=SUPERIMPOSE)
    dP.projectReferencesOntoSketch(sketch=dPP1, filter=COPLANAR_EDGES)
    dPP1.Line(point1=(c1/2, 0.0), point2=(c2/2, b1))
    dPPF1 = f.getSequenceFromMask(mask=('[#4 ]', ), )
    dP.PartitionFaceBySketch(sketchUpEdge=e[7], faces=dPPF1, 
        sketchOrientation=BOTTOM, sketch=dPP1)
    dPP1.unsetPrimaryObject()
    # B.1.2) create x-y partition
    t = dP.MakeSketchTransform(sketchPlane=f[5], sketchUpEdge=e[2], 
        sketchPlaneSide=SIDE1, sketchOrientation=BOTTOM, origin=(0.0, 0.0, 
        0.0))
    dPP2 = mdb.models[modelName].ConstrainedSketch(name='__profile__', 
        sheetSize=20.1, gridSpacing=0.5, transform=t)
    g = dPP2.geometry
    dPP2.setPrimaryObject(option=SUPERIMPOSE)
    dP.projectReferencesOntoSketch(sketch=dPP2, filter=COPLANAR_EDGES)
    dPP2.Line(point1=(c1/2, 0.0), point2=(c1/2, -a1))
    dPPF2 = f.getSequenceFromMask(mask=('[#20 ]', ), )
    dP.PartitionFaceBySketch(sketchUpEdge=e[2], faces=dPPF2, 
        sketchOrientation=BOTTOM, sketch=dPP2)
    dPP2.unsetPrimaryObject()
    # B.1.3) partition between two closest points
    dP.DatumAxisByTwoPoint(point1=v[0], point2=v[4])
    dP.DatumPointByCoordinate(coords=(c1/2, a1/2, c1/2))
    d = dP.datums
    dP.DatumPlaneByLinePoint(line=d[4], point=d[5])
    # B.2) CREATE TWO PARTS
    parts = ['1', '2']
    flip = [OFF, ON]
    partEdges = [2, 2]
    for i, j, k in zip(parts, flip, partEdges):
    # B.2.1) copy and cut duraMaterPart
        duraP = mdb.models[modelName].Part(name=duraMaterPart+i, 
            objectToCopy=mdb.models[modelName].parts[duraMaterPart])
        t = duraP.MakeSketchTransform(sketchPlane=d[6], sketchUpEdge=e[12], 
            sketchPlaneSide=SIDE1, sketchOrientation=BOTTOM, origin=(c1/2, 0, 0))
        duraPS = mdb.models[modelName].ConstrainedSketch(name='__profile__', 
            sheetSize=23.52, gridSpacing=0.58, transform=t)
        duraPS.setPrimaryObject(option=SUPERIMPOSE)
        duraP.projectReferencesOntoSketch(sketch=duraPS, filter=COPLANAR_EDGES)
    # B.2.2) rectangular cut
        duraPS.rectangle(point1=(40, 40), point2=(-40, -40))
        duraP.CutExtrude(sketchPlane=d[6], sketchUpEdge=e[12], sketchPlaneSide=SIDE1, 
            sketchOrientation=BOTTOM, sketch=duraPS, flipExtrudeDirection=j)
        duraPS.unsetPrimaryObject()
    # B.2.3) create cut partitions
        fDuraP, eDuraP = duraP.faces, duraP.edges
        t = duraP.MakeSketchTransform(sketchPlane=fDuraP[0], sketchUpEdge=eDuraP[k], 
            sketchPlaneSide=SIDE1, sketchOrientation=BOTTOM, origin=(0, 
            0.00, 0.00))
        duraPSPartition = mdb.models[modelName].ConstrainedSketch(name='__profile__', 
            sheetSize=10, gridSpacing=0.25, transform=t)
        duraPSPartition.setPrimaryObject(option=SUPERIMPOSE)
        duraP.projectReferencesOntoSketch(sketch=duraPSPartition, filter=COPLANAR_EDGES)
        if i == '1':
            # duraPSPartition.rectangle(point1=(-cutSize1, 0.0), point2=(-cutSize2, -b1))
            duraPSPartition.rectangle(point1=(0.0, -cutSize1), point2=(b1, -cutSize2))
        elif i == '2':
            # duraPSPartition.rectangle(point1=(cutSize1, 0.0), point2=(cutSize2, -b1))
            duraPSPartition.rectangle(point1=(0.0, -cutSize1), point2=(-b1, -cutSize2))
        duraPSPartitionFace = fDuraP.getSequenceFromMask(mask=('[#1 ]', ), )
        duraP.PartitionFaceBySketch(sketchUpEdge=eDuraP[k], faces=duraPSPartitionFace, 
            sketchOrientation=BOTTOM, sketch=duraPSPartition)
        duraPSPartition.unsetPrimaryObject()
        
    # C) CREATE RIGID PART FOR CONSTRAINED GEOMETRY
    rS = mdb.models[modelName].ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    g = rS.geometry
    rS.setPrimaryObject(option=STANDALONE)
    # C.1) Define geometry
    rS.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
    rS.FixedConstraint(entity=g[2])
    rS.EllipseByCenterPerimeter(center=(0.0, 0.0), axisPoint1=(0.0, a2), 
        axisPoint2=(b2, 0.0))
    rS.EllipseByCenterPerimeter(center=(0.0, 0.0), axisPoint1=(0.0, a2-0.01), 
        axisPoint2=(b2-0.01, 0.0))
    rS.Line(point1=(0.0, 0.0), point2=(0.0, a2))
    rS.Line(point1=(0.0, 0.0), point2=(b2, 0.0))
    rS.autoTrimCurve(curve1=g[3], point1=(-a2/3.5, a2/2))
    rS.autoTrimCurve(curve1=g[5], point1=(-a2/3.5, a2/2))
    rS.autoTrimCurve(curve1=g[9], point1=(a2/3.5, -a2/2))
    rS.autoTrimCurve(curve1=g[10], point1=(a2/3.5, -a2/2))
    rS.autoTrimCurve(curve1=g[7], point1=(0.05, a2/1.25))
    rS.autoTrimCurve(curve1=g[8], point1=(b2/3.75, -0.10))
    # C.1) Create part
    rP = mdb.models[modelName].Part(name=rigidPart, dimensionality=THREE_D, 
        type=DISCRETE_RIGID_SURFACE)
    rP.BaseSolidRevolve(sketch=rS, angle=90.0, flipRevolveDirection=OFF)
    rS.unsetPrimaryObject()
    rC, rE = rP.cells, rP.edges
    rP.ReferencePoint(point=rP.InterestingPoint(edge=rE[2], rule=MIDDLE))
    rP.RemoveCells(cellList = rC[0:1])
    rP.checkGeometry()

    # D) Delete duraMaterPart
    del mdb.models[modelName].parts[duraMaterPart]

#///// Material definition
def material(modelName, duraMaterPart):
    # A) Define material
    mdb.models[modelName].Material(name=duraMaterPart)
    # A.1) Hyperelastic
    mdb.models[modelName].materials[duraMaterPart].Hyperelastic(
        materialType=ISOTROPIC, testData=OFF, type=OGDEN, n=3, 
        volumetricResponse=VOLUMETRIC_DATA, table=((0.0066832903832431, 
        34.9785182961848, 0.0116066769387916, 34.980053927089, 
        0.0125350793232072, 34.9786299090531, 0.0, 0.0, 0.0), ))

#///// Section definition
def section(modelName, duraMaterPart):
    parts = ['1', '2']
    for i in parts:
    # A) Indentify cells
        dP = mdb.models[modelName].parts[duraMaterPart+i]
        dC = dP.cells
    # B) Create sections
        mdb.models[modelName].HomogeneousSolidSection(name=duraMaterPart+i,
            material=duraMaterPart, thickness=None)
        duraMater = dP.Set(cells=dC, name=duraMaterPart+i)
        dP.SectionAssignment(region=duraMater, sectionName=duraMaterPart+i, 
            offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', 
            thicknessAssignment=FROM_SECTION)

#///// Transverse model sets
def setsTransCut(modelName, duraMaterPart):
    # A) Dura mater sets for part 1
    dP = mdb.models[modelName].parts[duraMaterPart+'1']
    dF, dE = dP.faces, dP.edges
    # A.1) Dura mater faces
    xyAxisFace_1 = dF.getSequenceFromMask(mask=('[#10 ]', ), )
    yzAxisFace_1 = dF.getSequenceFromMask(mask=('[#40 ]', ), )
    xzCutFace_1  = dF.getSequenceFromMask(mask=('[#4 ]', ), )
    xzCutFaceAbove_1 = dF.getSequenceFromMask(mask=('[#2 ]', ), )
    xzCutFaceBelow_1 = dF.getSequenceFromMask(mask=('[#1 ]', ), )
    interiorDuralFace_1 = dF.getSequenceFromMask(mask=('[#8 ]', ), )
    # A.2) Dura mater sets
    dP.Set(faces=xyAxisFace_1, name='xyFace_1')
    dP.Set(faces=yzAxisFace_1, name='yzFace_1')
    dP.Set(faces=xzCutFace_1, name='xzCutFace_1')
    dP.Set(faces=xzCutFaceAbove_1, name='xzCutFaceAbove_1')
    dP.Set(faces=xzCutFaceBelow_1, name='xzCutFaceBelow_1')
    # A.3) Dura mater surfaces
    dP.Surface(side1Faces=xzCutFaceAbove_1, name='xzCutSurfaceAbove_1')
    dP.Surface(side1Faces=xzCutFaceBelow_1, name='xzCutSurfaceBelow_1')
    dP.Surface(side1Faces=xzCutFace_1, name='xzCutSurface_1')
    dP.Surface(side1Faces=interiorDuralFace_1, name='interiorDuralSurface_1')

    # B) Dura mater sets for part 2
    dP = mdb.models[modelName].parts[duraMaterPart+'2']
    dF, dE = dP.faces, dP.edges
    # B.1) Dura mater faces
    xyAxisFace_2 = dF.getSequenceFromMask(mask=('[#80 ]', ), )
    yzAxisFace_2 = dF.getSequenceFromMask(mask=('[#20 ]', ), )
    xzAxisFace_2 = dF.getSequenceFromMask(mask=('[#40 ]', ), )
    xzCutFace_2  = dF.getSequenceFromMask(mask=('[#4 ]', ), )
    xzCutFaceAbove_2 = dF.getSequenceFromMask(mask=('[#1 ]', ), )
    xzCutFaceBelow_2 = dF.getSequenceFromMask(mask=('[#2 ]', ), )
    interiorDuralFace_2 = dF.getSequenceFromMask(mask=('[#8 ]', ), )
    # B.2) Dura mater sets
    dP.Set(faces=xyAxisFace_2, name='xyFace_2')
    dP.Set(faces=yzAxisFace_2, name='yzFace_2')
    dP.Set(faces=xzAxisFace_2, name='xzFace_2')
    dP.Set(faces=xzCutFace_2, name='xzCutFace_2')
    dP.Set(faces=xzCutFaceAbove_2, name='xzCutFaceAbove_2')
    dP.Set(faces=xzCutFaceBelow_2, name='xzCutFaceBelow_2')
    # B.3) Dura mater surfaces
    dP.Surface(side1Faces=xzCutFaceAbove_2, name='xzCutSurfaceAbove_2')
    dP.Surface(side1Faces=xzCutFaceBelow_2, name='xzCutSurfaceBelow_2')
    dP.Surface(side1Faces=xzCutFace_2, name='xzCutSurface_2')
    dP.Surface(side1Faces=interiorDuralFace_2, name='interiorDuralSurface_2')

    # C) Rigid sets
    rP = mdb.models[modelName].parts[rigidPart]
    rF = rP.faces 
    # C.1) Rigid faces
    exteriorRigidFace = rF.getSequenceFromMask(mask=('[#2 ]', ), )
    # C.2) Rigid surfaces
    rP.Surface(side1Faces=exteriorRigidFace, name='exteriorRigidSurface')

#///// Longitudinal model sets
def setsLongCut(modelName, duraMaterPart):
    # A) Dura mater sets for part 1
    dP = mdb.models[modelName].parts[duraMaterPart+'1']
    dF, dE = dP.faces, dP.edges
    # A.1) Dura mater faces
    yzAxisFace_1 = dF.getSequenceFromMask(mask=('[#80 ]', ), )
    xyAxisFace_1 = dF.getSequenceFromMask(mask=('[#20 ]', ), )
    xzAxisFace_1 = dF.getSequenceFromMask(mask=('[#8 ]', ), )
    yzCutFace_1  = dF.getSequenceFromMask(mask=('[#2 ]', ), )
    yzCutFaceAbove_1 = dF.getSequenceFromMask(mask=('[#4 ]', ), )
    yzCutFaceBelow_1 = dF.getSequenceFromMask(mask=('[#1 ]', ), )
    interiorDuralFace_1 = dF.getSequenceFromMask(mask=('[#40 ]', ), )
    # A.2) Dura mater sets
    dP.Set(faces=xyAxisFace_1, name='xyFace_1')
    dP.Set(faces=xzAxisFace_1, name='xzFace_1')
    dP.Set(faces=yzAxisFace_1, name='yzFace_1')
    dP.Set(faces=yzCutFace_1, name='yzCutFace_1')
    dP.Set(faces=yzCutFaceAbove_1, name='yzCutFaceAbove_1')
    dP.Set(faces=yzCutFaceBelow_1, name='yzCutFaceBelow_1')
    # A.3) Dura mater surfaces
    dP.Surface(side1Faces=yzCutFaceAbove_1, name='yzCutSurfaceAbove_1')
    dP.Surface(side1Faces=yzCutFaceBelow_1, name='yzCutSurfaceBelow_1')
    dP.Surface(side1Faces=yzCutFace_1, name='yzCutSurface_1')
    dP.Surface(side1Faces=interiorDuralFace_1, name='interiorDuralSurface_1')

    # B) Dura mater sets for part 2
    dP = mdb.models[modelName].parts[duraMaterPart+'2']
    dF, dE = dP.faces, dP.edges
    # B.1) Dura mater faces
    xyAxisFace_2 = dF.getSequenceFromMask(mask=('[#100 ]', ), )
    xzAxisFace_2 = dF.getSequenceFromMask(mask=('[#40 ]', ), )
    yzCutFace_2  = dF.getSequenceFromMask(mask=('[#4 ]', ), )
    yzCutFaceAbove_2 = dF.getSequenceFromMask(mask=('[#2 ]', ), )
    yzCutFaceBelow_2 = dF.getSequenceFromMask(mask=('[#1 ]', ), )
    interiorDuralFace_2 = dF.getSequenceFromMask(mask=('[#10 ]', ), ) 
    # B.2) Dura mater sets
    dP.Set(faces=xyAxisFace_2, name='xyFace_2')
    dP.Set(faces=xzAxisFace_2, name='xzFace_2')
    dP.Set(faces=yzCutFace_2, name='yzCutFace_2')
    dP.Set(faces=yzCutFaceAbove_2, name='yzCutFaceAbove_2')
    dP.Set(faces=yzCutFaceBelow_2, name='yzCutFaceBelow_2')
    # B.3) Dura mater surfaces
    dP.Surface(side1Faces=yzCutFaceAbove_2, name='yzCutSurfaceAbove_2')
    dP.Surface(side1Faces=yzCutFaceBelow_2, name='yzCutSurfaceBelow_2')
    dP.Surface(side1Faces=yzCutFace_2, name='yzCutSurface_2')
    dP.Surface(side1Faces=interiorDuralFace_2, name='interiorDuralSurface_2')

    # C) Rigid sets
    rP = mdb.models[modelName].parts[rigidPart]
    rF = rP.faces 
    # C.1) Rigid faces
    exteriorRigidFace = rF.getSequenceFromMask(mask=('[#2 ]', ), )
    # C.2) Rigid surfaces
    rP.Surface(side1Faces=exteriorRigidFace, name='exteriorRigidSurface')

#///// Assembly
def assembly(modelName, duraMaterPart, duraMaterInstance, rigidPart, rigidInstance):
    # A) Dura mater
    duraP1 = mdb.models[modelName].parts[duraMaterPart+'1']
    duraP2 = mdb.models[modelName].parts[duraMaterPart+'2']
    # B) Rigid part
    rP = mdb.models[modelName].parts[rigidPart]
    # C) Assembly
    a = mdb.models[modelName].rootAssembly
    a.DatumCsysByDefault(CARTESIAN)
    instance = a.Instance(name=duraMaterInstance+'_1', part=duraP1, dependent=OFF)
    instance = a.Instance(name=duraMaterInstance+'_2', part=duraP2, dependent=OFF)
    instance = a.Instance(name=rigidInstance, part=rP, dependent=OFF)

#///// Mesh controls
def mesh(modelName, duraMaterInstance, rigidInstance, duraMeshSize, rigidMeshSize, duraMaterEdges_1, duraMaterEdges_2):
    # A) Set mesh controls
    a = mdb.models[modelName].rootAssembly  
    # A.1) Load the two dura mater part edges
    edP1 = a.instances[duraMaterInstance+'_1'].edges
    edP2 = a.instances[duraMaterInstance+'_2'].edges
    partInstances =(a.instances[duraMaterInstance+'_1'], 
        a.instances[duraMaterInstance+'_2'], )
    # A.2) Locally seed dura mater parts (do re-seed of first part just to run for loop)
    for i, j in zip(duraMaterEdges_1, duraMaterEdges_2):
        part1Edge = edP1.getSequenceFromMask(mask=(i, ), )
        a.seedEdgeBySize(edges=part1Edge, size=duraMeshSize, deviationFactor=duraMeshSize, 
        constraint=FINER)
        part2Edge = edP2.getSequenceFromMask(mask=(j, ), )
        a.seedEdgeBySize(edges=part2Edge, size=duraMeshSize, deviationFactor=duraMeshSize, 
        constraint=FINER)
    # A.3) Element types and mesh controls applied to part instances
    elemType1 = ElemType(elemCode=C3D8H, elemLibrary=STANDARD)
    elemType2 = ElemType(elemCode=C3D6H, elemLibrary=STANDARD)
    elemType3 = ElemType(elemCode=C3D4H, elemLibrary=STANDARD, 
        secondOrderAccuracy=OFF, distortionControl=DEFAULT)
    for i in partInstances:
        c = i.cells
        partCells = c.getSequenceFromMask(mask=('[#1 ]', ), )
        region =(partCells, )
        a.setElementType(regions=region, elemTypes=(elemType1, elemType2, 
        elemType3))
        a.setMeshControls(regions=partCells, elemShape=TET, technique=FREE)
    # A.4) Generate mesh
    a.generateMesh(regions=partInstances)

    # B) Rigid part
    rP =(a.instances[rigidInstance], )
    a.seedPartInstance(regions=rP, size=rigidMeshSize, deviationFactor=0.1, minSizeFactor=0.1)
    rC = a.instances[rigidInstance].cells
    a.setMeshControls(regions=rC, elemShape=TET, technique=FREE)
    a.generateMesh(regions=rP)

#///// Steps
def steps(modelName, step, timePeriod):
    # A) Define steps
    mdb.models[modelName].StaticStep(name=step, previous='Initial', 
                timePeriod=timePeriod, nlgeom=ON, stabilizationMethod=NONE, 
                timeIncrementationMethod=AUTOMATIC, maxNumInc=100, minInc=1e-25)
    mdb.models[modelName].steps[step].control.setValues(allowPropagation=OFF,
            resetDefaultValues=OFF, timeIncrementation=(4.0, 8.0, 200.0, 200.0, 
            200.0, 200.0, 200.0, 200.0, 200.0, 10.0, 200.0))
    mdb.models[modelName].steps[step].control.setValues(
        discontinuous=ON, displacementField=(0.005, 0.01, 0.0, 0.0, 0.02, 
        10, 1e-06, 1e-08, 1.0, 10, 1e-08), rotationField=(0.005, 0.01, 0.0, 0.0, 0.02, 
        10, 1e-06, 1e-08, 1.0, 1e-15))

#///// Boundary conditions
def boundaryConditions(modelName, duraMaterInstance, rigidInstance, step):
    a = mdb.models[modelName].rootAssembly
    # A) Dura mater sides fixed
    region = a.instances[duraMaterInstance+'_1'].sets['xyFace_1']
    mdb.models[modelName].ZsymmBC(name='xyFaceZFixed_1', createStepName=step,
        region=region, localCsys=None)
    region = a.instances[duraMaterInstance+'_2'].sets['xyFace_2']
    mdb.models[modelName].ZsymmBC(name='xyFaceZFixed_2', createStepName=step,
        region=region, localCsys=None)
    region = a.instances[duraMaterInstance+'_2'].sets['xzFace_2']
    mdb.models[modelName].YsymmBC(name='xzFaceYFixed_2', createStepName=step,
        region=region, localCsys=None)
    region = a.instances[duraMaterInstance+'_1'].sets['yzFace_1']
    mdb.models[modelName].XsymmBC(name='yzFaceXFixed_1', createStepName=step,
        region=region, localCsys=None)
    if modelName == 'neonateICRMouseTrans' or modelName == 'adultICRMouseTrans':
        region = a.instances[duraMaterInstance+'_2'].sets['yzFace_2']
        mdb.models[modelName].XsymmBC(name='yzFaceXFixed_2', createStepName=step,
            region=region, localCsys=None)
    elif modelName == 'neonateICRMouseLong' or modelName == 'adultICRMouseLong':
        region = a.instances[duraMaterInstance+'_1'].sets['xzFace_1']
        mdb.models[modelName].YsymmBC(name='xzFaceYFixed_1', createStepName=step,
            region=region, localCsys=None)

    # B) Rigid part: fixed structure
    r1 = a.instances[rigidInstance].referencePoints
    refPoints1 = (r1[2], )
    region = regionToolset.Region(referencePoints=refPoints1)
    mdb.models[modelName].EncastreBC(name='rigidFixedRefPoint', createStepName=step, 
        region=region, localCsys=None)

#///// Interaction with rigid part
def interaction(modelName, duraMaterInstance, rigidInstance, step, duraMaterCuts_1, duraMaterCuts_2):
    # A) Create contact property between dura mater and rigid part
    a = mdb.models[modelName].rootAssembly
    mdb.models[modelName].ContactProperty('contact')
    # A.1) tangential behavior
    mdb.models[modelName].interactionProperties['contact'].TangentialBehavior(
        formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, 
        pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, 
        table=((0.25, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, 
        fraction=0.005, elasticSlipStiffness=None)
    # A.2) normal behavior
    mdb.models[modelName].interactionProperties['contact'].NormalBehavior(
        pressureOverclosure=HARD, allowSeparation=ON, contactStiffness=DEFAULT, 
        contactStiffnessScaleFactor=1.0, clearanceAtZeroContactPressure=0.0, 
        stiffnessBehavior=NONLINEAR, stiffnessRatio=0.01, 
        upperQuadraticFactor=0.03, lowerQuadraticRatio=0.33333, 
        constraintEnforcementMethod=PENALTY)
    # A.3) geometric behavior
    mdb.models[modelName].interactionProperties['contact'].GeometricProperties(
        contactArea=1.0, padThickness=1.0)
    
    # B) Surface contact between dural parts and rigid.
    parts =['_1', '_2']
    for i in parts:
        region2 = a.instances[duraMaterInstance+i].surfaces['interiorDuralSurface'+i]
    # B.1) Rigid part interaction
        region1 = a.instances[rigidInstance].surfaces['exteriorRigidSurface']
        mdb.models[modelName].SurfaceToSurfaceContactStd(
            name='contactWhileDisplacing'+i, createStepName=step, 
            main=region1, secondary=region2, sliding=SMALL, thickness=ON, 
            interactionProperty='contact', adjustMethod=NONE, 
            initialClearance=OMIT, datumAxis=None, clearanceRegion=None)

    # C) Create tie contact property between cut surfaces
    tieNum = ['1', '2', '3']
    for i, j, k in zip(duraMaterCuts_1, duraMaterCuts_2, tieNum):
        region1 = a.instances[duraMaterInstance+'_1'].surfaces[i]
        region2 = a.instances[duraMaterInstance+'_2'].surfaces[j]
        mdb.models[modelName].Tie(name='tie'+k, main=region1, 
            secondary=region2, positionToleranceMethod=COMPUTED, adjust=ON, 
            tieRotations=ON, constraintEnforcement=SURFACE_TO_SURFACE, 
            thickness=ON)

#///// Job
def job(modelName, jobName):
    mdb.Job(name=jobName, model=modelName, description='', 
        type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, 
        memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=8, 
        numDomains=8, numGPUs=8)

#///// Initialization of models
class neonateICRMouseTransCut:
   def __init__(self):
    # Parts
    self.duraMeshSize = 0.050
    self.rigidMeshSize = 0.10
    # Geometric dimensions
    self.thickness = 0.008
    self.a1 = 4.337
    self.b1 = 4.278
    self.c1 = 4.278
    self.a2 = self.a1-self.thickness
    self.b2 = self.b1-self.thickness
    self.c2 = self.c1-self.thickness
    # Cut along z midline
    self.cutSize1 = (self.b1/2) - (self.b1/8.086)
    self.cutSize2 = (self.b1/2) + (self.b1/8.086)
    # Edge lists
    self.duraMaterEdges_1 = ['[#3804 ]','[#f0 ]','[#f ]','[#311 ]','[#5440 ]','[#5440 ]']
    self.duraMaterEdges_2 = ['[#28440 ]', '[#f0 ]', '[#f ]', '[#311 ]', '[#34800 ]','[#13004 ]']
    # Tie lists
    self.duraMaterCuts_1 = ['xzCutSurface_1', 'xzCutSurfaceAbove_1', 'xzCutSurfaceBelow_1']
    self.duraMaterCuts_2 = ['xzCutSurface_2', 'xzCutSurfaceAbove_2', 'xzCutSurfaceBelow_2']
    # Model name
    self.modelName = 'neonateICRMouseTrans'
    # Job properties
    self.jobName = 'neonateICRMouseDuraMaterTranCut'

class neonateICRMouseLongCut:
   def __init__(self):
    # Parts
    self.duraMeshSize = 0.050
    self.rigidMeshSize = 0.10
    # Geometric dimensions
    self.thickness = 0.008
    self.a1 = 4.337
    self.b1 = 4.278
    self.c1 = 4.278
    self.a2 = self.a1-self.thickness
    self.b2 = self.b1-self.thickness
    self.c2 = self.c1-self.thickness
    # Cut along y midline
    self.cutSize1 = (self.a1/2) - (self.a1/4.758)
    self.cutSize2 = (self.a1/2) + (self.a1/4.758)
    # Edge Lists
    self.duraMaterEdges_1 = ['[#6500 ]', '[#29004 ]','[#390 ]', '[#f ]', '[#71 ]','[#34800 ]']
    self.duraMaterEdges_2 = ['[#15040 ]', '[#e800 ]', '[#f0 ]', '[#f ]', '[#f ]', '[#311 ]']
    # Tie lists
    self.duraMaterCuts_1 = ['yzCutSurface_1', 'yzCutSurfaceAbove_1', 'yzCutSurfaceBelow_1']
    self.duraMaterCuts_2 = ['yzCutSurface_2', 'yzCutSurfaceAbove_2', 'yzCutSurfaceBelow_2']
    # Model name
    self.modelName = 'neonateICRMouseLong'
    # Job properties
    self.jobName = 'neonateICRMouseDuraMaterLongCut'

class adultICRMouseTransCut:
   def __init__(self):
    # Parts
    self.duraMeshSize = 0.050
    self.rigidMeshSize = 0.10
    # Geometric dimensions
    self.thickness = 0.0138
    self.a1 = 9.315
    self.b1 = 6.506
    self.c1 = 6.506
    self.a2 = self.a1-self.thickness
    self.b2 = self.b1-self.thickness
    self.c2 = self.c1-self.thickness
    # Cut along z midline
    self.cutSize1 = (self.b1/2) - (self.b1/7.423)
    self.cutSize2 = (self.b1/2) + (self.b1/7.423)
    # Edge Lists
    self.duraMaterEdges_1 = ['[#3804 ]','[#f0 ]','[#f ]','[#311 ]','[#5440 ]','[#5440 ]']
    self.duraMaterEdges_2 = ['[#28440 ]', '[#f0 ]', '[#f ]', '[#311 ]', '[#34800 ]','[#13004 ]']
    # Tie lists
    self.duraMaterCuts_1 = ['xzCutSurface_1', 'xzCutSurfaceAbove_1', 'xzCutSurfaceBelow_1']
    self.duraMaterCuts_2 = ['xzCutSurface_2', 'xzCutSurfaceAbove_2', 'xzCutSurfaceBelow_2']
    # Model name
    self.modelName = 'adultICRMouseTrans'
    # Job properties
    self.jobName = 'adultICRMouseDuraMaterTranCut'

class adultICRMouseLongCut:
   def __init__(self):
    # Parts
    self.duraMeshSize = 0.05
    self.rigidMeshSize = 0.10
    # Geometric dimensions
    self.thickness = 0.0138
    self.a1 = 9.315
    self.b1 = 6.506
    self.c1 = 6.506
    self.a2 = self.a1-self.thickness
    self.b2 = self.b1-self.thickness
    self.c2 = self.c1-self.thickness
    # Cut along y midline
    self.cutSize1 = (self.a1/2) - (self.a1/7.286)
    self.cutSize2 = (self.a1/2) + (self.a1/7.286)
    # Edge Lists
    self.duraMaterEdges_1 = ['[#6500 ]', '[#29004 ]','[#390 ]', '[#f ]', '[#71 ]','[#34800 ]']
    self.duraMaterEdges_2 = ['[#15040 ]', '[#e800 ]', '[#f0 ]', '[#f ]', '[#f ]', '[#311 ]']
    # Tie lists
    self.duraMaterCuts_1 = ['yzCutSurface_1', 'yzCutSurfaceAbove_1', 'yzCutSurfaceBelow_1']
    self.duraMaterCuts_2 = ['yzCutSurface_2', 'yzCutSurfaceAbove_2', 'yzCutSurfaceBelow_2']
    # Model name
    self.modelName = 'adultICRMouseLong'
    # Job properties
    self.jobName = 'adultICRMouseDuraMaterLongCut'

def jobs(m, modelProperties):
    if m.modelName == 'neonateICRMouseTrans' or m.modelName == 'adultICRMouseTrans':
        ICRMouseGeometryTransCut(m.modelName, duraMaterPart, rigidPart, m.a1, m.b1, m.c1, m.a2, m.b2, m.c2, m.cutSize1, m.cutSize2)
        setsTransCut(m.modelName, duraMaterPart)
    elif m.modelName == 'neonateICRMouseLong' or m.modelName == 'adultICRMouseLong':
        ICRMouseGeometryLongCut(m.modelName, duraMaterPart, rigidPart, m.a1, m.b1, m.c1, m.a2, m.b2, m.c2, m.cutSize1, m.cutSize2)
        setsLongCut(m.modelName, duraMaterPart)
    material(m.modelName, duraMaterPart)
    section(m.modelName, duraMaterPart)
    assembly(m.modelName, duraMaterPart, duraMaterInstance, rigidPart, rigidInstance)
    mesh(m.modelName, duraMaterInstance, rigidInstance, m.duraMeshSize, m.rigidMeshSize, m.duraMaterEdges_1, m.duraMaterEdges_2)
    steps(m.modelName, step, timePeriod)
    boundaryConditions(m.modelName, duraMaterInstance, rigidInstance, step)
    interaction(m.modelName, duraMaterInstance, rigidInstance, step, m.duraMaterCuts_1, m.duraMaterCuts_2)
    job(m.modelName, m.jobName)

# Part information
duraMaterPart = 'duraMater'
duraMaterInstance = duraMaterPart+'Instance'
rigidPart = 'rigidPart'
rigidInstance = 'rigidInstance'
# Step information
step = 'transverselyIsotropicStretching'
timePeriod = 1.0
modelProperties = [duraMaterPart, duraMaterInstance, rigidPart, rigidInstance, step, timePeriod]
# Run jobs
jobs(neonateICRMouseTransCut(), modelProperties)
jobs(neonateICRMouseLongCut(), modelProperties)
jobs(adultICRMouseTransCut(), modelProperties)
jobs(adultICRMouseLongCut(), modelProperties)


