""" Last updated: 2023/05/03

To be run in Abaqus standard/explicit

To run within Linux Terminal, Type: abaqus python SCRIPTNAME.py
        NOTE: To run, you must first initial abaqus through:
                    module load abaqus/2022
To run within ABAQUS/CAE, type (into Command Line): execfile('SCRIPTNAME.py')
        NOTE: To run through ABAQUS/CAE Kernal, you must add:
                    from abaqus import *
                    from abaqusConstants import *

Script for extracting incision opening ratio from numerical simulations included
in Consolini et al. 2023.
"""

from odbAccess import *
import numpy as np
import math
import csv

#///// Initialize odb
class odb:
   def __init__(self, odbName, part1, part2, stepSelection, elementType, elementSet1, 
    elementSet2, nodeBounds, frames):
        self.odbName = odbName
        self.odb = openOdb(odbName +'.odb', readOnly=False)
        self.part1 = self.odb.rootAssembly.instances[part1]
        self.part2 = self.odb.rootAssembly.instances[part2]
        self.stepSelection = self.odb.steps[stepSelection]
        self.elementType = elementType
        self.elementSet1 = self.part1.elementSets[elementSet1]
        self.elementSet2 = self.part2.elementSets[elementSet2]
        self.setName1 = elementSet1
        self.setName2 = elementSet2
        self.nodeBounds = nodeBounds
        self.frames = frames

#///// Extract v and a information, and define elements for stretch measurements
def cutLengthAndWidths_and_stretches(self):
    # A) Get nodes and elements in cut area
    cutNodeSet1 = self.part1.nodeSets[self.setName1]
    cutNodeSet2 = self.part2.nodeSets[self.setName2]
    cutElementSet1 = self.part1.elementSets[self.setName1]
    cutElementSet2 = self.part2.elementSets[self.setName2]
    index = [0, 1, 2]
    
    # B) Cut length and widths
    cutFrameList = []
    measurementLocationsList = []
    measurementNodes1List = []
    measurementNodes2List = []
    lengthList = []
    widthsList = []
    wOverlHalf = []
    lambdaPList = []
    for frame in range(len(self.frames)):
        print(frame)
    # B.1) Coordinates and node labels (Trans: x = aT , y = vT, Long: x = vL, y = aL)
        cutNodeCoordsDeformed1 = self.stepSelection.frames[frame].fieldOutputs['COORD'].getSubset(region=cutNodeSet1).values
        cutNodeCoordsDeformed2 = self.stepSelection.frames[frame].fieldOutputs['COORD'].getSubset(region=cutNodeSet2).values
    # B.2) Extracting a (length of cut), then 0.25*a, 0.5*a, 0.75*a
        lengths = []
        for node, coords in zip(cutNodeSet1.nodes, cutNodeCoordsDeformed1):
    # B.3) only extract nodes at tips of cut and compute half cut length (a/2)
            if node.label >= (self.nodeBounds[0]-0.5) and node.label <= (self.nodeBounds[1]+0.5):
                l = [coords.data[0], coords.data[1], coords.data[2]]
                lengths.append(l)
        lD = np.sqrt((lengths[1][0]-lengths[0][0])**2+(lengths[1][1]-lengths[0][1])**2+(lengths[1][2]-lengths[0][2])**2)
        l = 0.50*lD
    # B.4) Extracting v (width of cut)
        nodeHalf1 = []
        nodeHalf2 = []
        wCoords1_1 = []
        wCoords2_1 = []
        wCoords3_1 = []
        wCoords1_2 = []
        wCoords2_2 = []
        wCoords3_2 = []
        for node1, coords1, node2, coords2 in zip(cutNodeSet1.nodes, cutNodeCoordsDeformed1, cutNodeSet2.nodes, cutNodeCoordsDeformed2):
            l1 = np.sqrt((coords1.data[0]-lengths[0][0])**2+(coords1.data[1]-lengths[0][1])**2+(coords1.data[2]-lengths[0][2])**2)
            l2 = np.sqrt((coords2.data[0]-lengths[0][0])**2+(coords2.data[1]-lengths[0][1])**2+(coords2.data[2]-lengths[0][2])**2)
            if l1 >= (l-0.025) and l1 <= (l+0.025):
                nodeHalf1.append(node1.label)
                wCoords1_1.append(coords1.data[0])
                wCoords2_1.append(coords1.data[1])
                wCoords3_1.append(coords1.data[2])
            if l2 >= (l-0.025) and l2 <= (l+0.025):
                nodeHalf2.append(node2.label)
                wCoords1_2.append(coords2.data[0])
                wCoords2_2.append(coords2.data[1])
                wCoords3_2.append(coords2.data[2])
        wCoord1Ave_1 = np.average(wCoords1_1)
        wCoord2Ave_1 = np.average(wCoords2_1)
        wCoord3Ave_1 = np.average(wCoords3_1)
        wCoord1Ave_2 = np.average(wCoords1_2)
        wCoord2Ave_2 = np.average(wCoords2_2)
        wCoord3Ave_2 = np.average(wCoords3_2)
        wHalf = np.sqrt((wCoord1Ave_2-wCoord1Ave_1)**2+(wCoord2Ave_2-wCoord2Ave_1)**2+(wCoord3Ave_2-wCoord3Ave_1)**2)
        w = wHalf
        cutFrameList.append(frame)
        measurementLocationsList.append('w(l(half))')
        measurementNodes1List.append(nodeHalf1)
        measurementNodes2List.append(nodeHalf2)
        lengthList.append(l)
        widthsList.append(w)
        wOverlHalf.append(wHalf/lD)
    # B.7) Get lambda p for each frame
    for frame in range(len(self.frames)):
        print(frame)
        lambdaP1 = []
        for element1 in cutElementSet1.elements:
            for j in range(len(measurementNodes1List[frame])):
                if (element1.connectivity[0]==measurementNodes1List[frame][j]) or (element1.connectivity[1]==measurementNodes1List[frame][j]) or (element1.connectivity[2]==measurementNodes1List[frame][j]) or (element1.connectivity[3]==measurementNodes1List[frame][j]):
                    lambdaP = self.stepSelection.frames[frame].fieldOutputs['SDV1'].getSubset(region=element1, elementType=self.elementType).values[0].data
                    lambdaP1.append(lambdaP)
        aveLambdaP_1 = np.average(lambdaP1)
        lambdaPList.append(aveLambdaP_1)
    # B.8) Export results
    global cutResults
    cutResults = [cutFrameList, measurementLocationsList, measurementNodes1List, measurementNodes2List, lambdaPList, wOverlHalf, lengthList, widthsList]

#///// Results files
# Universal model properties
duraMaterInstances = ['DURAMATERINSTANCE_1', 'DURAMATERINSTANCE_2']
step = 'transverselyIsotropicStretching'
elementType = 'C3D4H'
modelProperties = [duraMaterInstances, step, elementType]
# Neonate and adult murine cranial dura mater models
class neonateTransverse:
    def __init__(self):
        self.odbName = 'neonateICRMouseDuraMaterTranCut'
        self.sets = ['XZCUTFACE_1','XZCUTFACE_2']
        self.nodeBounds = [9, 10, 7, 8]
        self.odb = openOdb(self.odbName +'.odb', readOnly=False)
        self.frames = self.odb.steps[step].frames
class neonateLongitudinal:
    def __init__(self):
        self.odbName = 'neonateICRMouseDuraMaterLongCut'
        self.sets = ['YZCUTFACE_1','YZCUTFACE_2']
        self.nodeBounds = [11, 12, 9, 10]
        self.odb = openOdb(self.odbName+'.odb', readOnly=False)
        self.frames = self.odb.steps[step].frames
class adultTransverse:
    def __init__(self):
        self.odbName = 'adultICRMouseDuraMaterTranCut'
        self.sets = ['XZCUTFACE_1','XZCUTFACE_2']
        self.nodeBounds = [9, 10, 7, 8]
        self.odb = openOdb(self.odbName+'.odb', readOnly=False)
        self.frames = self.odb.steps[step].frames
class adultLongitudinal:
    def __init__(self):
        self.odbName = 'adultICRMouseDuraMaterLongCut'
        self.sets = ['YZCUTFACE_1','YZCUTFACE_2']
        self.nodeBounds = [11, 12, 7, 9]
        self.odb = openOdb(self.odbName+'.odb', readOnly=False)
        self.frames = self.odb.steps[step].frames
# Effect of initial cut length vs. incision opening ratio
class adultTran32ndAxis:
    def __init__(self):
        self.odbName = 'adultICRMouseDuraMaterTranCut32ndAxis'
        self.sets = ['XZCUTFACE_1','XZCUTFACE_2']
        self.nodeBounds = [9, 10, 7, 8]
        self.odb = openOdb(self.odbName+'.odb', readOnly=False)
        self.frames = self.odb.steps[step].frames
class adultTran16thAxis:
    def __init__(self):
        self.odbName = 'adultICRMouseDuraMaterTranCut16thAxis'
        self.sets = ['XZCUTFACE_1','XZCUTFACE_2']
        self.nodeBounds = [9, 10, 7, 8]
        self.odb = openOdb(self.odbName+'.odb', readOnly=False)
        self.frames = self.odb.steps[step].frames
class adultTran8thAxis:
    def __init__(self):
        self.odbName = 'adultICRMouseDuraMaterTranCut8thAxis'
        self.sets = ['XZCUTFACE_1','XZCUTFACE_2']
        self.nodeBounds = [9, 10, 7, 8]
        self.odb = openOdb(self.odbName+'.odb', readOnly=False)
        self.frames = self.odb.steps[step].frames
class adultTran4thAxis:
    def __init__(self):
        self.odbName = 'adultICRMouseDuraMaterTranCut4thAxis'
        self.sets = ['XZCUTFACE_1','XZCUTFACE_2']
        self.nodeBounds = [9, 10, 7, 8]
        self.odb = openOdb(self.odbName+'.odb', readOnly=False)
        self.frames = self.odb.steps[step].frames
class adultTranHalfAxis:
    def __init__(self):
        self.odbName = 'adultICRMouseDuraMaterTranCutHalfAxis'
        self.sets = ['XZCUTFACE_1','XZCUTFACE_2']
        self.nodeBounds = [9, 10, 7, 8]
        self.odb = openOdb(self.odbName+'.odb', readOnly=False)
        self.frames = self.odb.steps[step].frames
class adultLong32ndAxis:
    def __init__(self):
        self.odbName = 'adultICRMouseDuraMaterLongCut32ndAxis'
        self.sets = ['YZCUTFACE_1','YZCUTFACE_2']
        self.nodeBounds = [11, 12, 7, 9]
        self.odb = openOdb(self.odbName+'.odb', readOnly=False)
        self.frames = self.odb.steps[step].frames
class adultLong16thAxis:
    def __init__(self):
        self.odbName = 'adultICRMouseDuraMaterLongCut16thAxis'
        self.sets = ['YZCUTFACE_1','YZCUTFACE_2']
        self.nodeBounds = [11, 12, 7, 9]
        self.odb = openOdb(self.odbName+'.odb', readOnly=False)
        self.frames = self.odb.steps[step].frames
class adultLong8thAxis:
    def __init__(self):
        self.odbName = 'adultICRMouseDuraMaterLongCut8thAxis'
        self.sets = ['YZCUTFACE_1','YZCUTFACE_2']
        self.nodeBounds = [11, 12, 7, 9]
        self.odb = openOdb(self.odbName+'.odb', readOnly=False)
        self.frames = self.odb.steps[step].frames
class adultLong4thAxis:
    def __init__(self):
        self.odbName = 'adultICRMouseDuraMaterLongCut4thAxis'
        self.sets = ['YZCUTFACE_1','YZCUTFACE_2']
        self.nodeBounds = [11, 12, 7, 9]
        self.odb = openOdb(self.odbName+'.odb', readOnly=False)
        self.frames = self.odb.steps[step].frames
class adultLongHalfAxis:
    def __init__(self):
        self.odbName = 'adultICRMouseDuraMaterLongCutHalfAxis'
        self.sets = ['YZCUTFACE_1','YZCUTFACE_2']
        self.nodeBounds = [11, 12, 7, 9]
        self.odb = openOdb(self.odbName+'.odb', readOnly=False)
        self.frames = self.odb.steps[step].frames
# Effect of curvature vs. incision opening ratio
class a2Tran:
    def __init__(self):
        self.odbName = 'a=2MouseDuraMaterTranCut'
        self.sets = ['XZCUTFACE_1','XZCUTFACE_2']
        self.nodeBounds = [9, 10, 7, 8]
        self.odb = openOdb(self.odbName+'.odb', readOnly=False)
        self.frames = self.odb.steps[step].frames
class a4Tran:
    def __init__(self):
        self.odbName = 'a=4MouseDuraMaterTranCut'
        self.sets = ['XZCUTFACE_1','XZCUTFACE_2']
        self.nodeBounds = [9, 10, 7, 8]
        self.odb = openOdb(self.odbName+'.odb', readOnly=False)
        self.frames = self.odb.steps[step].frames
class a6Tran:
    def __init__(self):
        self.odbName = 'a=6MouseDuraMaterTranCut'
        self.sets = ['XZCUTFACE_1','XZCUTFACE_2']
        self.nodeBounds = [9, 10, 7, 8]
        self.odb = openOdb(self.odbName+'.odb', readOnly=False)
        self.frames = self.odb.steps[step].frames
class a8Tran:
    def __init__(self):
        self.odbName = 'a=8MouseDuraMaterTranCut'
        self.sets = ['XZCUTFACE_1','XZCUTFACE_2']
        self.nodeBounds = [9, 10, 7, 8]
        self.odb = openOdb(self.odbName+'.odb', readOnly=False)
        self.frames = self.odb.steps[step].frames
class a10Tran:
    def __init__(self):
        self.odbName = 'a=10MouseDuraMaterTranCut'
        self.sets = ['XZCUTFACE_1','XZCUTFACE_2']
        self.nodeBounds = [9, 10, 7, 8]
        self.odb = openOdb(self.odbName+'.odb', readOnly=False)
        self.frames = self.odb.steps[step].frames
class a2Long:
    def __init__(self):
        self.odbName = 'a=2MouseDuraMaterLongCut'
        self.sets = ['YZCUTFACE_1','YZCUTFACE_2']
        self.nodeBounds = [11, 12, 7, 9]
        self.odb = openOdb(self.odbName+'.odb', readOnly=False)
        self.frames = self.odb.steps[step].frames
class a4Long:
    def __init__(self):
        self.odbName = 'a=4MouseDuraMaterLongCut'
        self.sets = ['YZCUTFACE_1','YZCUTFACE_2']
        self.nodeBounds = [11, 12, 7, 9]
        self.odb = openOdb(self.odbName+'.odb', readOnly=False)
        self.frames = self.odb.steps[step].frames
class a6Long:
    def __init__(self):
        self.odbName = 'a=6MouseDuraMaterLongCut'
        self.sets = ['YZCUTFACE_1','YZCUTFACE_2']
        self.nodeBounds = [11, 12, 7, 9]
        self.odb = openOdb(self.odbName+'.odb', readOnly=False)
        self.frames = self.odb.steps[step].frames
class a8Long:
    def __init__(self):
        self.odbName = 'a=8MouseDuraMaterLongCut'
        self.sets = ['YZCUTFACE_1','YZCUTFACE_2']
        self.nodeBounds = [11, 12, 7, 9]
        self.odb = openOdb(self.odbName+'.odb', readOnly=False)
        self.frames = self.odb.steps[step].frames
class a10Long:
    def __init__(self):
        self.odbName = 'a=10MouseDuraMaterLongCut'
        self.sets = ['YZCUTFACE_1','YZCUTFACE_2']
        self.nodeBounds = [11, 12, 7, 9]
        self.odb = openOdb(self.odbName+'.odb', readOnly=False)
        self.frames = self.odb.steps[step].frames

#///// Extraction of data
def solve(sim, modelProperties):
    analysis = odb(sim.odbName, duraMaterInstances[0], duraMaterInstances[1], step, elementType, sim.sets[0], sim.sets[1], sim.nodeBounds, sim.frames)
    cutLengthAndWidths_and_stretches(analysis)
    # Cut retraction dimensions 
    outFile = open(sim.odbName +'Results.csv', 'w+')
    header = ['Cut Frame','Location','Nodes1','Nodes2','lambda_p', 'wbar', 'l', 'w']
    output = zip(cutResults[0], cutResults[1], cutResults[2], cutResults[3], cutResults[4], cutResults[5], cutResults[6], cutResults[7])
    writer = csv.writer(outFile)
    writer.writerow(header)
    writer.writerows(output)
    outFile.close()

#///// Solve
# Neonate and adult murine cranial dura mater models
print('Neonate transverse cut')
solve(neonateTransverse(), modelProperties)
print('Neonate longitudinal cut')
solve(neonateLongitudinal(), modelProperties)
print('Adult transverse cut')
solve(adultTransverse(), modelProperties)
print('Adult longitudinal cut')
solve(adultLongitudinal(), modelProperties)
# Effect of initial cut length vs. incision opening ratio
print('Adult transverse 32nd of minor axis')
solve(adultTran32ndAxis(), modelProperties)
print('Adult transverse 16th of minor axis')
solve(adultTran16thAxis(), modelProperties)
print('Adult transverse 8th of minor axis')
solve(adultTran8thAxis(), modelProperties)
print('Adult transverse 4th of minor axis')
solve(adultTran4thAxis(), modelProperties)
print('Adult transverse half of minor axis')
solve(adultTranHalfAxis(), modelProperties)
print('Adult longitudinal 32nd of minor axis')
solve(adultLong32ndAxis(), modelProperties)
print('Adult longitudinal 16th of minor axis')
solve(adultLong16thAxis(), modelProperties)
print('Adult longitudinal 8th of minor axis')
solve(adultLong8thAxis(), modelProperties)
print('Adult longitudinal 4th of minor axis')
solve(adultLong4thAxis(), modelProperties)
print('Adult longitudinal half of minor axis')
solve(adultLongHalfAxis(), modelProperties)
# Effect of curvature vs. incision opening ratio
print('a=2 transverse cut')
solve(a2Tran(), modelProperties)
print('a=4 transverse cut')
solve(a4Tran(), modelProperties)
print('a=6 transverse cut')
solve(a6Tran(), modelProperties)
print('a=8 transverse cut')
solve(a8Tran(), modelProperties)
print('a=10 transverse cut')
solve(a10Tran(), modelProperties)
print('a=2 longitudinal cut')
solve(a2Long(), modelProperties)
print('a=4 longitudinal cut')
solve(a4Long(), modelProperties)
print('a=6 longitudinal cut')
solve(a6Long(), modelProperties)
print('a=8 longitudinal cut')
solve(a8Long(), modelProperties)
print('a=10 longitudinal cut')
solve(a10Long(), modelProperties)