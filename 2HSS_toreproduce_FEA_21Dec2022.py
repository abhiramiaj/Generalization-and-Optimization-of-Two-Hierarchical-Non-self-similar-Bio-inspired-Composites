import random
from abaqus import *
from abaqusConstants import *
import __main__
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
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
Mdb()
##########1H geometry details############
rho=18.00 #Aspect ratio
n = 5 #number of platelets in a period
phi = 0.8 #platelet volume fraction
############2H geometry starts here#############
rho2=20
n2=3
phi2=0.9
#############2H geometry ends here##########
b=1.#width of platelet

Ep = 7000 #Young's modulus of platelet
EpGm = 1000.
Gm = Ep/EpGm
nup = 0.22
num = 0.49
Em = 2*Gm*(1+num)
vg = b/phi - b #matrix thickness measured in y-dirn, (vertical matrix)
hg=0.25*vg#make small slits
w=b/2.0
lb=hg/2

ms=vg/2##mesh size




#Gm=Em/(2*(1+num))
Lp=rho*b
Wrve=n*(vg+b)
Lrve=(2*lb)+Lp
elsizeVert = vg/4. #element size
elsizeHor = elsizeVert*3

path_odb3='C:/Temp/SWAR%dVF%dG12.odb'%(rho,phi*100)
job_name3='SWAR%dVF%dG12'%(rho,phi*100) 
P=-10#load applied for finding out G12
flag=1
path_name='C:/Temp/SWAR%dVF%d'%(rho,phi*100)
path_odb1='C:/Temp/SWAR%dVF%dQ11.odb'%(rho,phi*100)
path_odb2='C:/Temp/SWAR%dVF%dQ22.odb'%(rho,phi*100)
job_name1='SWAR%dVF%dQ11'%(rho,phi*100)
job_name2='SWAR%dVF%dQ22'%(rho,phi*100)
Lp=rho*b



#totalwidth = b*rho + vg (length of rve)
#Wrve = (b+vg)*n
elsizeHor = (Lrve-vg) / int((Lrve-vg)/elsizeHor)
#
# main part
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(0.0, 0.0), 
    point2=(Lrve, Wrve))
mdb.models['Model-1'].Part(dimensionality=TWO_D_PLANAR, name='Part-1', type=
    DEFORMABLE_BODY)
mdb.models['Model-1'].parts['Part-1'].BaseShell(sketch=mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']
#
# cutout part
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(-hg/2., -b/2.), 
    point2=(-hg/2., b/2.))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(-hg/2., b/2.), 
    point2=(0., b/2.+vg))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.0, b/2.+vg), 
    point2=(hg/2., b/2.))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(hg/2., b/2.), 
    point2=(hg/2., -b/2.))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(hg/2., -b/2.), 
    point2=(0.0, -b/2.-vg))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.0, -b/2.-vg), 
    point2=(-hg/2., -b/2.))
cutoutpart = mdb.models['Model-1'].Part(dimensionality=TWO_D_PLANAR, name='cutout', 
    type=DEFORMABLE_BODY)
cutoutpart.BaseShell(sketch=mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']
#
# Assembly
#
mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
asm = mdb.models['Model-1'].rootAssembly
asm.Instance(dependent=ON, name='Part-1-1', 
    part=mdb.models['Model-1'].parts['Part-1'])
asm.Instance(dependent=ON, name='cutoutInst', part=cutoutpart)
newInstance = asm.InstanceFromBooleanCut(cuttingInstances=(asm.instances['cutoutInst'], ), 
    instanceToBeCut=asm.instances['Part-1-1'], 
    name='Part-2', originalInstances=DELETE)
del mdb.models['Model-1'].parts['Part-1']
asm.Instance(dependent=ON, name='cutoutInst', 
    part=cutoutpart)
asm.translate(instanceList=('cutoutInst', ), 
    vector=(Lrve, 0.0, 0.0))
newInstance = asm.InstanceFromBooleanCut(cuttingInstances=(asm.instances['cutoutInst'], ), 
    instanceToBeCut=asm.instances['Part-2-1'], 
    name='Part-3', originalInstances=DELETE)
del mdb.models['Model-1'].parts['Part-2']
#
xcoord = 0.
ycoord = 0.
currentInstance = 3
for stag in range(n-1):
    asm.Instance(dependent=ON, name='cutoutInst', 
        part=cutoutpart)
    ycoord += b+vg
    xcoord += Lrve/n
    asm.translate(instanceList=('cutoutInst', ), 
        vector=(xcoord, ycoord, 0.0))
    newPartname = "Part-%d"%(currentInstance+1)
    asm.InstanceFromBooleanCut(cuttingInstances=(asm.instances['cutoutInst'], ), 
        instanceToBeCut=asm.instances['Part-%d-1'%currentInstance], 
        name=newPartname, originalInstances=DELETE)
    del mdb.models['Model-1'].parts['Part-%d'%currentInstance]
    currentInstance += 1
#
asm.Instance(dependent=ON, name='cutoutInst', 
    part=cutoutpart)
asm.translate(instanceList=('cutoutInst', ), 
    vector=(0.0, Wrve, 0.0))
newPartname = "Part-%d"%(currentInstance+1)
asm.InstanceFromBooleanCut(cuttingInstances=(asm.instances['cutoutInst'], ), 
    instanceToBeCut=asm.instances['Part-%d-1'%currentInstance], 
    name=newPartname, originalInstances=DELETE)
del mdb.models['Model-1'].parts['Part-%d'%currentInstance]
currentInstance += 1
#
asm.Instance(dependent=ON, name='cutoutInst', 
    part=cutoutpart)
asm.translate(instanceList=('cutoutInst', ), 
    vector=(Lrve, Wrve, 0.0))
newPartname = "Part-%d"%(currentInstance+1)
asm.InstanceFromBooleanCut(cuttingInstances=(asm.instances['cutoutInst'], ), 
    instanceToBeCut=asm.instances['Part-%d-1'%currentInstance], 
    name=newPartname, originalInstances=DELETE)
del mdb.models['Model-1'].parts['Part-%d'%currentInstance]
currentInstance += 1

finalPart = mdb.models['Model-1'].parts['Part-%d'%currentInstance]
finalInstance = asm.instances['Part-%d-1'%currentInstance]
#
# Partitioning and set definitions
#
ycoord = b/2.
for i in range(n):
    mdb.models['Model-1'].ConstrainedSketch(gridSpacing=0.5, name='__profile__', 
        sheetSize=Lrve, 
        transform=finalPart.MakeSketchTransform(finalPart.faces[0], 
        sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0)))
    finalPart.projectReferencesOntoSketch(filter=COPLANAR_EDGES, 
        sketch=mdb.models['Model-1'].sketches['__profile__'])
    mdb.models['Model-1'].sketches['__profile__'].Line(
        point1=(0.0, ycoord), 
        point2=(Lrve, ycoord))
    finalPart.PartitionFaceBySketch(faces=finalPart.faces.findAt((hg*2, Wrve-b/4, 0.0),  ),
        sketch=mdb.models['Model-1'].sketches['__profile__'])
    del mdb.models['Model-1'].sketches['__profile__']
    #
    if i==0:
        plateletSet = finalPart.Set(faces=finalPart.faces.findAt(( (hg*2, ycoord-vg/2., 0.0), ), ),
            name='Platelet')
    else:
        addSet = finalPart.Set(faces=finalPart.faces.findAt( ((hg*2, ycoord-vg/2., 0.0), (0.0, 0.0, 1.0)), 
            ((Lrve-hg*2, ycoord-vg/2., 0.0), (0.0, 0.0, 1.0)), ), name='dummy')
        plateletSet = finalPart.SetByBoolean(name='Platelet', sets=(plateletSet, addSet))
    #
    ycoord += vg
    mdb.models['Model-1'].ConstrainedSketch(gridSpacing=0.5, name='__profile__', 
        sheetSize=Lrve, transform=finalPart.MakeSketchTransform(
        finalPart.faces[0], 
        sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0)))
    finalPart.projectReferencesOntoSketch(filter=COPLANAR_EDGES, 
        sketch=mdb.models['Model-1'].sketches['__profile__'])
    mdb.models['Model-1'].sketches['__profile__'].Line(
        point1=(0.0, ycoord), 
        point2=(Lrve, ycoord))
    finalPart.PartitionFaceBySketch(faces=finalPart.faces.findAt((hg*2, Wrve-b/4, 0.0),  ),
        sketch=mdb.models['Model-1'].sketches['__profile__'])
    del mdb.models['Model-1'].sketches['__profile__']
    #
    if i==0:
        softvgSet = finalPart.Set(faces=finalPart.faces.findAt(( (hg*2, ycoord-vg/2., 0.0), ), ),
            name='Softvg')
        addSet = finalPart.Set(faces=finalPart.faces.findAt( 
            ((Lrve-hg*2, ycoord-vg/2., 0.0), (0.0, 0.0, 1.0)), ), name='dummy')
    else:
        addSet = finalPart.Set(faces=finalPart.faces.findAt( 
            ((hg*2, ycoord-vg/2., 0.0), (0.0, 0.0, 1.0)), 
            ((Lrve*(i+0.5)/n, ycoord-vg/2., 0.0), (0.0, 0.0, 1.0)), 
            ((Lrve-hg*2, ycoord-vg/2., 0.0), (0.0, 0.0, 1.0)), ), 
            name='dummy')
    softvgSet = finalPart.SetByBoolean(name='Softvg', sets=(softvgSet, addSet))    
    #
    ycoord += b
addSet = finalPart.Set(faces=finalPart.faces.findAt( ((hg*2, Wrve-b/4., 0.0), (0.0, 0.0, 1.0)), ),
    name='dummy')
plateletSet = finalPart.SetByBoolean(name='Platelet', sets=(plateletSet, addSet))
#
asm.regenerate()
#
# Materials
#
mdb.models['Model-1'].Material(name='Platelet')
mdb.models['Model-1'].materials['Platelet'].Elastic(table=((Ep, nup), ))
mdb.models['Model-1'].Material(name='Softvg')
mdb.models['Model-1'].materials['Softvg'].Elastic(table=((Em, num), ))
mdb.models['Model-1'].HomogeneousSolidSection(material='Platelet', name='Platelet', thickness=None)
mdb.models['Model-1'].HomogeneousSolidSection(material='Softvg', name='Softvg', thickness=None)
finalPart.SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=finalPart.sets['Platelet'], 
    sectionName='Platelet', thicknessAssignment=FROM_SECTION)
finalPart.SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=finalPart.sets['Softvg'], 
    sectionName='Softvg', thicknessAssignment=FROM_SECTION)
#
# Step and Boundary conditions


mdb.models['Model-1'].rootAssembly.features.changeKey(fromName='Part-%d-1'%currentInstance,
	toName='Composite-1')
a1 = mdb.models['Model-1'].rootAssembly
a1.makeIndependent(instances=(a1.instances['Composite-1'], ))
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
a = mdb.models['Model-1'].rootAssembly.instances['Composite-1'].faces
rfm=[0]*len(a)
for i in range(0,len(a)):
    rfm[i]= a[i]
rfmr=tuple(rfm)                                         
a = mdb.models['Model-1'].rootAssembly
partInstances =(a.instances['Composite-1'], )
elemType1 = mesh.ElemType(elemCode=CPE4R, elemLibrary=STANDARD, 
        secondOrderAccuracy=OFF, hourglassControl=DEFAULT, 
        distortionControl=DEFAULT)
elemType2 = mesh.ElemType(elemCode=CPE3, elemLibrary=STANDARD)
a.setMeshControls(regions=rfmr, elemShape=QUAD,technique=STRUCTURED)
a.setElementType(regions=rfmr, elemTypes=(elemType1, elemType2))#uncomment to activate plane strain element
a.seedPartInstance(regions=partInstances, size=ms, deviationFactor=0.1, 
    minSizeFactor=0.1)
a.generateMesh(regions=partInstances)

session.viewports['Viewport: 1'].setValues(displayedObject=a)


#








mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial')
        



finalInstance=mdb.models['Model-1'].rootAssembly.instances['Composite-1']
#finalPart= mdb.models['Model-1'].parts['Composite']
a = mdb.models['Model-1'].rootAssembly
a.Set(edges=finalInstance.edges.findAt(((Lrve/2., Wrve, 0.0),), ), name="top")
a.Set(edges=finalInstance.edges.getByBoundingBox(xMin=Lrve-0.), name="right")
a.Set(edges=finalInstance.edges.getByBoundingBox(xMax=hg/5.), name="left")
a.Set(edges=finalInstance.edges.getByBoundingBox(yMax=vg), name="bottom")


# fixnode=mdb.models['Model-1'].rootAssembly.instances['Composite-1'].nodes.getClosest(coordinates=((Lrve/2., Wrve/2., 0.0),))[0].label
# a.Set(nodes=mdb.models['Model-1'].rootAssembly.instances['Composite-1'].nodes[fixnode-1:fixnode], name="fixNode")


# mdb.models['Model-1'].PinnedBC(createStepName='Step-1', localCsys=None, name='allfix',region=a.sets["fixNode"])
# generate the equations manually

nodesets = {}
for set,i in zip(["top","bottom","right","left"],[1,1,2,2]):
    nodesets[set]=[]
    for node in a.sets[set].nodes:
        nodesets[set].append((node.label,node.coordinates[0],node.coordinates[1]))
    nodesets[set].sort(key=lambda value: value[i])
mdb.models['Model-1'].keywordBlock.synchVersions(True)
keywordblock = mdb.models['Model-1'].keywordBlock
for i,kw in enumerate(keywordblock.sieBlocks):
    if kw.startswith("*End Part"):
        endpart = i
    elif kw.startswith("*End Instance"):
        endinstance = i
    elif kw.startswith("*Elset"):
        endassembly = i
    elif kw.startswith("*Static"):
        static = i
    elif kw.startswith("*Output, history"):
        history = i
        break


partstring="*Part, name=dummy-LR\n"
partstring+="*End Part\n"
partstring+="**\n"
partstring+="*Part, name=dummy-TB\n"
partstring+="*End Part\n"
partstring+="**\n"
keywordblock.insert(position=endpart+1,text=partstring)
endinstancestring="*Instance, name=dummy-LR-1, part=dummy-LR\n"
endinstancestring+="*Node\n"
endinstancestring+="100000, -10., 10., 0.\n"
endinstancestring+="**This dummy node can be arbitrary\n"
endinstancestring+="*Nset, nset=dummy-LR-1-RefPt_, internal\n"
endinstancestring+="100000,\n"
endinstancestring+="*End Instance\n"
endinstancestring+="*Instance, name=dummy-TB-1, part=dummy-TB\n"
endinstancestring+="*Node\n"
endinstancestring+="200000, 10., 0., 0.\n"
endinstancestring+="**This dummy node can be arbitrary\n"
endinstancestring+="*Nset, nset=dummy-TB-1-RefPt_, internal\n"
endinstancestring+="200000,\n"
endinstancestring+="*End Instance\n"
endinstancestring+="** Define nset Set-dummy-LR and Set-dummy-TB for the two dummy nodes\n"
endinstancestring+="**\n"
endinstancestring+="*Nset, nset=Set-dummy-LR, instance=dummy-LR-1\n"
endinstancestring+="100000,\n"
endinstancestring+="*Nset, nset=Set-dummy-TB, instance=dummy-TB-1\n"
endinstancestring+="200000,\n"
keywordblock.insert(position=endinstance+1,text=endinstancestring)

equationstring = "*Equation\n"
for ntop,nbot in zip(nodesets["top"],nodesets["bottom"]):
    equationstring += "3\n"
    equationstring += "Composite-1.%d, 1, 1.0, Composite-1.%d, 1, -1.0, %s, 1, -1.0\n" % (ntop[0],nbot[0],'Set-dummy-TB')
    equationstring += "3\n"
    equationstring += "Composite-1.%d, 2, 1.0, Composite-1.%d, 2, -1.0, %s, 2, -1.0\n" % (ntop[0],nbot[0],'Set-dummy-TB')
for nright,nleft in zip(nodesets["right"],nodesets["left"]):
    equationstring += "3\n"
    equationstring += "Composite-1.%d, 1, 1.0, Composite-1.%d, 1, -1.0, %s, 1, -1.0\n" % (nright[0],nleft[0],'Set-dummy-LR')
    equationstring += "3\n"
    equationstring += "Composite-1.%d, 2, 1.0, Composite-1.%d, 2, -1.0, %s, 2, -1.0\n" % (nright[0],nleft[0],'Set-dummy-LR')
keywordblock.insert(position=endassembly,text=equationstring)
bcstring = "**\n*Boundary\n"
bcstring += "Set-dummy-TB,1,1, 0.0\n" 
bcstring += "**\n*Boundary\n"
bcstring += "Set-dummy-LR,1,1, 0.01\n"
bcstring += "**\n*Boundary\n"
bcstring += "Set-dummy-LR,2,2, 0.0\n" 
bcstring += "**\n*Boundary\n"
bcstring += "Set-dummy-TB,2,2, 0.0\n"
keywordblock.insert(position=static+3,text=bcstring)
# ###################'Set-dummy-TB'##########################

job=mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
    multiprocessingMode=DEFAULT, name=job_name1, nodalOutputPrecision=SINGLE, 
    numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB,parallelizationMethodExplicit=DOMAIN, 
    numDomains=1,scratch='', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)

historystring = "*Output, history\n*Node Output, Nset=Set-dummy-LR\n RF\n" 
historystring += "*Node Output, Nset=Set-dummy-TB\n RF\n" 
keywordblock.replace(position=history+4,text=historystring)
mdb.jobs[job_name1].submit(consistencyChecking=OFF)
mdb.jobs[job_name1].waitForCompletion()


mdb.models['Model-1'].keywordBlock.setValues(edited = 0)


nodesets = {}
for set,i in zip(["top","bottom","right","left"],[1,1,2,2]):
    nodesets[set]=[]
    for node in a.sets[set].nodes:
        nodesets[set].append((node.label,node.coordinates[0],node.coordinates[1]))
    nodesets[set].sort(key=lambda value: value[i])
mdb.models['Model-1'].keywordBlock.synchVersions(True)
keywordblock = mdb.models['Model-1'].keywordBlock
for i,kw in enumerate(keywordblock.sieBlocks):
    if kw.startswith("*End Part"):
        endpart = i
    elif kw.startswith("*End Instance"):
        endinstance = i
    elif kw.startswith("*Elset"):
        endassembly = i
    elif kw.startswith("*Static"):
        static = i
    elif kw.startswith("*Output, history"):
        history = i
        break


partstring="*Part, name=dummy-LR\n"
partstring+="*End Part\n"
partstring+="**\n"
partstring+="*Part, name=dummy-TB\n"
partstring+="*End Part\n"
partstring+="**\n"
keywordblock.insert(position=endpart+1,text=partstring)
endinstancestring="*Instance, name=dummy-LR-1, part=dummy-LR\n"
endinstancestring+="*Node\n"
endinstancestring+="100000, -10., 10., 0.\n"
endinstancestring+="**This dummy node can be arbitrary\n"
endinstancestring+="*Nset, nset=dummy-LR-1-RefPt_, internal\n"
endinstancestring+="100000,\n"
endinstancestring+="*End Instance\n"
endinstancestring+="*Instance, name=dummy-TB-1, part=dummy-TB\n"
endinstancestring+="*Node\n"
endinstancestring+="200000, 10., 0., 0.\n"
endinstancestring+="**This dummy node can be arbitrary\n"
endinstancestring+="*Nset, nset=dummy-TB-1-RefPt_, internal\n"
endinstancestring+="200000,\n"
endinstancestring+="*End Instance\n"
endinstancestring+="** Define nset Set-dummy-LR and Set-dummy-TB for the two dummy nodes\n"
endinstancestring+="**\n"
endinstancestring+="*Nset, nset=Set-dummy-LR, instance=dummy-LR-1\n"
endinstancestring+="100000,\n"
endinstancestring+="*Nset, nset=Set-dummy-TB, instance=dummy-TB-1\n"
endinstancestring+="200000,\n"

keywordblock.insert(position=endinstance+1,text=endinstancestring)

equationstring = "*Equation\n"
for ntop,nbot in zip(nodesets["top"],nodesets["bottom"]):
    equationstring += "3\n"
    equationstring += "Composite-1.%d, 1, 1.0, Composite-1.%d, 1, -1.0, %s, 1, -1.0\n" % (ntop[0],nbot[0],'Set-dummy-TB')
    equationstring += "3\n"
    equationstring += "Composite-1.%d, 2, 1.0, Composite-1.%d, 2, -1.0, %s, 2, -1.0\n" % (ntop[0],nbot[0],'Set-dummy-TB')
for nright,nleft in zip(nodesets["right"],nodesets["left"]):
    equationstring += "3\n"
    equationstring += "Composite-1.%d, 1, 1.0, Composite-1.%d, 1, -1.0, %s, 1, -1.0\n" % (nright[0],nleft[0],'Set-dummy-LR')
    equationstring += "3\n"
    equationstring += "Composite-1.%d, 2, 1.0, Composite-1.%d, 2, -1.0, %s, 2, -1.0\n" % (nright[0],nleft[0],'Set-dummy-LR')
keywordblock.insert(position=endassembly,text=equationstring)
bcstring = "**\n*Boundary\n"
bcstring += "Set-dummy-TB,2,2, 0.01\n" 
bcstring += "**\n*Boundary\n"
bcstring += "Set-dummy-TB,1,1, 0.0\n"
bcstring += "**\n*Boundary\n"
bcstring += "Set-dummy-LR,2,2, 0.0\n" 
bcstring += "**\n*Boundary\n"
bcstring += "Set-dummy-LR,1,1, 0.0\n" 
keywordblock.insert(position=static+3,text=bcstring)


mdb.Job(name=job_name2, model='Model-1', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB,parallelizationMethodExplicit=DOMAIN, 
    numDomains=1,multiprocessingMode=DEFAULT, numCpus=1, 
    numGPUs=0)
historystring = "*Output, history\n*Node Output, Nset=Set-dummy-LR\n RF\n" 
historystring += "*Node Output, Nset=Set-dummy-TB\n RF\n"
keywordblock.replace(position=history+4,text=historystring)
mdb.jobs[job_name2].submit(consistencyChecking=OFF)
mdb.jobs[job_name2].waitForCompletion()


mesh_nodes=mdb.models['Model-1'].rootAssembly.instances['Composite-1'].nodes
a.Set(nodes=mesh_nodes.getByBoundingBox(xMin=Lrve-hg/2, yMin=vg/4, yMax=Wrve - vg/2.), name="BC")
a.Set(nodes=mesh_nodes.getByBoundingBox(xMax=hg,yMin=vg/4, yMax=Wrve - vg/2.), name="AD")
a.Set(nodes=mesh_nodes.getByBoundingBox(xMin=hg, xMax=Lrve - hg, yMax=vg/4.), name="AB")
a.Set(nodes=mesh_nodes.getByBoundingBox(yMin=Wrve - hg, xMin=hg, xMax=Lrve - hg), name="CD")


a.Set(nodes=mesh_nodes.getByBoundingBox(yMin=0, xMin=lb, xMax=lb, yMax=0.), name="A")
a.Set(nodes=mesh_nodes.getByBoundingBox(yMin=0, xMin=Lp+lb, xMax=Lp+lb, yMax=0.), name="B")
a.Set(nodes=mesh_nodes.getByBoundingBox(yMin=Wrve, xMin=Lp+lb, xMax=Lp+lb, yMax=Wrve), name="C") 
a.Set(nodes=mesh_nodes.getByBoundingBox(yMin=Wrve, xMin=lb, xMax=lb, yMax=Wrve), name="D") 



session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')
#mdb.models['Model-1'].boundaryConditions['allfix'].suppress()
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF, loads=ON, 
    bcs=ON, predefinedFields=ON, connectors=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
        meshTechnique=OFF)
mdb.models['Model-1'].keywordBlock.synchVersions(storeNodesAndElements=False)
mdb.models['Model-1'].keywordBlock.setValues(edited = 0)
mdb.models['Model-1'].keywordBlock.synchVersions(storeNodesAndElements=False)
a = mdb.models['Model-1'].rootAssembly
region = a.sets['A']
mdb.models['Model-1'].DisplacementBC(name='Au2zer', createStepName='Step-1', 
    region=region, u1=UNSET, u2=0.0, ur3=UNSET, amplitude=UNSET, fixed=OFF, 
                distributionType=UNIFORM, fieldName='', localCsys=None)
region = a.sets['B']
mdb.models['Model-1'].DisplacementBC(name='Bu2zer', createStepName='Step-1', 
    region=region, u1=UNSET, u2=0.0, ur3=UNSET, amplitude=UNSET, fixed=OFF, 
    distributionType=UNIFORM, fieldName='', localCsys=None)
region = a.sets['AB']
mdb.models['Model-1'].DisplacementBC(name='ABu1zer', createStepName='Step-1', 
    region=region, u1=0.0, u2=UNSET, ur3=UNSET, amplitude=UNSET, fixed=OFF, 
    distributionType=UNIFORM, fieldName='', localCsys=None)
region = a.sets['BC']
mdb.models['Model-1'].DisplacementBC(name='BCu2zer', createStepName='Step-1', 
    region=region, u1=UNSET, u2=0.0, ur3=UNSET, amplitude=UNSET, fixed=OFF, 
    distributionType=UNIFORM, fieldName='', localCsys=None)
region = a.sets['AD']
mdb.models['Model-1'].DisplacementBC(name='ADu2zer', createStepName='Step-1', 
    region=region, u1=UNSET, u2=0.0, ur3=UNSET, amplitude=UNSET, fixed=OFF, 
    distributionType=UNIFORM, fieldName='', localCsys=None)
region = a.sets['C']
mdb.models['Model-1'].ConcentratedForce(name='C', createStepName='Step-1', 
    region=region, cf1=P, distributionType=UNIFORM, field='', 
    localCsys=None)
region = a.sets['D']
mdb.models['Model-1'].ConcentratedForce(name='D', createStepName='Step-1', 
    region=region, cf1=P, distributionType=UNIFORM, field='', 
    localCsys=None)

# mdb.models['Model-1'].keywordBlock.setValues(edited = 0)


nodesets = {}
for set,i in zip(["top","bottom","right","left"],[1,1,2,2]):
    nodesets[set]=[]
    for node in a.sets[set].nodes:
        nodesets[set].append((node.label,node.coordinates[0],node.coordinates[1]))
    nodesets[set].sort(key=lambda value: value[i])
mdb.models['Model-1'].keywordBlock.synchVersions(True)
keywordblock = mdb.models['Model-1'].keywordBlock
for i,kw in enumerate(keywordblock.sieBlocks):
    if kw.startswith("*End Part"):
        endpart = i
    elif kw.startswith("*End Instance"):
        endinstance = i
    elif kw.startswith("*Elset"):
        endassembly = i
    elif kw.startswith("*Static"):
        static = i
    elif kw.startswith("*Output, history"):
        history = i
        break


equationstring = "*Equation\n"
for i in range(0,(len(nodesets["top"])-1)):
    equationstring += "2\n"
    equationstring += "Composite-1.%d, 1, 1.0, Composite-1.%d, 1, -1.0\n" % (nodesets["top"][i][0],nodesets["top"][i+1][0])
keywordblock.insert(position=endassembly,text=equationstring)



mdb.Job(name=job_name3, model='Model-1', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB,parallelizationMethodExplicit=DOMAIN, 
    numDomains=1,multiprocessingMode=DEFAULT, numCpus=1, 
    numGPUs=0)
historystring = "*Output, history\n" 
historystring += "*Node Output, Nset=D\n U\n"
keywordblock.replace(position=history,text=historystring)
mdb.jobs[job_name3].submit(consistencyChecking=OFF)
mdb.jobs[job_name3].waitForCompletion()


o3 = session.openOdb(name=path_odb1)
step = o3.steps["Step-1"]
horiRegion1 = step.historyRegions["Node DUMMY-LR-1.100000"]
vertRegion1 = step.historyRegions["Node DUMMY-TB-1.200000"]
fxx1 = horiRegion1.historyOutputs["RF1"].data[-1][1]
fxy1 = horiRegion1.historyOutputs["RF2"].data[-1][1]
fyx1 = vertRegion1.historyOutputs["RF1"].data[-1][1]
fyy1 = vertRegion1.historyOutputs["RF2"].data[-1][1]
exx=0.01/Lrve
Q11=fxx1/(Wrve*exx)
Q12=(fyy1)/(Lrve*exx)

with open("sw.txt","a") as outfile:
	outfile.write("----------------Results of 1st Hierarchy----------------------\n")
	outfile.write("------Simulation_using_RF,plane strain elements-----\n")
	outfile.write("n=%d, aspect=%d, vf=%.2f, Ep/Gm=%.1f\n" % (n,rho,phi, EpGm))
	outfile.write("---------------Horizontal displacement-----------------------\n")
	outfile.write("Q11=%2f\n" %(Q11))
	outfile.write("Q12=%2f\n" %(Q12))
	outfile.write("--------------------------------------\n")
o4 = session.openOdb(name=path_odb2)
step = o4.steps["Step-1"]
horiRegion2 = step.historyRegions["Node DUMMY-LR-1.100000"]
vertRegion2 = step.historyRegions["Node DUMMY-TB-1.200000"]
fxx2 = horiRegion2.historyOutputs["RF1"].data[-1][1]
fxy2 = horiRegion2.historyOutputs["RF2"].data[-1][1]
fyx2 = vertRegion2.historyOutputs["RF1"].data[-1][1]
fyy2 = vertRegion2.historyOutputs["RF2"].data[-1][1]
eyy=0.01/Wrve
Q12_vs=fxx2/(Wrve*eyy)
Q22=fyy2/(Lrve*eyy)
deno=(Q11*Q22) - (Q12*Q12)
S11=Q22/deno
S12=-Q12/deno
S22=Q11/deno
E11=1/S11
E22=1/S22
nu12=-S12/S11


with open("sw.txt","a") as outfile:
    outfile.write("---------------Vertical displacement-----------------------\n")
    outfile.write("Q12_vs=%2f\n" %(Q12_vs))
    outfile.write("Q22=%2f\n" %(Q22))
    outfile.write("E11=%2f\n" %(E11))
    outfile.write("E22=%2f\n" %(E22))
    outfile.write("nu12=%2f\n" %(nu12))
    outfile.write("--------------------------------------\n")


o3 = session.openOdb(name=path_odb3)
step = o3.steps["Step-1"]
U1reg = step.historyRegions["Node COMPOSITE-1.%d"%(nodesets["top"][0][0])]
delta = U1reg.historyOutputs["U1"].data[-1][1]
G12=2*P*Wrve/(Lrve*delta)

with open("sw.txt","a") as outfile:
    outfile.write("G12=%2f\n" %(2*P*Wrve/(Lrve*delta)))
    outfile.write("----------------End of simulation for 1st Hierarchy----------------------\n")


#######################################2H-SECOND HIERARCHY STARTS HERE###############

import random
Mdb()
rho=rho2 #Aspect ratio @ 2H
b=1.#width of platelet
n = n2 #number of platelets in a period, at second level of hierarchy
phi = phi2 #platelet volume fraction


num = 0.49

vg = b/phi - b #matrix thickness measured in y-dirn, (vertical matrix)
hg=0.25*vg#make small slits
w=b/2.0
lb=hg/2
ms=vg/2##mesh size


Wrve=n*(vg+b)
Lrve=(2*lb)+Lp
elsizeVert = vg/4. #element size
elsizeHor = elsizeVert*3

path_odb3='C:/Temp/SWAR%dVF%dn%dG12.odb'%(rho,phi*100,n)
job_name3='SWAR%dVF%dn%dG12'%(rho,phi*100,n) 
P=-10#load applied for finding out G12
flag=1
path_name='C:/Temp/SWAR%dVF%dn%d'%(rho,phi*100,n)
path_odb1='C:/Temp/SWAR%dVF%dn%dQ11.odb'%(rho,phi*100,n)
path_odb2='C:/Temp/SWAR%dVF%dn%dQ22.odb'%(rho,phi*100,n)
job_name1='SWAR%dVF%dn%dQ11'%(rho,phi*100,n)
job_name2='SWAR%dVF%dn%dQ22'%(rho,phi*100,n) 
Lp=rho*b


#totalwidth = b*rho + vg (length of rve)
#Wrve = (b+vg)*n
elsizeHor = (Lrve-vg) / int((Lrve-vg)/elsizeHor)
#
# main part
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(0.0, 0.0), 
    point2=(Lrve, Wrve))
mdb.models['Model-1'].Part(dimensionality=TWO_D_PLANAR, name='Part-1', type=
    DEFORMABLE_BODY)
mdb.models['Model-1'].parts['Part-1'].BaseShell(sketch=mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']
#
# cutout part
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(-hg/2., -b/2.), 
    point2=(-hg/2., b/2.))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(-hg/2., b/2.), 
    point2=(0., b/2.+vg))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.0, b/2.+vg), 
    point2=(hg/2., b/2.))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(hg/2., b/2.), 
    point2=(hg/2., -b/2.))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(hg/2., -b/2.), 
    point2=(0.0, -b/2.-vg))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.0, -b/2.-vg), 
    point2=(-hg/2., -b/2.))
cutoutpart = mdb.models['Model-1'].Part(dimensionality=TWO_D_PLANAR, name='cutout', 
    type=DEFORMABLE_BODY)
cutoutpart.BaseShell(sketch=mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']
#
# Assembly
#
mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
asm = mdb.models['Model-1'].rootAssembly
asm.Instance(dependent=ON, name='Part-1-1', 
    part=mdb.models['Model-1'].parts['Part-1'])
asm.Instance(dependent=ON, name='cutoutInst', part=cutoutpart)
newInstance = asm.InstanceFromBooleanCut(cuttingInstances=(asm.instances['cutoutInst'], ), 
    instanceToBeCut=asm.instances['Part-1-1'], 
    name='Part-2', originalInstances=DELETE)
del mdb.models['Model-1'].parts['Part-1']
asm.Instance(dependent=ON, name='cutoutInst', 
    part=cutoutpart)
asm.translate(instanceList=('cutoutInst', ), 
    vector=(Lrve, 0.0, 0.0))
newInstance = asm.InstanceFromBooleanCut(cuttingInstances=(asm.instances['cutoutInst'], ), 
    instanceToBeCut=asm.instances['Part-2-1'], 
    name='Part-3', originalInstances=DELETE)
del mdb.models['Model-1'].parts['Part-2']
#
xcoord = 0.
ycoord = 0.
currentInstance = 3
for stag in range(n-1):
    asm.Instance(dependent=ON, name='cutoutInst', 
        part=cutoutpart)
    ycoord += b+vg
    xcoord += Lrve/n
    asm.translate(instanceList=('cutoutInst', ), 
        vector=(xcoord, ycoord, 0.0))
    newPartname = "Part-%d"%(currentInstance+1)
    asm.InstanceFromBooleanCut(cuttingInstances=(asm.instances['cutoutInst'], ), 
        instanceToBeCut=asm.instances['Part-%d-1'%currentInstance], 
        name=newPartname, originalInstances=DELETE)
    del mdb.models['Model-1'].parts['Part-%d'%currentInstance]
    currentInstance += 1
#
asm.Instance(dependent=ON, name='cutoutInst', 
    part=cutoutpart)
asm.translate(instanceList=('cutoutInst', ), 
    vector=(0.0, Wrve, 0.0))
newPartname = "Part-%d"%(currentInstance+1)
asm.InstanceFromBooleanCut(cuttingInstances=(asm.instances['cutoutInst'], ), 
    instanceToBeCut=asm.instances['Part-%d-1'%currentInstance], 
    name=newPartname, originalInstances=DELETE)
del mdb.models['Model-1'].parts['Part-%d'%currentInstance]
currentInstance += 1
#
asm.Instance(dependent=ON, name='cutoutInst', 
    part=cutoutpart)
asm.translate(instanceList=('cutoutInst', ), 
    vector=(Lrve, Wrve, 0.0))
newPartname = "Part-%d"%(currentInstance+1)
asm.InstanceFromBooleanCut(cuttingInstances=(asm.instances['cutoutInst'], ), 
    instanceToBeCut=asm.instances['Part-%d-1'%currentInstance], 
    name=newPartname, originalInstances=DELETE)
del mdb.models['Model-1'].parts['Part-%d'%currentInstance]
currentInstance += 1

finalPart = mdb.models['Model-1'].parts['Part-%d'%currentInstance]
finalInstance = asm.instances['Part-%d-1'%currentInstance]
#
# Partitioning and set definitions
#
ycoord = b/2.
for i in range(n):
    mdb.models['Model-1'].ConstrainedSketch(gridSpacing=0.5, name='__profile__', 
        sheetSize=Lrve, 
        transform=finalPart.MakeSketchTransform(finalPart.faces[0], 
        sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0)))
    finalPart.projectReferencesOntoSketch(filter=COPLANAR_EDGES, 
        sketch=mdb.models['Model-1'].sketches['__profile__'])
    mdb.models['Model-1'].sketches['__profile__'].Line(
        point1=(0.0, ycoord), 
        point2=(Lrve, ycoord))
    finalPart.PartitionFaceBySketch(faces=finalPart.faces.findAt((hg*2, Wrve-b/4, 0.0),  ),
        sketch=mdb.models['Model-1'].sketches['__profile__'])
    del mdb.models['Model-1'].sketches['__profile__']
    #
    if i==0:
        plateletSet = finalPart.Set(faces=finalPart.faces.findAt(( (hg*2, ycoord-vg/2., 0.0), ), ),
            name='Platelet')
    else:
        addSet = finalPart.Set(faces=finalPart.faces.findAt( ((hg*2, ycoord-vg/2., 0.0), (0.0, 0.0, 1.0)), 
            ((Lrve-hg*2, ycoord-vg/2., 0.0), (0.0, 0.0, 1.0)), ), name='dummy')
        plateletSet = finalPart.SetByBoolean(name='Platelet', sets=(plateletSet, addSet))
    #
    ycoord += vg
    mdb.models['Model-1'].ConstrainedSketch(gridSpacing=0.5, name='__profile__', 
        sheetSize=Lrve, transform=finalPart.MakeSketchTransform(
        finalPart.faces[0], 
        sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0)))
    finalPart.projectReferencesOntoSketch(filter=COPLANAR_EDGES, 
        sketch=mdb.models['Model-1'].sketches['__profile__'])
    mdb.models['Model-1'].sketches['__profile__'].Line(
        point1=(0.0, ycoord), 
        point2=(Lrve, ycoord))
    finalPart.PartitionFaceBySketch(faces=finalPart.faces.findAt((hg*2, Wrve-b/4, 0.0),  ),
        sketch=mdb.models['Model-1'].sketches['__profile__'])
    del mdb.models['Model-1'].sketches['__profile__']
    #
    if i==0:
        softvgSet = finalPart.Set(faces=finalPart.faces.findAt(( (hg*2, ycoord-vg/2., 0.0), ), ),
            name='Softvg')
        addSet = finalPart.Set(faces=finalPart.faces.findAt( 
            ((Lrve-hg*2, ycoord-vg/2., 0.0), (0.0, 0.0, 1.0)), ), name='dummy')
    else:
        addSet = finalPart.Set(faces=finalPart.faces.findAt( 
            ((hg*2, ycoord-vg/2., 0.0), (0.0, 0.0, 1.0)), 
            ((Lrve*(i+0.5)/n, ycoord-vg/2., 0.0), (0.0, 0.0, 1.0)), 
            ((Lrve-hg*2, ycoord-vg/2., 0.0), (0.0, 0.0, 1.0)), ), 
            name='dummy')
    softvgSet = finalPart.SetByBoolean(name='Softvg', sets=(softvgSet, addSet))    
    #
    ycoord += b
addSet = finalPart.Set(faces=finalPart.faces.findAt( ((hg*2, Wrve-b/4., 0.0), (0.0, 0.0, 1.0)), ),
    name='dummy')
plateletSet = finalPart.SetByBoolean(name='Platelet', sets=(plateletSet, addSet))
#
asm.regenerate()
#
# Materials
#
mdb.models['Model-1'].Material(name='Platelet')
mdb.models['Model-1'].materials['Platelet'].Elastic(type=ORTHOTROPIC, table=((
        Q11, Q12, Q22, 1e-10, 1e-10, 1e-10, G12, 1e-10, 1e-10), )) #orthotropic
#mdb.models['Model-1'].materials['Platelet'].Elastic(table=((Ep, nup), ))
mdb.models['Model-1'].Material(name='Softvg')
mdb.models['Model-1'].materials['Softvg'].Elastic(table=((Em, num), ))
mdb.models['Model-1'].HomogeneousSolidSection(material='Platelet', name='Platelet', thickness=None)
mdb.models['Model-1'].HomogeneousSolidSection(material='Softvg', name='Softvg', thickness=None)
finalPart.SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=finalPart.sets['Platelet'], 
    sectionName='Platelet', thicknessAssignment=FROM_SECTION)

orientation=None

finalPart.MaterialOrientation(region=finalPart.sets['Platelet'], 
    orientationType=GLOBAL, axis=AXIS_3, 
    additionalRotationType=ROTATION_NONE, localCsys=None, fieldName='', 
    stackDirection=STACK_3)
finalPart.SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=finalPart.sets['Softvg'], 
    sectionName='Softvg', thicknessAssignment=FROM_SECTION)
#
# Step and Boundary conditions


mdb.models['Model-1'].rootAssembly.features.changeKey(fromName='Part-%d-1'%currentInstance,
	toName='Composite-1')
a1 = mdb.models['Model-1'].rootAssembly
a1.makeIndependent(instances=(a1.instances['Composite-1'], ))
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
a = mdb.models['Model-1'].rootAssembly.instances['Composite-1'].faces
rfm=[0]*len(a)
for i in range(0,len(a)):
    rfm[i]= a[i]
rfmr=tuple(rfm)                                         
a = mdb.models['Model-1'].rootAssembly
partInstances =(a.instances['Composite-1'], )
elemType1 = mesh.ElemType(elemCode=CPE4R, elemLibrary=STANDARD, 
        secondOrderAccuracy=OFF, hourglassControl=DEFAULT, 
        distortionControl=DEFAULT)
elemType2 = mesh.ElemType(elemCode=CPE3, elemLibrary=STANDARD)
a.setMeshControls(regions=rfmr, elemShape=QUAD,technique=STRUCTURED)
a.setElementType(regions=rfmr, elemTypes=(elemType1, elemType2))#uncomment to activate plane strain element
a.seedPartInstance(regions=partInstances, size=ms, deviationFactor=0.1, 
    minSizeFactor=0.1)
a.generateMesh(regions=partInstances)

session.viewports['Viewport: 1'].setValues(displayedObject=a)

mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial')
        



finalInstance=mdb.models['Model-1'].rootAssembly.instances['Composite-1']
#finalPart= mdb.models['Model-1'].parts['Composite']
a = mdb.models['Model-1'].rootAssembly
a.Set(edges=finalInstance.edges.findAt(((Lrve/2., Wrve, 0.0),), ), name="top")
a.Set(edges=finalInstance.edges.getByBoundingBox(xMin=Lrve-0.), name="right")
a.Set(edges=finalInstance.edges.getByBoundingBox(xMax=hg/5.), name="left")
a.Set(edges=finalInstance.edges.getByBoundingBox(yMax=vg), name="bottom")


# fixnode=mdb.models['Model-1'].rootAssembly.instances['Composite-1'].nodes.getClosest(coordinates=((Lrve/2., Wrve/2., 0.0),))[0].label
# a.Set(nodes=mdb.models['Model-1'].rootAssembly.instances['Composite-1'].nodes[fixnode-1:fixnode], name="fixNode")


# mdb.models['Model-1'].PinnedBC(createStepName='Step-1', localCsys=None, name='allfix',region=a.sets["fixNode"])
# generate the equations manually

nodesets = {}
for set,i in zip(["top","bottom","right","left"],[1,1,2,2]):
    nodesets[set]=[]
    for node in a.sets[set].nodes:
        nodesets[set].append((node.label,node.coordinates[0],node.coordinates[1]))
    nodesets[set].sort(key=lambda value: value[i])
mdb.models['Model-1'].keywordBlock.synchVersions(True)
keywordblock = mdb.models['Model-1'].keywordBlock
for i,kw in enumerate(keywordblock.sieBlocks):
    if kw.startswith("*End Part"):
        endpart = i
    elif kw.startswith("*End Instance"):
        endinstance = i
    elif kw.startswith("*Elset"):
        endassembly = i
    elif kw.startswith("*Static"):
        static = i
    elif kw.startswith("*Output, history"):
        history = i
        break


partstring="*Part, name=dummy-LR\n"
partstring+="*End Part\n"
partstring+="**\n"
partstring+="*Part, name=dummy-TB\n"
partstring+="*End Part\n"
partstring+="**\n"
keywordblock.insert(position=endpart+1,text=partstring)
endinstancestring="*Instance, name=dummy-LR-1, part=dummy-LR\n"
endinstancestring+="*Node\n"
endinstancestring+="100000, -10., 10., 0.\n"
endinstancestring+="**This dummy node can be arbitrary\n"
endinstancestring+="*Nset, nset=dummy-LR-1-RefPt_, internal\n"
endinstancestring+="100000,\n"
endinstancestring+="*End Instance\n"
endinstancestring+="*Instance, name=dummy-TB-1, part=dummy-TB\n"
endinstancestring+="*Node\n"
endinstancestring+="200000, 10., 0., 0.\n"
endinstancestring+="**This dummy node can be arbitrary\n"
endinstancestring+="*Nset, nset=dummy-TB-1-RefPt_, internal\n"
endinstancestring+="200000,\n"
endinstancestring+="*End Instance\n"
endinstancestring+="** Define nset Set-dummy-LR and Set-dummy-TB for the two dummy nodes\n"
endinstancestring+="**\n"
endinstancestring+="*Nset, nset=Set-dummy-LR, instance=dummy-LR-1\n"
endinstancestring+="100000,\n"
endinstancestring+="*Nset, nset=Set-dummy-TB, instance=dummy-TB-1\n"
endinstancestring+="200000,\n"
keywordblock.insert(position=endinstance+1,text=endinstancestring)

equationstring = "*Equation\n"
for ntop,nbot in zip(nodesets["top"],nodesets["bottom"]):
    equationstring += "3\n"
    equationstring += "Composite-1.%d, 1, 1.0, Composite-1.%d, 1, -1.0, %s, 1, -1.0\n" % (ntop[0],nbot[0],'Set-dummy-TB')
    equationstring += "3\n"
    equationstring += "Composite-1.%d, 2, 1.0, Composite-1.%d, 2, -1.0, %s, 2, -1.0\n" % (ntop[0],nbot[0],'Set-dummy-TB')
for nright,nleft in zip(nodesets["right"],nodesets["left"]):
    equationstring += "3\n"
    equationstring += "Composite-1.%d, 1, 1.0, Composite-1.%d, 1, -1.0, %s, 1, -1.0\n" % (nright[0],nleft[0],'Set-dummy-LR')
    equationstring += "3\n"
    equationstring += "Composite-1.%d, 2, 1.0, Composite-1.%d, 2, -1.0, %s, 2, -1.0\n" % (nright[0],nleft[0],'Set-dummy-LR')
keywordblock.insert(position=endassembly,text=equationstring)
bcstring = "**\n*Boundary\n"
bcstring += "Set-dummy-TB,1,1, 0.0\n" 
bcstring += "**\n*Boundary\n"
bcstring += "Set-dummy-LR,1,1, 0.01\n"
bcstring += "**\n*Boundary\n"
bcstring += "Set-dummy-LR,2,2, 0.0\n" 
bcstring += "**\n*Boundary\n"
bcstring += "Set-dummy-TB,2,2, 0.0\n"
keywordblock.insert(position=static+3,text=bcstring)
# ###################'Set-dummy-TB'##########################

job=mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
    multiprocessingMode=DEFAULT, name=job_name1, nodalOutputPrecision=SINGLE, 
    numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB,parallelizationMethodExplicit=DOMAIN, 
    numDomains=1,scratch='', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)

historystring = "*Output, history\n*Node Output, Nset=Set-dummy-LR\n RF\n" 
historystring += "*Node Output, Nset=Set-dummy-TB\n RF\n" 
keywordblock.replace(position=history+4,text=historystring)
mdb.jobs[job_name1].submit(consistencyChecking=OFF)
mdb.jobs[job_name1].waitForCompletion()


mdb.models['Model-1'].keywordBlock.setValues(edited = 0)


nodesets = {}
for set,i in zip(["top","bottom","right","left"],[1,1,2,2]):
    nodesets[set]=[]
    for node in a.sets[set].nodes:
        nodesets[set].append((node.label,node.coordinates[0],node.coordinates[1]))
    nodesets[set].sort(key=lambda value: value[i])
mdb.models['Model-1'].keywordBlock.synchVersions(True)
keywordblock = mdb.models['Model-1'].keywordBlock
for i,kw in enumerate(keywordblock.sieBlocks):
    if kw.startswith("*End Part"):
        endpart = i
    elif kw.startswith("*End Instance"):
        endinstance = i
    elif kw.startswith("*Elset"):
        endassembly = i
    elif kw.startswith("*Static"):
        static = i
    elif kw.startswith("*Output, history"):
        history = i
        break


partstring="*Part, name=dummy-LR\n"
partstring+="*End Part\n"
partstring+="**\n"
partstring+="*Part, name=dummy-TB\n"
partstring+="*End Part\n"
partstring+="**\n"
keywordblock.insert(position=endpart+1,text=partstring)
endinstancestring="*Instance, name=dummy-LR-1, part=dummy-LR\n"
endinstancestring+="*Node\n"
endinstancestring+="100000, -10., 10., 0.\n"
endinstancestring+="**This dummy node can be arbitrary\n"
endinstancestring+="*Nset, nset=dummy-LR-1-RefPt_, internal\n"
endinstancestring+="100000,\n"
endinstancestring+="*End Instance\n"
endinstancestring+="*Instance, name=dummy-TB-1, part=dummy-TB\n"
endinstancestring+="*Node\n"
endinstancestring+="200000, 10., 0., 0.\n"
endinstancestring+="**This dummy node can be arbitrary\n"
endinstancestring+="*Nset, nset=dummy-TB-1-RefPt_, internal\n"
endinstancestring+="200000,\n"
endinstancestring+="*End Instance\n"
endinstancestring+="** Define nset Set-dummy-LR and Set-dummy-TB for the two dummy nodes\n"
endinstancestring+="**\n"
endinstancestring+="*Nset, nset=Set-dummy-LR, instance=dummy-LR-1\n"
endinstancestring+="100000,\n"
endinstancestring+="*Nset, nset=Set-dummy-TB, instance=dummy-TB-1\n"
endinstancestring+="200000,\n"

keywordblock.insert(position=endinstance+1,text=endinstancestring)

equationstring = "*Equation\n"
for ntop,nbot in zip(nodesets["top"],nodesets["bottom"]):
    equationstring += "3\n"
    equationstring += "Composite-1.%d, 1, 1.0, Composite-1.%d, 1, -1.0, %s, 1, -1.0\n" % (ntop[0],nbot[0],'Set-dummy-TB')
    equationstring += "3\n"
    equationstring += "Composite-1.%d, 2, 1.0, Composite-1.%d, 2, -1.0, %s, 2, -1.0\n" % (ntop[0],nbot[0],'Set-dummy-TB')
for nright,nleft in zip(nodesets["right"],nodesets["left"]):
    equationstring += "3\n"
    equationstring += "Composite-1.%d, 1, 1.0, Composite-1.%d, 1, -1.0, %s, 1, -1.0\n" % (nright[0],nleft[0],'Set-dummy-LR')
    equationstring += "3\n"
    equationstring += "Composite-1.%d, 2, 1.0, Composite-1.%d, 2, -1.0, %s, 2, -1.0\n" % (nright[0],nleft[0],'Set-dummy-LR')
keywordblock.insert(position=endassembly,text=equationstring)
bcstring = "**\n*Boundary\n"
bcstring += "Set-dummy-TB,2,2, 0.01\n" 
bcstring += "**\n*Boundary\n"
bcstring += "Set-dummy-TB,1,1, 0.0\n"
bcstring += "**\n*Boundary\n"
bcstring += "Set-dummy-LR,2,2, 0.0\n" 
bcstring += "**\n*Boundary\n"
bcstring += "Set-dummy-LR,1,1, 0.0\n" 
keywordblock.insert(position=static+3,text=bcstring)


mdb.Job(name=job_name2, model='Model-1', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB,parallelizationMethodExplicit=DOMAIN, 
    numDomains=1,multiprocessingMode=DEFAULT, numCpus=1, 
    numGPUs=0)
historystring = "*Output, history\n*Node Output, Nset=Set-dummy-LR\n RF\n" 
historystring += "*Node Output, Nset=Set-dummy-TB\n RF\n"
keywordblock.replace(position=history+4,text=historystring)
mdb.jobs[job_name2].submit(consistencyChecking=OFF)
mdb.jobs[job_name2].waitForCompletion()





o3 = session.openOdb(name=path_odb1)
step = o3.steps["Step-1"]
horiRegion1 = step.historyRegions["Node DUMMY-LR-1.100000"]
vertRegion1 = step.historyRegions["Node DUMMY-TB-1.200000"]
fxx1 = horiRegion1.historyOutputs["RF1"].data[-1][1]
fxy1 = horiRegion1.historyOutputs["RF2"].data[-1][1]
fyx1 = vertRegion1.historyOutputs["RF1"].data[-1][1]
fyy1 = vertRegion1.historyOutputs["RF2"].data[-1][1]
exx=0.01/Lrve
Q11=fxx1/(Wrve*exx)
Q12=(fyy1)/(Lrve*exx)

with open("sw.txt","a") as outfile:
	outfile.write("----------------Results of 2ND Hierarchy----------------------\n")
	outfile.write("------Simulation_using_RF,plane strain elements-----\n")
	outfile.write("n=%d, aspect=%d, vf=%.2f, Ep/Gm=%.1f\n" % (n,rho,phi, EpGm))
	outfile.write("---------------Horizontal displacement-----------------------\n")
	outfile.write("Q11=%2f\n" %(Q11))
	outfile.write("Q12=%2f\n" %(Q12))
	outfile.write("--------------------------------------\n")


o4 = session.openOdb(name=path_odb2)
step = o4.steps["Step-1"]
horiRegion2 = step.historyRegions["Node DUMMY-LR-1.100000"]
vertRegion2 = step.historyRegions["Node DUMMY-TB-1.200000"]
fxx2 = horiRegion2.historyOutputs["RF1"].data[-1][1]
fxy2 = horiRegion2.historyOutputs["RF2"].data[-1][1]
fyx2 = vertRegion2.historyOutputs["RF1"].data[-1][1]
fyy2 = vertRegion2.historyOutputs["RF2"].data[-1][1]
eyy=0.01/Wrve
Q12_vs=fxx2/(Wrve*eyy)
Q22=fyy2/(Lrve*eyy)
deno=(Q11*Q22) - (Q12*Q12)
S11=Q22/deno
S12=-Q12/deno
S22=Q11/deno
E11=1/S11
E22=1/S22
nu12=-S12/S11


with open("sw.txt","a") as outfile:
    outfile.write("---------------Vertical displacement-----------------------\n")
    outfile.write("Q12_vs=%2f\n" %(Q12_vs))
    outfile.write("Q22=%2f\n" %(Q22))
    outfile.write("E11=%2f\n" %(E11))
    outfile.write("E22=%2f\n" %(E22))
    outfile.write("nu12=%2f\n" %(nu12))
    outfile.write("----------------End of simulation for 2nd Hierarchy----------------------\n")
    outfile.write("%######%############%##########%#############%\n")


# #######################################2H-SECOND HIERARCHY ENDS HERE###############







