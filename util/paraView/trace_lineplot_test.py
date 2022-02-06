# trace generated using paraview version 5.10.0
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'OpenFOAMReader'
openfoam = OpenFOAMReader(registrationName='open.foam', FileName='/home/andy/OpenFOAM/andy-v1706/run/rarefiedPlume/cases/1d/open.foam')
openfoam.MeshRegions = ['internalMesh', 'lagrangian/dsmc']
openfoam.CellArrays = ['Ma_Ar', 'SOFP_Ar', 'Ttra_Ar', 'U_Ar', 'cellLevel', 'dsmcNMean_Ar', 'dsmcN_Ar', 'dsmcSigmaTcRMax', 'fD_Ar', 'mctToDt_Ar', 'mct_Ar', 'mfpToDx_Ar', 'mfp_Ar', 'p_Ar', 'rhoM_Ar', 'rhoN_Ar', 'wallHeatFlux_Ar', 'wallShearStress_Ar']
openfoam.LagrangianArrays = ['U', 'classification', 'newParcel', 'origId', 'origProcId', 'typeId']

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# set active source
SetActiveSource(openfoam)

UpdatePipeline(time=0.01, proxy=openfoam)

animationScene1.GoToLast()

# set active source
SetActiveSource(openfoam)

# Properties modified on openfoam
openfoam.MeshRegions = ['internalMesh']

UpdatePipeline(time=0.2, proxy=openfoam)

# create a new 'Plot Over Line'
plotOverLine1 = PlotOverLine(registrationName='PlotOverLine1', Input=openfoam)
plotOverLine1.Point1 = [0.9975510239601135, 0.0, 0.0]
plotOverLine1.Point2 = [20.0, 1.3988569974899292, 0.0978228747844696]

UpdatePipeline(time=0.2, proxy=plotOverLine1)

UpdatePipeline(time=0.2, proxy=plotOverLine1)

# Maybe this will work? 
# kind of skeptical bc points are added after????
# SaveData('test2.csv', plotOverLine1)