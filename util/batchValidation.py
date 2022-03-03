import os
from paraview.simple import *


# Find path for cases
curr_dir_path = os.path.dirname(os.path.realpath(__file__))
cases = os.listdir(curr_dir_path + '/Cases')
pop = cases.index('baseCase')
cases.pop(pop)

# Initialise and run solver
for case_parent in cases:
    for case_child in os.listdir(curr_dir_path + '/Cases/' + case_parent):
        curr_case = curr_dir_path + '/Cases/' + case_parent + '/' + case_child
        # print(curr_case)
        openfoam = OpenFOAMReader(
            registrationName='open.foam',
            FileName = curr_case + '/open.foam')

        # # create a new 'OpenFOAMReader'
        # openfoam.MeshRegions = [
        #         'internalMesh', 'lagrangian/dsmc']

        # openfoam.CellArrays = [
        #         'Ma_Ar', 'SOFP_Ar', 'Ttra_Ar',
        #         'U_Ar', 'cellLevel', 'dsmcNMean_Ar',
        #         'dsmcN_Ar', 'dsmcSigmaTcRMax', 'fD_Ar',
        #         'mctToDt_Ar', 'mct_Ar', 'mfpToDx_Ar',
        #         'mfp_Ar', 'p_Ar', 'rhoM_Ar', 'rhoN_Ar',
        #         'wallHeatFlux_Ar', 'wallShearStress_Ar']

        # openfoam.LagrangianArrays = ['U', 'classification',
        #         'newParcel', 'origId', 'origProcId', 'typeId']

        # # get animation scene
        # animationScene1 = GetAnimationScene()

        # # update animation scene based on data timesteps
        # animationScene1.UpdateAnimationUsingDataTimeSteps()

        # # set active source
        # SetActiveSource(openfoam)

        # # Properties modified on openfoam
        # openfoam.MeshRegions = ['internalMesh']

        # UpdatePipeline(time=0.2, proxy=openfoam)

        # # create a new 'Plot Over Line'
        # plotOverLine1 = PlotOverLine(
        #         registrationName='PlotOverLine1',
        #         Input=openfoam)
        # plotOverLine1.Point1 = [0.9975, 0.0, 0.0]
        # plotOverLine1.Point2 = [20.0, 1.398, 0.0]

        # UpdatePipeline(time=0.2, proxy=plotOverLine1)
        #         # kind

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
        plotOverLine1.Point1 = [1.0, 0.0, 0.0]
        plotOverLine1.Point2 = [20.0, 0.748, 0.0]


        # lineChartView1 = GetActiveViewOrCreate('XYChartView')
        # plotOverLine1Display = GetDisplayProperties(plotOverLine1, view=lineChartView1)
        # plotOverLine1Display.SeriesVisibility = ['cellLevel', 'dsmcN_Ar', 'dsmcNMean_Ar', 'dsmcSigmaTcRMax', 'fD_Ar_Magnitude', 'Ma_Ar', 'mct_Ar', 'mctToDt_Ar', 'mfp_Ar', 'mfpToDx_Ar', 'p_Ar', 'Points_Magnitude', 'rhoM_Ar', 'rhoN_Ar', 'SOFP_Ar', 'Ttra_Ar', 'U_Ar_Magnitude', 'wallHeatFlux_Ar', 'wallShearStress_Ar']


        UpdatePipeline(time=0.2, proxy=plotOverLine1)

        UpdatePipeline(time=0.2, proxy=plotOverLine1)


        spreadSheetView1 = CreateView('SpreadSheetView')
        spreadSheetView1.ColumnToSort = ''
        spreadSheetView1.BlockSize = 1024

        # get active source.
        plotOverLine1 = GetActiveSource()

        # show data in view
        plotOverLine1Display = Show(plotOverLine1, spreadSheetView1, 'SpreadSheetRepresentation')

        # assign view to a particular cell in the layout
        layout1 = GetLayoutByName("Layout #1")

        AssignViewToLayout(view=spreadSheetView1, layout=layout1, hint=6)

        # Properties modified on plotOverLine1Display
        plotOverLine1Display.Assembly = ''

        # resize frame
        layout1.SetSplitFraction(2, 0.14847161572052403)

        # export view
        ExportView(curr_case + '/results.csv', view=spreadSheetView1)
        # SaveData(curr_case + '/results.csv', plotOverLine1)