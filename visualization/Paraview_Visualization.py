# state file generated using paraview version 5.10.0-RC1

# uncomment the following three lines to ensure this script works in future versions
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [965, 368]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [0.0012874962794600098, 0.010300671532124994, 0.6982609236440099]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [-6.2917245962670165, 1.950585388838533, 1.7271000754001522]
renderView1.CameraFocalPoint = [0.001287496279460042, 0.010300671532124961, 0.6982609236440099]
renderView1.CameraViewUp = [0.14259629996396966, -0.061335336233023036, 0.9878786725938402]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 1.7250870083692138
renderView1.CameraParallelProjection = 1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(965, 368)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Legacy VTK Reader'
equilibrium_90vtk = LegacyVTKReader(registrationName='equilibrium_90.vtk', FileNames=['/home/BU/dredwan1/membrane_dynamics_results/membrane_dynamics_data/Rp_0.3/forced_wrapping/outside/rho_001/chi90/equilibrium.vtk'])

# create a new 'Sphere'
sphere1 = Sphere(registrationName='Sphere1')
sphere1.Center = [0.0, 0.0, 1.28853]
sphere1.Radius = 0.3
sphere1.ThetaResolution = 20
sphere1.PhiResolution = 20

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from equilibrium_90vtk
equilibrium_90vtkDisplay = Show(equilibrium_90vtk, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
equilibrium_90vtkDisplay.Representation = 'Surface With Edges'
equilibrium_90vtkDisplay.AmbientColor = [0.611764705882353, 0.611764705882353, 0.611764705882353]
equilibrium_90vtkDisplay.ColorArrayName = [None, '']
equilibrium_90vtkDisplay.DiffuseColor = [0.611764705882353, 0.611764705882353, 0.611764705882353]
equilibrium_90vtkDisplay.Opacity = 0.3
equilibrium_90vtkDisplay.SelectTCoordArray = 'None'
equilibrium_90vtkDisplay.SelectNormalArray = 'None'
equilibrium_90vtkDisplay.SelectTangentArray = 'None'
equilibrium_90vtkDisplay.EdgeColor = [0.0, 0.0, 0.0]
equilibrium_90vtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
equilibrium_90vtkDisplay.SelectOrientationVectors = 'None'
equilibrium_90vtkDisplay.ScaleFactor = 0.208971747171132
equilibrium_90vtkDisplay.SelectScaleArray = 'None'
equilibrium_90vtkDisplay.GlyphType = 'Arrow'
equilibrium_90vtkDisplay.GlyphTableIndexArray = 'None'
equilibrium_90vtkDisplay.GaussianRadius = 0.0104485873585566
equilibrium_90vtkDisplay.SetScaleArray = [None, '']
equilibrium_90vtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
equilibrium_90vtkDisplay.OpacityArray = [None, '']
equilibrium_90vtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
equilibrium_90vtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
equilibrium_90vtkDisplay.PolarAxes = 'PolarAxesRepresentation'
equilibrium_90vtkDisplay.ScalarOpacityUnitDistance = 0.19492739506988482
equilibrium_90vtkDisplay.OpacityArrayName = [None, '']

# show data from sphere1
sphere1Display = Show(sphere1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
sphere1Display.Representation = 'Surface'
sphere1Display.AmbientColor = [0.6666666666666666, 0.3333333333333333, 0.0]
sphere1Display.ColorArrayName = [None, '']
sphere1Display.DiffuseColor = [0.6666666666666666, 0.3333333333333333, 0.0]
sphere1Display.Specular = 0.8
sphere1Display.SelectTCoordArray = 'None'
sphere1Display.SelectNormalArray = 'Normals'
sphere1Display.SelectTangentArray = 'None'
sphere1Display.OSPRayScaleArray = 'Normals'
sphere1Display.OSPRayScaleFunction = 'PiecewiseFunction'
sphere1Display.SelectOrientationVectors = 'None'
sphere1Display.ScaleFactor = 0.059999996423721315
sphere1Display.SelectScaleArray = 'None'
sphere1Display.GlyphType = 'Arrow'
sphere1Display.GlyphTableIndexArray = 'None'
sphere1Display.GaussianRadius = 0.002999999821186066
sphere1Display.SetScaleArray = ['POINTS', 'Normals']
sphere1Display.ScaleTransferFunction = 'PiecewiseFunction'
sphere1Display.OpacityArray = ['POINTS', 'Normals']
sphere1Display.OpacityTransferFunction = 'PiecewiseFunction'
sphere1Display.DataAxesGrid = 'GridAxesRepresentation'
sphere1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
sphere1Display.ScaleTransferFunction.Points = [-0.9965844750404358, 0.0, 0.5, 0.0, 0.9965844750404358, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
sphere1Display.OpacityTransferFunction.Points = [-0.9965844750404358, 0.0, 0.5, 0.0, 0.9965844750404358, 1.0, 0.5, 0.0]

# ----------------------------------------------------------------
# restore active source
SetActiveSource(sphere1)
# ----------------------------------------------------------------


if __name__ == '__main__':
    # generate extracts
    SaveExtracts(ExtractsOutputDirectory='extracts')
