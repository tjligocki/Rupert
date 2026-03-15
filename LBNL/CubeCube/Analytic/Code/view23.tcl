#
# This example creates a polygonal model of a mace made of a sphere
# and a set of cones adjusted on its surface using glyphing. 
#
# The sphere is rendered to the screen through the usual VTK render window
# and interactions is performed using vtkRenderWindowInteractor.
# The basic setup of source -> mapper -> actor -> renderer ->
# renderwindow is typical of most VTK programs.  
#

#
# First we include the VTK Tcl packages which will make available 
# all of the vtk commands to Tcl
#
package require vtk
package require vtkinteraction

#
# Next we create an instance of vtkSphereSource and set some of its 
# properties
#
vtkOutlineSource cube

vtkTubeFilter cubeTube
    cubeTube SetInput [cube GetOutput]
    cubeTube SetRadius 0.05
    cubeTube SetNumberOfSides 12
#    cubeTube CappingOn

#
# We create an instance of vtkPolyDataMapper to map the polygonal data 
# into graphics primitives. We connect the output of the cube source
# to the input of this mapper 
#
vtkPolyDataMapper cubeMapper
    cubeMapper SetInput [cubeTube GetOutput]

#
# Create an actor to represent the cube. The actor coordinates rendering of
# the graphics primitives for a mapper. We set this actor's mapper to be
# the mapper which we created above.
#
vtkActor cubeActor
    cubeActor SetMapper cubeMapper

vtkSphereSource sphere
    sphere SetThetaResolution 16 
    sphere SetPhiResolution 16
    sphere SetRadius [expr 1.5 * [cubeTube GetRadius]]

vtkGlyph3D cubeCorners
    cubeCorners SetInput [cube GetOutput]
    cubeCorners SetSource [sphere GetOutput]

vtkPolyDataMapper cubeCornersMapper
    cubeCornersMapper SetInput [cubeCorners GetOutput]

vtkActor cubeCornersActor
    cubeCornersActor SetMapper cubeCornersMapper
    [cubeCornersActor GetProperty] SetColor 0 0 1

vtkCubeSource cubeFaces
    cubeFaces SetXLength 2
    cubeFaces SetYLength 2
    cubeFaces SetZLength 2

vtkPolyDataMapper cubeFacesMapper
    cubeFacesMapper SetInput [cubeFaces GetOutput]

vtkActor cubeFacesActor
    cubeFacesActor SetMapper cubeFacesMapper
    [cubeFacesActor GetProperty] SetOpacity 0.3

vtkPlaneSource square
    square SetOrigin -0.644405 -0.775910 -1.000000
    square SetPoint1 1.000000 -1.000000 0.131506
    square SetPoint2 -1.000000 1.000000 -0.131506

vtkPolyDataMapper squareMapper
    squareMapper SetInput [square GetOutput]

vtkActor squareActor
    squareActor SetMapper squareMapper
    [squareActor GetProperty] SetColor 0 1 0

vtkGlyph3D squareCorners
    squareCorners SetInput [square GetOutput]
    squareCorners SetSource [sphere GetOutput]

vtkPolyDataMapper squareCornersMapper
    squareCornersMapper SetInput [squareCorners GetOutput]

vtkActor squareCornersActor
    squareCornersActor SetMapper squareCornersMapper
    [squareCornersActor GetProperty] SetColor 1 1 0

#
# Create the Renderer and assign actors to it. A renderer is like a
# viewport. It is part or all of a window on the screen and it is responsible
# for drawing the actors it has. We also set the background color here.
#
vtkRenderer renderer
    renderer AddActor cubeActor
    renderer AddActor cubeCornersActor
    renderer AddActor cubeFacesActor
    renderer AddActor squareActor
    renderer AddActor squareCornersActor
#    renderer SetBackground 1 1 1

[renderer GetActiveCamera] ParallelProjectionOn

#
# We create the render window which will show up on the screen
# We put our renderer into the render window using AddRenderer. We also
# set the size to be 600 pixels by 600
#
vtkRenderWindow renWin
    renWin AddRenderer renderer
    renWin SetSize 600 600

#
# Finally we create the render window interactor handling user
# interactions. vtkRenderWindowInteractor provides a
# platform-independent interaction mechanism for mouse/key/time
# events. vtkRenderWindowInteractor also provides controls for
# picking, rendering frame rate, and headlights. It is associated
# to a render window.
#
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

#
# vtkRenderWindowInteractor provides default key bindings.  The 'u'
# key will trigger its "user method", provided that it has been
# defined. Similarly the 'e' or 'q' key will trigger its "exit
# method". The lines below set these methods through the AddObserver
# method with the events "UserEvent" and "ExitEvent". The corresponding 
# "user-method" Tcl code will bring up the .vtkInteract widget and 
# allow the user to evaluate any Tcl code and get access to all 
# previously-created VTK objects. The
# "exit-method" Tcl code will exit (do not try to free up any objects
# we created using 'vtkCommand DeleteAllObjects' because you are right
# inside a VTK object.
#
iren AddObserver UserEvent {wm deiconify .vtkInteract}
iren AddObserver ExitEvent {exit}

#
# Render the image
#
renWin Render 

# 
# Hide the default . widget
#
wm withdraw .

#
# You only need this line if you run this script from a Tcl shell
# (tclsh) instead of a Tk shell (wish) 
#
tkwait window .
