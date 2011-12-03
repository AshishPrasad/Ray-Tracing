################################################################################
DATA STRUCTURES USED:
################################################################################
class Ray
--------------------------------------------------------------------------------
-Defined in "Ray.h"
-Contains information about the source, direction, 
 nearest intersection distance,intersected object and the color it assumes.
--------------------------------------------------------------------------------
class Vector4 and class Matrix4
--------------------------------------------------------------------------------
-Defined in Vector4.h and Matrix4.h
-Supports various vector and matrix operations.

################################################################################
ASSUMPTIONS:
################################################################################
-Camera position is [0,0,0].
-Screen Position along z-axis: -1.000000
-Ambient Light Coefficient = Diffused Light Coefficient.
-The depth of  the recursive mode for ray tracing is 3 by default.
 One can set the depth of the recursive ray tracing by changing the default 
 value of "#define MAX_DEPTH 3" in the file ImageDataStructure.h
################################################################################
