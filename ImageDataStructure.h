/* 
 * File:   ImageDataStructure.h
 * Author: ashish
 *
 * Created on 30 October, 2011, 1:55 AM
 */

#ifndef IMAGEDATASTRUCTURE_H
#define	IMAGEDATASTRUCTURE_H

#define MAX_TRIANGLES 1500
#define MAX_SPHERES 10
#define MAX_LIGHTS 10

#define WHITE_LIGHT 255
#define NONE -5
#define SPHERE -15
#define TRIANGLE -20

// Edit this macro to specify the depth of recursive ray tracing
#define MAX_DEPTH 3

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

typedef struct _Triangle
{
  struct Vertex v[3];
} Triangle;

typedef struct _Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
} Sphere;

typedef struct _Light
{
  double position[3];
  double color[3];
} Light;

extern Triangle triangles[MAX_TRIANGLES];
extern Sphere spheres[MAX_SPHERES];
extern Light lights[MAX_LIGHTS];
extern double ambient_light[3];

#endif	/* IMAGEDATASTRUCTURE_H */
