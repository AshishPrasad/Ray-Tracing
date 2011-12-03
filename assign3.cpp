//============================================================================
// Name        : Assignment3.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

/*
Author: Frank Pfenning
*/

#include <cassert>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*********************Files to be Included*************************************/
#include "ImageDataStructure.h"
#include "pic.h"
#include "Ray.h"
#include "Vector4.h"
/******************************************************************************/

using namespace std;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

char *filename=0;
int mode = MODE_DISPLAY;

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles;
int num_spheres;
int num_lights;

unsigned char buffer[HEIGHT][WIDTH][3];

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

// The field of view of the camera
#define fov 60.0

/***************************Additional Parameters******************************/

// The place where camera is placed
const double cameraZ = 0.0;

// Screen plane
const double screenZ = -1.0;
const double screenY = tan((fov/2)*(M_PI/180))*(cameraZ-screenZ);
const double screenX = screenY*((double)WIDTH/(double)HEIGHT);

//// For Debugging Purposes
//int intersectTriangle[3] = {0,0,0};

/***************************Solving Equations**********************************/
// (A)x^2 + (B)x + C = 0
void SolveQuadraticEquation(double A, double B, double C, double* soln1, double* soln2){
    A=-A;B=-B;C=-C;
    *soln1 = NONE; *soln2 = NONE;
    double discriminantSquare = pow(B,2) - 4*A*C;
    
    if(discriminantSquare >= 0){
        if (A == 0){
            *soln1 = -C/B;
            *soln2 = -C/B;
        }        
        else {
            double discriminant = sqrt(discriminantSquare);
            *soln1 = (-1*B - discriminant)/(2*A);
            *soln2 = (-1*B + discriminant)/(2*A);
        }
    }
}


// A1x+B1y = C1
// A2x+B2y = C2
void SolveLinearEquation(double A1, double B1, double C1, double A2, double B2, double C2,
        double* x, double* y){
    
    C1=-C1; C2=-C2;
    
    if (A1*B2!=B1*A2){
        *x = (B1*C2 - C1*B2)/(A1*B2 - B1*A2);
        *y = (C1*A2 - A1*C2)/(A1*B2 - B1*A2);
    }
    else {
        *x = NONE;
        *y = NONE;
    }
}

/*******************Interpolating Triangle Properties**************************/

// Obtain intersection point ratios: alpha and beta
// V' = V0 + alpha(V0P) = V1 + beta(V1V2) on the line V1V2 
void ObtainIntermediateTriangleIntersection(int index, Vector4* intersectPoint, 
        double* alpha, double* beta){
    
    // (V0P)alpha -(V1V2)beta = V0V1
    Vector4* V0 = new Vector4(triangles[index].v[0].position);
    Vector4* V1 = new Vector4(triangles[index].v[1].position);
    Vector4* V2 = new Vector4(triangles[index].v[2].position);

    Vector4* V0P = new Vector4(intersectPoint);
    V0P->Subtraction(V0);

    Vector4* V1V2 = new Vector4(V2);
    V1V2->Subtraction(V1);

    Vector4* V0V1 = new Vector4(V1);
    V0V1->Subtraction(V0);

    double A1 = V0P->vector[0];
    double B1 = -1 * V1V2->vector[0];
    double C1 = V0V1->vector[0];

    double A2 = V0P->vector[1];
    double B2 = -1 * V1V2->vector[1];
    double C2 = V0V1->vector[1];
    
    double A3 = V0P->vector[2];
    double B3 = -1 * V1V2->vector[2];
    double C3 = V0V1->vector[2];

    SolveLinearEquation(A1, B1, C1, A2, B2, C2, alpha, beta);

    if (*alpha == NONE) {
        SolveLinearEquation(A2, B2, C2, A3, B3, C3, alpha, beta);
    }
    
    if (*alpha == NONE) {
        SolveLinearEquation(A1, B1, C1, A3, B3, C3, alpha, beta);
    }
    
//    // For Debugging Purposes:
//    printf("\n (%f,%f)", *alpha, *beta);
//    // intermediatePoint = V0 + alphaV0P
//    Vector4* intermediatePoint = new Vector4(V0P->vector);
//    intermediatePoint->Scale(*alpha);
//    intermediatePoint->Addition(V0);
//    delete(intermediatePoint);
    
    assert((*alpha>=1 || *alpha==0)&& *beta<=1 && *beta>=0);
    
    delete(V0V1);
    delete(V1V2);
    delete(V0P);
    delete(V2);
    delete(V1);
    delete(V0);
}

void InterpolateTriangleProperty(double V0Property[3], double V1Property[3], double V2Property[3], 
        double alpha, double beta, Vector4* property, bool normalize = false){
    
    assert((alpha>=1 || alpha==0)&& beta<=1 && beta>=0);
    
    if (alpha!=0){    
        
        // _V = betaV2.....(1-beta)V1        
        Vector4* _VPropertyHelper = new Vector4(V2Property);
        _VPropertyHelper->Scale(beta);

        Vector4* _VProperty = new Vector4(V1Property);
        _VProperty->Scale(1 - beta);
        _VProperty->Addition(_VPropertyHelper);
            delete(_VPropertyHelper);
        
        if (normalize) {
            _VProperty->ConvertToUnitVector();
        }    
            
        // property = ((alpha-1)*V0Property + _VProperty)/alpha    
        property->ReInitialize(V0Property);
        property->Scale(alpha-1);
        property->Addition(_VProperty);
        property->Scale(1/alpha);    
        
        if (normalize) {
            property->ConvertToUnitVector();
        }
        
        delete(_VProperty);
    }
    else{
        property->ReInitialize(V0Property);
        
        if (normalize){
            property->ConvertToUnitVector();
        }
    }
}

/**********************Ray Intersection Functions*******************************/

double FindNearestTIntersectionWithSphere(Ray* ray, int* index, int noCheckID = NONE) {
    
    assert(noCheckID==NONE || (noCheckID>=0 && noCheckID<num_spheres));
    
    double nearestTIntersect = INFINITY;
    *index = NONE;
    
    /* (x_d^2 + y_d^2 + z_d^2)t^2 + 2((x_o-x_c)x_d +(y_o-y_c)y_d + (z_o-z_c)z_d)t
     *  + (x_o-x_c)^2 + (y_o-y_c)^2 + (z_o-z_c)^2 - r^2 = 0
     */
    unsigned int i;
    for (i = 0; i < num_spheres; i++) {
        
        if (i==noCheckID){
            continue;
        }
        
        Vector4* sphereCenter = new Vector4(spheres[i].position);
        
        Vector4* sphereRay = new Vector4(ray->source);
        sphereRay->Subtraction(sphereCenter);
        
        double A = pow(ray->dir->Mod(), 2);
        double B = 2 * (ray->dir->DotProduct(sphereRay));
        double C = pow(sphereRay->Mod(), 2) - pow(spheres[i].radius, 2);
        
        double soln1=0.0,soln2=0.0;
        SolveQuadraticEquation(A,B,C,&soln1,&soln2);       
        
        if (soln1 >= 0 && soln2 >= 0){
//            //For Debugging Purposes
//            printf("(%f,%f)\n",soln1,soln2);
            
            double candidateT = NONE;
            if(soln1<soln2){
                candidateT = soln1;
            }
            else{
                candidateT = soln2;
            }
                
            //assert(candidateT>=soln1 && candidateT>=soln2);
            if(nearestTIntersect > candidateT /*&& candidateT>(cameraZ-screenZ)*/) {
                nearestTIntersect = candidateT;            
                *index = i;           
            }
        }
        
        delete(sphereCenter);
        delete(sphereRay);
    }       
    
    return nearestTIntersect;
}

double FindNearestTIntersectionWithTriangle(Ray* ray, int* index, int noCheckID = NONE){
    
    assert(noCheckID==NONE || (noCheckID>=0 && noCheckID<num_triangles));
    
    double nearestTIntersect = INFINITY;
    *index = NONE;
    
    unsigned int i;    
    for (i = 0; i < num_triangles; i++) {
        
        if (i==noCheckID){
            continue;
        }
        // Vertices V0, V1 and V2
        Vector4* V0 = new Vector4(triangles[i].v[0].position);
        Vector4* V1 = new Vector4(triangles[i].v[1].position);
        Vector4* V2 = new Vector4(triangles[i].v[2].position);
        
        // Edges: VOV1, V0V2
        Vector4* V0V1 = new Vector4(V1);
        V0V1->Subtraction(V0);
        Vector4* V0V2 = new Vector4(V2);
        V0V2->Subtraction(V0);
        
        // Normal: N = V0V1 X V0V2
        Vector4* N = new Vector4(V0V1);
        N->CrossProductPost(V0V2);
        N->ConvertToUnitVector();
                
        // ray and normal should be opposite
        if (N->DotProduct(ray->dir)>0){
            N->Scale(-1);
        }
        
        // N = [A B C]
        // AX0 + BY0 + CZ0 + D = 0
        double D = -1*N->DotProduct(V0);                       
        
        if (N->DotProduct(ray->dir) != (double) 0) {
            // t = -(N.Ro + D)/ (N.Rd)
            double candidateT = -1 * (N->DotProduct(ray->source) + D) / N->DotProduct(ray->dir);            
            
            if (candidateT > 0){
                
                // P = Ro + tRd
                Vector4* point = new Vector4(ray->dir);
                point->Scale(candidateT);
                point->Addition(ray->source);

                // alpha*V0V1 + beta*V0V2 = V0P
                Vector4* V0P = new Vector4(point);
                V0P->Subtraction(V0);

                double alpha = 0.0;
                double beta = 0.0;

                //(x01)alpha + (x02)beta = x0p;
                //(y01)alpha + (y02)beta = y0p;
                //(z01)alpha + (z02)beta = z0p;
                double A1 = V0V1->vector[0];
                double B1 = V0V2->vector[0];
                double C1 = V0P->vector[0];

                double A2 = V0V1->vector[1];
                double B2 = V0V2->vector[1];
                double C2 = V0P->vector[1];
                
                double A3 = V0V1->vector[2];
                double B3 = V0V2->vector[2];
                double C3 = V0P->vector[2];

                SolveLinearEquation(A1, B1, C1, A2, B2, C2, &alpha, &beta);

                if (alpha == NONE) {
                    SolveLinearEquation(A2, B2, C2, A3, B3, C3, &alpha, &beta);
                }
                
                if (alpha == NONE) {
                    SolveLinearEquation(A1, B1, C1, A3, B3, C3, &alpha, &beta);
                }
                
//                // For Debugging Purposes
//                printf("\nTriangle: %d ==> CandidateT: %f ==> (%f,%f)",i,candidateT,alpha,beta);
                
                if (alpha > 0 && beta > 0 && (alpha + beta) < 1) {
                    if (nearestTIntersect > candidateT) {
                        nearestTIntersect = candidateT;
                        *index = i;
                    }
                }

                delete(V0P);
                delete(point); 
            }
        }
        
        delete(N);
        delete(V0V2);
        delete(V0V1);
        delete(V2);
        delete(V1);
        delete(V0);
    }
    
    return nearestTIntersect;
}

double FindNearestTIntersection(Ray* ray, int noCheckObject = NONE, int noCheckID = NONE){
    
    assert((noCheckObject==NONE && noCheckID==NONE) 
            || (noCheckObject==SPHERE && noCheckID>=0 && noCheckID<num_spheres)
            || (noCheckObject==TRIANGLE && noCheckID>=0 && noCheckID<num_triangles));
            
    
    int dontCheckSphereID = NONE;
    int dontCheckTriangleID = NONE;
    
    if(noCheckObject==SPHERE){
        dontCheckSphereID = noCheckID;
    }
    if(noCheckObject==TRIANGLE){
        dontCheckTriangleID = noCheckID;
    }
    
    int index1 = NONE;
    double nearestSphereIntersect = FindNearestTIntersectionWithSphere(ray,&index1, dontCheckSphereID);    
        
    int index2 = NONE;    
    double nearestTriangleIntersect = FindNearestTIntersectionWithTriangle(ray,&index2, dontCheckTriangleID);    
    
//    // For Debugging Purposes
//    printf("[%f,%f]\t",nearestSphereIntersect,nearestTriangleIntersect);
    
    int index = NONE;
    double nearestTIntersect = INFINITY;
    
    if (nearestSphereIntersect < nearestTriangleIntersect){
        nearestTIntersect = nearestSphereIntersect;
        index = index1;
    }
    else{
        nearestTIntersect = nearestTriangleIntersect;
        index = index2;
        
//        // For Debugging Purposes
//        intersectTriangle[index] = 1;
    }
    
    if (index != NONE) {

        ray->intersectID = index;
        ray->tIntersect = nearestTIntersect;

        // Obtaining Intersection Point = P0 + t*dir
        ray->intersectPoint->ReInitialize(ray->dir->vector);
        ray->intersectPoint->Scale(nearestTIntersect);
        ray->intersectPoint->Addition(ray->source);
        
        ray->color->ReInitialize(ambient_light);             
                
        if (index == index1){
            
            ray->intersectedObject = SPHERE;            

            // Obtaining normal at Intersection Point = intersectPoint - center
            ray->normal->ReInitialize(spheres[index].position);
            ray->normal->Scale(-1);
            ray->normal->Addition(ray->intersectPoint);
            ray->normal->ConvertToUnitVector();

            // Assuming ambient coefficient = diffused coefficient, add ambient light
            ray->color->Scale(spheres[index].color_diffuse[0], spheres[index].color_diffuse[1],
                    spheres[index].color_diffuse[2]);                
        }
        else {
            
            ray->intersectedObject = TRIANGLE;
            
            // Obtaining normal at Intersection Point
            double alpha = 0.0, beta = 0.0;
            ObtainIntermediateTriangleIntersection(index, ray->intersectPoint, &alpha, &beta);
            //printf("\n(%f,%f):%f\n",alpha,beta,ray->tIntersect);
            
            InterpolateTriangleProperty(triangles[index].v[0].normal, triangles[index].v[1].normal,
                    triangles[index].v[2].normal, alpha, beta, ray->normal, true);

            // Assuming ambient coefficient = diffused coefficient, add ambient light
            Vector4* diffusedCoeff = new Vector4();
            InterpolateTriangleProperty(triangles[index].v[0].color_diffuse, triangles[index].v[1].color_diffuse, triangles[index].v[2].color_diffuse,
                    alpha, beta, diffusedCoeff);
                      
            assert(diffusedCoeff->HasNonNegativeEntries());
            ray->color->Scale(diffusedCoeff->vector[0], diffusedCoeff->vector[1], diffusedCoeff->vector[2]);
            delete(diffusedCoeff);
        }
    }    
    
    return nearestTIntersect;
}

double FindDistanceOfIntersectedPointWithLight(Ray* ray, int index){
    
    Vector4* P2L = new Vector4(lights[index].position);
    P2L->Subtraction(ray->intersectPoint);
    double P2Ldist =  P2L->Mod();
    delete(P2L);
    return P2Ldist;
}


/***********************Tracing Functions**************************************/

// Forward Declaration
void RT_trace(Ray*,int,int,int);

// For Shading Purposes
void RT_shade(Ray* ray, int depth) {
    
    if (depth <= MAX_DEPTH){                
        
        // Color Coefficients
        Vector4* diffuseCoeff = new Vector4();
        Vector4* specularCoeff = new Vector4();
        double shininess;
        
        // Diffuse and Specular Terms
        double diffuseTerm[3] = {0.0, 0.0, 0.0};
        double specularTerm[3] = {0.0, 0.0, 0.0};
        
        if (ray->intersectedObject == SPHERE) {

            diffuseCoeff->ReInitialize(spheres[ray->intersectID].color_diffuse);
            specularCoeff->ReInitialize(spheres[ray->intersectID].color_specular);
            shininess = spheres[ray->intersectID].shininess;
        }
        else if (ray->intersectedObject == TRIANGLE) {

            // Interpolate diffuse and specular coefficients
            double alpha, beta;
            ObtainIntermediateTriangleIntersection(ray->intersectID, ray->intersectPoint, &alpha, &beta);

            diffuseCoeff->ReInitialize();
            InterpolateTriangleProperty(triangles[ray->intersectID].v[0].color_diffuse,
                    triangles[ray->intersectID].v[1].color_diffuse,
                    triangles[ray->intersectID].v[2].color_diffuse,
                    alpha, beta, diffuseCoeff, false);

            specularCoeff->ReInitialize();
            InterpolateTriangleProperty(triangles[ray->intersectID].v[0].color_specular,
                    triangles[ray->intersectID].v[1].color_specular,
                    triangles[ray->intersectID].v[2].color_specular,
                    alpha, beta, specularCoeff, false);

            // Interpolate shininess
            if (alpha > 0) {

                double intermediateShininess = beta * triangles[ray->intersectID].v[2].shininess
                        + (1 - beta) * triangles[ray->intersectID].v[1].shininess;

                shininess = ((alpha - 1) * triangles[ray->intersectID].v[0].shininess
                        + intermediateShininess) / alpha;
            } 
            else {
                shininess = triangles[ray->intersectID].v[0].shininess;
            }
        }

        // Boundary Checking        
        assert(diffuseCoeff->HasNonNegativeEntries());
        assert(!diffuseCoeff->HasGreaterThanOneEntries());
        assert(specularCoeff->HasNonNegativeEntries());
        assert(!specularCoeff->HasGreaterThanOneEntries());
        assert(shininess >= 0);
        
        // Normal Vector: N
        Vector4* N = new Vector4(ray->normal);
        N->ConvertToUnitVector();
        
        // Viewer Direction: V = ray->intersectedPoint - ray->source = -ray->dir
        Vector4* V = new Vector4(ray->dir);
        V->Scale(-1);
        V->ConvertToUnitVector();
        
        int dontCheckObject = NONE;
        int dontCheckObjectID = NONE;
        
        if (ray->intersectedObject == TRIANGLE) {
            dontCheckObject = ray->intersectedObject;
            dontCheckObjectID = ray->intersectID;
        }
        
                
        unsigned int i;
        for (i = 0; i < num_lights; i++) {            
            
            // Shoot a ray in the direction of the light
            Vector4* lightPos = new Vector4(lights[i].position);
            Ray* checkShadow = new Ray(ray->intersectPoint, lightPos);
            bool shadow = false;            
            
            double nearestTIntersect = FindNearestTIntersection(checkShadow,dontCheckObject,dontCheckObjectID);
            double point2LightDistance = FindDistanceOfIntersectedPointWithLight(ray,i);
            if (nearestTIntersect != INFINITY 
                && nearestTIntersect< point2LightDistance)
            {
                shadow = true;    
            }

            // Point to Light Vector: L = light_pos - intersectedPoint
            Vector4* L = new Vector4(lights[i].position);
            L->Subtraction(ray->intersectPoint);
            L->ConvertToUnitVector();            

            // Reflected Ray: R = 2(N.L)N - L
            Vector4* R = new Vector4(ray->normal);
            R->Scale(2 * N->DotProduct(L));
            R->Subtraction(L);
            R->ConvertToUnitVector();

            // Dot Product of N and L: NL
            double NL = N->DotProduct(L);
            if (NL < 0) {
                NL = 0;
            }

            // Dot Product of R and V: RV
            double RV = R->DotProduct(V);
            if (RV < 0) {
                RV = 0;
            }

            assert(NL <= 1 && RV <= 1);

            double attenuationTerm =  depth==1?1:1; 
            
            if (!shadow) {                
                unsigned int j;
                for (j = 0; j < 3; j++) {
                    diffuseTerm[j] += attenuationTerm * lights[i].color[j] * diffuseCoeff->vector[j] * NL;
                    specularTerm[j] += attenuationTerm * lights[i].color[j] * specularCoeff->vector[j] * pow(RV, shininess);
                }
            }

            
            
            delete(L);
            delete(R);

            delete(lightPos);
            delete(checkShadow);
        }

        // Adding diffused light
        Vector4* diffusedColor = new Vector4(diffuseTerm);
        //          //For Debugging Purposes
        //          if(ray->intersectedObject == SPHERE){                
        //            ray->color->Display();
        //            diffusedColor->Display();
        //            cout<<"................................................."<<endl;
        //          }
        ray->color->Addition(diffusedColor);
        delete(diffusedColor);

        // Adding specular light
        Vector4* specularColor = new Vector4(specularTerm);
        ray->color->Addition(specularColor);        
        delete(specularColor);

        //R = 2*(NI)N - I ; N = N; I = V = ray->source - ray->IntersectPoint
        Vector4* R_ = new Vector4(ray->normal);
        Vector4* I = new Vector4(ray->dir);
        I->Scale(-1);
        R_->Scale(2*N->DotProduct(I));
        R_->Subtraction(I);
        Vector4* end = new Vector4(ray->intersectPoint);
        end->Addition(R_);       
        
        
        Ray* reflectedRay = new Ray(ray->intersectPoint, end);
        delete(R_);
        delete(I);
        delete(end);
        assert(dontCheckObject!=SPHERE);
        RT_trace(reflectedRay, depth+1, dontCheckObject, dontCheckObjectID);

        //For Debugging
//        if(ray->intersectedObject==TRIANGLE/* && ray->intersectPoint->vector[1]<-2*screenX/4
//                && ray->intersectPoint->vector[1]>-3*screenX/4*/){
//            cout<<ray->intersectID;            
//            cout<<"---->";
//            ray->normal->Display();
//            cout<<"------>";
//            I->Display();
//            cout<<I->DotProduct(ray->normal)<<"\t";
//            cout<<reflectedRay->dir->DotProduct(ray->normal);
//            reflectedRay->dir->Display();
//            cout<<": "<<reflectedRay->tIntersect<<"\t";            
//            reflectedRay->color->Display();
//            cout<<endl;
//        }
        
        if(reflectedRay->tIntersect!=INFINITY){          
            
            reflectedRay->color->Scale(specularCoeff->vector[0], specularCoeff->vector[1], specularCoeff->vector[2]);                       
            ray->color->Addition(reflectedRay->color);
        }
        
        delete(reflectedRay);

        delete(N);
        delete(V);
        delete(specularCoeff);
        delete(diffuseCoeff);
    }    
}

// Tracing ray light
void RT_trace(Ray* ray, int depth,int dontCheckObject=NONE,int dontCheckObjectID=NONE){
    ray->color->ReInitialize(ambient_light);   
    
    if (FindNearestTIntersection(ray,dontCheckObject,dontCheckObjectID)!= (double)INFINITY){        
        RT_shade(ray, depth);       
    }        
}
/************************Draw Scene********************************************/

//MODIFY THIS FUNCTION
void draw_scene()
{ 
    // Mapping screen coordinates
    double x2i = screenX/((double)WIDTH/2.0);
    double y2j = screenY/((double)HEIGHT/2.0);
    
    printf("\n Camera Position along z-axis: %f",cameraZ);
    printf("\n Screen Position along z-axis: %f",screenZ);
    printf("\n Number of Triangles: %d",num_triangles);
    printf("\n Number of spheres: %d",num_spheres);
    printf("\n Number of lights: %d\n\n",num_lights);
    printf("\n Ambient Light: (%f,%f,%f)\n\n",ambient_light[0],ambient_light[1],ambient_light[2]);
    
    int x,y;
    for (x = -WIDTH / 2; x < WIDTH / 2; x++) 
    {
        glBegin(GL_POINTS);
        for (y = -HEIGHT / 2; y < HEIGHT / 2; y++) 
        {
            //printf("\n(%d,%d):",x,y);
            Ray* ray = new Ray(0, 0, cameraZ, (double) x*x2i, (double) y*y2j, screenZ);
            
            RT_trace(ray, 1);
            
            if (ray->tIntersect != INFINITY){
                                
                if(ray->color->vector[0]>1){                    
                    //printf("\n(%d,%d):%f",x,y,ray->color->vector[0]);
                    ray->color->vector[0] = 1;                    
                }
                if(ray->color->vector[1]>1){                    
                    //printf("\n(%d,%d):%f",x,y,ray->color->vector[1]);
                    ray->color->vector[1] = 1;
                }
                if(ray->color->vector[2]>1){                    
                    //printf("\n(%d,%d):%f",x,y,ray->color->vector[1]);
                    ray->color->vector[2] = 1;
                }
                
               ray->color->Scale(WHITE_LIGHT);
               plot_pixel(x, y, (int)ray->color->vector[0], (int)ray->color->vector[1], (int)ray->color->vector[2]);                
            }
            else{
                plot_pixel(x, y, WHITE_LIGHT, WHITE_LIGHT, WHITE_LIGHT);
            }
            
//            // For Debugging Purposes              
//            if (double t = ray->CheckPointOnRay(triangles[0].v[0].position)>0){
//                printf("\n(%f,%f): %f",x,y,t);                
//            }

            delete(ray);
        }
        glEnd();
        glFlush();
    }
  
  printf("Done!\n");
  fflush(stdout);
// //For Debugging Purposes
//  printf("Triangles:(%d,%d,%d): ", intersectTriangle[0],intersectTriangle[1],intersectTriangle[2]);  
}

/*****************************Remaining Functions*********************************************************************/

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  glColor3f(((double)r)/256.f,((double)g)/256.f,((double)b)/256.f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  buffer[HEIGHT-y-1][x][0]=r;
  buffer[HEIGHT-y-1][x][1]=g;
  buffer[HEIGHT-y-1][x][2]=b;
}

void plot_pixel(int x,int y,unsigned char r,unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
      plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  Pic *in = NULL;

//  in = pic_alloc(640, 480, 3, NULL);
//  printf("Saving JPEG file: %s\n", filename);
//
//  memcpy(in->pix,buffer,3*WIDTH*HEIGHT);
//  if (jpeg_write(filename, in))
//    printf("File saved Successfully\n");
//  else
//    printf("Error in Saving\n");
//
//  pic_free(in);

}

void parse_check(char *expected,char *found)
{
  if(strcasecmp(expected,found))
    {
      char error[100];
      printf("Expected '%s ' found '%s '\n",expected,found);
      printf("Parse error, abnormal abortion\n");
      exit(0);
    }

}

void parse_doubles(FILE*file, char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE*file,double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE *file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  int i;
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i",&number_of_objects);

  printf("number of objects: %i\n",number_of_objects);
  char str[200];

  parse_doubles(file,"amb:",ambient_light);

  for(i=0;i < number_of_objects;i++)
    {
      fscanf(file,"%s\n",type);
      printf("%s\n",type);
      if(strcasecmp(type,"triangle")==0)
	{

	  printf("found triangle\n");
	  int j;

	  for(j=0;j < 3;j++)
	    {
	      parse_doubles(file,"pos:",t.v[j].position);
	      parse_doubles(file,"nor:",t.v[j].normal);
	      parse_doubles(file,"dif:",t.v[j].color_diffuse);
	      parse_doubles(file,"spe:",t.v[j].color_specular);
	      parse_shi(file,&t.v[j].shininess);
	    }

	  if(num_triangles == MAX_TRIANGLES)
	    {
	      printf("too many triangles, you should increase MAX_TRIANGLES!\n");
	      exit(0);
	    }
	  triangles[num_triangles++] = t;
	}
      else if(strcasecmp(type,"sphere")==0)
	{
	  printf("found sphere\n");

	  parse_doubles(file,"pos:",s.position);
	  parse_rad(file,&s.radius);
	  parse_doubles(file,"dif:",s.color_diffuse);
	  parse_doubles(file,"spe:",s.color_specular);
	  parse_shi(file,&s.shininess);

	  if(num_spheres == MAX_SPHERES)
	    {
	      printf("too many spheres, you should increase MAX_SPHERES!\n");
	      exit(0);
	    }
	  spheres[num_spheres++] = s;
	}
      else if(strcasecmp(type,"light")==0)
	{
	  printf("found light\n");
	  parse_doubles(file,"pos:",l.position);
	  parse_doubles(file,"col:",l.color);

	  if(num_lights == MAX_LIGHTS)
	    {
	      printf("too many lights, you should increase MAX_LIGHTS!\n");
	      exit(0);
	    }
	  lights[num_lights++] = l;
	}
      else
	{
	  printf("unknown type in scene description:\n%s\n",type);
	  exit(0);
	}
    }
  return 0;
}

void display()
{

}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
      draw_scene();
      if(mode == MODE_JPEG)
    	  save_jpg();
  }
  once=1;
}

/*****************Initialize Projection Properties*****************************/
void init()
{
    glMatrixMode(GL_PROJECTION);
    glOrtho(-WIDTH/2,WIDTH/2,-HEIGHT/2,HEIGHT/2,1,-1);    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glClearColor(0,0,0,0);
    glClear(GL_COLOR_BUFFER_BIT);
}

/********************************Main Function*********************************/

int main (int argc, char ** argv)
{
//    if (argc < 2 || argc > 3) {
//        printf("usage: %s <scenefile> [jpegname]\n", argv[0]);
//        exit(0);
//    }
//    if (argc == 3) {
//        mode = MODE_JPEG;
//        filename = argv[2];
//    } 
//    else if (argc == 2)
//        mode = MODE_DISPLAY;

	mode = MODE_DISPLAY;
        
        //char* scene_file = "test1.scene";
	//char* scene_file = "test2.scene";
        //char* scene_file = "spheres.scene";
        //char* scene_file = "table.scene";
        //char* scene_file = "SIGGRAPH.scene";
        //loadScene(scene_file);
    
	glutInit(&argc,argv);
	loadScene(argv[1]);
	

	glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
	glutInitWindowPosition(0,0);
	glutInitWindowSize(WIDTH,HEIGHT);
	int window = glutCreateWindow("Ray Tracer");
	glutDisplayFunc(display);
	glutIdleFunc(idle);
	init();
	glutMainLoop();    
}
