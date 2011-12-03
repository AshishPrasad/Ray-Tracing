/* 
 * File:   Pixel.cpp
 * Author: ashish
 * 
 * Created on 27 October, 2011, 4:07 PM
 */

#include <stdio.h>
#include <math.h>
#include <cassert>

#include "Ray.h"
#include "Vector4.h"
#include "ImageDataStructure.h"

Ray::Ray(){
    this->source = new Vector4(0,0,0);
    this->dir = new Vector4(0, 0, -1);  
    this->Initialize();
}

Ray::Ray(const double p0_x, const double p0_y, const double p0_z, 
        const double p1_x, const double p1_y, const double p1_z) {
    
    this->source = new Vector4(p0_x,p0_y,p0_z);
    
    // dir = end - source
    this->dir = new Vector4(p1_x, p1_y, p1_z);    
    this->dir->Subtraction(source);        
    this->dir->ConvertToUnitVector();
    
    this->Initialize();
}

Ray::Ray(const Vector4* source, const Vector4* end){
    
    this->source = new Vector4(source);
    
    // dir = end - source
    this->dir = new Vector4(end);
    this->dir->Subtraction(source);
    this->dir->ConvertToUnitVector();
    
    this->Initialize();
}

void Ray::Initialize(){
    this->tIntersect = INFINITY;    
    this->intersectPoint = new Vector4();    
    this->intersectedObject = NONE;
    this->intersectID = NONE;
    
    this->normal = new Vector4();    
    this->color = new Vector4();    
}

Ray::Ray(const Ray& orig) {
}

Ray::~Ray() {
}

