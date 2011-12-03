/* 
 * File:   Pixel.h
 * Author: ashish
 *
 * Created on 27 October, 2011, 4:07 PM
 */
#ifndef PIXEL_H
#define	PIXEL_H

class Vector4;

class Ray {    
public:
    Vector4* source;
    Vector4* dir;
    double tIntersect;
    Vector4* intersectPoint;
    int intersectedObject;
    int intersectID;   
    Vector4* normal;
    Vector4* color;
    
    Ray();    
    Ray(const double p0_x, const double p0_y, const double p0_z, 
        const double p1_x, const double p1_y, const double p1_z);
    Ray(const Vector4* source, const Vector4* end);
    
    Ray(const Ray& orig);    
    virtual ~Ray();
private:
    void Initialize();
};

#endif	/* PIXEL_H */