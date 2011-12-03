/* 
 * File:   Vector4.h
 * Author: ASHISH
 *
 * Created on September 20, 2011, 10:55 AM
 */

#ifndef VECTOR4_H
#define	VECTOR4_H

// Forward Declaration
class Matrix4;

class Vector4 {
public:
    double vector[4];
    
    // Constructors for the class
    Vector4();
    Vector4(const double,const double, const double);
    Vector4(const double[3]);
    Vector4(const Vector4*);
    
    // Reinitialize to the values passed as parameters
    void ReInitialize(const double, const double, const double);
    void ReInitialize(const double[3]);
    void ReInitialize();
    
    // Check whether the vector contains positive entries
    bool HasNonNegativeEntries();
    
    // Check whether the vector contains positive entries
    bool HasGreaterThanOneEntries();
    
    // Gives the mod of the x, y and z elements
    double Mod();
    
    // Adds the x, y and z components to the given vector
    void Addition(const Vector4*);
    
    // Subtracts the x, y and z components from the given vector
    void Subtraction(const Vector4*); 
    
    // Gives the Dot Product of the two vectors
    double DotProduct(const Vector4*);
    
    // Gives the Cross Product of the two vectors
    void CrossProductPost(const Vector4*);
    void CrossProductPre(const Vector4*);
    
    // Multiplies the vector with the Matrix
    void PreMultiplication(const Matrix4*); 
    
    // Scales the x, y and z componenets by new values
    void Scale(const double, const double, const double);
    void Scale(const double);
    
    // Converts the vector to unit vector if mod > 0 else returns false
    bool ConvertToUnitVector();
    
    // Homogenizes the vector if value of w!=0 else returns false
    bool Homogenize();
    
    // Displays the given vector
    void Display();       
    
    virtual ~Vector4();
private:
    void Rotate(double, char);
};

#endif	/* VECTOR4_H */