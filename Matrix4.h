/* 
 * File:   Matrix4.h
 * Author: ASHISH
 *
 * Created on September 20, 2011, 9:40 AM
 */

#ifndef MATRIX4_H
#define	MATRIX4_H

// Forward Declaration
class Vector4;

class Matrix4 {
public:
    double matrix[4][4];  
    
    // Constructor functions
    Matrix4();
    Matrix4(const Vector4*, const Vector4*, const Vector4*);
    Matrix4(const Matrix4*);    
    
    // Initialize to identity
    void InitializeMatrix();
    
    void Addition(const Matrix4*);
    void Subtraction(const Matrix4*);
    void ScalarMultiplication(const double);
    void PreMultiplication(const Matrix4*);
    void PostMultiplication(const Matrix4*);
    void TransposeMatrix();
    void InverseMatrix();    
    void SwapRows(const int i, const int j);
    
    void Display();
    void DisplayRow(const int);    
    
    virtual ~Matrix4();  
    
private:
    int FindRowHavingGreatestElement(const int, const int);    
};
#endif	/* MATRIX4_H */