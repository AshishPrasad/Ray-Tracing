/* 
 * File:   Vector4.cpp
 * Author: ASHISH
 * 
 * Created on September 20, 2011, 10:55 AM
 */
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <cassert>

#include "Matrix4.h"
#include "Vector4.h"

using namespace std;

Vector4::Vector4() 
{    
    this->ReInitialize(0,0,0);
}

Vector4::Vector4(const double x, const double y, const double z)
{   
    this->ReInitialize(x,y,z);
}

Vector4::Vector4(const double vect[3])
{   
    this->ReInitialize(vect[0],vect[1],vect[2]);
}

Vector4::Vector4(const Vector4* vect) 
{
    assert(vect);    
    this->ReInitialize(vect->vector);
}

void Vector4::ReInitialize(const double x, const double y, const double z)
{
    vector[0] = x;
    vector[1] = y;
    vector[2] = z;
    vector[3] = 1;
}

void Vector4::ReInitialize(const double v[3])
{
    this->ReInitialize(v[0], v[1], v[2]);
}

void Vector4::ReInitialize()
{
    this->ReInitialize(0,0,0);
}

bool Vector4::HasNonNegativeEntries(){
    return (this->vector[0]>=0 && this->vector[0]>=0 && this->vector[0]>=0);
}

bool Vector4::HasGreaterThanOneEntries(){
    return (this->vector[0]>=1 || this->vector[0]>=1 || this->vector[0]>=1);
}

double Vector4::Mod()
{
    double mod = 0.0;    
    int row;
    for(row=0; row<3; row++)
    {
        mod += pow(vector[row],2);
    }
    
    return sqrt(mod);
}

void Vector4::Display()
{    
    int row;
    cout<<"[";
    for(row=0; row<3; row++)
    {
        cout<<vector[row]<<",";
    }
    cout<<vector[3]<<"]";
}

void Vector4::Addition(const Vector4* vect)
{
    assert(vect);
    int row;
    for(row=0; row<3; row++)
    {
        vector[row] += vect->vector[row];
    }
}

void Vector4::Subtraction(const Vector4* vect)
{
    assert(vect);
    int row;
    for(row=0; row<3; row++)
    {
        vector[row] -= vect->vector[row];
    }
}

void Vector4::PreMultiplication(const Matrix4* matrix)
{   
    assert(matrix);    
    double temp[4] = {0.0, 0.0, 0.0, 0.0};    
    int row, col, k;        
    for(row = 0; row<4; row++)
    {           
        for(col = 0; col<4; col++)
        {            
            temp[row] += matrix->matrix[row][col]*vector[col];                        
        }        
    }
    
    for(row = 0; row<4; row++)
    {
        vector[row] = temp[row];
    }
}

double Vector4::DotProduct(const Vector4* vect)
{
    assert(vect);
    double val = 0.0;
    int row;
    for(row = 0; row<3; row++)
    {
        val += vector[row]*vect->vector[row];
    }    
    return val;
}

void Vector4::CrossProductPost(const Vector4* vect)
{
    assert(vect);
    double temp[4] = {0.0, 0.0, 0.0, 0.0};    
    int row;
    for (row = 0; row<3; row++)
    {
        temp[row] = vector[(row+1)%3]*vect->vector[(row+2)%3] - 
                vector[(row+2)%3]*vect->vector[(row+1)%3];
    }
    
    for (row = 0; row<3; row++)
    {
        vector[row] = temp[row];
    }    
}

void Vector4::CrossProductPre(const Vector4* vect)
{
    assert(vect);
    double temp[4] = {0.0, 0.0, 0.0, 0.0};    
    int row;
    for (row = 0; row<3; row++)
    {
        temp[row] = vect->vector[(row+1)%3]*vector[(row+2)%3] - 
                vect->vector[(row+2)%3]*vector[(row+1)%3];
    }
    
    for (row = 0; row<3; row++)
    {
        vector[row] = temp[row];
    }    
}

void Vector4::Scale(const double x, const double y, const double z)
{
    vector[0] *= x;
    vector[1] *= y;
    vector[2] *= z;
}

void Vector4::Scale(const double scale)
{
    Scale(scale,scale,scale);
}

bool Vector4::Homogenize()
{
    if(vector[3] == 0)
    {
        return false;
    }
    
    if(vector[3]!=1)
    {
        int row;
        for(row=0; row<4; row++)
        {
            vector[row] /= vector[3];
        }
    }
    
    return true;
}

bool Vector4::ConvertToUnitVector()
{
    bool convert = false;
    
    double mod = this->Mod();
    int row;
    
    if(mod>0)
    {               
        for(row=0;row<3;row++)
        {
            vector[row] = vector[row]/mod;
        }
        convert = true;
    }  
    
    return convert;
}

Vector4::~Vector4() {
 //   free(vector);
}