/* 
 * File:   Matrix4.cpp
 * Author: ASHISH
 * 
 * Created on September 20, 2011, 9:40 AM
 */

#include "Matrix4.h"
#include "Vector4.h"
#include <iostream>
#include <stdlib.h>
#include <cassert>
#include <math.h>

using namespace std;

Matrix4::Matrix4()
{   
    InitializeMatrix();
}

Matrix4::Matrix4(const Vector4* u, const Vector4* v, const Vector4* w)
{
    assert(u);
    assert(v);
    assert(w);
    
    InitializeMatrix();
    
    int col;
    for(col = 0; col<3; col++)
    {
        matrix[0][col] = u->vector[col];
        matrix[1][col] = v->vector[col];
        matrix[2][col] = w->vector[col];
    }   
}

Matrix4::Matrix4(const Matrix4* mat)
{
    assert(mat);
    int row, col;
    for(row = 0; row<4; row++)
    {
        for(col = 0; col<4; col++)
        {
            matrix[row][col] = mat->matrix[row][col];
        }
    }
}

void Matrix4::InitializeMatrix()
{
    int row, col;
    for(row = 0; row<4; row++)
    {
        for(col = 0; col<4; col++)
        {
            matrix[row][col] = 0;
        }
    }
    for(row = 0; row<4; row++)
    {
        matrix[row][row] = 1;
    }
}

void Matrix4::Addition(const Matrix4* mat)
{
    assert(mat);    
    int row, col;
    for(row = 0; row<4; row++)
    {
        for(col = 0; col<4; col++)
        {
            matrix[row][col] += mat->matrix[row][col];
        }
    }
}

void Matrix4::Subtraction(const Matrix4* mat)
{
    assert(mat);
    int row, col;
    for(row = 0; row<4; row++)
    {
        for(col = 0; col<4; col++)
        {
            matrix[row][col] -= mat->matrix[row][col];
        }
    }
}

void Matrix4::ScalarMultiplication(const double scale)
{
    int row, col;
    for(row = 0; row<4; row++)
    {
        for(col = 0; col<4; col++)
        {
            matrix[row][col] *= scale;
        }
    }
}

void Matrix4::PreMultiplication(const Matrix4* mat)
{
    assert(mat);
    double temp[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    int row, col, k;
    for(row = 0; row<4; row++)
    {
        for(col = 0; col<4; col++)
        {            
            for(k = 0; k<4; k++)
            {
                temp[row][col] += mat->matrix[row][k]*matrix[k][col];
            }            
        }
    }
    
    for(row = 0; row<4; row++)
    {
        for(col = 0; col<4; col++)
        {    
            matrix[row][col] = temp[row][col];
        }
    }    
}

void Matrix4::PostMultiplication(const Matrix4* mat)
{
    assert(mat);
    double temp[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    int row, col, k;
    for(row = 0; row<4; row++)
    {
        for(col = 0; col<4; col++)
        {            
            for(k = 0; k<4; k++)
            {
                temp[row][col] += matrix[row][k]*mat->matrix[k][col];
            }            
        }
    }
    
    for(row = 0; row<4; row++)
    {
        for(col = 0; col<4; col++)
        {    
            matrix[row][col] = temp[row][col];
        }
    }    
}

void Matrix4:: TransposeMatrix()
{
    int i, j;
    double temp;
    for (i = 0; i<4; i++)
    {
        for (j = i+1; j<4; j++)
        {
            temp =  matrix[i][j];
            matrix[i][j] = matrix[j][i];
            matrix[j][i] = temp;
        }
    }
}

void Matrix4::SwapRows(const int i, const int j)
{
    assert(i<4);
    assert(j<4);
   
    if (i!=j)
    {    
        double temp = 0.0;
        int col;
        for(col=0; col<4; col++)
        {
            temp = matrix[i][col];
            matrix[i][col] = matrix[j][col];
            matrix[j][col] = temp;
        }
    }
}

int Matrix4::FindRowHavingGreatestElement(const int start_row, const int col)
{
    int row;
    int max_row = start_row;
    double max = matrix[start_row][col];
    
    for(row = start_row+1; row<4; row++)
    {
        if(matrix[row][col]> max)
        {
            max = matrix[row][col];
            max_row = row;
        }
    }
    
    return max_row;
}

// Gauss-Jordan Method
void Matrix4:: InverseMatrix()
{   
    double temp[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    
    int row, col, i;
    for(row = 0; row<4; row++)
    {
        temp[row][row] = 1;
    }
    
    for (i = 0; i<4; i++)
    {        
//        int greatest = FindRowHavingGreatestElement(i, i);
//        SwapRows(i, greatest);
        double m = matrix[i][i];
        if (m!=0)
        {
            for(row = 0; row<4; row++)
            {
                if (row!=i)
                {
                    double n = matrix[row][i];                    
                    for (col = 0; col<4; col++)
                    {
                        matrix[row][col] -= (matrix[i][col]/m)*n;
                        temp[row][col] -= (temp[i][col]/m)*n;                        
                    }                
                }
            }
        }                
    }
    
    // Homogenizing
    for(row = 0; row<4; row++)
    {        
        for(col = 0; col<4; col++)
        {    
            temp[row][col] /= matrix[row][row];            
        }
        matrix[row][row] = 1;
    }
    
    for(row = 0; row<4; row++)
    {
        for(col = 0; col<4; col++)
        {    
            matrix[row][col] = temp[row][col];
        }
    }
}

void Matrix4::Display()
{    
    int row, col;    
    for (row = 0; row<4; row++)
    {
        for (col = 0; col<4; col++)
        {
            cout<<matrix[row][col]<<"\t";
        }    
        cout<<"\n";
    }    
    cout<<"\n";
}

void Matrix4::DisplayRow(const int row)
{
    int col;
    for (col = 0; col < 4; col++)
    {
        cout << matrix[row][col] << "\t";
    }
    
    cout << "\n";
}

Matrix4::~Matrix4() 
{    
    //free(matrix);
}

