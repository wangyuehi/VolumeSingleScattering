//
//  Vector3.h
//  RayTracer
//
//  Created by WangYue on 16/9/20.
//  Copyright © 2016年 WangYue. All rights reserved.
//

#ifndef Vector3_h
#define Vector3_h

#if defined(_WIN64) ||defined(_WIN32)
/* Microsoft Windows  */
#define M_PI 3.14159265358979323846264
#endif



class Vector3{
public:
    
    Vector3():x(.0),y(.0),z(.0){}
    
    Vector3(const double x,const double y,const double z){
        this->x = x;
        this->y = y;
        this->z = z;
    }
    
    Vector3(const Vector3 &v){
        this->x = v.x;
        this->y = v.y;
        this->z = v.z;
    }



    
    Vector3 operator+(const Vector3 &b) const { return Vector3(x + b.x, y + b.y, z + b.z); }
    
    Vector3 operator+(const double &b) const { return Vector3(x + b, y + b, z + b); }
    
    Vector3 operator+=(const double &b) const { return Vector3(x + b, y + b, z + b); }
    
    Vector3 operator+=(const Vector3 &b) const { return Vector3(x + b.x, y + b.y, z + b.z); }
    
    Vector3 operator-(const Vector3 &b) const { return Vector3(x - b.x, y - b.y, z - b.z); }
    
    Vector3 operator*(double b) const { return Vector3(x * b, y * b, z * b); }
    
    Vector3 operator*(Vector3 b) const { return Vector3(x * b.x, y * b.y, z * b.z); }
    
     Vector3 operator*=(double b) const { return Vector3(x * b, y * b, z * b); }
    
    Vector3 operator/(double b) const { return Vector3(x / b, y / b, z / b); }
    
       bool operator> (Vector3 &b)
    {
        if((x>b.x)&&(y>b.y)&&(z>b.z))
            return true;
        else
            return false;
    }
    bool operator< (Vector3 &b)
    {
        if((x<b.x)&&(y<b.y)&&(z<b.z))
            return true;
        else
            return false;
    }

    
    Vector3 cross(const Vector3 &b){return Vector3(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);}
   
    double dot(const Vector3 &b) const { return x * b.x + y * b.y + z * b.z; }
    Vector3 mult(const Vector3 &b) const { return Vector3(x * b.x, y * b.y, z * b.z); }
    
    double length(){ return sqrt(x * x + y * y + z * z);
    };
   
    Vector3 norm(){
        if(length()  < eps){
            std::cout<<"norm is zero"<<std::endl;
    }
        return  Vector3(x/length(), y / length(), z / length()); }

    double x;
    double y;
    double z;

        
};

#endif /* Vector3_h */

