//
//  Sphere.h
//  RayTracer
//
//  Created by WangYue on 16/9/20.
//  Copyright © 2016年 WangYue. All rights reserved.
//
#ifndef Sphere_h
#define Sphere_h

#include <stdio.h>
#include "Vector3.h"
#include "Object.h"

class Sphere : public Object{
protected:
    double r;
    Vector3 center;
    Color color;
    int exponent;
    Vector3 radiance;
    Material material;
public:
    
    
    Sphere(double r,Vector3 &center,  Color &color){
        this->r = r;
        this->center = center;
        this->color = color;
        this->exponent = 2;
        radiance = Vector3( 0,0,0);
        material = DIFFUSE;
    };
    
    Sphere(double r,Vector3 &center,  Color &color, int exponent){
        this->r = r;
        this->center = center;
        this->color = color;
        this->exponent = exponent;
        radiance = Vector3( 0,0,0);
        material = DIFFUSE;
    };
    Sphere(double r,Vector3 &center,  Color &color, int exponent,Vector3 radiance){
        this->r = r;
        this->center = center;
        this->color = color;
        this->exponent = exponent;
        this->radiance = radiance;
        material = DIFFUSE;
    };
    
    Sphere(double r,Vector3 &center,  Color &color, int exponent,Vector3 radiance,Material m){
        this->r = r;
        this->center = center;
        this->color = color;
        this->exponent = exponent;
        this->radiance = radiance;
        material = m;
    };
    
    Sphere(): r(1.0), center(Vector3(0.0,0.0,0.0)), color (Color(1.0,0.4,0.4)), exponent(25),radiance(0,0,0),material(DIFFUSE){};
    
    virtual Color getObjectColor(){
        //std::cout<<"call Sphere Color"<<std::endl;
        return color;}
    
    virtual Vector3 getSphereCenter(){return center;}
    
     virtual double getSphereRadius(){return r;}
   
    virtual int getExponent(){return exponent;}
    
    virtual Vector3 getNormalAt(Vector3 &point){
        // The normal vector direction is from the center of sphere, ad point out of the sphere
        //std::cout<<"Sphere normalAt"<<std::endl;
        Vector3 normal = point - center;
        normal = normal.norm();
        double dist = normal.dot(normal);
        if(dist - getSphereRadius() * getSphereRadius() > 0.01)
        {
            std::cout<<"Sphere::getNormalAt: not on the sphere!"<<std::endl;
        }
        
        return normal;
    
    }


    virtual double findIntersection(Ray ray){
       // std::cout<<"spere intersections"<<std::endl;
        Vector3 op = center - ray.getOrigin();

        double t;
        double b = op.dot(ray.getDirection());
        double det = b * b - op.dot(op) + r * r;

        if (det < eps) return -1;
        else det = sqrt(det);

        return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : -1);
 
    }
    
    virtual bool insideObject(Vector3 lightPos){
        
        
        Vector3 d = lightPos - getSphereCenter();
        double dist2 = d.dot(d);
        return dist2 < getSphereRadius() * getSphereRadius() ;
        
        
    }
   
    virtual Material getMaterial(){
        return material;
    }
    
    virtual Vector3 getRadiance(){
        
        return radiance;
    }
};

#endif /* Sphere_h */


