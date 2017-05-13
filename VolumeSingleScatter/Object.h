//
//  Object.h
//  RayTracer
//
//  Created by WangYue on 16/9/21.
//  Copyright © 2016年 WangYue. All rights reserved.
//

#ifndef Object_h
#define Object_h

#include <stdio.h>
#include "Vector3.h"
#include "Ray.h"
#include "Color.h"
enum Material {GLOSSY, DIFFUSE,SPECULAR};
class Object{
    
    int exponent;
    Color color;
    Material material;
public:
    
    Object();

    virtual Color getObjectColor(){
        std::cout<<"!!!call Object color"<<std::endl;
        return Color(0.0,0.0,0.0);
    };
    
    virtual double findIntersection(Ray ray){
        std::cout<<"obj intersection"<<std::endl;
        
        return -1;
    
    };
    
    virtual Vector3 getNormalAt(Vector3 &intersectionPoint){
        std::cout<<"!!!call Object Normal"<<std::endl;
        return Vector3(0,0,0);
    };
    virtual int getExponent(){
        std::cout<<"!!!call Object Exponent"<<std::endl;
        return exponent;
    }

    // Computer reflect direction
    virtual Vector3 reflect(const Vector3 &I, const Vector3 &N){
        // The direction of Wi is into the point
        //std::cout<<"calling the object reflect"<<std::endl;
        return N * (I.dot(N))* 2 - I;
    }
    virtual bool insideObject(Vector3 lightPos){
        std::cout<<"!!!call Object Inside"<<std::endl;        
        return false;
        
    }
    virtual Vector3 getRadiance(){
         std::cout<<"!!!call Object GetRadiance"<<std::endl;
        return Vector3(0,0,0);
    }
    virtual Material getMaterial (){
        std::cout<<"!!!call Object Material"<<std::endl;
        return DIFFUSE;
    }
    
};

Object::Object(){};


#endif /* Object_h */
