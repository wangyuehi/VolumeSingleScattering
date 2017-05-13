//
//  Light.h
//  LightTracer
//
//  Created by WangYue on 16/9/20.
//  Copyright © 2016年 WangYue. All rights reserved.
//

#ifndef Light_h
#define Light_h


#include "Vector3.h"
#include "Color.h"

class Light
{
    Vector3 position;
    Color color;
    double Radiance;
    
    
public:
    
    
    Light(const Vector3 &o, const Color &dir,double &radiance);
    Light();
    
    Vector3 getLightPosition(){return position;};
    Color getLightColor(){return color;};
    double getRadiance(){return Radiance;}
    
    
};

Light::Light(){
    position = Vector3(0,0,0);
    color = Color(1,1,1);
    Radiance = 10000;
    
}

Light::Light(const Vector3 &o, const Color &c,double &radiance){
    position = o;
    color = c;
    Radiance = radiance;
}

#endif /* Light_h */
