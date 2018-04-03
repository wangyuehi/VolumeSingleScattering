#if defined(_WIN32)
#include <Windows.h>
#endif
//
//  main.cpp
//  RayTracer
//
//  Created by WangYue on 16/9/20.
//  Copyright © 2016年 WangYue. All rights reserved.
//
//  g++ -c main.cpp

#define inf 1e20
#define eps 1e-20

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <math.h>
#include <time.h>
#include <algorithm> 
#include <stdlib.h>
#include <stdio.h>



#include "Vector3.h"
#include "Camera.h"
#include "Color.h"
#include "Light.h"
#include "Sphere.h"
#include "Plane.h"
#include "Triangle.h"
#include "MotionSphere.h"
#include "Sampler.h"
#include "Random.h"
#include "Volume.h"
#include "MCSampler1D.h"
#include "MCSampler2D.h"

using namespace std;


inline double clamp(double x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
inline int toDisplayValue(double x){ return int( pow(clamp(x), 1.0/2.2) * 255 + .5); }

void savePPM(const char *filename, int w, int h, Color *pixelColors){

#if defined(_WIN32)
    FILE *f;
    fopen_s(&f, filename, "w");
#else
    FILE *f = fopen(filename, "w");
#endif
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    
    for (int p = 0; p < w * h; p++) {
        fprintf(f,"%d %d %d ", toDisplayValue(pixelColors[p].x), toDisplayValue(pixelColors[p].y), toDisplayValue(pixelColors[p].z));
    }
    fclose(f);

}


bool collectIntersections(vector<double> &intersections,vector<Object*> &objects, Ray &cameraRay){
    // Find intersections of cameraRay with every object, if there's no intersection, return inf.
    // If there's no intersection, return false.
    bool IntersectionExists = false;
    for( int index = 0 ; index < objects.size(); index ++) {
        double t = objects[index]->findIntersection(cameraRay);
        if(t > 0){
            intersections.push_back(t);
            if(!IntersectionExists)
                IntersectionExists = true;
        }
        else
            intersections.push_back(inf);
    }
    
    return IntersectionExists;
  
}


bool ClosestIntersectionObject(vector<Object*> &objects, Ray &cameraRay, double &t, int &id){
    // Find the object that has the closest intersection point
    // Return true if there's intersection, then t is the first intersection point among all spheres
    
    vector<double> intersections;
    
    t = inf;
  
    if( collectIntersections(intersections, objects, cameraRay))
        for( int index = 0; index < intersections.size(); index ++ ) {
            // Find the smallest t;
            if(intersections.at(index)<t){
                t = intersections.at(index);
                id = index;
            }
        }
    
    return t < inf;
    
};

bool trace(
           Vector3 &orig,  Vector3 &dir,
           vector<Object*> &objects,
           double &tNear, int &index)
{
    tNear = inf;
    Ray ray (orig, dir);
    for (int k = 0; k < objects.size() - 1; ++k) {
        if (objects[k]->findIntersection(ray)  < tNear && objects[k]->findIntersection(ray) > 0  ) {
            
            tNear = objects[k]->findIntersection(ray) ;
            index = k;
            
        }
    }
    return (tNear<inf);
    
}

bool traceAll(
           Vector3 &orig,  Vector3 &dir,
            vector<Object*> &objects,
           double &tNear, int &index)
{
    tNear = inf;
    Ray ray (orig, dir);
    for (int k = 0; k < objects.size()  ; ++k) {

        if (objects[k]->findIntersection(ray)  < tNear && objects[k]->findIntersection(ray) > 0 ) {
           
            tNear = objects[k]->findIntersection(ray) ;
            index = k;
         
        }
    }
    return (tNear<inf);

}



Color getPixelColorPointLight(vector<Object*> &objects, int id, double &t, Ray &ray, Light &light){
    
    
    Color c;
    Object* obj;
    
    
    // Intersection Point
    Vector3 intersectionPoint = ray.getOrigin() + ray.getDirection() * t;
   
    Vector3 normal;
    if(id>=0){
         obj = objects.at(id);
    // Normal of intersection point, point out of the sphere
    normal = obj->getNormalAt(intersectionPoint);
    }
    else{
        obj = objects.at(objects.size() - 1);
        normal = Vector3 (0,0,1);
    }
    // If there are multiple lights, should be a for loop here
    
    // Light direction, from intersection point to the light position

    Vector3 lightDirection =  light.getLightPosition() - intersectionPoint ;
    
    // Square of Distance between light and point
    double distance2 = lightDirection.dot(lightDirection);
    
    // Light direction
    lightDirection = lightDirection.norm();
    double tNear = inf;
    int shadowObjIndex = 0;
    // Shadow check
    // Move the shadow point a little bit along normal
    double offset = 1e-9;
    Vector3 shadowPointOrig =
    (ray.getDirection().dot(normal) < 0) ?
    intersectionPoint + normal * offset :
    intersectionPoint - normal * offset ;
    
    bool inShadow = trace(shadowPointOrig, lightDirection, objects, tNear, shadowObjIndex ) &&
    tNear * tNear < distance2 ;
    
    double cosine;
    if(obj->insideObject(light.getLightPosition())){
        cosine = clamp(- lightDirection.dot(normal));
    }
    else{
        cosine = clamp(lightDirection.dot(normal));
    }
    
    
    double L_r;
    
    double fr = 1 / M_PI;
    if(id>=0){
        L_r =  clamp(light.getRadiance()  *  fr * cosine *  (1 - inShadow) / (distance2 * 4 * M_PI));
    }
    else{
        L_r = light.getRadiance()  *  fr *  (1 - inShadow) / ( distance2 * 4 * M_PI);
// cout<<distance2<<endl;
    }
    //task1
    //double L_r =  light.getRadiance() * fr * cosine  / (distance2 * 4 * M_PI);
//    cout<<distance2<<endl;
    c = obj->getObjectColor() * L_r  ;
    
    return c;
    
}

Color getPixelColorWithSphereLight(vector<Object*> &objects, int id, double &t_ray, Ray &ray, Vector3& hitPoint, const int N, const SampleMode sampleMode){
    
    
    
    Color c;
    
   
    
    // Intersection Point
    Vector3 intersectionPoint = ray.getOrigin() + ray.getDirection() * t_ray;
    // Get the light Obj
    // Here the last object in vector is light by default.
    Sphere* LightObj =(Sphere*)objects.at(objects.size() - 1);
    // Normal of intersection point
    Object* obj;
    Vector3 normal;
    if(id>=0){
        obj = objects.at(id);
        normal = obj->getNormalAt(intersectionPoint);
        if(obj->insideObject(LightObj->getSphereCenter())){
            normal = normal * -1;
        }
    }
    else{
        obj = objects.at(objects.size() - 1);
        normal = Vector3(0,0,1);
        
    }
  
    
    
   
    
    // vector3 from intersection point X to light center C
    Vector3 XC = LightObj->getSphereCenter() - intersectionPoint;
    
    // ThetaMax
    double sinThetaMax = LightObj->getSphereRadius() / XC.length();
    double cosThetaMax = sqrt(1 - sinThetaMax * sinThetaMax);
    
    double lightR = LightObj->getSphereRadius();

    
    if( id == objects.size() - 1){
        
        c = obj->getObjectColor();
    }
    
    else{
        
        // Monte Carlo Begins
        
        double L_r  = 0.0;
        
        for(int m = 0; m < N; m++){
            
            Vector3  sampleLightDirection;
            
            sampleLightDirection = getSample(*LightObj,intersectionPoint, sampleMode, normal, sinThetaMax).norm();
            
            Ray secondary (intersectionPoint, sampleLightDirection);
            double t = LightObj->findIntersection(secondary);
            
            if(t > eps)
                
            {
                
                Vector3 sampleLightPoint = intersectionPoint + sampleLightDirection * t;
                hitPoint = sampleLightPoint;
                
                // Check if there are intersections between light sample and intersection point.
                
                // Light direction, from intersection point to the light position
                Vector3 lightDirection = sampleLightPoint - intersectionPoint ;
                
                // Distance between light sample point and the intersection point
                double distance2 = lightDirection.dot(lightDirection);
                
                // Light direction
                lightDirection = lightDirection.norm();
                double tNear = inf;
                int shadowObjIndex = 0;
                // Shadow check
                // Move the shadow point a little bit along normal
                double offset = 1e-9;
                Vector3 shadowPointOrig =
                (ray.getDirection().dot(normal) < 0) ?
                intersectionPoint + normal * offset :
                intersectionPoint - normal * offset ;
                
                bool inShadow = trace(shadowPointOrig, lightDirection, objects, tNear, shadowObjIndex ) &&
                tNear * tNear < distance2 ;
                
                
                
                
                double L_i = 10.0;
                double fr = 0.0;
                Vector3 lightSampleNormal = LightObj->getNormalAt(sampleLightPoint);
                
                double cosThetaI = clamp(lightDirection.dot(normal));
                double cosThetaO = clamp(-lightSampleNormal.dot(lightDirection));
                
                if(id<0){
                    cosThetaI = 1;
                }
                
                
                if(obj->getMaterial() == DIFFUSE){
                    fr = 1 / M_PI;
                }
                else if(obj->getMaterial() == SPECULAR)// Blinn-phong
                {
                    Vector3 HalfWay = normal + lightDirection;
                    HalfWay = HalfWay.norm();
                    fr = pow( HalfWay.dot(normal), obj->getExponent()) * (obj->getExponent() + 2) / (2 * M_PI);
                }
                else {
                    cout<< "no model for such material!";
                }
                
                switch (sampleMode) {
                    case Importance_Solid_Angle:
                    {
                        double pdfConstant = 2 * M_PI * (1 - cosThetaMax);
                        L_r =  L_r + fr * L_i * cosThetaI *  (1 - inShadow) * pdfConstant ;
                        
                    }
                        
                        break;
                    case Importance_Surface_Area:
                    {   // pdf = distance * distance / (CosThetaO * Area)
                        double  Area =  2 * M_PI * (1 - sinThetaMax) * lightR * lightR / distance2;
                        
                        // Fr = (fr(x, Wi, Wr) * Li(r(x, Wik), -Wik)) * cosThetaI / pdf;
                        // L_r = fr * L_i * (1 - inShadow) * cosThetaI * cosThetaO * area / distance
                        
                        L_r += fr * L_i * (1 - inShadow) * cosThetaI * cosThetaO * Area ;
                        
                    }
                        
                        break;
                    case Uniform_Solid_Angle:
                    {
                        L_r += fr * L_i * (1 - inShadow) * cosThetaI * 4 * M_PI ;
                    }
                        break;
                    case Cosine_Weighted:
                    {
                        L_r += fr * L_i * (1 - inShadow)  * M_PI;
                    }
                        break;
                    case Uniform_Hemisphere_Solid_Angle:
                    {
                        L_r += fr * L_i * (1 - inShadow) * cosThetaI * 2 * M_PI ;
                    }
                        break;
                    default:
                        cout<<"inavailable sample mode"<<endl;
                        
                        
                }
                
            }
            else // No intersection for sample light ray with sphere light
            {
                
                // cout<<"miss out"<<endl;
            }
            
        }
        L_r = clamp(L_r/N);
        
        c = obj->getObjectColor() * L_r  ;
    }
    
    return c;
    
}


Color getPixelColorWithSphereLightDBG(vector<Object*> &objects, int id, double &t, Ray &ray, Light &light, const int N, const SampleMode sampleMode){
    
    
    
    Color c;
    
    Object* obj = objects.at(id);
    
    // Intersection Point
    Vector3 intersectionPoint = ray.getOrigin() + ray.getDirection() * t;
    
    // Normal of intersection point
    Vector3 normal = obj->getNormalAt(intersectionPoint);
    if(obj->insideObject(light.getLightPosition())){
        normal = normal * -1;
    }
    
    
    // Get the light Obj
    // Here the last object in vector is light by default.
    Sphere* LightObj =(Sphere*)objects.at(objects.size() - 1);
    
    
    // vector3 from intersection point X to light center C
    Vector3 XC = LightObj->getSphereCenter() - intersectionPoint;
    
    // ThetaMax
    double sinThetaMax = LightObj->getSphereRadius() / XC.length();
    double cosThetaMax = sqrt(1 - sinThetaMax * sinThetaMax);
    
    double lightR = LightObj->getSphereRadius();
    
    
    if( id == objects.size() - 1){
        
        c = obj->getObjectColor();
    }
    
    else{
        
        // Monte Carlo Begins
        
        Vector3 L_r (0.0,0.0,0.0);
        
        for(int m = 0; m < N; m++){
            
            Vector3  sampleLightDirection;
            
            sampleLightDirection = getSample(*LightObj,intersectionPoint, sampleMode, normal, sinThetaMax).norm();
            //sampleLightDirection = CosWeightedRandomHemisphereDirection2(normal);
            Ray secondary (intersectionPoint, sampleLightDirection);
            double t = LightObj->findIntersection(secondary);
            
            if(t > eps)
                
            {
                //cout<<"no miss "<<endl;
                Vector3 sampleLightPoint = intersectionPoint + sampleLightDirection * t;
                
                
                // Check if there are intersections between light sample and intersection point.
                
                // Light direction, from intersection point to the light position
                Vector3 lightDirection = sampleLightPoint - intersectionPoint ;
                
                // Distance between light sample point and the intersection point
                double distance2 = lightDirection.dot(lightDirection);
                
                // Light direction
                lightDirection = lightDirection.norm();
                double tNear = inf;
                int shadowObjIndex = 0;
                // Shadow check
                // Move the shadow point a little bit along normal
                double offset = 1e-9;
                Vector3 shadowPointOrig =
                (ray.getDirection().dot(normal) < 0) ?
                intersectionPoint + normal * offset :
                intersectionPoint - normal * offset ;
                
                bool inShadow = trace(shadowPointOrig, lightDirection, objects, tNear, shadowObjIndex ) &&
                tNear * tNear < distance2 ;
                
                
                
                
                Vector3 L_i (10.0,10.0,10.0);
                
                Vector3 lightSampleNormal = LightObj->getNormalAt(sampleLightPoint);
                
                double cosThetaI = clamp(lightDirection.dot(normal));
                double cosThetaO = clamp(-lightSampleNormal.dot(lightDirection));
                
                double brdf = 0.0;
                if(obj->getMaterial() == DIFFUSE){
                    brdf = 1 / M_PI;
                }
                else if(obj->getMaterial() == SPECULAR)
                {
                 
       
                    Vector3 outgoing  = (normal * (normal.dot(lightDirection))* 2 - lightDirection).norm() ;
                    //Vector3 Half = (outgoing + lightDirection).norm();
                    double temp =(outgoing.dot(ray.getDirection()));
                    brdf = pow( -1 * temp, obj->getExponent()) * (obj->getExponent() + 2) / (2 * M_PI);
                     //brdf = pow( Half.dot(normal), obj->getExponent()) * (obj->getExponent() + 2) / (2 * M_PI);
//                    if((normal.dot(lightDirection)<=0)||brdf <eps)
//                        cout<<"???";
                    cout<<"still but with phong cosine power";

                }
                else {
                    cout<< "no model for such material!";
                }
                
                switch (sampleMode) {
                    case Importance_Solid_Angle:
                    {
                        double pdfConstant = 2 * M_PI * (1 - cosThetaMax);
                        L_r =  L_r +  L_i * cosThetaI *  (1 - inShadow) * pdfConstant * brdf;
                     
                    }
                        
                        break;
                    case Importance_Surface_Area:
                    {   // pdf = distance * distance / (CosThetaO * Area)
                        double  Area =  2 * M_PI * (1 - sinThetaMax) * lightR * lightR / distance2;
                        
                        // Fr = (fr(x, Wi, Wr) * Li(r(x, Wik), -Wik)) * cosThetaI / pdf;
                        // L_r = fr * L_i * (1 - inShadow) * cosThetaI * cosThetaO * area / distance
                        
                        L_r += L_i * (1 - inShadow) * cosThetaI * cosThetaO * Area * brdf ;
                        
                    }
                        
                        break;
                    case Uniform_Solid_Angle:
                    {
                        L_r +=  L_i * (1 - inShadow) * cosThetaI * 4 * M_PI * brdf;
                    }
                        break;
                    case Importance_Phong:{
                        L_r += L_i * (1 - inShadow) * cosThetaI;
                    }
                        break;
                    case Cosine_Weighted:{
                        L_r += L_i * (1 - inShadow) * brdf;
                    }
                        break;
                    default:
                        cout<<"inavailable sample mode"<<endl;
                        
                        
                }
                
            }
//            else // No intersection for sample light ray with sphere light
//            {
//                
//                 cout<<"miss out"<<endl;
//            }
            
        }
        L_r = L_r/N;
        
        c = obj->getObjectColor() * L_r  ;
       // c =  L_r  ;
    }
    
    return c;
    
}

Color getPixelColorAmbientOcclusion(vector<Object*> &objects, int id, double &t, Ray &ray, Light &light, const int N, const SampleMode sampleMode){
    
    
    
    Color c;
    
    Object* obj = objects.at(id);
    
    // Intersection Point
    Vector3 intersectionPoint = ray.getOrigin() + ray.getDirection() * t;
    
    // Normal of intersection point, point out of the sphere
    Vector3 normal = obj->getNormalAt(intersectionPoint);
  
    double L_r  = 0.0;
    
        for(int m = 0; m < N; m++){
            
            
            Vector3  sampleLightDirection = getSample((Sphere&)*obj,intersectionPoint, sampleMode, normal, 0.0).norm() ;
            
            if(sampleLightDirection.dot(normal)<=0){
                //cout<<"hh"<<endl;
                sampleLightDirection = sampleLightDirection * (-1);
            }
            
            Ray secondary (intersectionPoint, sampleLightDirection);
         
            
            // Check if there are intersections in the direction.
            
            double tNear = inf;
            int shadowObjIndex = 0;
            // Shadow check
            // Move the shadow point a little bit along normal
            double offset = 1e-9;
            Vector3 shadowPointOrig =
            (ray.getDirection().dot(normal) < 0) ?
            intersectionPoint + normal * offset :
            intersectionPoint - normal * offset ;
            
            bool inShadow = traceAll(shadowPointOrig, sampleLightDirection, objects, tNear, shadowObjIndex ) ;

            double cosThetaI = clamp(sampleLightDirection.dot(normal));
     
            switch (sampleMode) {
                case Cosine_Weighted:
                {
                    double ro = 1.0;
                    L_r += ro * (1 - inShadow) / M_PI  ;
                }
                    break;
        
                default:// Uniform_Solid_Angle:
                
                   L_r += (1 - inShadow) * cosThetaI  * 2 / M_PI;
                   break;
            }
        }
        L_r = L_r/N;
        
        c = obj->getObjectColor() * L_r  ;
    
    
    return c;
    
}

Color getPixelColor_ExplicitPathTracing(vector<Object*> &objects, int id, double &t, Ray &ray, Light &light,  int N, const SampleMode sampleMode){
    // N is iteration times
    
    Color c (0,0,0), L_direct (0,0,0), L_indirect (0,0,0);

    Object* hitObj = objects.at(id);
    
    
    if(id == objects.size() - 1){
        
        // The hit point is on light
        c = hitObj->getRadiance();
        
    }
    else{
        //if(N < 8)
        {
            // The hit point is on another object
            
            // Ld += contribution from light;
            L_direct  = getPixelColorWithSphereLightDBG(objects, id, t, ray, light, 1, sampleMode);
            
            // ω’ = random direction in hemisphere above n;
            // Li += brdf * shade(trace(x, ω’)) * dot(n, ω’) / (p(ω’));
            Vector3 intersectionPoint = ray.getOrigin() + ray.getDirection() * t;
            
            Vector3 normal = hitObj->getNormalAt(intersectionPoint);
            if(hitObj->insideObject(ray.getOrigin())){
                normal = normal * -1;
            }
            double offset = 1e-9;
            Vector3 offsetIntersectionPoint =
            (ray.getDirection().dot(normal) < 0) ?
            intersectionPoint + normal * offset :
            intersectionPoint - normal * offset ;
            
            
            // russian roulette
            Color ObjColor = hitObj->getObjectColor();
            
            // Q is the largest among r g b
            double Q = max(ObjColor.x, max(ObjColor.y, ObjColor.z));
            
            // if random is smaller than Q
            // if objColor is very near to black, then Q should be very small, so the chance of randf() < Q should be small, so it have a better chance to not recurse again.
            if(randf() < Q){
                // Shoot a new ray
                Vector3 w ;
                Ray second_ray;
                int IndexOfClosestObject;
                double hit_t;
                bool missScene = true;
                
                // Avoid double count
                do{
                    w = getSampleUniformSolidAngleHemisphere(intersectionPoint, normal);
                    second_ray  = Ray (offsetIntersectionPoint,w);
                    missScene = !ClosestIntersectionObject(objects, second_ray, hit_t, IndexOfClosestObject);
                }while(IndexOfClosestObject == objects.size() -1);
                
                //if(!missScene){
                    double brdf = 1 / M_PI;
                    double cosine = clamp(normal.dot(w));
                    double pdf = 1 / (2 * M_PI);
                    L_indirect = hitObj ->getObjectColor() * getPixelColor_ExplicitPathTracing(objects, IndexOfClosestObject, hit_t, second_ray, light, ++N, sampleMode) * brdf * cosine / (Q * pdf);
                //}
            }
        }
            
        c = (L_direct + L_indirect );
     
    }
    
    return c;
    
}

Color getPixelColor_ImplicitPathTracing(vector<Object*> &objects, int id, double &t, Ray &ray, Light &light,  int Depthi, const SampleMode sampleMode){
    // N is iteration times
    
    Color c (0,0,0), L_indirect (0,0,0);
   

    Object* hitObj = objects.at(id);

    if(id == objects.size() - 1){
       
        // The hit point is on light
      
        c = c + hitObj->getRadiance() ;
        
    }
    else{
    if(Depthi < 2)
    {
        // The hit point is on another object
        
        // ω’ = random direction in hemisphere above n;
        // Li += brdf * shade(trace(x, ω’)) * dot(n, ω’) / (p(ω’));
        Vector3 intersectionPoint = ray.getOrigin() + ray.getDirection() * t;
        Vector3 normal = hitObj->getNormalAt(intersectionPoint);
        if(hitObj->insideObject(ray.getOrigin())){
            normal = normal * -1;
        }
      
        double offset = 1e-9;
        Vector3 offsetIntersectionPoint =
        (ray.getDirection().dot(normal) < 0) ?
        intersectionPoint + normal * offset :
        intersectionPoint - normal * offset ;
        
        Vector3  w = getSample((Sphere&)objects, offsetIntersectionPoint, sampleMode, normal, 0.0);
        
        Ray second_ray (offsetIntersectionPoint,w);
        int IndexOfClosestObject;
        double hit_t;
        bool MissScene = !ClosestIntersectionObject(objects, second_ray, hit_t, IndexOfClosestObject);
        
        double throughput = 1;
        
        if(!MissScene){

            double cosine = clamp(normal.dot(w));
            double brdf = 1 / M_PI;
            double pdf = 1 / (2 * M_PI);
            switch (sampleMode) {
                case Uniform_Hemisphere_Solid_Angle:
                    
                    throughput *=  brdf * cosine / pdf;
                    
                    break;
                case Cosine_Weighted:
                    
                    pdf = 1 / M_PI;
                    throughput *= brdf / pdf;
                    
                    break;
                   
                    
                default:
                    
                    cout<< " no such sampling strategy! " << endl;
                    
                    break;
            }
            
            L_indirect = hitObj->getObjectColor() * getPixelColor_ImplicitPathTracing(objects, IndexOfClosestObject, hit_t, second_ray, light, ++Depthi, sampleMode) * throughput;
           
        }
        
        c = c + L_indirect;

    }
    else{
        Vector3 hitpoint ;
        c = c + getPixelColorWithSphereLight(objects, id, t, ray, hitpoint, 100, sampleMode);
    }
    }


    return c;
    
}

Color getPixelColor_VolumeSingleScatter(vector<Object*> &objects, int id, double &t_MAX, Ray &ray, int STEPS, Volume* volume, Light& light){
   
    Color pixelColor;
    Color Li(0,0,0);
  
    Sphere* hitObj = (Sphere*)objects.at(id);
    Color background = hitObj->getObjectColor();
    
    Vector3 hitPoint = ray.getOrigin() + ray.getDirection() * t_MAX;
    
    double transmittance = volume->getTransmittanceHomogeneos(ray.getOrigin(), hitPoint);
    Li = Li +  background * transmittance;
    
    double stepSize = t_MAX / STEPS;
    double t_march = 0;
    Vector3 x_march = ray.getOrigin();
    double transmittanceFromEye = 0.0;

    // ray marching in the volume before hit the surface
    for(int i = 0 ; i < STEPS - 1; i++){
        t_march += stepSize;
        x_march = ray.getOrigin() + ray.getDirection() * t_march;
        transmittanceFromEye = volume->getTransmittanceHomogeneos(ray.getOrigin(), x_march);
        double sigmaS = volume->getSigmaS();
        double pf = volume->getPhaseFuctionValue();
        
        //Vector3 lightHitPoint;
       // Color Le = getPixelColorWithSphereLight(objects, -1, t_march, ray, lightHitPoint, 1, Importance_Solid_Angle);
      
        Color Le = getPixelColorPointLight(objects, -1, t_march, ray, light);
        double transmittanceFromLight = volume->getTransmittanceHomogeneos(x_march, light.getLightPosition());
        
        Li = Li + Le * transmittanceFromEye * sigmaS * pf * transmittanceFromLight * stepSize;
    
    }

   return Li;
}
    
Color getPixelColor_VolumeSingleScatter_point(vector<Object*> &objects, int id, double &t_MAX, Ray &ray, int STEPS, Volume* volume, Light& light){
    
    Color pixelColor;
    Color Li(0,0,0);

    
   // Sphere* hitObj = (Sphere*)objects.at(id);
    
    Color background =getPixelColorPointLight(objects, id, t_MAX, ray, light);
   // Vector3 hitPoint = ray.getOrigin() + ray.getDirection() * t_MAX;

  
    double transmittance = volume->getTransmittanceHomogeneos(t_MAX);
    //double transmittance = volume->getTransmittanceHomogeneos(ray.getOrigin(), hitPoint);
    Li =  background * transmittance;


    double stepSize = t_MAX / STEPS;
    double t_march = 0;
    Vector3 x_march;
    double transmittanceFromEye = 0.0;
    
    // ray marching in the volume before hit the surface
    for(int i = 0 ; i < STEPS - 1; i++){
        t_march += stepSize;
        x_march = ray.getOrigin() + ray.getDirection() * t_march;
        //transmittanceFromEye = volume->getTransmittanceHomogeneos(ray.getOrigin(), x_march);
        transmittanceFromEye = volume->getTransmittanceHomogeneos(t_march);

        double sigmaS = volume->getSigmaS();
        double pf = volume->getPhaseFuctionValue();
        Color Le = getPixelColorPointLight(objects, -1, t_march, ray, light);
        double transmittanceFromLight = volume->getTransmittanceHomogeneos(x_march, light.getLightPosition());
        
        Li = Li + Le * transmittanceFromEye * sigmaS * pf * transmittanceFromLight * stepSize;
        
    }
    
    return Li;
}
Color getPixelColor_VolumeSingleScatter_point2(vector<Object*> &objects, int id, double &t_MAX, Ray &ray, int MCsamples, Volume* volume, Light& light){
    
    Color pixelColor;
    Color Li(0,0,0), Ls(0,0,0);
    double sigmaT = volume->getSigmaT();
    double sigmaS = volume->getSigmaS();
    Vector3 lightPos = light.getLightPosition();
//        if(id == objects.size() - 1){
//            cout<<endl;
//        }
    Color background =getPixelColorPointLight(objects, id, t_MAX, ray, light);
    double transmittance = volume->getTransmittanceHomogeneos(t_MAX);
    Li =  background * transmittance;

    sampleMode1D sampleMode = Distance;
    for(int i = 0 ; i < MCsamples; i++){
   
        double pdf;
        double tempDist;

        double t_sample = sample1DDBG(ray, lightPos, t_MAX, sigmaT, pdf, sampleMode,tempDist);
        Vector3 x_sample = ray.getOrigin() + ray.getDirection() * t_sample;
        double transmittanceFromEye = volume->getTransmittanceHomogeneos(t_sample);
        
        double pf = volume->getPhaseFuctionValue();
        Color Le = getPixelColorPointLight(objects, -1, t_sample, ray, light);
        double transmittanceFromLight = volume->getTransmittanceHomogeneos(x_sample, light.getLightPosition());
//   
//        cout<<tempDist<<endl;
//
//        cout<<"--"<<endl;

        Ls = Ls + Le * transmittanceFromEye * sigmaS * pf * transmittanceFromLight / pdf;
        
    }
    Li = Li + Ls / MCsamples;
    
    return Li;
}
Color getPixelColor_VolumeSingleScatter_sphere(vector<Object*> &objects, int id, double &t_MAX, Ray &ray, int STEPS, Volume* volume, Light& light){
    
    Color pixelColor;
    Color Li(0,0,0);
    

    if(id == objects.size() - 1){
        cout<<endl;
    }
    Vector3 hitPoint;
    Color background = getPixelColorWithSphereLight(objects, id, t_MAX, ray, hitPoint,1,Importance_Solid_Angle);

    //Vector3 hitPoint = ray.getOrigin() + ray.getDirection() * t_MAX;
    double transmittance = volume->getTransmittanceHomogeneos(t_MAX);
    //double transmittance = volume->getTransmittanceHomogeneos(ray.getOrigin(), hitPoint);
    Li =  background * transmittance;

       
        double stepSize = t_MAX / STEPS;
        double t_march = 0;
        Vector3 x_march;
        double transmittanceFromEye = 0.0;
        
        // ray marching in the volume before hit the surface
        for(int i = 0 ; i < STEPS - 1; i++){
            t_march += stepSize;
            x_march = ray.getOrigin() + ray.getDirection() * t_march;
            //transmittanceFromEye = volume->getTransmittanceHomogeneos(ray.getOrigin(), x_march);
            transmittanceFromEye = volume->getTransmittanceHomogeneos(t_march);
            
            double sigmaS = volume->getSigmaS();
            double pf = volume->getPhaseFuctionValue();
            Vector3 lightHitPoint;
            Color Le = getPixelColorWithSphereLight(objects, -1, t_march, ray, lightHitPoint, 1, Importance_Solid_Angle);
            double transmittanceFromLight = volume->getTransmittanceHomogeneos(x_march, lightHitPoint);
            
            Li = Li + Le * transmittanceFromEye * sigmaS * pf * transmittanceFromLight * stepSize;
            
        }
     
    return Li;
}
Color getPixelColor_VolumeSingleScatter_sphere2(vector<Object*> &objects, int id, double &t_MAX, Ray &ray, int MCsamples, Volume* volume, Light& light){
    
    Color pixelColor;
    Color Li(0,0,0), Ls(0,0,0);
    double sigmaT = volume->getSigmaT();
    double sigmaS = volume->getSigmaS();
    Vector3 lightPos = light.getLightPosition();
    Vector3 hitPoint;
//            if(id == objects.size() - 1){
//                cout<<endl;
//            }
    Color background = getPixelColorWithSphereLight(objects, id, t_MAX, ray, hitPoint,1,Importance_Solid_Angle);
    double transmittance = volume->getTransmittanceHomogeneos(t_MAX);
    Li =  background * transmittance;
    
    sampleMode1D sampleMode = EquiAngular;
    for(int i = 0 ; i < MCsamples; i++){
        
        double pdf;
        double tempDist;
        
        double t_sample = sample1DDBG(ray, lightPos, t_MAX, sigmaT, pdf, sampleMode,tempDist);
        Vector3 x_sample = ray.getOrigin() + ray.getDirection() * t_sample;
        double transmittanceFromEye = volume->getTransmittanceHomogeneos(t_sample);
        
        double pf = volume->getPhaseFuctionValue();
        Color Le = getPixelColorWithSphereLight(objects, -1, t_sample, ray, hitPoint,1,Importance_Solid_Angle);
        double transmittanceFromLight = volume->getTransmittanceHomogeneos(x_sample, light.getLightPosition());
        //
        //        cout<<tempDist<<endl;
        //
        //        cout<<"--"<<endl;
        
        Ls = Ls + Le * transmittanceFromEye * sigmaS * pf * transmittanceFromLight / pdf;
        
    }
    Li = Li + Ls / MCsamples;
    
    return Li;
}

//Color vPT(x, ω)
//tmax = nearestSurface(x, ω)
//t = -log(1 - rand()) / σt // Sample free path
//if t < tmax: // Volume interaction
//x += t * ω
//pdf_t = σt * exp(-σt * t)
//(ω’, pdf_ω’) = samplePF(ω)
//return σs * Tr(t) * PF(ω) / (pdf_t * pdf_ω’) * vPT(x,ω’)
//else: // Surface interaction
//x += tmax * ω
//Pr_tmax = exp(-σt * tmax)
//(ω’, pdf_ω’) = sampleBRDF(n, ω)
//return Tr(tmax) * BRDF(ω, ω’) / (Pr_tmax * pdf_ω’) * vPT(x, ω’)


void renderDBG(const int width, const int height, vector<Object*> objects, const char *filename, Camera& Cam, Light& light, Volume* volume,const int Samples, const SampleMode mode, int AASample, AASamplerStrategy aaSampleStrategy){
    
  
    int n = width * height;
    Color *pixels = new Color[n];

    
    int idx;
    Color pixelColor;
    
    // Anti-aliasing
    
    // sample of each pixel
    int AAdepth = AASample;
    int SampleSum = AAdepth * AAdepth;
    Sampler sampler(width, height, AAdepth);
    Color sumColor (0,0,0);
    
    #pragma omp parallel for schedule(dynamic, 1) private(pixelColor)
    
    
    for (int y = 0; y < height; y++){
        
        fprintf(stderr,"\r%5.2f%%",100.*y/(height-1));
 
        for(int x = 0; x < width; x++){
            
            idx = (height - y - 1) * width + x;
            
            sumColor = Color(0,0,0);
//            if(x == 255 && y == 289){
//                
//                cout<<"dd"<<endl;
//            }
            

            double sampleX = double(x), sampleY = double(y);
            
            Sample s(sampleX,sampleY);
            vector<Sample> samples;
            sampler.getSample(&s, samples, aaSampleStrategy);
        
            for (int current = 0; current < samples.size(); current ++) {
                Vector3 camRayOrigin = Cam.getCamPosition();
                
                sampleX = samples.at(current).getX();
                sampleY = samples.at(current).getY();
                Vector3 camRayDirection = Cam.getCamUp() * ( sampleX /width - .5) +  Cam.getCamRight() * ( sampleY /height - .5) + Cam.getCamDirection();

                camRayDirection = camRayDirection.norm();
               
                double motionTime = randf();
                Ray cameraRay(camRayOrigin,camRayDirection,motionTime);
                 //Ray cameraRay(camRayOrigin,camRayDirection);
                
                double t;
                int IndexOfClosestObject;
            
                if(ClosestIntersectionObject(objects, cameraRay, t, IndexOfClosestObject)){
                    Vector3 hitPoint;
                   // pixelColor = getPixelColor_ImplicitPathTracing(objects,IndexOfClosestObject, t,cameraRay, light, 0, mode);
                    pixelColor = getPixelColorPointLight(objects, IndexOfClosestObject, t, cameraRay, light);
                    //pixelColor = getPixelColor_VolumeSingleScatter_point2(objects, IndexOfClosestObject, t, cameraRay, 8, volume,light);
                    //pixelColor = getPixelColorWithSphereLight(objects,IndexOfClosestObject, t,cameraRay, hitPoint, Samples, mode);
                    //pixelColor = getPixelColorWithSphereLightDBG(objects, IndexOfClosestObject, t, cameraRay, light, Samples, mode);
                    //pixelColor = getPixelColorAmbientOcclusion(objects,IndexOfClosestObject, t,cameraRay, light, Samples, mode);
                }
                else{
                    pixelColor = Color(1,1,1);
                }
                
                sumColor = sumColor + pixelColor;

            }
           
        pixels[idx] = sumColor / (SampleSum) ;
    }
}

    
    
    savePPM(filename, width, height, pixels);
    delete []pixels;



}


int main(){
    cout << "rendering..." << endl;
    

    int width = 512, height  = 384;
    double AspectRatio = double(width)/double(height);
    double TanHalfFOV = 0.57735 ;
    
    Vector3 CamPos (50, 50, 275.0);
    //Vector3 CamPos (0, 0, 200.0);
    Vector3 LookAt (0, -0.05, -1);
    //Vector3 LookAt (0, 0, -1);
    Vector3 CamDir  = (LookAt).norm();
    Vector3 CamUp ( AspectRatio * TanHalfFOV, 0., 0.);
    Vector3 CamRight  = (CamUp.cross(CamDir)).norm() * TanHalfFOV;
    
    Camera Cam (CamPos,CamDir,CamUp,CamRight);

    
    Vector3 light_position(50, 68.6 - .27, 81.6);
    
    double radiance = 40000.0;
    Light light(light_position, Color(1,1,1),radiance);
    
    Volume homogeneousVolume(0.0001,0.003);
    
    Vector3 sphereO0  =  light_position;
    Vector3 sphereO1(1e5 + 1, 40.8, 81.6);
    Vector3 sphereO2(-1e5 + 99, 40.8, 81.6);
    Vector3 sphereO3(50, 40.8, 1e5);
    Vector3 sphereO4(50, 1e5, 81.6);
    Vector3 sphereO5(50, -1e5 + 81.6, 81.6);
    Vector3 sphereO6(27, 16.5, 47);
    Vector3 sphereO7(73, 16.5, 78);
    Vector3 sphereO8(50, 68.6 - .27, 81.6);
    Vector3 sphereO9(50, -1e5, 81.6);
    
    
    Vector3 bigSphereVolume((1e3)/2, 40.8, 81.6);
    
    
    Color c1 (.75, .25, .25);
    Color c2 (.25, .25, .75);
    Color c3(.75, .75, .75);
    Color c4(.999, .999, .999);
    Color c5(1,1,1);
    Color black(0,0,0);
    Vector3 dir1 (5,0,0);
    Vector3 radiance1 (10,10,10);
    Vector3 radiance2 (0,0,0);
    
    double exp = 2;
    Sphere spheres[] = {
        
//        
//        Sphere(1e5,  sphereO1,  c1,exp,radiance2, DIFFUSE),
//        Sphere(1e5,  sphereO2,  c2,exp,radiance2, DIFFUSE),
//        Sphere(1e5,   sphereO3,  c3,exp,radiance2, DIFFUSE),
//        Sphere(1e5,   sphereO4,  c3,exp,radiance2, DIFFUSE),
//        Sphere(1e5,  sphereO5,  c3,exp,radiance2, DIFFUSE),
        Sphere(16.5, sphereO6,  c4, exp,radiance2, DIFFUSE),
       Sphere(16.5, sphereO7,  c4, exp,radiance2, DIFFUSE),
        
        Sphere(1e3, bigSphereVolume,  black, exp, radiance1, DIFFUSE),
        Sphere(2,   sphereO8,  c5, exp, radiance1, DIFFUSE),
        
        
    };
    MotionSphere mSphere[2]
       = {
        MotionSphere(16.5, sphereO6,  c4, dir1, 0.0),
        MotionSphere(16.5, sphereO7,  c4, dir1, 0.0),
      };

//    bool temp1 = spheres[0].insideObject(CamPos);
//    bool temp2 = spheres[0].insideObject(spheres[1].getSphereCenter());
      vector<Object*> objects;
//    for (int i = 0; i< int(sizeof(mSphere)/sizeof(MotionSphere)); i++){
//       objects.push_back(&mSphere[i]);
//    }
    for (int i = 0; i< int(sizeof(spheres)/sizeof(Sphere)); i++){
        objects.push_back((Object*)&spheres[i]);
    }
    renderDBG(width,height, objects, "D:/Yue/VolumeSingleScattering/dist-step8.ppm", Cam, light, &homogeneousVolume,1, Importance_Solid_Angle,4, Jittered);

    return 0;

  
}
