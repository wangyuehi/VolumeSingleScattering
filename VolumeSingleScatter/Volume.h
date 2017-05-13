//
//  Volume.h
//  VolumeSingleScatter
//
//  Created by WangYue on 2016-11-22.
//  Copyright © 2016 WangYue. All rights reserved.
//

#ifndef Volume_h
#define Volume_h

#include <stdio.h>
#include "Vector3.h"

//enum VolumeProperty{Isotropic, Anisotropic};

class Volume {
    double sigmaT;
    double sigmaA;
    double sigmaS;
    double pf;// phase function
   // VolumeProperty mode;
    
public:
    Volume(double sigmaa, double sigmas){
        sigmaA = sigmaa;
        sigmaS = sigmas;
        sigmaT = sigmaA + sigmaS;
        pf = 1 / M_PI;// Isotropic
    }
    
    double getSigmaS(){
        return sigmaS;
    }
    double getSigmaT(){
        return sigmaT;
    }
    double getTransmittanceHomogeneos(Vector3 x, Vector3 y){
        
        Vector3 temp = x - y;
        
        double distance = sqrt(temp.dot(temp));
        
        double transmittance = exp(- distance * sigmaT);
        
        return transmittance;
    }
    double getTransmittanceHomogeneos(double t){
        
        
       
        
        double transmittance = exp( - t * sigmaT);
        
        return transmittance;
    }
    double getPhaseFuctionValue(){
        return pf;
    }
};

#endif /* Volume_h */
