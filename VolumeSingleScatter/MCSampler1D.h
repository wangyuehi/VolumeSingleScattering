//
//  MCSampler1D.h
//  VolumeSingleScatter
//
//  Created by WangYue on 2016-11-30.
//  Copyright © 2016 WangYue. All rights reserved.
//

#ifndef MCSampler1D_h
#define MCSampler1D_h

enum sampleMode1D{
    Distance,
    EquiAngular
};

double EquiAngularSampling(Ray& ray, Vector3& lightPos, double t_MAX,double sigmaT, double& pdf,double& temp ){
    double t_sample;

    
    // PDF: p(t) = D / ((theta_b − theta_a)(D2 +t2));
    // t(x) = D * tan(1− randf()) * theta_a + randf() * theta_b;
    // theta_x = atan(x/D);
    
    Vector3 rayOrigin = ray.getOrigin();
    Vector3 rayDir = ray.getDirection();
    // get coord of closest point to light along (infinite) ray
    double delta = (lightPos - rayOrigin).dot(rayDir);
    
    // get distance this point is from light
    double D = (rayOrigin + rayDir * delta - lightPos).length();
    
    // get angle of endpoints
    double thetaA = - atan(delta / D);
    double thetaB = atan((t_MAX - delta) / D);
    
    // take sample
    double u = randf();
    double t = D * tan(thetaA *(1 - u)+ thetaB * u);
    
//    // debug
    temp = (D*D + t*t);
    
    t_sample = delta + t;
    pdf = D/((thetaB - thetaA)*(D*D + t*t));
    
    return t_sample;
};

double distanceSampling(double t_MAX,double sigmaT, double& pdf){
    double t_sample;
    
    // PDF: p(t) = (sigmaT / (1 - exp(sigmaT * t_MAX))) * exp(-sigmaT * t);
    // CDF: P(t) = (1 - exp(-sigmaT * t))/(1 - exp(-sigmaT * t_MAX));
    // INVERSE: t = - ln(1 - randf() * (1 - exp(sigmaT * t_MAX))) / sigmaT;
    
    t_sample = - log(1 - randf() * (1 - exp(- sigmaT * t_MAX))) / sigmaT;
    pdf = sigmaT * exp( - sigmaT * t_sample) / (1 - exp( - sigmaT * t_MAX));

    return t_sample;
};

//double sample1D(Ray& ray, Vector3& lightPos, double t_MAX,double sigmaT, double& pdf,sampleMode1D sampleMode){
//    double t_sample;
//    
//    switch (sampleMode) {
//        case EquiAngular:
//                t_sample = EquiAngularSampling(ray, lightPos, t_MAX, sigmaT, pdf);
//            break;
//            
//        default:
//            t_sample = distanceSampling(t_MAX, sigmaT,  pdf);
//    }
//    
//
//    
//    return t_sample;
//    
//};
double sample1DDBG(Ray& ray, Vector3& lightPos, double t_MAX,double sigmaT, double& pdf,sampleMode1D sampleMode, double& temp){
    double t_sample;
    
    switch (sampleMode) {
        case EquiAngular:
            t_sample = EquiAngularSampling(ray, lightPos, t_MAX, sigmaT, pdf, temp);
            break;
            
        default:
            t_sample = distanceSampling(t_MAX, sigmaT,  pdf);
    }
    
    
    
    return t_sample;
    
};
#endif /* MCSampler1D_h */
