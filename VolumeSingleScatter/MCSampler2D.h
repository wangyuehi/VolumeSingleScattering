//
//  MCSampler2D.h
//  VolumeSingleScatter
//
//  Created by WangYue on 2016-11-30.
//  Copyright © 2016 WangYue. All rights reserved.
//

#ifndef MCSampler2D_h
#define MCSampler2D_h

enum SampleMode {
    Importance_Solid_Angle,
    Importance_Surface_Area,
    Uniform_Solid_Angle,
    Uniform_Area,
    Cosine_Weighted,
    Importance_Phong,
    Uniform_Hemisphere_Solid_Angle
};


void createCoordinateSystem(Vector3 &N, Vector3 &Nt, Vector3 &Nb)
{
    if (abs(N.x) > abs(N.y))
        Nt = Vector3(N.z, 0, -N.x) / sqrt(N.x * N.x + N.z * N.z);
    else
        Nt = Vector3(0, -N.z, N.y) / sqrt(N.y * N.y + N.z * N.z);
    Nb = N.cross(Nt);
    //    Nt = Vector3(0,1,0).cross(N).norm();
    //    Nb = N.cross(Nt).norm();
    
    
    
}

Vector3 getSampleImportanceSurfaceArea(Sphere &sphere,const Vector3 X){
    Vector3 Y;
    // Rejection sampling
    // Test Y is in the area
    Vector3 CY ;
    Vector3 YX ;
    double eta1, eta2,theta, phi;
    double r = sphere.getSphereRadius();
    Vector3 center = sphere.getSphereCenter();
    // When the cosine is bigger than 0, it is in the area
    // if not, generate a new point.
    
    
    do{
        eta1 = randf();
        eta2 = randf();
        
        theta = acos(1 - 2 * eta1);
        phi = 2 * M_PI * eta2;
        //double x =  2 * cos(2 * M_PI * eta2) * sqrt(eta1 * (1 - eta1));
        //double y =  2 * sin(2 * M_PI * eta2) * sqrt(eta1 * (1 - eta1));
        //double z =   (1 - 2 * eta1);
        
        // Y = Vector3(x, y,z) * r + center;
        //Y = Vector3( r * cos(theta) * sin(phi), r * sin(theta) * sin(phi), r * cos(phi)) ;
        Y = Vector3( r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta));
        // Now Y is on the sphere.
        
        // Test Y is in the area
        Y = Y + center;
        CY = sphere.getNormalAt(Y) ;
        
        YX = X - Y;
        
    }while(CY.dot(YX)<=eps);
    Vector3 direction = Y - X;
    return direction;
    
}
Vector3 CosWeightedRandomHemisphereDirection2( Vector3 n )
{
    double Xi1 = randf();
    double Xi2 = randf();
    
    double  theta = acos(sqrt(1.0-Xi1));
    double  phi = 2.0 * 3.1415926535897932384626433832795 * Xi2;
    
    double xs = sinf(theta) * cosf(phi);
    double ys = cosf(theta);
    double zs = sinf(theta) * sinf(phi);
    
    Vector3 y(n.x, n.y, n.z);
    Vector3 h = y;
    if (fabs(h.x)<=fabs(h.y) && fabs(h.x)<=fabs(h.z))
        h.x= 1.0;
    else if (fabs(h.y)<=fabs(h.x) && fabs(h.y)<=fabs(h.z))
        h.y= 1.0;
    else
        h.z= 1.0;
    
    
    Vector3 x = (h.cross(y)).norm();
    Vector3 z = (x.cross(y)).norm();
    
    Vector3 direction = x * xs+ y * ys + z * zs;
    return direction.norm();
}

Vector3 getSampleCosineWeighted(Sphere &sphere,const Vector3 X, Vector3 &normal,double sinThataMax){
    // Generate the direction according to cosine
    
    double eta1, eta2; // For angles Theta and Fi.
    double Fi;
    eta1 = randf(); // For theta
    eta2 = randf(); // For Fi
    
    Vector3 U,V,W;
    Fi = 2 * M_PI * eta2;
    
    // from reference
    double x = cos(Fi) * sqrt(eta1);
    double y = sin(Fi) * sqrt(eta1);
    double z = sqrt(1 - eta1);
    Vector3 sample (x,y,z);
    
    
    
    
    createCoordinateSystem(normal, W, V);
    // sample = V * sample.x + W * sample.y + normal * sample.z;
    sample = Vector3 (
                      sample.x * V.x + sample.y * W.x + sample.z * normal.x,
                      sample.x * V.y + sample.y * W.y + sample.z * normal.y,
                      sample.x * V.z + sample.y * W.z + sample.z * normal.z
                      ).norm();
    //    if(sample.dot(normal)<=0){
    //        //cout<<"hh"<<endl;
    //        sample = sample * (-1);
    //    }
    
    return sample;
}
Vector3 getSampleUniformSolidAngle(const Vector3 X, Vector3 &normal){
    // Sample the sphere generate the direction uniformly
    // Return the direction
    
    double eta1, eta2; // For angles Theta and Fi.
    
    Vector3 CY ;
    Vector3 YX ;
    Vector3 W, V;
    
    Vector3  SampleDirection,sample;
    
    eta1 = randf(); // For theta
    eta2 = randf(); // For Fi
    
    //double r = sphere.getSphereRadius();
    
    double x =  2 * cos(2 * M_PI * eta2) * sqrt(eta1 * (1 - eta1));
    double y =  2 * sin(2 * M_PI * eta2) * sqrt(eta1 * (1 - eta1));
    double z =   (1 - 2 * eta1);
    
    sample = Vector3 (x,y,z);
    
    createCoordinateSystem(normal, V, W);
    
    sample = Vector3 (
                      sample.x * V.x + sample.y * W.x + sample.z * normal.x,
                      sample.x * V.y + sample.y * W.y + sample.z * normal.y,
                      sample.x * V.z + sample.y * W.z + sample.z * normal.z
                      );
    //    sample = Vector3 (
    //                      sample.x * W.x + sample.y * V.x + sample.z * normal.x,
    //                      sample.x * W.y + sample.y * V.y + sample.z * normal.y,
    //                      sample.x * W.z + sample.y * V.z + sample.z * normal.z
    //                      );
    
    
    return sample.norm();
}
Vector3 getSampleUniformSolidAngleHemisphere(const Vector3 X, Vector3 &normal){
    
    Vector3 sample = getSampleUniformSolidAngle(X, normal);
    
    if(sample.dot(normal) < 0){
        sample = sample * -1;
    }
    return sample;
}

Vector3 getSampleImportanceSolidAngle(Sphere &sphere,const Vector3 X, Vector3 &normal, double sinThetaMax){
    // Randomly generate the direction from intersection point to the sphere area.
    double radius = sphere.getSphereRadius();
    Vector3 XC =sphere.getSphereCenter() - X;
    double distance2XC = XC.dot(XC);
    
    Vector3  SampleDirection, SampleDirectionTemp;
    Vector3 W = XC.norm();
    // Vector3 V,W;
    //createCoordinateSystem(normal, V, W);
    Vector3 V = W.cross(normal).norm();
    Vector3 U = V.cross(W).norm();
    
    double eta1 = randf(); // For theta
    double eta2 = randf(); // For Fi
    
    double Theta = acos(1 - eta1 + eta1 * sqrt(1 - (radius * radius / distance2XC)));
    double Fi = 2 * M_PI * eta2;
    SampleDirection =  U * (cos(Fi) * sin(Theta)) + V * sin(Fi) * sin(Theta) + W * cos(Theta);
    
    return SampleDirection.norm();
}

Vector3 getSampleImportancePhong(Sphere &sphere,const Vector3 X, Vector3 &normal){
    Vector3 V, W;
    double eta1 = randf(), eta2 = randf();
    int n = sphere.getExponent();
    double Theta = acos(pow(1.0-eta1, 1.0 / (n + 1.0)));//acos(sqrt(eta1));
    double Phi = 2 * M_PI * eta2;
    //    double phi = 2.0f * PI * xi.x;
    //    double theta = acos(pow(1.0f - xi.y, 1.0f/(n+1.0f)));
    //    return vec2(phi, theta);
    
    double x = cos(Phi) * sin(Theta);
    double y = sin(Phi) * sin(Theta);
    double z = cos(Theta);
    
    Vector3 sample = Vector3 (x,y,z);
    
    // need rotate??
    
    createCoordinateSystem(normal, V, W);
    
    sample = Vector3 (
                      sample.x * V.x + sample.y * W.x + sample.z * normal.x,
                      sample.x * V.y + sample.y * W.y + sample.z * normal.y,
                      sample.x * V.z + sample.y * W.z + sample.z * normal.z
                      );
    
    return sample.norm();
}



Vector3 getSample(Sphere &sphere,const Vector3 X, SampleMode mode, Vector3 &normal, double sinThetaMax){
    switch (mode) {
        case Cosine_Weighted:
            // if use the cosine weighted solid angle, sin(ThetaMax) need to be provided
            return getSampleCosineWeighted(sphere,X, normal,sinThetaMax);
            break;
        case Importance_Surface_Area:
            return getSampleImportanceSurfaceArea(sphere, X);
            break;
        case Uniform_Solid_Angle:
            return getSampleUniformSolidAngle( X, normal);
            break;
        case Uniform_Hemisphere_Solid_Angle:
            return getSampleUniformSolidAngleHemisphere(X, normal);
            break;
        case Importance_Phong:
            return getSampleImportancePhong(sphere, X, normal);
        default:// Importance solid angle
            return getSampleImportanceSolidAngle(sphere,X, normal, sinThetaMax);
            break;
            
            
    }
}

#endif /* MCSampler2D_h */
