//
//  Sampler.h
//  RayTracer
//
//  Created by WangYue on 2016-10-30.
//  Copyright © 2016 WangYue. All rights reserved.
//

#ifndef Sampler_h
#define Sampler_h
#include "Random.h"
#include "Sample.h"
enum AASamplerStrategy{
    Jittered, Random, Halton, Poisson, Hammersley
};

class Sampler {

    
    int xPixelEnd;
    int yPixelEnd;
    int current;
    int Adepth;
    bool first = true;
    std::vector<Sample> sampleBuf;
    
public:
    Sampler (int x, int y)
    {
        xPixelEnd = x;
        yPixelEnd = y;
        current = 0;
    };
    
    Sampler (int x, int y, int sample)
    {
        xPixelEnd = x;
        yPixelEnd = y;
        current = 0;
        Adepth = sample;
    };

    
    
    inline double RadicalInverse(int n, int base){
        double val = 0;
        double invBase = 1. / base, invBi = invBase;
        int d_i = 0;
        while (n > 0) {

            d_i = (n % base);
            val += d_i * invBi;
            n *= invBase;
            invBi *= invBase;
        }
        return val;
    }
    
    double Lerp(double t, double a, double b){
        return (1-t)*a + t*b;
    }
    
    void HaltonSample(Sample* sample, std::vector<Sample> &samples){
        
        samples.clear();
        int k = (sample->getY() * yPixelEnd + sample->getX()) * Adepth * Adepth;
        
        if(first){
            first = false;
            for(int ax = 0; ax < Adepth; ax++ ){
                for (int ay = 0; ay < Adepth; ay++){
                    
                    // (k/n, phi_2(k))
                    double xOffset = RadicalInverse(k, 2);
                    double yOffset = RadicalInverse(k++, 3);
                    
                    //std::cout<<xOffset <<" "<<yOffset <<std::endl;
                    
                    Sample s(sample->getX() + xOffset, sample->getY() + yOffset);
                    Sample offset(xOffset, yOffset);
                    samples.push_back(s);
                    sampleBuf.push_back(offset);
                    
                }
            }
        }
        else{
        
            for(int ax = 0; ax < Adepth; ax++ ){
                for (int ay = 0; ay < Adepth; ay++){
                    
                    Sample s(sample->getX() + sampleBuf.at(k).getX(), sample->getY() + sampleBuf.at(k++).getY());
                    
                    samples.push_back(s);
                    
                }
            }
        }
        
    }
    
//    void PoissonSample(Sample* sample, std::vector<Sample> &samples){
//        
//        samples.clear();
//        int k = 0;
//        for(int ax = 0; ax < Adepth; ax++ ){
//            for (int ay = 0; ay < Adepth; ay++){
//                
//                // (k/n, phi_2(k))
//                double xOffset = RadicalInverse(k, 2);
//                double yOffset = RadicalInverse(k++, 3);
//                
//                // std::cout<<xOffset <<" "<<yOffset <<std::endl;
//                
//                Sample s(sample->getX() + xOffset, sample->getY() + yOffset);
//                samples.push_back(s);
//                
//            }
//        }
//        
//    }
    void HammersleySample(Sample* sample, std::vector<Sample> &samples){
        
        samples.clear();
        int N = Adepth * Adepth;
        int k = 0;
        for(int ax = 0; ax < Adepth; ax++ ){
            for (int ay = 0; ay < Adepth; ay++){
                
                // (k/n, phi_2(k))
                double xOffset = double(k) / N;
                double yOffset = RadicalInverse(k++, 2);
                
                //std::cout<<xOffset <<" "<<yOffset <<std::endl;
                
                Sample s(sample->getX() + xOffset, sample->getY() + yOffset);
                samples.push_back(s);
                
            }
        }
        
    }
    
    void JitteredSample(Sample* sample, std::vector<Sample> &samples){
      
        samples.clear();
        
        for(int ax = 0; ax < Adepth; ax++ ){
            for (int ay = 0; ay < Adepth; ay++){
                
                double xOffset = double(ax) / Adepth;
                double yOffset = double(ay) / Adepth;
           
                // Jittered
                double xOffsetRandom = randf() * xOffset;
                double yOffsetRandom = randf() * yOffset;

                Sample s(sample->getX() + xOffset + xOffsetRandom , sample->getY() + yOffset + yOffsetRandom);
                samples.push_back(s);
                
            }
        }
        
    }
    void getRandomSample(Sample* sample, std::vector<Sample> &samples){
        
        samples.clear();
        
        for(int ax = 0; ax < Adepth; ax++ ){
            for (int ay = 0; ay < Adepth; ay++){
                
                double xOffset = randf();
                double yOffset = randf();
                
                Sample s(sample->getX() + xOffset  , sample->getY() + yOffset);
                samples.push_back(s);
                
            }
        }
        
    }
    
    void getSample(Sample* sample, std::vector<Sample> &samples, AASamplerStrategy strategy){
        
        switch (strategy){
            case Jittered:
                JitteredSample(sample, samples);
                break;
            case Hammersley:
                HammersleySample(sample, samples);
                break;
            case Halton:
                HaltonSample(sample, samples);
                break;
            default:// Random
                getRandomSample(sample, samples);
        }
    }
    
};
#endif /* Sampler_h */
