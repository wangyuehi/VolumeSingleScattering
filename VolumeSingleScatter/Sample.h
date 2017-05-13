//
//  Sample.h
//  RayTracer
//
//  Created by WangYue on 2016-10-31.
//  Copyright © 2016 WangYue. All rights reserved.
//

#ifndef Sample_h
#define Sample_h

class Sample {
    double imageX, imageY;
public:
    Sample(double x, double y){
        imageX = x;
        imageY = y;
    }
    double getX(){
        return imageX;
    }
    double getY(){
        return imageY;
    }
};
#endif /* Sample_h */
