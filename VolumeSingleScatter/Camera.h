//
//  Camera.h
//  CameraTracer
//
//  Created by WangYue on 16/9/20.
//  Copyright © 2016年 WangYue. All rights reserved.
//

#ifndef Camera_h
#define Camera_h


#include "Vector3.h"

class Camera
{
    Vector3 camera_pos;
    Vector3 camera_dir;
    Vector3 camera_right;
    Vector3 camera_up;
    
public:
    
    
    Camera(const Vector3 &pos, const Vector3 &dir, const Vector3 &right, const Vector3 &up);
    Camera();
    
    Vector3 getCamPosition(){return camera_pos;};
    Vector3 getCamDirection(){return camera_dir;};
    Vector3 getCamRight(){return camera_right;};
    Vector3 getCamUp(){return camera_up;};
    
    
};

Camera::Camera(){
    camera_pos = Vector3(0,0,0);
    camera_dir = Vector3(0,0,1);
    camera_right = Vector3(0,0,0);
    camera_up = Vector3(0,0,0);
    
    
    
}

Camera::Camera(const Vector3 &p, const Vector3 &d,const Vector3 &up,const Vector3 &right){
    camera_pos = p;
    camera_dir = d;
    camera_up = up;
    camera_right = right;
}

#endif /* Camera_h */
