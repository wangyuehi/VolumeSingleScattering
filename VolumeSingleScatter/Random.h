//
//  Random.h
//  RayTracer
//
//  Created by WangYue on 2016-11-02.
//  Copyright © 2016 WangYue. All rights reserved.
//

#ifndef Random_h
#define Random_h
#include <boost/random.hpp>
#include <boost/format.hpp>
double randf(void)
{
    boost::mt19937 rng(67);
    static boost::uniform_01<boost::mt19937> zeroone(rng);
    return zeroone();
}

#endif /* Random_h */
