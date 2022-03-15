//
//  main.cpp
//  Minkowski
//
//  Created by Jean-Marie Mirebeau on 02/09/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#include <iostream>
#include "Meissner_3_Test.h"
#include "Meissner_3.h"
#include "Meissner_2_Test.h"
#include "Meissner_2.h"
#include "ConstraintsProduct.h"
//#include "Minkowski_3_Test.h"
//#include "Minkowski_3.h"
#include "Minkowski_2_Test.h"
#include "Minkowski_2.h"
#include "QuotientedNewton.h"

int main(int argc, const char * argv[]) {
    Minkowski_2_Test::Test2(); 
//    Minkowski_3_Test::Test1();
    
//    Minkowski_3_Test::MinkowskiInterpolate({&Minkowski_3_Test::octahedronEqns, &Minkowski_3_Test::tetraEqns},30);
    
//    Meissner_2_Test::Test1(200,0.05);
    
/*
    if(argc>1){
        const int i= atoi(argv[1]);
        std::cout << i << std::endl;
        
        double stopVal = 0.4197;
        if(argc>2) stopVal = atof(argv[2]);
        
        if(2<=i && i<=10)
            Meissner_3_Test::Test1(i,stopVal);
    }
 */
	
//    for(int i=5; i<10; ++i)
//        Meissner_3_Test::Test1(i);
    
    std::cout << "Hello, World!\n";
    return 0;
}
