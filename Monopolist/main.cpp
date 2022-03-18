

#include <string>
#include <iostream>

//#include "ConvexityConstraint_Test.hxx"


//#include "Regular_Test.h"

//#include "DGO_CallTest.h"
/*
#include "NewtonSolvers_Test.h"
#include "PrincipalAgent_Test.h"
#include "Geometry_2.h"
*/

/*
#include "LipschitzConstraint.h"
#include "ConvexityConstraint.h"
#include "PrincipalAgent.h"
#include "PrincipalAgent_Test.h"
#include "PrincipalAgent_NLOpt.h"
 */


#include "Geometry_3.h"
#include "PrincipalAgent_3.h"
#include "ConvexityConstraint_3.h"
#include "PrincipalAgent_3_Test.h"



/*
#include "Minkowski_3_Traits.h"
#include "Minkowski_3_Test.h"
#include "MinkowskiConstraint_3.h"
*/

/*#include "MinkowskiConstraint_2.h"
#include "Minkowski_2_Test.h"
#include "MinkowskiEnergy_2.h"
 */

//#include "Refinement_2_Traits.h"

int main(){
    
//    Minkowski_2_Test::Test1();
//    Refinement_2_Test::Test2();
    
//    Minkowski_2_Test::Test0();
    
//    Minkowski_3::Test3();
    
//    PrincipalAgent_3_Test::PA0(20,1e-4);
    
    
//    PrincipalAgent_3_Test::Octa0();
//        PrincipalAgent_3_Test::Tetra0();
    
//    PrincipalAgent_Test::PA0(30 , PrincipalAgent_Test::ShapeType::Triangle);
    
/*    PrincipalAgent_Test::PA1(PrincipalAgent_Test::MakeShape(5 , PrincipalAgent_Test::ShapeType::Square),
                             8, "PA1_Triangle_2");*/

    
/*    PrincipalAgent_Test::PALinear(PrincipalAgent_Test::MakeShape(10,PrincipalAgent_Test::ShapeType::Square),
                             1, "PARef_Square_Linear");*/

/*    PrincipalAgent_Test::PALinear2(PrincipalAgent_Test::MakeShape(10,PrincipalAgent_Test::ShapeType::Square),
                                  1, "PARef2_Square_Linear");*/
    
//    PrincipalAgent_Test::PA_NLOpt_0(10);
//    PrincipalAgent_Test::PA0();


    
    
//    PrincipalAgent_Test::Ref0(6);
//    PrincipalAgent_Test::Ref0(3,3);
    
//    PrincipalAgent_Test::Opt0();
  
//    CGalTest_Manual();
//    CGal_MyTest();
//    ConvexityConstraint({}, {});
//    ConvexityConstraint_Test::Pos0();
//    ConvexityConstraint_Test::Cvx0();
    
/*    std::pair<int,int> myPair{1,2};
    std::cout << std::get<1>(myPair) << std::endl;*/
    
//    ConvexityConstraint_Test::Pos0();
    
//    ConvexityConstraint_Test::Cvx0();
    /*
    const std::string dgo_mps = "/Users/mirebeau/mosek/7/tools/examples/data/dgo.mps";
    const std::string dgo_f = "/Users/mirebeau/mosek/7/tools/examples/data/dgo.f";
    const int argc=3;
    const char *argv[argc]={nullptr,dgo_mps.c_str(),dgo_f.c_str()};
    dgo_main(argc, const_cast<char**>(argv));
    */
    return EXIT_SUCCESS;
}

/*
 // Basic printing is not ok ...
 
 Geometry_3::CGT::Vector v(1,2,3);
 std::cout << v << "\n";
 
 Geometry_3::NS::VectorType w(2);
 w[0]=1; w[1]=2;
 std::cout << w << "\n";
 */

