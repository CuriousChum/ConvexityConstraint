//
//  DGO_Traits.h
//  CGalTest
//
//  Created by Jean-Marie Mirebeau on 29/01/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef CGalTest_DGO_Traits_h
#define CGalTest_DGO_Traits_h

/*
 Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.
 
 File     : dgopt.h
 */

#include "mosek.h"

typedef void *dgohand_t;

MSKrescodee
MSK_dgoread(MSKtask_t task,
            char      *nldatafile,
            MSKintt   *numvar,       /* numterms in primal */
            MSKintt   *numcon,       /* numvar in primal */
            MSKintt   *t,            /* number of constraints in primal*/
            double    **v,           /* coiefients for terms in primal*/
            MSKintt   **p            /* corresponds to number of terms
                                      in each constraint in the
                                      primal */
);

MSKrescodee
MSK_dgosetup(MSKtask_t task,
             MSKintt   numvar,
             MSKintt   numcon,
             MSKintt   t,
             double    *v,
             MSKintt   *p,
             dgohand_t *nlh);

MSKrescodee
MSK_freedgo(MSKtask_t task,
            dgohand_t  *nlh);


#endif
