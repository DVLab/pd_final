#include "GlobalPlacer.h"
#include "ExampleFunction.h"
#include "NumericalOptimizer.h"
#include <stdlib.h>

GlobalPlacer::GlobalPlacer(Placement &placement,LayerMgr &layer)
	:_placement(placement),_layer(layer)
{

}


void GlobalPlacer::place()
{
	///////////////////////////////////////////////////////////////////
	// The following example is only for analytical methods.
	// if you use other methods, you can skip and delete it directly.
	//////////////////////////////////////////////////////////////////

	ExampleFunction ef(_placement,_layer); // require to define the object function and gradient function

    vector<double> x(2); // solution vector, size: num_blocks*2 
                         // each 2 variables represent the X and Y dimensions of a block
    x[0] = 100; // initialize the solution vector
    x[1] = 100;

    for(unsigned i = 0;i<_placement.numModules();++i)
      _layer.addModule(0,&_placement.module(i));

    NumericalOptimizer no(ef);
    no.setX(x); // set initial solution
    no.setNumIteration(35); // user-specified parameter
    no.setStepSizeBound(5); // user-specified parameter
    no.solve(); // Conjugate Gradient solver

    cout << "Current solution:" << endl;
    for (unsigned i = 0; i < no.dimension(); i++) {
        cout << "x[" << i << "] = " << no.x(i) << endl;
    }
    cout << "Objective: " << no.objective() << endl;
	////////////////////////////////////////////////////////////////


}
