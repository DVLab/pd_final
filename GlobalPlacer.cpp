#include "GlobalPlacer.h"
#include "ExampleFunction.h"
#include "NumericalOptimizer.h"
#include <stdlib.h>
#include "Util.h"
#include <queue>
#include <fstream>


GlobalPlacer::GlobalPlacer(Placement &placement, LayerMgr &layer)
	:_placement(placement),_layer(layer)
{
}


void GlobalPlacer::place()
{
	///////////////////////////////////////////////////////////////////
	// The following example is only for analytical methods.
	// if you use other methods, you can skip and delete it directly.
	//////////////////////////////////////////////////////////////////
    double wirelength,wirelength_2=0.0,wirelength_3=0.0;
    priority_queue<double> Q_x;
    priority_queue<double> Q_y;
    double c_X = (_placement.boundryRight() + _placement.boundryLeft())/2;
    double c_Y = (_placement.boundryTop() + _placement.boundryBottom())/2;
    double d_X = _placement.boundryRight() - _placement.boundryLeft();
    double d_Y = _placement.boundryTop() - _placement.boundryBottom();
    double num = _placement.numModules();

    ofstream outfile_x("output_x");
    ofstream outfile_y("output_y");

    ExampleFunction ef(_placement); // require to define the object function and gradient function

    vector<double> x(2*num,0); // solution vector, size: num_blocks*2 
    vector<double> d(2*num,0);
    vector<double> P_x;
    vector<double> P_y;

    for(unsigned i = 0;i<num;++i){
        if(_placement.module(i).isFixed()){
            x[i] = _placement.module(i).x();
            x[i+num] = _placement.module(i).y();
        }
        else{
            x[i] = c_X;
            x[i+num] = c_Y;
        }
        _placement.module(i).setCenterPosition(x[i],x[i+num]);
    }
    wirelength = _placement.computeHpwl();

    for(unsigned i = 0;i<_placement.numNets();++i){
        double up=-10000,down=10000,right=-10000,left=10000;
        for(unsigned j=0;j<_placement.net(i).numPins();++j){
            double _x = _placement.net(i).pin(j).x();
            double _y = _placement.net(i).pin(j).y();
            if(_x > right)
                right = _x;
            if(_x < left)
                left = _x;
            if(_y > up)
                up = _y;
            if(_y < down)
                down = _y;
        }
        wirelength_2 += up-down+right-left;
    }

    NumericalOptimizer no(ef);

    no.setX(x); // set initial solution
    no.setNumIteration(1000/*1.5*ite/(num/st)*/); // user-specified parameter
    no.setStepSizeBound(10/*step*(num/st)*/); // user-specified parameter
    no.solve(); // Conjugate Gradient solver

    for(unsigned i = 0;i<num;++i){
        _placement.module(i).setCenterPosition(no.x(i),no.x(i+num));
        x[i] = no.x(i);
        x[i+num] = no.x(i+num);
    }

    ef.reset();
    no.setX(x); 
    no.setNumIteration(50000/*ite/(num/st)*/); 
    no.setStepSizeBound(100/*step*(num/st)*/); 
    no.solve();
    

    for(unsigned i = 0;i<num;++i){
        _placement.module(i).setCenterPosition(no.x(i),no.x(i+num));
        outfile_x<<no.x(i)<<endl;
        outfile_y<<no.x(i+num)<<endl;
    }

    outfile_x.close();
    outfile_y.close();


    cout << "Objective: " << no.objective() << endl;
    printf( "\nLast HPWL: %.0f\n",wirelength);
    printf( "\nLast HPWL_2: %.0f\n",wirelength_2);
    printf( "\nNew HPWL: %.0f\n",_placement.computeHpwl());

    for(unsigned i = 0;i<_placement.numNets();++i){
        double up=-1000000,down=1000000,right=-1000000,left=1000000;
        for(unsigned j=0;j<_placement.net(i).numPins();++j){
            double _x = _placement.net(i).pin(j).x();
            double _y = _placement.net(i).pin(j).y();
            if(_x > right)
                right = _x;
            if(_x < left)
                left = _x;
            if(_y > up)
                up = _y;
            if(_y < down)
                down = _y;
        }
        wirelength_3 += up-down+right-left;
    }

    printf( "\nNew HPWL_2: %.0f\n",wirelength_3);
	////////////////////////////////////////////////////////////////
}
