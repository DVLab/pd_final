#include "GlobalPlacer.h"
#include "ExampleFunction.h"
#include "NumericalOptimizer.h"
#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>
#include <cassert>
#define TSV_SIZE 50

vector<Placement> pLayer;

GlobalPlacer::GlobalPlacer(Placement &placement,LayerMgr &layer)
	:_placement(placement),_layer(layer)
{

}

string int2str(int i) {
  string s;
  stringstream ss(s);
  ss << i;
  return ss.str();
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
	_placement.setBoundary(
	   _placement.boundryLeft()/_layer.getLayerCount(),
	   _placement.boundryBottom()/_layer.getLayerCount(),
	   _placement.boundryRight()/_layer.getLayerCount(),
	   _placement.boundryTop()/_layer.getLayerCount()); 

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
	pLayer.resize(_layer.getLayerCount());
//	vector<unsigned> moduleHasAdd;
    // iterator to vector element:
	map<Module*,Module*> moduleMap;
	map<Module*,unsigned> moduleIDMap;	
	std::map<unsigned,unsigned>::iterator it;
	for(unsigned i=0; i<pLayer.size();i++){
		pLayer[i].setBoundary(
	   _placement.boundryLeft()/_layer.getLayerCount(),
	   _placement.boundryBottom()/_layer.getLayerCount(),
	   _placement.boundryRight()/_layer.getLayerCount(),
	   _placement.boundryTop()/_layer.getLayerCount()); 

	}
	for(unsigned i=0; i<_placement.numNets();i++){
		Net& n0= _placement.net(i);
		LayerMgr tempLayer(_layer.getLayerCount());
//		cout<<"layer_count:"<<_layer.getLayerCount();
		int min_z=_layer.getLayerCount();
		int max_z=0;
		for(unsigned j=0; j<n0.numPins(); j++){
			Pin & p=n0.pin(j);
			Module & m0 = _placement.module(p.moduleId());
			tempLayer.addModule(_layer.getModuleLayer(&m0),&m0);
			cout<<"module_layer:"<<_layer.getModuleLayer(&m0)<<endl;
			max_z=max(max_z,_layer.getModuleLayer(&m0));
			min_z=min(min_z,_layer.getModuleLayer(&m0));
		}
		for(unsigned j=0;j<tempLayer.getLayerCount();j++ ){
			Net* n = new Net();
			pLayer[j].addNet(*n);
			double sum_x=0;
			double sum_y=0;
			Module* m=0;				
			for(unsigned k=0;k<tempLayer.getLayerSize(j);k++){
				Module* m0=tempLayer.getModule(j,k);
				map<Module*,Module*>::iterator it=moduleMap.find(m0);
				sum_x+=m0->centerX();
				sum_y+=m0->centerY();
				m=0;
				if(it==moduleMap.end()){
					m=new Module(m0->name(),m0->width(),m0->height(),m0->isFixed());
					pLayer[j].addModule(*m);
					m->setCenterPosition(m0->centerX(),m0->centerY());
					m->setOrient(m0->orient());
					moduleMap[m0]=m;
					moduleIDMap[m]=pLayer[j].numModules();
				}
				else{
					m=it->second;	
				}
				assert(m);
				Pin * p=new Pin();
				pLayer[j].addPin(*p);
				Pin( moduleIDMap[m]  ,pLayer[j].numNets(),0 , 0);
				n->addPin(p);
				m->addPin(p);
			}
			cout<<"max_z:"<<max_z<<" min_z:"<<min_z<<endl;
			if(max_z-min_z>0){
				m= new Module("tsv_"+int2str((int) i )+"_"+int2str((int) j ),TSV_SIZE,TSV_SIZE,false);
				pLayer[j].addModule(*m);
				m->setCenterPosition(sum_x/tempLayer.getLayerSize(j),sum_y/tempLayer.getLayerSize(j));
				Pin * p=new Pin();
				pLayer[j].addPin(*p);
				Pin( pLayer[j].numModules() ,pLayer[j].numNets(),0 , 0);
				n->addPin(p);
				m->addPin(p);
			}
		}
	}

/*
	for(unsigned i=0; i<_layer.getLayerCount();i++){
		for(unsigned j=0; j<_layer.getLayerSize(i);j++){
			Module * m= _layer.getModule(i,j);
			pLayer[i].addModule(*m);
		}
	}
*/
	/*
	for(unsigned i=0; i<_layer.getLayerCount();i++){
		for(unsigned j=0;j<pLayer[i].numModules();j++){
			cout<<pLayer[i].module(j).name()<<endl;
		}
	}*/
}
