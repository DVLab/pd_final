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

    double CenterX = (_placement.boundryRight() + _placement.boundryLeft())/2;
    double CenterY = (_placement.boundryTop() + _placement.boundryBottom())/2;
    double CenterZ = 1 + (_layer.getLayerCount())/2;
    //double placement_width = _placement.boundryRight() - _placement.boundryLeft();
    //double placement_height = _placement.boundryTop() - _placement.boundryBottom();
    //double placement_z = _layer.getLayerCount();
    double num = _placement.numModules();
    //double max_x = -10000000,max_y = -10000000,max_z = -10000000;
    //double min_x = 10000000,min_y = 10000000,min_z = 10000000;
		_placement.setBoundary(
	   _placement.boundryLeft()/_layer.getLayerCount(),
	   _placement.boundryBottom()/_layer.getLayerCount(),
	   _placement.boundryRight()/_layer.getLayerCount(),
	   _placement.boundryTop()/_layer.getLayerCount()); 


    ExampleFunction ef(_placement,_layer); // require to define the object function and gradient function

    vector<double> x(2*num,0);

    for(size_t i = 0 ; i < num ; i++){
        int xx = rand()%200;
        if(xx%2 == 0){
            x[i] = CenterX + (double)xx;
            x[i+num] = CenterY + (double)xx;
        }
        else{
            x[i] = CenterX - (double)xx;
            x[i+num] = CenterY - (double)xx;
        }
    }

    for(size_t i = 0 ; i < num ; i++){
        _placement.module(i).setCenterPosition(x[i], x[i+num]);
    }

/*
	cout<<"********0_start**********"<<endl;
	for(size_t i = 0 ; i < num ; i++){
        cout<<_placement.module(i).name()<<endl;
    }
	cout<<"********0_end**********"<<endl;*/
    double minBX = 10000.0;
    double maxBX = -10000.0;
    double minBY = 10000.0;
    double maxBY = -10000.0;

    for(size_t i = 0 ; i < num ; i++){
        if(x[i] > maxBX)
            maxBX = x[i];
        if(x[i] < minBX)
            minBX = x[i];
        if(x[i+num] > maxBY)
            maxBY = x[i+num];
        if(x[i+num] < minBY)
            minBY = x[i+num];
    }
    //cout<<minBX<<" "<<maxBX<<" "<<minBY<<" "<<maxBY<<endl;
		
    srand(time(NULL));
    //double r = (double)rand()/RAND_MAX;

    NumericalOptimizer no(ef);
    no.setX(x); // set initial solution
    no.setNumIteration(5000); // user-specified parameter
    no.setStepSizeBound(50); // user-specified parameter
    no.solve(); // Conjugate Gradient solver
	
    for(size_t i = 0 ; i < num ; i++){
        _placement.module(i).setCenterPosition(no.x(i), no.x(i+num));
    }

    cout<<"before:  "<<minBX<<" "<<maxBX<<" "<<minBY<<" "<<maxBY<<endl;
    for(size_t i = 0 ; i < num ; i++){
        if(no.x(i) > maxBX)
            maxBX = no.x(i);
        if(no.x(i) < minBX)
            minBX = no.x(i);
        if(no.x(i+num) > maxBY)
            maxBY = no.x(i+num);
        if(no.x(i+num) < minBY)
            minBY = no.x(i+num);
    }
    cout<<"after:   "<<minBX<<" "<<maxBX<<" "<<minBY<<" "<<maxBY<<endl;
	cout<<CenterX<<"  "<<CenterY<<"  "<<CenterZ<<endl;
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

	for(unsigned j=0;j<_placement.numModules();j++){
		_layer.moveModule(&_placement.module(j),0,j%5);
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
			cout<<"x:"<<m0.x()<<" y:"<<m0.y()<<" z:"<<_layer.getModuleLayer(&m0)<<endl;
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
			//cout<<"max_z:"<<max_z<<" min_z:"<<min_z<<endl;
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
	cout<<"********1_start**********"<<endl;
	for(unsigned i=0; i<_layer.getLayerCount();i++){
		for(unsigned j=0;j<pLayer[i].numModules();j++){
			cout<<pLayer[i].module(j).name()<<endl;
		}
	}
	cout<<"********1_end**********"<<endl;
*/
}
