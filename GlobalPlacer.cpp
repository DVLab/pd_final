#include "GlobalPlacer.h"
#include "ExampleFunction.h"
#include "NumericalOptimizer.h"
#include <stdlib.h>
#include <iostream>
#include <string>
#include <cassert>
#include <math.h>

#include "Util.h"
#define TSV_SIZE 50
#include <queue>
#include <fstream>

//extern Placement pLayer;
extern vector<Placement*> pLayer;
extern string int2str(int i);


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
		_placement.setBoundary(
	   _placement.boundryTop(),///sqrt((double)_layer.getLayerCount()),
	   _placement.boundryLeft(),//sqrt((double)_layer.getLayerCount()),
	   _placement.boundryBottom(),//sqrt((double)_layer.getLayerCount()),
	   _placement.boundryRight());//sqrt((double)_layer.getLayerCount()) ); 
	
    double wirelength,wirelength_2=0.0,wirelength_3=0.0;
    priority_queue<double> Q_x;
    priority_queue<double> Q_y;
    double c_X = (_placement.boundryRight() + _placement.boundryLeft())/2;
    double c_Y = (_placement.boundryTop() + _placement.boundryBottom())/2;
    double d_X = _placement.boundryRight() - _placement.boundryLeft();
    double d_Y = _placement.boundryTop() - _placement.boundryBottom();
    double num = _placement.numModules();
    //double max_x = -10000000,max_y = -10000000,max_z = -10000000;
    //double min_x = 10000000,min_y = 10000000,min_z = 10000000;

		

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


	unsigned limit=_placement.numNets();
	for(unsigned i=0;i<_placement.numRows();i++){
		

	}
/*	for(unsigned j=0;j<_placement.numModules();j++){
		Module& m=_placement.module(j);
		cout<<m.name()<<" "<<m.width()<<" "<<m.height()<<endl;
	}*/
/*
cout<<"********0_start**********"<<endl;
	for(unsigned j=0;j<_placement.numModules();j++){
		Module& m=_placement.module(j);
		cout<<m.name()<<endl;
		cout<<m.name()<<"_";
		for(unsigned k=0;k<m.numPins();k++){
			Net & n0=_placement.net(m.pin(k).netId());
			for(unsigned z2=0; z2<n0.numPins(); z2++){
				Pin & p=n0.pin(z2);
				cout<<_placement.module(p.moduleId()).name()<<"_";
			}
		}
		cout<<endl;
	}
	for(unsigned z1=0; z1<limit;z1++){
		cout<<"net_";
		Net & n0=_placement.net(z1);
		for(unsigned z2=0; z2<n0.numPins(); z2++){
			Pin & p=n0.pin(z2);
			cout<<_placement.module(p.moduleId()).name()<<"_";
		}
		cout<<endl;
	}
cout<<"********0_end**********"<<endl;*/

	///////////////////////////////////////////////////////////////
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
    no.setNumIteration(3000/*1.5*ite/(num/st)*/); // user-specified parameter
    no.setStepSizeBound(10/*step*(num/st)*/); // user-specified parameter
    no.solve(); // Conjugate Gradient solver

    for(unsigned i = 0;i<num;++i){
        _placement.module(i).setCenterPosition(no.x(i),no.x(i+num));
        x[i] = no.x(i);
        x[i+num] = no.x(i+num);
    }

    ef.reset();
    no.setX(x); 
    no.setNumIteration(30000/*ite/(num/st)*/); 
    no.setStepSizeBound(100/*step*(num/st)*/); 
    no.solve();

    for(unsigned i = 0;i<num;++i){
        _placement.module(i).setCenterPosition(no.x(i),no.x(i+num));
        outfile_x<<no.x(i)<<endl;
        outfile_y<<no.x(i+num)<<endl;
    }

    outfile_x.close();
    outfile_y.close();

    vector< vector<unsigned> > recordModule(4);
    vector< vector<unsigned> > recordNumTSV(4);
    for(unsigned i = 0 ; i < _placement.numNets() ; i++){
        double numModuleI = 0;
	double numModuleII = 0;
	double numModuleIII = 0;
	double numModuleIV = 0;
	bool hasTSV12 = false;  // 12 -> has TSV between I & II
	bool hasTSV23 = false;
	bool hasTSV34 = false;
	bool hasTSV41 = false;
	for(unsigned j = 0 ; j < _placement.net(i).numPins() ; j++){
	    if(_placement.module(_placement.net(i).pin(j).moduleId()).centerX() > c_X){
		if(_placement.module(_placement.net(i).pin(j).moduleId()).centerY() > c_Y)
		    numModuleI += 1;
		else
		    numModuleIV += 1;
	    }
	    else{
		if(_placement.module(_placement.net(i).pin(j).moduleId()).centerY() > c_Y)
		    numModuleII += 1;
		else
		    numModuleIII += 1;
	    }

	    if(numModuleI != 0 && numModuleII != 0){ // need TSV between I & II	
		hasTSV12 = true;
	    }
	    if(numModuleII != 0 && numModuleIII != 0){ // need TSV between II & III	
		hasTSV23 = true;
	    }
	    if(numModuleIII != 0 && numModuleIV != 0){ // need TSV between III & IV	
		hasTSV34 = true;
	    }
	    if(numModuleIV != 0 && numModuleI != 0){ // need TSV between IV & I		
		hasTSV41 = true;
	    }
	    if(hasTSV12 && hasTSV23 && hasTSV34 && hasTSV41)
		break;
	    }
	    if(hasTSV12)
		recordNumTSV[0].push_back(i);
    	    if(hasTSV23)
		recordNumTSV[1].push_back(i);
            if(hasTSV34)
		recordNumTSV[2].push_back(i);
	    if(hasTSV41)
		recordNumTSV[3].push_back(i);
        }

        for(unsigned i = 0 ; i < num ; i++){
	    if(_placement.module(i).centerX() > c_X){
		if(_placement.module(i).centerY() > c_Y)
		    recordModule[0].push_back(i);
		else
	   	    recordModule[3].push_back(i);
	    }
	    else{
		if(_placement.module(i).centerY() > c_Y)
	   	    recordModule[1].push_back(i);
		else
		    recordModule[2].push_back(i);
	    } 	
	}

	for(size_t i = 0 ; i < recordNumTSV.size() ; i++){
	    cout<<"recordNumTSV["<<i<<"] "<<recordNumTSV[i].size()<<endl;
	}
	for(size_t i = 0 ; i < recordModule.size() ; i++){
	    cout<<"recordModule["<<i<<"] "<<recordModule[i].size()<<endl;
	}
	cout<<"aaa "<<endl;
	for(size_t i = 0 ; i < recordModule.size() ; i++){
	    for(size_t j = 0 ; j < recordModule[i].size() ; j++){
                    double _xx = _placement.module(recordModule[i].at(j)).x();
                    double _yy = _placement.module(recordModule[i].at(j)).y();
	

		    if(i == 0){
		        _layer.addModule(3, &_placement.module(recordModule[i].at(j)));
                        _placement.module(recordModule[i].at(j)).setCenterPosition(_xx,_yy);
                    }
		    else if(i == 1){
			_layer.addModule(1, &_placement.module(recordModule[i].at(j)));
                        _placement.module(recordModule[i].at(j)).setCenterPosition(c_X-_xx,_yy);
                    }
		    else if(i == 2){
			_layer.addModule(2, &_placement.module(recordModule[i].at(j)));
                        _placement.module(recordModule[i].at(j)).setCenterPosition(c_X-_xx,c_Y-_yy);
                    }
		    else if(i == 3){
			_layer.addModule(0, &_placement.module(recordModule[i].at(j)));
                        _placement.module(recordModule[i].at(j)).setCenterPosition(_xx,c_Y-_yy);
                    }
	        
	    }
	}
	vector<int> numTSV(3,0);
	for(unsigned i = 0;i<_placement.numNets();++i){
		int z_top = -10;
		int z_bottom = 10;
		for(unsigned j=0;j<_placement.net(i).numPins();++j){
			unsigned id = _placement.net(i).pin(j).moduleId();
			Module& mm = _placement.module(id);
			if(_layer.getModuleLayer(&mm)>z_top)
				z_top = _layer.getModuleLayer(&mm);
			if(_layer.getModuleLayer(&mm)<z_bottom)
				z_bottom = _layer.getModuleLayer(&mm);
		}
		for(unsigned j=0;j<numTSV.size();++j){
			if(z_top>=j+1 && z_bottom<=j)
				numTSV[j] += 1;
		}cout<<i<<"  "<<z_top<<" "<<z_bottom<<" "<<numTSV[0]<<" "<<numTSV[1]<<" "<<numTSV[2]<<endl;
	}
	cout<<"1  "<<numTSV[0]<<endl;
	cout<<"2  "<<numTSV[1]<<endl;
	cout<<"3  "<<numTSV[2]<<endl;
	//cout<<"bbb"<<endl;
	//for(size_t i = 0 ; i < num ; i++)
	//	cout<<_layer.getModuleLayer(&_placement.module(i))<<" ";
	//int count = 0 ;
 	//for(size_t i = 0 ; i < recordModule.size() ; i++)
	//	count += recordModule[i].size();
	//cout<<"count "<<count<<endl;

    cout << "Objective: " << no.objective() << endl;
    printf( "\nLast HPWL: %.0f\n",wirelength);
    printf( "\nLast HPWL_2: %.0f\n",wirelength_2);
    printf( "\nNew HPWL: %.0f\n",_placement.computeHpwl());

    for(unsigned i = 0;i<_placement.numNets();++i){
        double up=-1000000,down=1000000,right=-1000000,left=1000000;
        for(unsigned j=	0;j<_placement.net(i).numPins();++j){
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
/*    cout<<"after:   "<<minBX<<" "<<maxBX<<" "<<minBY<<" "<<maxBY<<endl;
	cout<<CenterX<<"  "<<CenterY<<"  "<<CenterZ<<endl;
    cout << "Objective: " << no.objective() << endlu;*/
////////////////////////////////////////////////////////////////
	pLayer.resize(_layer.getLayerCount());
//	vector<unsigned> moduleHasAdd;
    // iterator to vector element:
	map<Module*,Module*> moduleMap;
	map<Module*,unsigned> moduleIDMap;	
	std::map<unsigned,unsigned>::iterator it;
	for(unsigned i=0; i<pLayer.size();i++){
		pLayer[i]=new Placement();
		pLayer[i]->setBoundary(
	   _placement.boundryTop(),
	   _placement.boundryLeft(),//_layer.getLayerCount(),
	   _placement.boundryBottom(),//_layer.getLayerCount(),
	   _placement.boundryRight());//_layer.getLayerCount()); 
	   pLayer[i]->setRowHeight(_placement.getRowHeight());
	   pLayer[i]->setRectangle(_placement.rectangleChip());
	   pLayer[i]->setName(_placement.name());
	   pLayer[i]->setPlName(_placement.plname());
	   for(unsigned j=0;j<_placement.numRows();j++){//TODO
			pLayer[i]->addRow(_placement.row(j));
	   }
	}

	for(unsigned j=0;j<_placement.numModules();j++){
		_layer.moveModule(&_placement.module(j),0,j%_layer.getLayerCount());
	}



//	int limit=_placement.numNets();
	for(unsigned i=0; i<limit;i++){
		Net& n0= _placement.net(i);
		LayerMgr tempLayer(_layer.getLayerCount());
//		cout<<"layer_count:"<<_layer.getLayerCount();
		int min_z=_layer.getLayerCount();
		int max_z=0;
		double sum_x=0;
		double sum_y=0;
		for(unsigned j=0; j<n0.numPins(); j++){
			Pin & p=n0.pin(j);
			Module & m0 = _placement.module(p.moduleId());
			tempLayer.addModule(_layer.getModuleLayer(&m0),&m0);
			sum_x+=m0.centerX();
			sum_y+=m0.centerY();

		//	cout<<"x:"<<m0.x()<<" y:"<<m0.y()<<" z:"<<_layer.getModuleLayer(&m0)<<endl;
			max_z=max(max_z,_layer.getModuleLayer(&m0));
			min_z=min(min_z,_layer.getModuleLayer(&m0));
		}
		for(unsigned j=0;j<tempLayer.getLayerCount();j++ ){
			if(tempLayer.getLayerSize(j)!=0){
				Net* np = new Net();
				pLayer[j]->addNet(*np);
				Net& n=pLayer[j]->net(pLayer[j]->numNets()-1);
		//		double sum_x=0;
		//		double sum_y=0;
				Module* mp=0;//new Module();
			//	mp=0;
			//	delete mp;
			//	cout<<"n NetNum 0:"<<n.numPins()<<endl;
			//	cout<<"pLayer[i] NetNum 0:"<<pLayer[i]->net(pLayer[i]->numNets()-1).numPins()<<endl;
				for(unsigned k=0;k<tempLayer.getLayerSize(j);k++){
					Module* m0=tempLayer.getModule(j,k);
					map<Module*,Module*>::iterator it=moduleMap.find(m0);
					//	cout<<m0->name()<<endl;	
					if(it==moduleMap.end()){
				//		cout<<"type1"<<endl;
						mp=new Module(m0->name(),m0->width(),m0->height(),m0->isFixed());
						pLayer[j]->addModule(*mp);

						Module& m=pLayer[j]->module(pLayer[j]->numModules()-1);
						m.setCenterPosition(m0->centerX(),m0->centerY());
						m.setOrient(m0->orient());
				//		cout<<"m0:"<<m0->name()<<" m:"<<m.name()<<" id:"<<pLayer[j]->numModules()-1<<endl;
						moduleMap[m0]=&m;
						moduleIDMap[&m]=pLayer[j]->numModules()-1;

						Pin * p=new Pin( pLayer[j]->numModules()-1  ,pLayer[j]->numNets()-1,0 , 0);
						pLayer[j]->addPin(*p);
						n.addPin(p);
						m.addPin(p);
					}
					else{
				//		cout<<"type2"<<endl;
						Module& m=pLayer[j]->module(moduleIDMap[it->second]);
				//		cout<<"m0:"<<it->first->name()<<" m:"<<m.name()<<" id:"<< moduleIDMap[it->second] <<endl;
						Pin * p=new Pin( moduleIDMap[it->second]  ,pLayer[j]->numNets()-1,0 , 0);
						pLayer[j]->addPin(*p);
						n.addPin(p);
						m.addPin(p);
					}
				}
				//cout<<"max_z:"<<max_z<<" min_z:"<<min_z<<endl;
				if(max_z-min_z>0){
					//string tsv_name="tsv_"+int2str((int) i )+"_"+int2str((int) j );
					string tsv_name="tsv";//_"+int2str((int) i )+"_"+int2str((int) j );
					mp=new Module(tsv_name ,TSV_SIZE,TSV_SIZE,true);
					pLayer[j]->addModule(*mp);
					Module& m=pLayer[j]->module(pLayer[j]->numModules()-1);
					m.setCenterPosition(sum_x/tempLayer.getTotalModuleNum(),sum_y/tempLayer.getTotalModuleNum());
					Pin * p=new Pin( pLayer[j]->numModules()-1 ,pLayer[j]->numNets()-1,0 , 0);
					pLayer[j]->addPin(*p);				
					n.addPin(p);
					m.addPin(p);
				}
			//	cout<<"n NetNum 1:"<<n.numPins()<<endl;
			//	cout<<"pLayer[j] NetNum 1:"<<pLayer[j]->net(pLayer[j]->numNets()-1).numPins()<<endl;
			}
		}
	}

	for(unsigned i=0; i<pLayer.size();i++){
		pLayer[i]->connectPinsWithModulesAndNets();
		pLayer[i]->updateDesignStatistics_public();
	}
	/*
	cout<<"********1_start**********"<<endl;
	for(unsigned i=0; i<_layer.getLayerCount();i++){
		for(unsigned j=0;j<pLayer[i]->numModules();j++){
			Module& m=pLayer[i]->module(j);
			cout<<m.name()<<endl;
			cout<<m.name()<<"_";
			for(unsigned k=0;k<m.numPins();k++){
				Net & n0=pLayer[i]->net(m.pin(k).netId());
				for(unsigned z2=0; z2<n0.numPins(); z2++){
					Pin & p=n0.pin(z2);
					cout<<pLayer[i]->module(p.moduleId()).name()<<"_";
				}
			}
			cout<<endl;
		}
		for(unsigned z1=0; z1<pLayer[i]->numNets();z1++){
			cout<<"net_";
			Net & n0=pLayer[i]->net(z1);
			for(unsigned z2=0; z2<n0.numPins(); z2++){
				Pin & p=n0.pin(z2);
				cout<<pLayer[i]->module(p.moduleId()).name()<<"_";
			}
			cout<<endl;
		}
	}
	cout<<"********1_end**********"<<endl;*/

}
