#ifndef GLOBALPLACER_H
#define GLOBALPLACER_H

#include "Placement.h"
#include <map>

//int Layer::_totalId=0;
typedef vector<Module*> Layer; 

class LayerMgr
{
public:
	//LayerMgr(){}
	
	LayerMgr( unsigned size)
	//	:_placement(placement)
	{
		resizeLayer(size);
	}

	//LayerMgr(unsigned size){
	//	resizeLayer(size); 
	//}
	~LayerMgr(){}

	void clearAll(){
		_layerVec.clear();
		_moduleMap.clear();
	}

	void resizeLayer(unsigned size){
//		clearAll();
		_layerVec.resize(size);
	}
	void addLayer(unsigned count=1){
		_layerVec.resize(_layerVec.size()+count);
	}
	unsigned getLayerSize(unsigned layerId){
		return _layerVec[layerId].size();
	}
	unsigned getTotalModuleNum(){
		unsigned s=0;
		for(unsigned i=0;i<_layerVec.size();i++){
			s+=getLayerSize(i);
		}
		return s;
	}
	unsigned getLayerCount(){
		return _layerVec.size();
	}
	Module* removeModule(int layerId,unsigned i){
		if((layerId<0)||(layerId>(int)_layerVec.size())){return 0;}
		return removeModuleFromLayer( _layerVec[(unsigned)layerId],i);
	}
	bool removeModule(int layerId,Module* m){
		if((layerId<0)||(layerId>(int)_layerVec.size())){return 0;}
		return removeModuleFromLayer( _layerVec[(unsigned)layerId],m);
	}
	bool addModule(int layerId,Module* m){
		if((layerId<0)||(layerId>(int)_layerVec.size())){return 0;}
		addModuleToLayer( _layerVec[(unsigned)layerId],m,layerId);
		return 1;
	}
	Module* getModule(int layerId,unsigned i){
		if((layerId<0)||(layerId>(int)_layerVec.size())){return 0;}
		return getModuleFromLayer( _layerVec[(unsigned)layerId],i);
	}
	bool moveModule(Module *m,int layer_1,int layer_2){
		if(m==0 || layer_1<0 || layer_2>(int)_layerVec.size()){
			cout<<"Wrong!!!"<<endl;
			cout<<"m: "<<m<<"  layer_1: "<<layer_1<<"  layer_2: "<<layer_2<<endl;
			return false;
		}
		removeModule(layer_1,m);
		addModule(layer_2,m);
		return true;
   }
	const  int getModuleLayer( Module* m){
		map<Module*,int>::const_iterator it=_moduleMap.find(m);
		if(it==_moduleMap.end()){
			return -1;
		}
		else{
			return it->second;
		}
	}
/*
	unsigned getLayerTSVNum(unsigned layer){
		unsigned num = 0;
		for(unsigned i = 0 ; i < _placement.numNets() ; i++){
			bool hasModuleAbove = false;
			bool hasModuleBelow = false;
			for(unsigned j = 0 ; j < _placement.net(i).numPins() ; j++){
				unsigned tmp = getModuleLayer(&_placement.module(_placement.net(i).pin(j).moduleId()));
				if(tmp > layer)
					hasModuleAbove = true;
				if(tmp < layer)
					hasModuleBelow = true;
				if(hasModuleBelow && hasModuleAbove){
					num += 1;
					break;
				}
			}
		}
		return num;
	}
*/


private:
	Module* getModuleFromLayer(Layer& layer,unsigned i){
		if(i>=layer.size()){
			return 0;
		}
		else{
			return layer[i];
		}
	}
	Module* removeModuleFromLayer(Layer& layer,unsigned i){
		if(layer.size()==0){
			return 0;
		}
		i = (i <= layer.size()) ? i : layer.size();
		vector<Module* >::iterator it =layer.begin();
		for(unsigned j =0 ;j<i;j++){
			it++;
		}
		Module* m2=(*it);
		layer.erase(it);
		_moduleMap[m2]=-1;
		return m2;
	}
	bool removeModuleFromLayer(Layer& layer,Module* m){
		if(m==0){
			return false;
		}
		vector<Module* >::iterator it =layer.begin();
		while(it!=layer.end()){
			if((*it)==m){
				layer.erase(it);
				_moduleMap[m]=-1;
				break;
			}
			else{
				it++;
			}
		}
		return true;
	}
	void addModuleToLayer(Layer& layer,Module *m,int layerId){
		if(m!=0){
			layer.push_back(m);
			_moduleMap[m]=layerId;
		}
	}
	vector<Layer> _layerVec;
	map<Module*,int> _moduleMap;
//	Placement& _placement;
};


class GlobalPlacer 
{
public:
    GlobalPlacer(Placement &placement,LayerMgr &layer);
	void place();

private:
    Placement& _placement;
    LayerMgr&  _layer;
	
};

#endif // GLOBALPLACER_H
