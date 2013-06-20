#ifndef GLOBALPLACER_H
#define GLOBALPLACER_H

#include "Placement.h"
#include <map>

//int Layer::_totalId=0;
typedef vector<Module*> Layer; 

class LayerMgr
{
public:
	LayerMgr(){}

	LayerMgr(unsigned size){
		resizeLayer(size); 
	}
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
	const  int getModuleLayer( Module* m){
		map<Module*,int>::const_iterator it=_moduleMap.find(m);
		if(it==_moduleMap.end()){
			return -1;
		}
		else{
			return it->second;
		}
	}
private:
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
};


class GlobalPlacer 
{
public:
    GlobalPlacer(Placement &placement);
	void place();

private:
    Placement& _placement;
	
};

#endif // GLOBALPLACER_H
