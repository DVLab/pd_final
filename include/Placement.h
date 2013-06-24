#ifndef PLACEMENT_H
#define PLACEMENT_H

#include <vector>
#include <string>
using namespace std;

#include "Module.h"
#include "Net.h"
#include "Pin.h"
#include "Row.h"
#include "Rectangle.h"

class Placement
{
public:
    Placement();

    /////////////////////////////////////////////
    // input/output
    /////////////////////////////////////////////
    void readBookshelfFormat(string filePathName, string plFileName);
    void outputBookshelfFormat(string filePathName); // output file function

    /////////////////////////////////////////////
    // get
    /////////////////////////////////////////////
    string name() {return _name;}
    string plname() {return _loadplname;}
    double computeHpwl();
    double computeTotalNetLength(int cellid);
    Rectangle rectangleChip() {return _rectangleChip;}  //Chip rectangle
    double boundryTop() {return _boundryTop;}
    double boundryLeft() {return _boundryLeft;}
    double boundryBottom() {return _boundryBottom;}
    double boundryRight() {return _boundryRight;}

    /////////////////////////////////////////////
    // operation
    /////////////////////////////////////////////
    void moveDesignCenter(double xOffset, double yOffset);

    /////////////////////////////////////////////
    // get design objects/properties
    /////////////////////////////////////////////
    Module& module(unsigned moduleId) {return _modules[moduleId];}
    Net& net(unsigned netId) {return _nets[netId];}
    Pin& pin(unsigned pinId) {return _pins[pinId];}
    Row& row(unsigned rowId) {return _rows[rowId];}

    double getRowHeight() {return _rowHeight;}

    unsigned numModules() {return _modules.size();}
    unsigned numNets() {return _nets.size();}
    unsigned numPins() {return _pins.size();}
    unsigned numRows() {return _rows.size();}

    /////////////////////////////////////////////
    // methods for design (hypergraph) construction
    /////////////////////////////////////////////
    void addModule(const Module &module) {_modules.push_back(module);}
    void addPin(const Pin &pin) {_pins.push_back(pin),_pins.back().setPinId(_pins.size());}
    void addRow(const Row &row) {_rows.push_back(row);}
    void addNet(const Net &net) {_nets.push_back(net);} //ADD




    void setNumModules(unsigned size) {_modules.resize(size);}
    void setNumNets(unsigned size) {_nets.resize(size);}
    void setNumPins(unsigned size) {_pins.resize(size);}
    void setNumRows(unsigned size) {_rows.resize(size);}

    void clearModules() {_modules.clear();}
    void clearNets() {_nets.clear();}
    void clearPins() {_pins.clear();}
    void clearRows() {_rows.clear();}

    // initialize pins for modules and nets (construct hypergraph)
    void connectPinsWithModulesAndNets();

    ////////////////////
    vector<Row> m_sites; // for Legalization and Detailplace
    vector<Module> modules_bak; //for Detailplace

	
	void setBoundary(double t, double l ,double b ,double r ){//ADD
		 _boundryTop=t;
		 _boundryLeft=l;
		 _boundryBottom=b;
		 _boundryRight=r;
	}
	void setRectangle(Rectangle r){ //ADD
    	_rectangleChip=r;
	}
	void setRowHeight(double h){_rowHeight=h;} //ADD
	void setName(string n){_name=n;}//ADD
	void setPlName(string n){_loadplname=n;}//ADD
	void updateDesignStatistics_public(){
		updateDesignStatistics();
	}
private:
    /////////////////////////////////////////////
    // properties
    /////////////////////////////////////////////
    string _name;
    string _loadplname;

    /////////////////////////////////////////////
    // design data
    /////////////////////////////////////////////
    vector<Module> _modules;
    vector<Net> _nets;
    vector<Pin> _pins;
    vector<Row> _rows;

	
    /////////////////////////////////////////////
    // design statistics
    /////////////////////////////////////////////
    void updateDesignStatistics();
    Rectangle _rectangleChip;
    double _rowHeight;
    double _boundryTop;
    double _boundryLeft;
    double _boundryBottom;
    double _boundryRight;
};

#endif // PLACEMENT_H
