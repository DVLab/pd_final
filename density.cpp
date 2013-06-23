#include "ExampleFunction.h"
#include <math.h>
#include <cmath>
#include <iostream>

ExampleFunction::ExampleFunction(Placement &placement,LayerMgr &layer):_placement(placement),_layer(layer)
{
}

void ExampleFunction::evaluateFG(const vector<double> &x, double &f, vector<double> &g)
{
   f = 0;
   unsigned num = _placement.numModules();
   double dx,dy,dz;
   double Px,Py,Pz;
   double dPx,dPy,dPz; 
   double D = 0;

	double a_x, b_x;
	double a_y, b_y;
	double a_z, b_z;

	double binNum = 1.0;
	double binWidth = (_placement.boundryRight() - _placement.boundryLeft()) / sqrt(binNum);
	//cout<<"binWidth "<<binWidth<<endl;
	double binHeight = (_placement.boundryTop() - _placement.boundryBottom()) / sqrt(binNum);
	//cout<<"binHeight "<<binHeight<<endl;
	double M = 0*(binWidth*binHeight - 0); // assume t_density = 1.0 , area of preplaced block = 0
	//cout<<"M "<<M<<endl;

	vector<double> binCenterX,binCenterY,binCenterZ;
	for(int i = 0 ; i < sqrt(binNum) ; i++){
		binCenterX.push_back(_placement.boundryLeft() + (2*i+1)*((_placement.boundryRight() - _placement.boundryLeft())/ (sqrt(binNum)*2)));
		binCenterY.push_back(_placement.boundryTop() - (2*i+1)*((_placement.boundryTop() - _placement.boundryBottom())/ (sqrt(binNum)*2)));
		//binCenterZ.push_back(0 + (2*i+1)*0.5);
	}

	// bin[x][y][z]
//	for(unsigned z = 0 ; z < sqrt(binNum) ; z++){
		for(unsigned y = 0 ; y < sqrt(binNum) ; y++){
			for(unsigned i = 0 ; i < sqrt(binNum) ; i++){
				D = 0;
				for(unsigned j = 0 ; j < num ; j++){
					// x coordinate
					a_x = 4/((_placement.module(j).width() + 2*binWidth)*(_placement.module(j).width() + 4*binWidth));
					b_x = 2/(binWidth*(_placement.module(j).width() + 4*binWidth));
					dx = abs(x[j] - binCenterX[i]);
					if(dx >= 0 && dx <= _placement.module(j).width()/2 + binWidth){
						Px = 1 - a_x*dx*dx;
						dPx = -2*a_x*dx;
					}
					else if (dx >= (_placement.module(j).width()/2 + 2*binWidth)){
						Px = 0;
						dPx = 0;
					}
					else{
						Px = b_x*(dx - _placement.module(j).width()/2 - 2*binWidth)*(dx - _placement.module(j).width()/2 - 2*binWidth);
						dPx = 2*b_x*(dx - _placement.module(j).width()/2 - 2*binWidth);
					}
					// y coordinate
					a_y = 4/((_placement.module(j).height() + 2*binHeight)*(_placement.module(j).height() + 4*binHeight));
					b_y = 2/(binHeight*(_placement.module(j).height() + 4*binHeight));
					dy = abs(x[j+num] - binCenterY[y]);
					if(dy >= 0 && dy <= _placement.module(j).height()/2 + binHeight){
						Py = 1 - a_y*dy*dy;
						dPy = -2*a_y*dy;
					}
					else if (dy >= (_placement.module(j).height()/2 + 2*binHeight)){
						Py = 0;
						dPy = 0;
					}
					else{
						Py = b_y*(dy - _placement.module(j).height()/2 - 2*binHeight)*(dy - _placement.module(j).height()/2 - 2*binHeight);
						dPy = 2*b_y*(dy - _placement.module(j).height()/2 - 2*binHeight);
					}/*
					// z coordinate
					a_z = 4/15;
					b_z = 2/5;
					dz = abs(x[j+2*num] - binCenterZ[z]);
					if(dz >= 0 && dz <= 1.5){
						Pz = 1 - a_z*dz*dz;
						dPz = -2*a_z*dz;
					}
					else if (dz >= 2.5){
						Pz = 0;
						dPz = 1;//0;
					}
					else{
						Pz = b_z*(dz - 2.5)*(dz - 2.5);
						dPz = 2*b_z*(dz - 2.5);
					}*/
					D += Px*Py;//*Pz;
					g[j] = 2*(D-M)*Py*dPx;
					g[j+num] = 2*(D-M)*Px*dPy;
					//cout<<D<<" "<<dx<<"  "<<dy<<"  "<<g[j]<<" "<<g[j+num]<<"  "<<dPx<<"  "<<dPy<<"  "<<endl;
					//g[j] = 2*(D - M)*Py*Pz*dPx;
					//g[j+num] = 2*(D-M)*Px*Pz*dPy;
					//g[j+2*num] = 2*(D-M)*Px*Py*dPz;
				}
				f += (D - M)*(D - M);
			}
		}
//	}
}

void ExampleFunction::evaluateF(const vector<double> &x, double &f)
{
   f = 0;
   unsigned num = _placement.numModules();
   double dx,dy,dz;
   double Px,Py,Pz;
   //double dPx,dPy,dPz; 
   double D = 0;

	double a_x, b_x;
	double a_y, b_y;
	double a_z, b_z;

	double binNum = 1.0;
	double binWidth = (_placement.boundryRight() - _placement.boundryLeft()) / sqrt(binNum);
	//cout<<"binWidth "<<binWidth<<endl;
	double binHeight = (_placement.boundryTop() - _placement.boundryBottom()) / sqrt(binNum);
	//cout<<"binHeight "<<binHeight<<endl;
	double M = 0*(binWidth*binHeight - 0); // assume t_density = 1.0 , area of preplaced block = 0
	//cout<<"M "<<M<<endl;

	vector<double> binCenterX,binCenterY,binCenterZ;
	for(int i = 0 ; i < sqrt(binNum) ; i++){
		binCenterX.push_back(_placement.boundryLeft() + (2*i+1)*((_placement.boundryRight() - _placement.boundryLeft())/ (sqrt(binNum)*2)));
		binCenterY.push_back(_placement.boundryTop() - (2*i+1)*((_placement.boundryTop() - _placement.boundryBottom())/ (sqrt(binNum)*2)));
		//binCenterZ.push_back(0 + (2*i+1)*0.5);
	}

	// bin[x][y][z]
//	for(unsigned z = 0 ; z < sqrt(binNum) ; z++){
		for(unsigned y = 0 ; y < sqrt(binNum) ; y++){
			for(unsigned i = 0 ; i < sqrt(binNum) ; i++){
				D = 0;
				for(unsigned j = 0 ; j < num ; j++){
					// x coordinate
					a_x = 4/((_placement.module(j).width() + 2*binWidth)*(_placement.module(j).width() + 4*binWidth));
					b_x = 2/(binWidth*(_placement.module(j).width() + 4*binWidth));
					dx = abs(x[j] - binCenterX[i]);
					if(dx >= 0 && dx <= _placement.module(j).width()/2 + binWidth){
						Px = 1 - a_x*dx*dx;
						//dPx = -2*a_x*dx;
					}
					else if (dx >= (_placement.module(j).width()/2 + 2*binWidth)){
						Px = 0;
						//dPx = 0;
					}
					else{
						Px = b_x*(dx - _placement.module(j).width()/2 - 2*binWidth)*(dx - _placement.module(j).width()/2 - 2*binWidth);
						//dPx = 2*b_x*(dx - _placement.module(j).width()/2 - 2*binWidth);
					}
					// y coordinate
					a_y = 4/((_placement.module(j).height() + 2*binHeight)*(_placement.module(j).height() + 4*binHeight));
					b_y = 2/(binHeight*(_placement.module(j).height() + 4*binHeight));
					dy = abs(x[j+num] - binCenterY[y]);
					if(dy >= 0 && dy <= _placement.module(j).height()/2 + binHeight){
						Py = 1 - a_y*dy*dy;
						//dPy = -2*a_y*dy;
					}
					else if (dy >= (_placement.module(j).height()/2 + 2*binHeight)){
						Py = 0;
						//dPy = 0;
					}
					else{
						Py = b_y*(dy - _placement.module(j).height()/2 - 2*binHeight)*(dy - _placement.module(j).height()/2 - 2*binHeight);
						//dPy = 2*b_y*(dy - _placement.module(j).height()/2 - 2*binHeight);
					}/*
					// z coordinate
					a_z = 4/15;
					b_z = 2/5;
					dz = abs(x[j+2*num] - binCenterZ[z]);
					if(dz >= 0 && dz <= 1.5){
						Pz = 1 - a_z*dz*dz;
						//dPz = -2*a_z*dz;
					}
					else if (dz >= 2.5){
						Pz = 0;
						//dPz = 0;
					}
					else{
						Pz = b_z*(dz - 2.5)*(dz - 2.5);
						//dPz = 2*b_z*(dz - 2.5);
					}*/
					D += Px*Py;//*Pz;
				}
				f += (D - M)*(D - M);
			}
		}
//	}
	/*
	// Definition of Z(z)
	double numTSV = 0 ;
	for(size_t i = 0 ; i < _placement.numNets() ; i++){
		double tmpMax = -1000.0;
		double tmpMin =  1000.0;
		for(size_t j = 0 ; j < _placement.net(i).numPins() ; j++){
			if(x[_placement.net(i).pin(j).moduleId() + 2*num] > tmpMax)
				tmpMax = x[_placement.net(i).pin(j).moduleId() + 2*num];
			if(x[_placement.net(i).pin(j).moduleId() + 2*num] < tmpMin)
				tmpMin = x[_placement.net(i).pin(j).moduleId() + 2*num];
		}
		numTSV += (tmpMax - tmpMin);
	}
	f += numTSV;
	// End of definition of Z(z)
	*/

}

unsigned ExampleFunction::dimension()
{
    return  _placement.numModules()*2;
    // each two dimension represent the X and Y dimensions of each block
}
