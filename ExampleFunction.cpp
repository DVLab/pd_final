#include "ExampleFunction.h"
#include <math.h>
#include <cmath>
#include <iostream>
// minimize 3*x^2 + 2*x*y + 2*y^2 + 7

ExampleFunction::ExampleFunction(Placement &placement)
	:_placement(placement),count(1)
{
}

void ExampleFunction::evaluateFG(const vector<double> &x, double &f, vector<double> &g)
{

    if(count == -1){
   f = 0;
   for(unsigned i = 0;i<2*_placement.numModules();++i)
      g[i] = 0;
   unsigned num = _placement.numModules();
   double dx,dy,absdx,absdy;
   double Px,Py;
   double dPx,dPy;
   double D = 0;

   double a_x, b_x;
   double a_y, b_y;

   double binNum = 25.0;
   double binWidth = (_placement.boundryRight() - _placement.boundryLeft()) / sqrt(binNum);
   double binHeight = (_placement.boundryTop() - _placement.boundryBottom()) / sqrt(binNum);
   double M = 800;//0*(binWidth*binHeight - 0);

   vector<double> binCenterX,binCenterY,binCenterZ;
   vector<double> dP(2*num,0);
   vector<double> P(2*num,0);
   for(int i = 0 ; i < sqrt(binNum) ; i++){
      binCenterX.push_back(_placement.boundryLeft() + (2*i+1)*((_placement.boundryRight() - _placement.boundryLeft())/ (sqrt(binNum)*2)));
      binCenterY.push_back(_placement.boundryTop() - (2*i+1)*((_placement.boundryTop() - _placement.boundryBottom())/ (sqrt(binNum)*2)));
   }


   for(unsigned y = 0 ; y < sqrt(binNum) ; y++){
      for(unsigned i = 0 ; i < sqrt(binNum) ; i++){
         D = 0;
         for(unsigned j = 0 ; j < num ; j++){
            // x coordinate
            a_x = 4/((_placement.module(j).width() + 2*binWidth)*(_placement.module(j).width() + 4*binWidth));
            b_x = 2/(binWidth*(_placement.module(j).width() + 4*binWidth));
            dx = x[j] - binCenterX[i];
            absdx = abs(dx);
            if(absdx >= 0 && absdx <= (_placement.module(j).width()/2 + binWidth) ){
               Px = 1 - a_x*dx*dx;
               dPx = -2*a_x*dx;
            }
            else if (absdx >= (_placement.module(j).width()/2 + 2*binWidth)){
               Px = 0;
               dPx = 0;
            }
            else{
               Px = b_x*(absdx - _placement.module(j).width()/2 - 2*binWidth)*(absdx - _placement.module(j).width()/2 - 2*binWidth);
               if(dx>0)
                  dPx = 2*b_x*(dx - _placement.module(j).width()/2 - 2*binWidth);
               else
                  dPx = 2*b_x*(dx + _placement.module(j).width()/2 + 2*binWidth);
            }
            // y coordinate
            a_y = 4/((_placement.module(j).height() + 2*binHeight)*(_placement.module(j).height() + 4*binHeight));
            b_y = 2/(binHeight*(_placement.module(j).height() + 4*binHeight));
            dy = x[j+num] - binCenterY[y];
            absdy = abs(dy);
            if(absdy >= 0 && absdy <= (_placement.module(j).height()/2 + binHeight) ){
               Py = 1 - a_y*dy*dy;
               dPy = -2*a_y*dy;
            }
            else if (absdy >= (_placement.module(j).height()/2 + 2*binHeight)){
               Py = 0;
               dPy = 0;
            }
            else{
               Py = b_y*(absdy - _placement.module(j).height()/2 - 2*binHeight)*(absdy - _placement.module(j).height()/2 - 2*binHeight);
               if(dy>0)
                  dPy = 2*b_y*(dy - _placement.module(j).height()/2 - 2*binHeight);
               else
                  dPy = 2*b_y*(dy + _placement.module(j).height()/2 + 2*binHeight);
            }
            dP[j] = dPx;
            dP[j+num] = dPy;
            P[j] = Px;
            P[j+num] = Py;
            D += Px*Py;
         }
         f += (D - M)*(D - M);
         for(unsigned j = 0 ; j < num ; j++){
            g[j] += 2*(D-M)*dP[j]*P[j+num];
            g[j+num] += 2*(D-M)*dP[j+num]*P[j];
            if( x[j]+2*_placement.module(j).width() < _placement.boundryLeft() || x[j]-2*_placement.module(j).width() > _placement.boundryRight() )
               g[j] = 0;
            if( x[j+num]+2*_placement.module(j).height() < _placement.boundryBottom() || x[j+num]-2*_placement.module(j).height() > _placement.boundryTop() )
               g[j+num] = 0;
         }
      }
   }
    }
///////////////////////////////////////////////////////////////////////
    else{
        f = 0;
        unsigned num = _placement.numModules();
        double r = 0.01*(_placement.boundryRight() - _placement.boundryLeft());

        for(unsigned i = 0;i<2*num;++i)
            g[i] = 0;


        for(unsigned i = 0 ; i<_placement.numNets();++i){
            double f_1=0.0,f_2=0.0,f_3=0.0,f_4=0.0;
            vector<unsigned> moduleID;
            Net& net = _placement.net(i);
            for(unsigned j = 0; j<net.numPins();++j){
                moduleID.push_back(net.pin(j).moduleId());
            }

            for(unsigned k = 0; k<moduleID.size();++k){
                if(_placement.module(moduleID[k]).isFixed()){
                    Module& module = _placement.module(moduleID[k]);
                    f_1 = f_1 + exp( (module.centerX()+net.pin(k).xOffset()) /r);
                    f_2 = f_2 + exp((-1*(module.centerX()+net.pin(k).xOffset()))/r);
                    f_3 = f_3 + exp((module.centerY()+net.pin(k).yOffset())/r);
                    f_4 = f_4 + exp((-1*(module.centerY()+net.pin(k).yOffset()))/r);
                }    
                else{
                    f_1 = f_1 + exp((x[moduleID[k]]+net.pin(k).xOffset())/r);
                    f_2 = f_2 + exp((-1*(x[moduleID[k]]+net.pin(k).xOffset()))/r);
                    f_3 = f_3 + exp((x[moduleID[k]+num]+net.pin(k).yOffset())/r);
                    f_4 = f_4 + exp((-1*(x[moduleID[k]+num])+net.pin(k).xOffset())/r);
                }

            }
            f = f + log(f_1) + log(f_2) + log(f_3) + log(f_4) ;

            for(unsigned k=0; k<moduleID.size();++k){
                if(_placement.module(moduleID[k]).isFixed()==false){
                    g[moduleID[k]]     += ( (exp((x[moduleID[k]]+net.pin(k).xOffset())/r))*(1/f_1) + (exp((-1*(x[moduleID[k]]+net.pin(k).xOffset()))/r)*(-1))*(1/f_2) ); 
                    g[moduleID[k]+num] += ( (exp((x[moduleID[k]+num]+net.pin(k).yOffset())/r))*(1/f_3) + (exp((-1*(x[moduleID[k]+num])+net.pin(k).xOffset())/r)*(-1))*(1/f_4) );
                }
            }

        }
        f = f * r;
    }
}

void ExampleFunction::evaluateF(const vector<double> &x, double &f)
{
    if(count == -1){
   f = 0;
   unsigned num = _placement.numModules();
   double dx,dy,absdx,absdy;
   double Px,Py;
   double D = 0;

   double a_x, b_x;
   double a_y, b_y;

   double binNum = 25.0;
   double binWidth = (_placement.boundryRight() - _placement.boundryLeft()) / sqrt(binNum);
   double binHeight = (_placement.boundryTop() - _placement.boundryBottom()) / sqrt(binNum);
   double M = 800;//0*(binWidth*binHeight - 0);

   vector<double> binCenterX,binCenterY,binCenterZ;
   for(int i = 0 ; i < sqrt(binNum) ; i++){
      binCenterX.push_back(_placement.boundryLeft() + (2*i+1)*((_placement.boundryRight() - _placement.boundryLeft())/ (sqrt(binNum)*2)));
      binCenterY.push_back(_placement.boundryTop() - (2*i+1)*((_placement.boundryTop() - _placement.boundryBottom())/ (sqrt(binNum)*2)));
   }


   for(unsigned y = 0 ; y < sqrt(binNum) ; y++){
      for(unsigned i = 0 ; i < sqrt(binNum) ; i++){
         D = 0;
         for(unsigned j = 0 ; j < num ; j++){
            // x coordinate
            a_x = 4/((_placement.module(j).width() + 2*binWidth)*(_placement.module(j).width() + 4*binWidth));
            b_x = 2/(binWidth*(_placement.module(j).width() + 4*binWidth));
            dx = x[j] - binCenterX[i];
            absdx = abs(dx);
            if(absdx >= 0 && absdx <= (_placement.module(j).width()/2 + binWidth) ){
               Px = 1 - a_x*dx*dx;
            }
            else if (absdx >= (_placement.module(j).width()/2 + 2*binWidth)){
               Px = 0;
            }
            else{
               Px = b_x*(absdx - _placement.module(j).width()/2 - 2*binWidth)*(absdx - _placement.module(j).width()/2 - 2*binWidth);
            }
            // y coordinate
            a_y = 4/((_placement.module(j).height() + 2*binHeight)*(_placement.module(j).height() + 4*binHeight));
            b_y = 2/(binHeight*(_placement.module(j).height() + 4*binHeight));
            dy = x[j+num] - binCenterY[y];
            absdy = abs(dy);
            if(absdy >= 0 && absdy <= (_placement.module(j).height()/2 + binHeight)){
               Py = 1 - a_y*dy*dy;
            }
            else if (absdy >= (_placement.module(j).height()/2 + 2*binHeight)){
               Py = 0;
            }
            else{
               Py = b_y*(absdy - _placement.module(j).height()/2 - 2*binHeight)*(absdy - _placement.module(j).height()/2 - 2*binHeight);
            }
            D += Px*Py;
         }
         f += (D - M)*(D - M);
      }
    }}
/////////////////////////////////////
    else{
        f = 0;
        unsigned num = _placement.numModules();
        double r = 0.01*(_placement.boundryRight() - _placement.boundryLeft());

        for(unsigned i = 0 ; i<_placement.numNets();++i){
            double f_1=0.0,f_2=0.0,f_3=0.0,f_4=0.0;
            vector<unsigned> moduleID;
            Net& net = _placement.net(i);
            for(unsigned j = 0; j<net.numPins();++j){
                moduleID.push_back(net.pin(j).moduleId());
            }

            for(unsigned k = 0; k<moduleID.size();++k){
                if(_placement.module(moduleID[k]).isFixed()){
                    Module& module = _placement.module(moduleID[k]);
                    f_1 = f_1 + exp((module.centerX()+net.pin(k).xOffset()) /r);
                    f_2 = f_2 + exp((-1*(module.centerX()+net.pin(k).xOffset()))/r);
                    f_3 = f_3 + exp((module.centerY()+net.pin(k).yOffset())/r);
                    f_4 = f_4 + exp((-1*(module.centerY()+net.pin(k).yOffset()))/r);
                }    
                else{
                    f_1 = f_1 + exp((x[moduleID[k]]+net.pin(k).xOffset())/r);
                    f_2 = f_2 + exp((-1*(x[moduleID[k]]+net.pin(k).xOffset()))/r);
                    f_3 = f_3 + exp((x[moduleID[k]+num]+net.pin(k).yOffset())/r);
                    f_4 = f_4 + exp((-1*(x[moduleID[k]+num])+net.pin(k).xOffset())/r);
                }
            }
            f = f + log(f_1) + log(f_2) + log(f_3) + log(f_4) ;
        }
        f = f * r;

    }
}

unsigned ExampleFunction::dimension()
{
    return 2*_placement.numModules(); // num_blocks*2 
    // each two dimension represent the X and Y dimensions of each block
}

/*double ExampleFunction::abs(double a,bool &inv){
    if(a>0){inv = false;
        return a;}
    else{inv = true;
        return (-1)*a;}
}*/

void ExampleFunction::reset(){
    count = -count;
}
