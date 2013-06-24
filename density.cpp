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
   for(unsigned i = 0;i<2*_placement.numModules();++i)
      g[i] = 0;
   unsigned num = _placement.numModules();
   double dx,dy,absdx,absdy;
   double Px,Py;
   double dPx,dPy;
   double D = 0;

   double a_x, b_x;
   double a_y, b_y;

   double binNum = 1.0;
   double binWidth = (_placement.boundryRight() - _placement.boundryLeft()) / sqrt(binNum);
   double binHeight = (_placement.boundryTop() - _placement.boundryBottom()) / sqrt(binNum);
   double M = 0*(binWidth*binHeight - 0);

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
         }
      }
   }
}
void ExampleFunction::evaluateF(const vector<double> &x, double &f)
{
   f = 0;
   unsigned num = _placement.numModules();
   double dx,dy,absdx,absdy;
   double Px,Py;
   double D = 0;

   double a_x, b_x;
   double a_y, b_y;

   double binNum = 1.0;
   double binWidth = (_placement.boundryRight() - _placement.boundryLeft()) / sqrt(binNum);
   double binHeight = (_placement.boundryTop() - _placement.boundryBottom()) / sqrt(binNum);
   double M = 0*(binWidth*binHeight - 0); 

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
   }
}

unsigned ExampleFunction::dimension()
{
    return  _placement.numModules()*2;
    // each two dimension represent the X and Y dimensions of each block
}
