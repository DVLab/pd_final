#include "ExampleFunction.h"
#include <math.h>
#include <iostream>

ExampleFunction::ExampleFunction(Placement &placement,LayerMgr &layer):_placement(placement),_layer(layer)
{
}

void ExampleFunction::evaluateFG(const vector<double> &x, double &f, vector<double> &g)
{
   f = 0;
   unsigned num = _placement.numModules();
   double r = 0.01*(_placement.boundryRight() - _placement.boundryLeft());

   for(unsigned i = 0;i<3*num;++i)
      g[i] = 0;


   for(unsigned i = 0 ; i<_placement.numNets();++i){
      double f_x_1=0.0,f_x_2=0.0,f_y_1=0.0,f_y_2=0.0,f_z_1=0.0,f_z_2=0.0;
      vector<unsigned> moduleID;
      Net& net = _placement.net(i);
      for(unsigned j = 0; j<net.numPins();++j){
         moduleID.push_back(net.pin(j).moduleId());
      }

      for(unsigned k = 0; k<moduleID.size();++k){
         double xi,yi,zi;
         if(_placement.module(moduleID[k]).isFixed()){
            Module& module = _placement.module(moduleID[k]);
            xi = (module.centerX()+net.pin(k).xOffset());
            yi = (module.centerY()+net.pin(k).yOffset());
            zi = _layer.getModuleLayer(&module);
         }    
         else{
            xi = (x[moduleID[k]]+net.pin(k).xOffset());
            yi = (x[moduleID[k]+num]+net.pin(k).yOffset());
            zi = (x[moduleID[k]+2*num]);
         }

         f_x_1 = f_x_1 + exp( xi /r);     //exp( (module.centerX()+net.pin(k).xOffset()) /r)
         f_x_2 = f_x_2 + exp((-1*xi)/r);  //exp((-1*(module.centerX()+net.pin(k).xOffset()))/r);
         f_y_1 = f_y_1 + exp( yi /r);     //exp((module.centerY()+net.pin(k).yOffset())/r);
         f_y_2 = f_y_2 + exp((-1*yi)/r);  //exp((-1*(module.centerY()+net.pin(k).yOffset()))/r);
         f_z_1 = f_z_1 + exp( zi /r);
         f_z_2 = f_z_2 + exp((-1*zi)/r);

      }
      f = f + log(f_x_1) + log(f_x_2) + log(f_y_1) + log(f_y_2) + log(f_z_1) + log(f_z_2);

      for(unsigned k=0; k<moduleID.size();++k){
         if(_placement.module(moduleID[k]).isFixed()==false){
            g[moduleID[k]]     += ( (exp((x[moduleID[k]]+net.pin(k).xOffset())/r))*(1/f_x_1) + (exp((-1*(x[moduleID[k]]+net.pin(k).xOffset()))/r)*(-1))*(1/f_x_2) ); 
            g[moduleID[k]+num] += ( (exp((x[moduleID[k]+num]+net.pin(k).yOffset())/r))*(1/f_y_1) + (exp((-1*(x[moduleID[k]+num])+net.pin(k).xOffset())/r)*(-1))*(1/f_y_2) );
            g[moduleID[k]+2*num]+=( (exp((x[moduleID[k]+2*num])/r))*(1/f_z_1) + (exp((-1*(x[moduleID[k]+2*num]))/r)*(-1))*(1/f_z_2) );
         }
      }

   }
   f = f * r;
}

void ExampleFunction::evaluateF(const vector<double> &x, double &f)
{
   f = 0;
   unsigned num = _placement.numModules();
   double r = 0.01*(_placement.boundryRight() - _placement.boundryLeft());

   for(unsigned i = 0 ; i<_placement.numNets();++i){
      double f_x_1=0.0,f_x_2=0.0,f_y_1=0.0,f_y_2=0.0,f_z_1=0.0,f_z_2=0.0;
      vector<unsigned> moduleID;
      Net& net = _placement.net(i);
      for(unsigned j = 0; j<net.numPins();++j){
         moduleID.push_back(net.pin(j).moduleId());
      }

      for(unsigned k = 0; k<moduleID.size();++k){
         double xi,yi,zi;
         if(_placement.module(moduleID[k]).isFixed()){
            Module& module = _placement.module(moduleID[k]);
            xi = (module.centerX()+net.pin(k).xOffset());
            yi = (module.centerY()+net.pin(k).yOffset());
            zi = _layer.getModuleLayer(&module);
         }    
         else{
            xi = (x[moduleID[k]]+net.pin(k).xOffset());
            yi = (x[moduleID[k]+num]+net.pin(k).yOffset());
            zi = (x[moduleID[k]+2*num]);
         }

         f_x_1 = f_x_1 + exp( xi /r);     //exp( (module.centerX()+net.pin(k).xOffset()) /r)
         f_x_2 = f_x_2 + exp((-1*xi)/r);  //exp((-1*(module.centerX()+net.pin(k).xOffset()))/r);
         f_y_1 = f_y_1 + exp( yi /r);     //exp((module.centerY()+net.pin(k).yOffset())/r);
         f_y_2 = f_y_2 + exp((-1*yi)/r);  //exp((-1*(module.centerY()+net.pin(k).yOffset()))/r);
         f_z_1 = f_z_1 + exp( zi /r);
         f_z_2 = f_z_2 + exp((-1*zi)/r);

      }
      f = f + log(f_x_1) + log(f_x_2) + log(f_y_1) + log(f_y_2) + log(f_z_1) + log(f_z_2);
   }
   f = f * r;
}

unsigned ExampleFunction::dimension()
{
    return  _placement.numModules()*3;
    // each two dimension represent the X and Y dimensions of each block
}
