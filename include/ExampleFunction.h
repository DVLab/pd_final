#ifndef EXAMPLEFUNCTION_H
#define EXAMPLEFUNCTION_H

#include "NumericalOptimizerInterface.h"
#include "Placement.h"
#include "GlobalPlacer.h"
class ExampleFunction : public NumericalOptimizerInterface
{
public:
    ExampleFunction(Placement &placement,LayerMgr &layer);

    void evaluateFG(const vector<double> &x, double &f, vector<double> &g);
    void evaluateF(const vector<double> &x, double &f);
    unsigned dimension();

private:
    Placement& _placement;
    LayerMgr& _layer;
};
#endif // EXAMPLEFUNCTION_H
