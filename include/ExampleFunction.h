#ifndef EXAMPLEFUNCTION_H
#define EXAMPLEFUNCTION_H

#include "NumericalOptimizerInterface.h"
#include "Placement.h"

class ExampleFunction : public NumericalOptimizerInterface
{
public:
    ExampleFunction(Placement &placement);

    void evaluateFG(const vector<double> &x, double &f, vector<double> &g);
    void evaluateF(const vector<double> &x, double &f);
    void density(const vector<double> &x, double &f);
    unsigned dimension();
    double abs(double a,bool& inv);
    void reset();

private:
    Placement& _placement;
    double _last;
    int count;
};
#endif // EXAMPLEFUNCTION_H
