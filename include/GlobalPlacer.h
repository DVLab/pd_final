#ifndef GLOBALPLACER_H
#define GLOBALPLACER_H

#include "Placement.h"

class GlobalPlacer 
{
public:
    GlobalPlacer(Placement &placement);
	void place();

private:
    Placement& _placement;
	
};

#endif // GLOBALPLACER_H
