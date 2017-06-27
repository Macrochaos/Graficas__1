#ifndef UAD_PRIMITIVEMANAGER_H
#define UAD_PRIMITIVEMANAGER_H

#include <vector>
#include <d3dx9math.h>
#include "PrimitiveBase.h"
#include "Cmatrix4D.h"

class PrimitiveManager {
public:
	void SetVP(Cmatrix4D *vp) {
		pVP = vp;
	}
	int  CreateTriangle();
	int	 CreateCube();

	void DrawPrimitives();
	void DestroyPrimitives();
	PrimitiveBase*	GetPrimitive(unsigned int);

	std::vector<PrimitiveBase*> primitives;

	Cmatrix4D *pVP;
};

#endif