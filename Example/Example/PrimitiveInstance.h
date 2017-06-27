#ifndef UAD_PRIMITIVE_INSTANCE_H
#define UAD_PRIMITIVE_INSTANCE_H

#include "PrimitiveBase.h"
#include <d3dx9math.h>
#include "Cmatrix4D.h"



class PrimitiveInst {
	public:
		void	CreateInstance(PrimitiveBase *pPrim, Cmatrix4D *pVP) {
			pBase = pPrim;
			pViewProj = pVP;

			Position = _Identity();
			Scale = _Identity();
			RotationX = _Identity();
			RotationY = _Identity();
			RotationZ = _Identity();
			//Cmatrix4D RotationZ = _Identity();
			Final = _Identity();
			/*
			D3DXMatrixIdentity(&Position);
			D3DXMatrixIdentity(&Scale);
			D3DXMatrixIdentity(&RotationX);
			D3DXMatrixIdentity(&RotationY);
			D3DXMatrixIdentity(&RotationZ);
			D3DXMatrixIdentity(&RotationZ);
			D3DXMatrixIdentity(&Final);	
			*/
		}

		void	TranslateAbsolute(float x, float y, float z);
		void	RotateXAbsolute(float ang);
		void	RotateYAbsolute(float ang);
		void	RotateZAbsolute(float ang);
		void	ScaleAbsolute(float sc);

		void	TranslateRelative(float x, float y, float z);
		void	RotateXRelative(float ang);
		void	RotateYRelative(float ang);
		void	RotateZRelative(float ang);
		void	ScaleRelative(float sc);

		void	Update();
		void	Draw();

		Cmatrix4D		Position;
		Cmatrix4D		Scale;
		Cmatrix4D		RotationX;
		Cmatrix4D		RotationY;
		Cmatrix4D		RotationZ;
		Cmatrix4D		Final;

		Cmatrix4D		*pViewProj;

		PrimitiveBase	*pBase;
};


#endif