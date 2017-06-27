#include "PrimitiveInstance.h"

void PrimitiveInst::TranslateAbsolute(float x, float y, float z){
	//D3DXMatrixTranslation(&Position, x, y, z);
	Position = Translation_RH(x, y, z);
}

void PrimitiveInst::RotateXAbsolute(float ang) {
	//D3DXMatrixRotationX(&RotationX,D3DXToRadian(ang));
	RotationX = RotationX_LH(ang);
}

void PrimitiveInst::RotateYAbsolute(float ang) {
	//D3DXMatrixRotationY(&RotationY, D3DXToRadian(ang));
	RotationY = RotationY_LH(ang);
}

void PrimitiveInst::RotateZAbsolute(float ang) {
	//D3DXMatrixRotationZ(&RotationZ, D3DXToRadian(ang));
	RotationZ = RotationZ_LH(ang);
}

void PrimitiveInst::ScaleAbsolute(float sc) {
	//D3DXMatrixScaling(&Scale,sc,sc,sc);
	Scale = Scaling(sc, sc, sc);
}

void PrimitiveInst::TranslateRelative(float x, float y, float z) {
	//D3DXMATRIX tmp;
	//D3DXMatrixTranslation(&tmp, x, y, z);
	Cmatrix4D tmp;
	tmp= Translation_RH(x, y, z);
	Position = Position * tmp;
}

void PrimitiveInst::RotateXRelative(float ang) {
	//D3DXMATRIX tmp;
	//D3DXMatrixRotationX(&tmp, D3DXToRadian(ang));
	Cmatrix4D tmp;
	tmp = RotationX_LH(ang);
	RotationX = RotationX * tmp;
}

void PrimitiveInst::RotateYRelative(float ang) {
	//D3DXMATRIX tmp;
	//D3DXMatrixRotationY(&tmp, D3DXToRadian(ang));
	Cmatrix4D tmp;
	tmp = RotationY_LH(ang);
	RotationY = RotationY * tmp;
}

void PrimitiveInst::RotateZRelative(float ang) {
	//D3DXMATRIX tmp;
	//D3DXMatrixRotationZ(&tmp, D3DXToRadian(ang));
	Cmatrix4D tmp;
	tmp = RotationZ_LH(ang);
	RotationZ = RotationZ * tmp;
}

void PrimitiveInst::ScaleRelative(float sc) {
	//D3DXMATRIX tmp;
	//D3DXMatrixScaling(&tmp, sc, sc, sc);
	Cmatrix4D tmp;
	tmp = Scaling(sc, sc, sc);
	Scale = Scale * tmp;
}

void PrimitiveInst::Update() {
	Final = Scale*RotationX*RotationY*RotationZ*Position;
}

void PrimitiveInst::Draw(){
	pBase->Draw(&Final.m[0][0],&(*pViewProj).m[0][0]);
}
