#pragma once
class CVector4D
{
public:
	union
	{
		struct { float x, y, z, w; };
		struct { float r, g, b, a; };
		struct { unsigned long ulx, uly, ulz, ulw; };
		float v[4];
		int i[4];

	};
};
class Cmatrix4D
{

public:
	Cmatrix4D operator=(float *t) {
		Cmatrix4D R;
		for (int i = 0; i < 16; i++) {
			R.v[i] = t[i];
		}
		return R;
	}
	union
	{
		struct
		{
			float m00, m01, m02, m03;
			float m10, m11, m12, m13;
			float m20, m21, m22, m23;
			float m30, m31, m32, m33;
		};
		float m[4][4];
		CVector4D vec[4];
		float v[16];
	};
	Cmatrix4D();
	~Cmatrix4D();
	Cmatrix4D(float *vp);
};
Cmatrix4D Zero_init();
Cmatrix4D operator*(Cmatrix4D& A, Cmatrix4D& B);

CVector4D operator*(const CVector4D& A, float f);
CVector4D operator*(CVector4D& V, Cmatrix4D& M);
CVector4D operator*(Cmatrix4D& M, CVector4D& V);
CVector4D operator*(CVector4D& A, CVector4D& B);
CVector4D operator-(const CVector4D& cola, const CVector4D& punta);
CVector4D operator+(const CVector4D& A, const CVector4D& B);
CVector4D operator/(const CVector4D& A, float B);
Cmatrix4D _Identity();
Cmatrix4D RotationX_RH(float theta);
Cmatrix4D RotationY_RH(float theta);
Cmatrix4D RotationZ_RH(float theta);
Cmatrix4D RotationX_LH(float theta);
Cmatrix4D RotationY_LH(float theta);
Cmatrix4D RotationZ_LH(float theta);
Cmatrix4D Transposed(Cmatrix4D& M);
Cmatrix4D Fast_Inverse(Cmatrix4D& M);
float Inverse(Cmatrix4D & M, Cmatrix4D & R);
Cmatrix4D Translation_RH(float dx, float dy, float dz);
Cmatrix4D Scaling(float sx, float sy, float sz);
float Dot(const CVector4D& A, CVector4D& B);
float Mag(const CVector4D& A);
CVector4D BuildPlane(CVector4D vertex0, CVector4D vertex1);
CVector4D Cross3(const CVector4D& A, const CVector4D& B);
CVector4D Normalize(const CVector4D& A);
