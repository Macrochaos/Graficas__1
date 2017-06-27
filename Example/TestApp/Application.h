#include <Core.h>
#include <PrimitiveManager.h>
#include <PrimitiveInstance.h>

#include <d3dx9math.h>
#include <Timer.h>

class TestApp : public AppBase {
public:
	TestApp() : AppBase() {}
	void InitVars();
	void CreateAssets();
	void DestroyAssets();

	void OnUpdate();
	void OnDraw();
	void OnInput();

	void OnPause();
	void OnResume();

	void OnReset();

	PrimitiveManager PrimitiveMgr;
	PrimitiveInst	Cubes[10];
	PrimitiveInst	Triangle[10];

	CVector4D		Position;
	CVector4D		Orientation;
	CVector4D		Scaling;
	Cmatrix4D		WorldTransform;

	Cmatrix4D		View;
	Cmatrix4D		Projection;
	Cmatrix4D		VP;

	Timer			DtTimer;
};