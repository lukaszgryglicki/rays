#ifndef _VECTOR_H__MDBMA_
#define _VECTOR_H__MDBMA_
#include "Main.h"

typedef struct Vector3
{			
	Vector3() {}
	Vector3 (float new_x, float new_y, float new_z) {x = new_x; y = new_y; z = new_z;}
	Vector3 operator+(Vector3 vVector) {return Vector3(vVector.x+x, vVector.y+y, vVector.z+z);}
	Vector3 operator-(Vector3 vVector) {return Vector3(x-vVector.x, y-vVector.y, z-vVector.z);}
	Vector3 operator*(float number)	 {return Vector3(x*number, y*number, z*number);}
	Vector3 operator/(float number)	 {return Vector3(x/number, y/number, z/number);}
	float x, y, z;
} Vector3;				
struct Vector2 
{
	Vector2() {}
	Vector2 (float new_x, float new_y) {x = new_x; y = new_y;}
	Vector2 operator+(Vector2 vVector) {return Vector2(vVector.x+x, vVector.y+y);}
	Vector2 operator-(Vector2 vVector) {return Vector2(x-vVector.x, y-vVector.y);}
	Vector2 operator*(float number)	 {return Vector2(x*number, y*number);}
	Vector2 operator/(float number)	 {return Vector2(x/number, y/number);}
	float x, y;
};

Vector3 Vector(Vector3 vPoint1, Vector3 vPoint2);
Vector3 AddVector(Vector3 vVector1, Vector3 vVector2);
Vector3 DivideVectorByScaler(Vector3 vVector1, float Scaler);
Vector3 Cross(Vector3 vVector1, Vector3 vVector2);
float Dot(Vector3 vVector1, Vector3 vVector2);
Vector3 Normalize(Vector3 vNormal);
Vector3 Normal(Vector3 vTriangle[]);
float Magnitude(Vector3 vNormal);
float Absolute(float num);


#endif

