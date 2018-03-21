#include "Vector.h"

Vector3 Vector(Vector3 vPoint1, Vector3 vPoint2)
{
	//printf("Vector from points.\n");
	Vector3 vVector;							
	vVector.x = vPoint1.x - vPoint2.x;			
	vVector.y = vPoint1.y - vPoint2.y;			
	vVector.z = vPoint1.z - vPoint2.z;			
	return vVector;								
}

Vector3 AddVector(Vector3 vVector1, Vector3 vVector2)
{
	//printf("Add vector.\n");
	Vector3 vResult;							
	vResult.x = vVector2.x + vVector1.x;		
	vResult.y = vVector2.y + vVector1.y;		
	vResult.z = vVector2.z + vVector1.z;		
	return vResult;								
}

Vector3 DivideVectorByScaler(Vector3 vVector1, float Scaler)
{
	//printf("Divide vector by scaler.\n");
	Vector3 vResult;							
	vResult.x = vVector1.x / Scaler;			
	vResult.y = vVector1.y / Scaler;			
	vResult.z = vVector1.z / Scaler;			
	return vResult;								
}

Vector3 Cross(Vector3 vVector1, Vector3 vVector2)
{
	//printf("Cross product.\n");
	Vector3 vNormal;									
	vNormal.x = ((vVector1.y * vVector2.z) - (vVector1.z * vVector2.y));
	vNormal.y = ((vVector1.z * vVector2.x) - (vVector1.x * vVector2.z));
	vNormal.z = ((vVector1.x * vVector2.y) - (vVector1.y * vVector2.x));
	return vNormal;										
}

float Dot(Vector3 vVector1, Vector3 vVector2) 
{
 //printf("Dot product.\n");
 return ( (vVector1.x * vVector2.x) + (vVector1.y * vVector2.y) + (vVector1.z * vVector2.z) );
}

Vector3 Normalize(Vector3 vNormal)
{
	//printf("Normalize.\n");
	float magnitude = Magnitude(vNormal);				
	vNormal.x /= magnitude;								
	vNormal.y /= magnitude;								
	vNormal.z /= magnitude;								
	return vNormal;							
}

Vector3 Normal(Vector3 vTriangle[])					
{	
	//printf("Normal.\n");
	Vector3 vVector1 = vTriangle[2] - vTriangle[0];
	Vector3 vVector2 = vTriangle[1] - vTriangle[0];
	Vector3 vNormal = Cross(vVector1, vVector2);		
	vNormal = Normalize(vNormal);						
	return vNormal;										
}

float Magnitude(Vector3 vNormal)
{
 //printf("Magnitude.\n");
 return (float)sqrt( (vNormal.x * vNormal.x) + (vNormal.y * vNormal.y) + (vNormal.z * vNormal.z) );
}

float Absolute(float num)
{
 printf("Absolute.\n");
	if(num < 0) return (0 - num);
	return num;	
}

