#ifndef _H_JPEG_LOADER_H__MDBMA_
#define _H_JPEG_LOADER_H__MDBMA_

typedef unsigned char                                   BYTE, uchar;
typedef unsigned short int                              WORD;
typedef unsigned long                                   DWORD;
typedef unsigned int                                    UINT;
//typedef unsigned __int64                                QWORD;

#define RGB16(r, g, b) ((unsigned short)(((r >> 3) << 11)+((g >> 2) << 5)+(b >> 3)))
#define RGB15(r, g, b) ((unsigned short)(((r >> 3) << 10)+((g >> 3) << 5)+(b >> 3)))
#define RGB32(r, g, b) ((unsigned long)((r << 16)+(g << 8)+b))

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define FLOAT_PI 3.141592654f
#define FLOAT_EPSILON 1.0e-6f
 
#define RADTODEG(RAD) (((float)RAD*180.0f)/FLOAT_PI)

#define DEGTORAD(DEG) (((float)DEG*FLOAT_PI)/180.0f)

// for matriser     CMatrix matM; matM(_Y_, _X_) = 1.0f;
#define _X_		0
#define _Y_		1
#define _Z_		2
#define _W_		3

class JpegLoader
{
	public:
		bool Load(char *filename, int bitsperpixel, void *where);
		void GetSize(char *filename, unsigned long &width, unsigned long &height);
};


#endif
