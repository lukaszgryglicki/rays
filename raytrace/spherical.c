/* Written by MorgothDBMA, morgothdbma@o2.pl, tel: +48693582014 */
/* Lukasz Gryglicki MiNI M1 CC */
/* License BSD */
#include <stdlib.h>
#include <stdio.h>
#include <GL/glut.h>
#include <jpeglib.h>
#include <setjmp.h>
#include <math.h>

#define PI 3.1415926F
#define ERR_CANNOTREAD 1
#define ERR_BADJPEG    2
#define ERR_GRAYJPEG   3
#define ERR_256CJPEG   4
#define  ERR_NOMEMORY  5

double angX, angY, angZ, obsD, traX, traY, traZ, tStep;	
int wx,wy;		
int cx,cy;
unsigned int tid;
int tri_mode;
int ndivs1, ndivs2;
double scale;
double angle_cut;

struct my_error_mgr
{
 struct jpeg_error_mgr pub;
 jmp_buf setjmp_buffer;
};

typedef struct my_error_mgr * my_error_ptr;
static my_error_ptr errptr = NULL;

static void my_error_exit(j_common_ptr cinfo)
{
 my_error_ptr myerr = (my_error_ptr) cinfo->err;
 (*cinfo->err->output_message) (cinfo);
 longjmp(myerr->setjmp_buffer, 1);
}


void help()
{
	printf("argument: file_in.jpeg --> spherical map\n");
	printf("keys\n\twsadez:\trotations\n\tk:\ttriangles/lines toggle\n");
	printf("\t1234:\tmanipulate ndivs1,ndivs2\n\t56:\tscaling world\n");
	printf("\t78:\tmanipulate observer distance\n");
	printf("\tikjlom:\ttranslations\n\t90:\tangle strip cut\n\th:\thelp\n");
	printf("\tn:\tinfo curr params\n");
}

void normals()
{
 angX = angZ = 0.;
 angY = 180.;
 tri_mode = 1;
 ndivs1 = 180;
 ndivs2 = 90;
 scale = 10.;
 obsD = 2.;
 traX = traY = traZ = 0.;
 angle_cut = 25.;
}

void info()
{
	printf("angle: (%f, %f, %f)\n", angX, angY, angZ);
	printf("translate: (%f, %f, %f), obsD=%f\n", traX, traY, traZ, obsD);
	printf("triangle_net: %d(%dx%d)\n", tri_mode, ndivs1, ndivs2);
	printf("scale = %f, angle_cut = %f\n", scale, angle_cut);

}

void check_scene()
{
	if (angX < 0.) angX += 360.;
	if (angY < 0.) angY += 360.;
	if (angZ < 0.) angZ += 360.;

	if (angX >= 360.) angX -= 360.;
	if (angY >= 360.) angY -= 360.;
	if (angZ >= 360.) angZ -= 360.;

	if (ndivs1 < 8) ndivs1 = 8;
	if (ndivs2 < 6) ndivs2 = 6;
	if (ndivs1 > 720) ndivs1 = 720;
	if (ndivs2 > 720) ndivs2 = 720;

	if (scale < 2.) scale = 2.;
	if (scale > 15000.) scale = 15000.;
	if (obsD < 0.5) obsD = 0.5;
	if (obsD > 20.) obsD = 20.;

	if (angle_cut < 0.) angle_cut = 0.;
	if (angle_cut > 50.) angle_cut = 50.;

	//printf("scale = %f\n", scale);
	//printf("angle_cut = %f\n", angle_cut);

}

void keyboard(unsigned char key, int x, int y)
{
 switch (key)
   {
        case 27: case 'q':  exit(0); break;
		case 'a': angY += 1.; break;
		case 'd': angY -= 1.; break;
		case 'w': angX += 1.; break;
		case 's': angX -= 1.; break;
		case 'e': angZ += 1.; break;
		case 'z': angZ -= 1.; break;
		case 't': tri_mode = !tri_mode; break;
		case '1': ndivs1 = (int)((double)ndivs1 / 1.2); break;
		case '2': ndivs1 = (int)((double)ndivs1 * 1.2); break;
		case '3': ndivs2 = (int)((double)ndivs2 / 1.2); break;
		case '4': ndivs2 = (int)((double)ndivs2 * 1.2); break;
		case '5': scale  /= 1.1; break;
		case '6': scale  *= 1.1; break;
		case '7': obsD  -= .03; break;
		case '8': obsD  += .03; break;
		case 'l': traX -= .1; break;
		case 'j': traX += .1; break;
		case 'i': traZ += .1; break;
		case 'k': traZ -= .1; break;
		case 'm': traY += .1; break;
		case 'o': traY -= .1; break;
		case '9': angle_cut -= .5; break;
		case '0': angle_cut += .5; break;
		case ' ': normals();  break;
		case 'h': help(); break;
		case 'n': info(); break;
   }
   //printf("Net: (%dx%d)\n", ndivs1, ndivs2);
   check_scene();
}

void resize_scene(int w, int h)
{
 wx = w;
 wy = h;
 glViewport(0, 0, (GLsizei)wx, (GLsizei)wy);
 glMatrixMode(GL_PROJECTION);
 glLoadIdentity();
 glFrustum(-1., 1., -1., 1., 1.1, 30000.);
 glMatrixMode(GL_MODELVIEW);
}

void renderSpherical()
{
	int i, j;
	double x[4],y[4],z[4], s[4], t[4], s1, s2, t1, t2;
	double angles1, angles2, anglet1, anglet2;
	double acut;

	if (tri_mode) glBegin(GL_TRIANGLES);
	else glBegin(GL_LINES);
	glColor3f(1., 1., 1.);

	acut = PI - ((angle_cut * PI) / 90.);

	for (i=0;i<ndivs1;i++)
	{
		s1 = (double)i / (double)ndivs1;
		s2 = (double)(i+1) / (double)ndivs1;
		angles1 = s1 * 2. * PI;
		angles2 = s2 * 2. * PI;

		for (j=0;j<ndivs2;j++)
		{
			t1 = (double)j / (double)ndivs2;
			t2 = (double)(j+1) / (double)ndivs2;
			anglet1 = (t1 * acut) - acut/2.;
			anglet2 = (t2 * acut) - acut/2.;

			s[0] = s1;
			t[0] = t1;
			s[1] = s1;
			t[1] = t2;
			s[2] = s2;
			t[2] = t2;
			s[3] = s2;
			t[3] = t1;

			x[0] = sin(angles1) * cos(anglet1);
			y[0] = sin(anglet1);
			z[0] = -cos(angles1) * cos(anglet1);
			x[1] = sin(angles1) * cos(anglet2);
			y[1] = sin(anglet2);
			z[1] = -cos(angles1) * cos(anglet2);
			x[2] = sin(angles2) * cos(anglet2);
			y[2] = sin(anglet2);
			z[2] = -cos(angles2) * cos(anglet2);
			x[3] = sin(angles2) * cos(anglet1);
			y[3] = sin(anglet1);
			z[3] = -cos(angles2) * cos(anglet1);

			glTexCoord2d(s[0], t[0]);
			glVertex3d(x[0], y[0], z[0]);
			glTexCoord2d(s[1], t[1]);
			glVertex3d(x[1], y[1], z[1]);
			glTexCoord2d(s[2], t[2]);
			glVertex3d(x[2], y[2], z[2]);

			glTexCoord2d(s[0], t[0]);
			glVertex3d(x[0], y[0], z[0]);
			glTexCoord2d(s[2], t[2]);
			glVertex3d(x[2], y[2], z[2]);
			glTexCoord2d(s[3], t[3]);
			glVertex3d(x[3], y[3], z[3]);
		}
	}
	glEnd();
}

void render_scene(void)
{
 glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
 glLoadIdentity();
 gluLookAt(0., 0., obsD, 0. , 0., 1., 0., 1., 0.);
 glTranslated(traX, traY, traZ);
 glScaled(scale, scale, scale);
 glRotated(angX, 1., 0., 0.);
 glRotated(angY, 0., 1., 0.);
 glRotated(angZ, 0., 0., 1.);
/* glTranslated(-(GLdouble)(wx/2), - (GLdouble)(wy/2), 0.);*/
/* glColor3d(0.,1.,0.);*/
 renderSpherical();
 glFlush();
 glutSwapBuffers();
}

void anim(void)
{
 glutPostRedisplay();
}

void visible(int vis)
{
 if (vis == GLUT_VISIBLE) glutIdleFunc(anim);
 else                     glutIdleFunc(NULL);
}

void set_rows(int cur_row, JSAMPARRAY pixel_data, int width, unsigned char** bits)
{
	JSAMPROW ptr;
	register int j;
	ptr = pixel_data[0];
	for (j=0;j<width;j++)
	{
		(*bits)[(cur_row*width+j)*3] 	= GETJSAMPLE(*ptr++);
		(*bits)[((cur_row*width+j)*3)+1] 	= GETJSAMPLE(*ptr++);
		(*bits)[((cur_row*width+j)*3)+2] 	= GETJSAMPLE(*ptr++);
	}
}


int load_jpeg_file(unsigned char** bits, int* x, int* y, char* filename)
{
    struct jpeg_decompress_struct cinfo;
    struct my_error_mgr jerr;
    FILE* infile;
    JSAMPARRAY buffer;
    int row_stride;
    int i;
    *x = *y = 0;
    *bits = NULL;
    if ((infile = fopen(filename, "rb")) == NULL) return ERR_CANNOTREAD;
    cinfo.err = jpeg_std_error(&jerr.pub);
    jerr.pub.error_exit = my_error_exit;
    errptr = &jerr;
    if (setjmp(jerr.setjmp_buffer))
    {
		jpeg_destroy_decompress(&cinfo);
		fclose(infile);
		return ERR_BADJPEG;
    }
    jpeg_create_decompress(&cinfo);
    jpeg_stdio_src(&cinfo, infile);
    jpeg_read_header(&cinfo, TRUE);
    if (cinfo.jpeg_color_space == JCS_GRAYSCALE)
    {
		jpeg_destroy_decompress(&cinfo);
		fclose(infile);
		return ERR_GRAYJPEG;
    }
    else cinfo.quantize_colors = FALSE;
    jpeg_start_decompress(&cinfo);
    if (cinfo.output_components == 1)
    {
		jpeg_destroy_decompress(&cinfo);
		fclose(infile);
		return ERR_256CJPEG;
    }
    *x = cinfo.output_width;
    *y = cinfo.output_height;
    *bits = calloc(3*(*y)*(*x), sizeof(unsigned char));
    if (!(*bits)) return ERR_NOMEMORY;
    row_stride = cinfo.output_width * cinfo.output_components;
    buffer = (*cinfo.mem->alloc_sarray)((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);
    i = (*y)-1;
    while (cinfo.output_scanline < cinfo.output_height)
       {
		jpeg_read_scanlines(&cinfo, buffer, 1);		/* FIXME howto read multiple lines at the time?? */
		set_rows(i,buffer,*x, bits);
		i--;
       }
    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    fclose(infile);
    return 0;
}

void TransformTexture(unsigned char** bits, int* cx, int* cy)
{
 int i,j;
 int x,y,xx,yy;
 unsigned char* nbits;
 i = j = 1;
 while (i<(*cx)) i*=2;
 while (j<(*cy)) j*=2;
 nbits = malloc(3*i*j);
 if (!nbits) return;
 for (y=0;y<j;y++)
 {
	yy = (y * (*cy)) / j;
	for (x=0;x<i;x++)
	{
		xx = (x * (*cx)) / i;
		nbits[3*(i*y+x)]     = (*bits)[3*(yy*(*cx)+xx)];
		nbits[(3*(i*y+x))+1] = (*bits)[(3*(yy*(*cx)+xx))+1];
		nbits[(3*(i*y+x))+2] = (*bits)[(3*(yy*(*cx)+xx))+2];
	}
 }
 free(*bits);
 *bits = nbits;
 *cx = i;
 *cy = j;
}


void LoadTexture(char *filename, GLuint* texture)
{
  int err;
  int cx,cy;
  unsigned char* bits;
  err = load_jpeg_file(&bits, &cx, &cy, filename);
  glGenTextures(1, texture);
  glBindTexture(GL_TEXTURE_2D, *texture);
  //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  if (err) 
    { 
     printf("Cannot load texture: %s... generating dynamic texture...\n",filename); 
     exit(1);
    }
  else TransformTexture(&bits, &cx, &cy);
  glTexImage2D(GL_TEXTURE_2D, 0, 3, cx, cy, 0, GL_RGB, GL_UNSIGNED_BYTE, bits);
  if (bits) free(bits);
}

void Init(char* fn)
{
 unsigned int tex;
 help();
 normals();

 glClearColor(0.,0.,0.,0.);
 glShadeModel(GL_SMOOTH);
 glEnable(GL_DEPTH_TEST);
 glEnable(GL_TEXTURE_2D);
 LoadTexture(fn, &tid);
}


int main(int lb, char** par)
{
 if (lb < 2) 
 {
		printf("Argument required: file_in.jpeg\n");
		return 1;
 }
 wx = 800;
 wy = 800;
 glutInit(&lb, par);
 glutInitDisplayMode(GLUT_DOUBLE);
 glutInitWindowSize(wx, wy);
 glutInitWindowPosition(100, 100);
 glutCreateWindow(par[0]);
 Init(par[1]);
 glutDisplayFunc(render_scene);
 glutReshapeFunc(resize_scene);
 glutKeyboardFunc(keyboard);
 glutVisibilityFunc(visible);
 glutMainLoop();
 return 0;
}

/* Written by MorgothDBMA, morgothdbma@o2.pl, tel: +48693582014 */
/* Lukasz Gryglicki MiNI M1 CC */
/* License BSD */
