#include "qtgui.hxx"
#include <qapplication.h>
#include <string.h>

#define MAXCMDS 512

int numargs;
int current_arg;
char** cmdline;

void setup_defaults()
{
 numargs = 1;
 cmdline = new char*[MAXCMDS];
 for (int i=0;i<MAXCMDS;i++) 
   {
    cmdline[i] = NULL;
   }
 current_arg = 0;
 cmdline[0] = new char[16];
 strcpy(cmdline[0], "QTraytracer");
 current_arg ++;
}

int main( int argc, char **argv)
{
    setup_defaults();	
    QApplication a(argc,argv);
    QtUI w;
    w.resize( 300, 300 );
    a.setMainWidget(&w);
    w.show();
    a.exec();
    return 0;
}
