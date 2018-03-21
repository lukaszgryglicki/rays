#include "qtgui.hxx"
#include "../rayslib.h"

extern int numargs;
extern int current_arg;
extern char** cmdline;

QtUI::QtUI( QWidget* parent, const char* name )
    : QWidget( parent, name )
{

    QPopupMenu *file = new QPopupMenu( this );
    file->insertItem( "&Go!",  this, SLOT(go()), CTRL+Key_G);
    file->insertItem( "&Help",  this, SLOT(help()), CTRL+Key_H);
    file->insertItem( "&Quit",  qApp, SLOT(quit()), CTRL+Key_Q);

    // Create a menu bar
    QMenuBar *m = new QMenuBar( this );
    m->setSeparator( QMenuBar::InWindowsStyle );
    m->insertItem("&File", file );
}

void QtUI::help()
{
 QMessageBox::about(this, 
	QString("Help"), 
	QString(
		"This is Help"
		));

}

void QtUI::go()
{
 run_RT(numargs, cmdline);	
}

#include "moc_qtgui.cxx"

