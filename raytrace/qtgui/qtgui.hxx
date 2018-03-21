#ifndef __QT_UI_MDBMA__
#define __QT_UI_MDBMA__

#include <qwidget.h>
#include <qfiledialog.h>
#include <qpushbutton.h>
#include <qslider.h>
#include <qlayout.h>
#include <qframe.h>
#include <qmenubar.h>
#include <qpopupmenu.h>
#include <qapplication.h>
#include <qkeycode.h>
#include <qmessagebox.h>

class QtUI : public QWidget
{
    Q_OBJECT
public slots:
    void help();
    void go();
public:
    QtUI(QWidget* parent = 0, const char* name = 0);
private:
};


#endif

