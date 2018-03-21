/****************************************************************************
** QtUI meta object code from reading C++ file 'qtgui.hxx'
**
** Created: Tue Feb 14 16:57:07 2006
**      by: The Qt MOC ($Id: qt/moc_yacc.cpp   3.3.4   edited Jan 21 18:14 $)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#undef QT_NO_COMPAT
#include "qtgui.hxx"
#include <qmetaobject.h>
#include <qapplication.h>

#include <private/qucomextra_p.h>
#if !defined(Q_MOC_OUTPUT_REVISION) || (Q_MOC_OUTPUT_REVISION != 26)
#error "This file was generated using the moc from 3.3.4. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

const char *QtUI::className() const
{
    return "QtUI";
}

QMetaObject *QtUI::metaObj = 0;
static QMetaObjectCleanUp cleanUp_QtUI( "QtUI", &QtUI::staticMetaObject );

#ifndef QT_NO_TRANSLATION
QString QtUI::tr( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "QtUI", s, c, QApplication::DefaultCodec );
    else
	return QString::fromLatin1( s );
}
#ifndef QT_NO_TRANSLATION_UTF8
QString QtUI::trUtf8( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "QtUI", s, c, QApplication::UnicodeUTF8 );
    else
	return QString::fromUtf8( s );
}
#endif // QT_NO_TRANSLATION_UTF8

#endif // QT_NO_TRANSLATION

QMetaObject* QtUI::staticMetaObject()
{
    if ( metaObj )
	return metaObj;
    QMetaObject* parentObject = QWidget::staticMetaObject();
    static const QUMethod slot_0 = {"help", 0, 0 };
    static const QUMethod slot_1 = {"go", 0, 0 };
    static const QMetaData slot_tbl[] = {
	{ "help()", &slot_0, QMetaData::Public },
	{ "go()", &slot_1, QMetaData::Public }
    };
    metaObj = QMetaObject::new_metaobject(
	"QtUI", parentObject,
	slot_tbl, 2,
	0, 0,
#ifndef QT_NO_PROPERTIES
	0, 0,
	0, 0,
#endif // QT_NO_PROPERTIES
	0, 0 );
    cleanUp_QtUI.setMetaObject( metaObj );
    return metaObj;
}

void* QtUI::qt_cast( const char* clname )
{
    if ( !qstrcmp( clname, "QtUI" ) )
	return this;
    return QWidget::qt_cast( clname );
}

bool QtUI::qt_invoke( int _id, QUObject* _o )
{
    switch ( _id - staticMetaObject()->slotOffset() ) {
    case 0: help(); break;
    case 1: go(); break;
    default:
	return QWidget::qt_invoke( _id, _o );
    }
    return TRUE;
}

bool QtUI::qt_emit( int _id, QUObject* _o )
{
    return QWidget::qt_emit(_id,_o);
}
#ifndef QT_NO_PROPERTIES

bool QtUI::qt_property( int id, int f, QVariant* v)
{
    return QWidget::qt_property( id, f, v);
}

bool QtUI::qt_static_property( QObject* , int , int , QVariant* ){ return FALSE; }
#endif // QT_NO_PROPERTIES
