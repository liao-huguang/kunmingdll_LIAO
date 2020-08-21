#-------------------------------------------------
#
# Project created by QtCreator 2020-08-20T10:23:29
#
#-------------------------------------------------

QT       -= gui

TARGET = kunmingdll
TEMPLATE = lib

DEFINES += KUNMINGDLL_LIBRARY

SOURCES += kunmingdll.cpp

HEADERS += kunmingdll.h \
    H/myq_dtranssuub_wz.h \
    H/math/mysum.h \
    H/myinitsub_m_wz.h

unix {
    target.path = /usr/lib
    INSTALLS += target
}
