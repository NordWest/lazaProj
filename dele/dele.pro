#-------------------------------------------------
#
# Project created by QtCreator 2014-01-14T12:56:07
#
#-------------------------------------------------

QT       -= gui

TARGET = dele
TEMPLATE = lib

DEFINES += DELE_LIBRARY

SOURCES += dele.cpp \
        ephem_util.cpp

HEADERS += dele.h\
        dele_global.h \
        ephem_types.h \
        ephem_util.h

unix:!symbian {
    maemo5 {
        target.path = /opt/usr/lib
    } else {
        target.path = /usr/lib
        header.path=/usr/local/include
        header.files= dele.h ephem_types.h ephem_util.h
    }

    INSTALLS+=header
    INSTALLS += target
}
