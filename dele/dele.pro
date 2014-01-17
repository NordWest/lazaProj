#-------------------------------------------------
#
# Project created by QtCreator 2014-01-15T11:01:22
#
#-------------------------------------------------

QT       -= gui

TARGET = dele
TEMPLATE = lib

DEFINES += DELE_LIBRARY

VERSION = 1.0.0

SOURCES += dele.cpp
HEADERS += dele.h \
        dele_global.h

LIBS += -lephem_read

unix:!symbian {
    maemo5 {
        target.path = /opt/usr/lib
    } else {
        target.path = /usr/local/lib
        header.path=/usr/local/include
    }
    header.files=dele.h
    INSTALLS += target
    INSTALLS+=header
}
