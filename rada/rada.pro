#-------------------------------------------------
#
# Project created by QtCreator 2014-01-14T12:45:48
#
#-------------------------------------------------

QT       -= gui

TARGET = rada
TEMPLATE = lib

LIBS += -lgomp
QMAKE_CXXFLAGS+=-fopenmp

DEFINES += RADA_LIBRARY

SOURCES += rada.cpp \
    force_ev.cpp

HEADERS += rada.h\
        rada_global.h

unix:!symbian {
    maemo5 {
        target.path = /opt/usr/lib
    } else {
        target.path = /usr/local/lib
        header.path = /usr/local/include

    }
    header.files=rada.h
    INSTALLS += target header
}

LIBS += -lgomp -ldele -lephem_read
QMAKE_CXXFLAGS+=-fopenmp
