win32 {
 TARGET = ./../../libs/win32/mb
}
unix {
TARGET = mb
}
TEMPLATE = lib
SOURCES += mb.cpp
HEADERS += mb.h
CONFIG += staticlib warn_off
header.path=/usr/local/include
header.files=mb.h
target.path +=/usr/local/lib
INSTALLS+=header
INSTALLS+=target
