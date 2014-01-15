#ifndef DELE_GLOBAL_H
#define DELE_GLOBAL_H

#include <QtCore/qglobal.h>

#if defined(DELE_LIBRARY)
#  define DELESHARED_EXPORT Q_DECL_EXPORT
#else
#  define DELESHARED_EXPORT Q_DECL_IMPORT
#endif

#endif // DELE_GLOBAL_H
