//--- Aurora/QtSupport/QtSupport.hpp -------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_QT_SUPPORT_QT_SUPPORT_HPP
#define AURORA_QT_SUPPORT_QT_SUPPORT_HPP

#include "AuroraLib/Image.hpp"

#include <QImage>

#if (AURORA_PLATFORM == AURORA_PLATFORM_WINDOWS)
# if (AURORA_QT_SUPPORT_BUILD == 1)
#   define AURORA_QT_SUPPORT_API __declspec(dllexport)
# else
#   define AURORA_QT_SUPPORT_API __declspec(dllimport)
# endif
#elif (AURORA_PLATFORM == AURORA_PLATFORM_LINUX)
# define AURORA_QT_SUPPORT_API
#endif

namespace Aurora
{

AURORA_QT_SUPPORT_API std::string saveFileDialog(QWidget* parent,
                                                 const QString& caption,
                                                 const QString& filter,
                                                 const QString& defaultSuffix,
                                                 QString& directory);

AURORA_QT_SUPPORT_API QImage imageToQImage(const Image& image);

inline QSize calculateSize(int cols, Real aspectRatio)
{
  if (aspectRatio >= 1)
    return QSize(cols, static_cast<int>(cols / aspectRatio));
  else
    return QSize(static_cast<int>(cols * aspectRatio), cols);
}

} // namespace Aurora

#endif
