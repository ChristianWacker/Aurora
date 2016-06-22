//--- Aurora/AuroraLib/Exception.cpp -------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/Exception.hpp"

namespace Aurora
{

Exception::Exception(const std::string& message,
                     const std::string& function,
                     const std::string& filename,
                     size_t line,
                     const std::string& type) noexcept :
  mMessage(message), mFunction(function), mFilename(filename), mLine(line),
  mType(type),
  mString(mFilename + "(" + std::to_string(mLine) + "): In function \"" +
          mFunction + "\" an exception of type \"" + mType +
          "\" has been thrown. Message: " + mMessage)
{ }

std::string Exception::message() const
{
  return mMessage;
}

std::string Exception::string() const
{
  return mString;
}

std::string Exception::type() const
{
  return mType;
}

std::string Exception::function() const
{
  return mFunction;
}

std::string Exception::filename() const
{
  return mFilename;
}

size_t Exception::line() const
{
  return mLine;
}

const char* Exception::what() const noexcept
{
  return mString.c_str();
}

EInvalidParameter::EInvalidParameter(const std::string& message,
                                     const std::string& function,
                                     const std::string& filename,
                                     size_t line) :
  Exception(message, function, filename, line, "EInvalidParameter")
{ }


ENotSupported::ENotSupported(const std::string& message,
                             const std::string& function,
                             const std::string& filename,
                             size_t line) :
  Exception(message, function, filename, line, "ENotSupported")
{ }

EUnknown::EUnknown(const std::string& message,
                   const std::string& function,
                   const std::string& filename,
                   size_t line) :
  Exception(message, function, filename, line, "EUnknown")
{ }

EInOut::EInOut(const std::string& message,
               const std::string& function,
               const std::string& filename,
               size_t line,
               const std::string& type) :
  Exception(message, function, filename, line, type)
{ }

EFileNotFound::EFileNotFound(const std::string& path,
                             const std::string& function,
                             const std::string& filename,
                             size_t line) :
  EInOut("The file \"" + path + "\" could not be found.",
         function, filename, line, "EFileNotFound")
{ }

EDirectoryNotFound::EDirectoryNotFound(const std::string& path,
                                       const std::string& function,
                                       const std::string& filename,
                                       size_t line) :
  EInOut("The file \"" + path + "\" could not be found.",
         function, filename, line, "EFileNotFound")
{ }

} // namespace Aurora
