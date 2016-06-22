//--- Aurora/AuroraLib/Exception.hpp -------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_EXCEPTION_HPP
#define AURORA_AURORA_LIB_EXCEPTION_HPP

#include "AuroraLib/Prerequisites.hpp"

#include <exception>
#include <string>

namespace Aurora
{

/// Base class for all exceptions that are thrown from within AuroraLib.
/// @remarks
///  As the constructor of this class is protected one has to use child
///  classes which clarify the type of the exception.
class AURORA_API Exception : public std::exception
{
protected:
  /// Constructs a new exception.
  /// @param message
  ///  A string that contains more information concerning the exception. This
  ///  string is passed through the constructors of the child classes.
  /// @param function
  ///  Name of the function that has thrown this exception
  /// @param filename
  ///  Name of the file that contains the throw-statement.
  /// @param line
  ///  Line number of the throw-statement.
  /// @param type
  ///  The name of the exception class as a string.
  Exception(const std::string& message,
            const std::string& function,
            const std::string& filename,
            size_t line,
            const std::string& type = "Exception") noexcept;

public:
  /// Returns the message passed to the constructor
  std::string message() const;

  /// Returns a string describing the exception. This string includes
  /// the type of the exception, the message, the function name, the filename
  /// and line.
  std::string string() const;

  /// Returns the type of the exception
  std::string type() const;

  /// Returns the name of the function that has thrown this exception
  std::string function() const;

  /// Returns the name of the file that contains the throw-statement.
  std::string filename() const;

  /// Returns line number of the throw-statement for this exception.
  size_t line() const;

  /// Inherited function from std::exception. Return the string from @ref
  /// string() as c-string.
  const char* what() const noexcept override;

private:
  const std::string mMessage;
  const std::string mFunction;
  const std::string mFilename;
  const size_t mLine;
  const std::string mType;
  const std::string mString;
};

/// Exception that is thrown if an invalid parameter is passed to a function.
class AURORA_API EInvalidParameter : public Exception
{
public:
  /// Construct a new exception.
  /// @param message
  ///  A string that describes the cause for this exception being thrown.
  /// @param function
  ///  Name of the function that has thrown this exception
  /// @param filename
  ///  Name of the file that contains the throw-statement.
  /// @param line
  ///  Line number of the throw-statement.
  EInvalidParameter(const std::string& message, const std::string& function,
                    const std::string& filename, size_t line);
};

/// Exception that is thrown if a not supported function is called
class AURORA_API ENotSupported : public Exception
{
public:
  /// @copydoc EInvalidParameter::EInvalidParameter
  ENotSupported(const std::string& message, const std::string& function,
                const std::string& filename, size_t line);
};

/// Exception with an unknown cause.
class AURORA_API EUnknown : public Exception
{
public:
  /// @copydoc EInvalidParameter::EInvalidParameter
  EUnknown(const std::string& message, const std::string& function,
           const std::string& filename, size_t line);
};

/// Base class for all IO-exceptions
class AURORA_API EInOut : public Exception
{
public:
  /// @copydoc Exception::Exception
  EInOut(const std::string& message, const std::string& function,
         const std::string& filename, size_t line,
         const std::string& type = "EInOut");
};

/// This exception is thrown if a file is not found.
class AURORA_API EFileNotFound : public EInOut
{
public:
  /// Constructs a new EFileNotFound exception.
  /// @param path
  ///  Path to the file that could not be found.
  /// @param function
  ///  Name of the function that has thrown this exception.
  /// @param filename
  ///  Name of the file that contains the throw-statement.
  /// @param line
  ///  Line number of the throw-statement.
  EFileNotFound(const std::string& path, const std::string& function,
                const std::string& filename, size_t line);
};

/// This exception is thrown if a directory is not found.
class AURORA_API EDirectoryNotFound : public EInOut
{
public:
  /// Constructs a new EDirectoryNotFound exception.
  /// @param path
  ///  Path to the file that could not be found.
  /// @param function
  ///  Name of the function that has thrown this exception.
  /// @param filename
  ///  Name of the file that contains the throw-statement.
  /// @param line
  ///  Line number of the throw-statement.
  EDirectoryNotFound(const std::string& path, const std::string& function,
                     const std::string& filename, size_t line);
};

/// Macro to throw Exceptions derived from Aurora::Exception. It provides the
/// constructor with the function name, filename and line number of the
/// throw-statement.
#define AURORA_THROW(Type, message) \
  throw Type(message, __FUNCTION__, __FILE__, __LINE__)

} // namespace Aurora

#endif
