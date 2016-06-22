//--- Aurora/AuroraLib/CompilerAbstraction.hpp ---------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------
//
/// @file
/// @brief Header to hide and abstract compiler properties.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_COMPILER_ABSTRACTION_HPP
#define AURORA_AURORA_LIB_COMPILER_ABSTRACTION_HPP

// This header is meant to be included by Prerequisites.hpp. So do not include
// Prerequisites.hpp here.

// Determine Platform
#define AURORA_PLATFORM_WINDOWS 1
#define AURORA_PLATFORM_LINUX 2

#ifdef _WIN32
# define AURORA_PLATFORM AURORA_PLATFORM_WINDOWS
# if (AURORA_DYNAMIC_LIBRARY == 1)
#   if (AURORA_NONCLIENT_BUILD == 1)
#     define AURORA_API __declspec(dllexport)
#   else
#     define AURORA_API __declspec(dllimport)
#   endif
# else
#   define AURORA_API
# endif
# define AURORA_BREAK __debugbreak();

#else
# define AURORA_PLATFORM AURORA_PLATFORM_LINUX
# define AURORA_API
# define AURORA_BREAK __builtin_trap();

#endif


// Determine Compiler
#define AURORA_COMPILER_MSVC  1
#define AURORA_COMPILER_INTEL 2
#define AURORA_COMPILER_GCC   3
#define AURORA_COMPILER_CLANG 4

#if defined(__INTEL_COMPILER)
# define AURORA_COMPILER AURORA_COMPILER_INTEL
# if (AURORA_PLATFORM == AURORA_PLATFORM_LINUX)
// On Linux the Intel compiler tries to be compatible with gcc
#   define AURORA_COMPILER_IS_COMPATIBLE_WITH AURORA_COMPILER_GCC
# elif (AURORA_PLATFORM == AURORA_PLATFORM_WINDOWS)
// On Windows the Intel compiler tries to compatible with the Microsoft
// compiler
#   define AURORA_COMPILER_IS_COMPATIBLE_WITH AURORA_COMPILER_MSVC
# endif

#elif defined(__clang__)
# define AURORA_COMPILER AURORA_COMPILER_CLANG
# define AURORA_COMPILER_IS_COMPATIBLE_WITH AURORA_COMPILER_GCC

#elif defined(__GNUG__)
# define AURORA_COMPILER AURORA_COMPILER_GCC
# define AURORA_COMPILER_IS_COMPATIBLE_WITH AURORA_COMPILER_GCC

#elif defined(_MSC_VER)
# define AURORA_COMPILER AURORA_COMPILER_MSVC
# define AURORA_COMPILER_IS_COMPATIBLE_WITH AURORA_COMPILER_MSVC

// disable: "<type> needs to have dll-interface to be used by clients'
// false error for all standard templates
# pragma warning (disable : 4251)

// disable: "non-DLL-interface classkey 'identifier' used as base for
// DLL-interface classkey 'identifier'"
# pragma warning (disable : 4275)

// disable: not all control paths return a value
# pragma warning (disable : 4715)

//# pragma warning (disable : 4800)

// Disable warnings about potentially unsafe function in the Standard C++
// Library (SCL)
# define _SCL_SECURE_NO_WARNINGS
# define _CRT_SECURE_NO_WARNINGS 1
# define _CRT_NONSTDC_NO_DEPRECATE

#else
# error "Unknown Compiler"

#endif

#if (AURORA_COMPILER_IS_COMPATIBLE_WITH == AURORA_COMPILER_GCC)
# define AURORA_ASSUME_INTERN(expression) \
   do { if (!(expression)) __builtin_unreachable(); } while (0)
# define AURORA_CONSTEXPR constexpr

#endif // AURORA_COMPILER_IS_COMPATIBLE_WITH == AURORA_COMPILER_GCC

#if (AURORA_COMPILER_IS_COMPATIBLE_WITH == AURORA_COMPILER_MSVC)
# define AURORA_ASSUME_INTERN(expression) __assume(expression)
# define AURORA_CONSTEXPR

#endif

#endif
