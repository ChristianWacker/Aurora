//--- Aurora/AuroraLib/CommandLine.cpp -----------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#include "AuroraLib/CommandLine.hpp"

#include "AuroraLib/Utils.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>

namespace Aurora
{

ECommandLineError::ECommandLineError(const std::string& message,
                                     const std::string& function,
                                     const std::string& filename,
                                     size_t line) :
  Exception(message, function, filename, line, "ECommandLineError")
{ }

namespace CL
{

OptionBase::Registry& OptionBase::registry()
{
  static Registry registry;
  return registry;
}

void OptionBase::help(const std::string& description)
{
  if (!description.empty())
  {
    std::cerr << description << '\n';
    std::cerr << '\n';
  }
  std::cerr << "USAGE\n";
  std::cerr << std::left;

  for (auto i : registry())
  {
    std::cerr << " -"  << std::setw(18) << i.first
              <<          std::setw(20) << i.second->string()
              <<                           i.second->description() << '\n';
  }
}

std::string OptionBase::configString()
{
  std::string result;
  for (auto i : registry())
    result += i.first + " " + i.second->string() + "\n";
  return result;
}

void parseConfigFile(const std::string& filename)
{
  std::ifstream inFile(filename);

  if (!inFile)
    AURORA_THROW(EFileNotFound, filename);

  int lineNumber = 0;

  // Read the file line by line
  while (inFile)
  {
    ++lineNumber;

    std::string line;
    getline(inFile, line);

    // look for an comment
    const size_t hashPosition = line.find_first_of('#');

    // remove the comment
    if (std::string::npos != hashPosition)
      line = line.substr(0, hashPosition);

    line = trim(line);

    // skip empty lines
    if (line.length() == 0)
      continue;

    size_t pos = line.find(' ');
    if (std::string::npos == pos)
    {
      AURORA_THROW(EInOut, "Config error, line: " + std::to_string(lineNumber) +
                           "): Space ' ' expected");
    }

    std::string name = line.substr(0, pos);
    std::string value = line.substr(pos + 1);
    OptionBase::parse(name, value);
  }
}

bool parse(int argc, char** argv, const std::string& description)
{
  try
  {
    // start with i = 1 since the first argument is the program name
    for (int i = 1; i < argc; ++i)
    {
      // convert c-string into std::string
      const std::string argument(argv[i]);

      if ('-' != argument[0])
        AURORA_THROW(ECommandLineError, "Options must begin with '-'");

      // remove the leading '-'
      std::string name = argument.substr(1, argument.length() - 1);

      if ("-help" == name)
      {
        OptionBase::help(description);
        return false;
      }

      ++i;
      if (i >= argc)
      {
        AURORA_THROW(ECommandLineError,
                     "No value specified for \"" + name + "\"");
      }

      const std::string value = argv[i];

      if ("configFile" == name)
      {
        parseConfigFile(value);
        continue;
      }

      OptionBase::parse(name, value);
    }
  }
  catch (const Exception& e)
  {
    OptionBase::help(description);
    throw e;
  }

  return true;
}

std::string configString()
{
  return OptionBase::configString();
}

} // namespace CommandLine

CL::Option<int64_t> verbosity("verbosity", "", 1);

} // namespace Aurora
