//--- Aurora/AuroraLib/CommandLine.hpp -----------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_COMMAND_LINE_HPP
#define AURORA_AURORA_LIB_COMMAND_LINE_HPP

#include "AuroraLib/Exception.hpp"

#include <map>
#include <sstream>
#include <vector>

namespace Aurora
{

/// Exception class that is thrown, if the command line cannot be parsed.
class AURORA_API ECommandLineError : public Exception
{
public:
  /// @copydoc EInvalidParameter::EInvalidParameter
  ECommandLineError(const std::string& message, const std::string& function,
                    const std::string& filename, size_t line);
};

/// namespace for the command line
namespace CL
{

/// General parser class. The standard implementation uses a stringstream for
/// for parsing.
template<class DataType>
class Parser
{
public:
  static DataType parse(const std::string& s)
  {
    std::stringstream stream(s);
    DataType value;
    stream >> value;
    return value;
  }

  static std::string toString(DataType value)
  {
    return std::to_string(value);
  }
};

/// Specialization of the Parser class for std::string. The values are simply
/// passed through.
template<>
class Parser<std::string>
{
public:
  static std::string parse(const std::string& s)
  {
    return s;
  }

  static std::string toString(const std::string& value)
  {
    return value;
  }
};

/// Concrete base class for all option classes
class AURORA_API OptionBase
{
public:
  OptionBase(const std::string& name, const std::string& description) :
    mName(name), mDescription(description)
  {
    registerOption(name, this);
  }

  virtual ~OptionBase() {}

  /// Write a help text to stderr
  static void help(const std::string& description);
  /// Return a string with all command line options and their values.
  static std::string configString();

  /// Parse a name, value pair
  static void parse(const std::string& name, const std::string& value)
  {
    auto it = registry().find(name);

    if (registry().end() == it)
      AURORA_THROW(EInvalidParameter, "Option \"" + name + "\" is unknown.");

    it->second->parse(value);
  }

  /// Set a decsription for the option
  void setDescription(const std::string& description)
  {
    mDescription = description;
  }

  /// Return the name of the option.
  std::string name() const
  {
    return mName;
  }

  /// Return the description of the option
  std::string description() const
  {
    return mDescription;
  }

  virtual void parse(const std::string& value) = 0;
  virtual std::string string() = 0;

private:
  std::string mName;
  std::string mDescription;

  typedef std::map<std::string, OptionBase*> Registry;
  static Registry& registry();

  static void registerOption(const std::string& name, OptionBase* option)
  {
    if (!registry().insert(std::make_pair(name, option)).second)
    {
      AURORA_THROW(EInvalidParameter,
                   "Option \"" + name + "\" has already been registered");
    }
  }
};

/// Storage for Options value. Non-class variant.
template<class DataType, bool isClass>
class OptionStorage
{
public:
  OptionStorage(const DataType& initialValue) :
    mValue(initialValue)
  { }

  operator DataType() const
  {
    return mValue;
  }

  operator DataType&()
  {
    return mValue;
  }

  void setValue(DataType value)
  {
    mValue = value;
  }

  DataType value()
  {
    return mValue;
  }

private:
  DataType mValue;
};

template<class DataType>
class OptionStorage<DataType, true> : public DataType
{
public:
  OptionStorage(const DataType& initialValue) :
    DataType(initialValue)
  { }

  void setValue(DataType value)
  {
    DataType::operator=(value);
  }

  DataType value()
  {
    return *this;
  }
};

/// OptionHelper for non enumeration types
template<class DataType, bool isEnum>
class OptionHelper :
  public OptionBase,
  public OptionStorage<DataType, std::is_class<DataType>::value>
{
public:
  OptionHelper(const std::string& name, const std::string& description,
      const DataType& initialValue = DataType()) :
    OptionBase(name, description),
    OptionStorage<DataType, std::is_class<DataType>::value>(initialValue)
  {}

  void parse(const std::string& value) override
  {
    this->setValue(Parser<DataType>::parse(value));
  }

  std::string string() override
  {
    return Parser<DataType>::toString(this->value());
  }
};

#define enumVal(value) #value, (value)

/// Specialization of OptionHelper to represent am enum
template<class DataType>
class OptionHelper<DataType, true> :
  public OptionBase, public OptionStorage<DataType, false>
{
public:
  template<class... Values>
  OptionHelper(const std::string& name, const std::string& description,
               const DataType& initialValue, Values... values) :
    OptionBase(name, description),
    OptionStorage<DataType, false>(initialValue)
  {
    addValue(values...);
  }

  void parse(const std::string& value) override
  {
    for (auto i : mPossibleValues)
      if (value == i.first)
      {
        this->setValue(i.second);
        return;
      }

    AURORA_THROW(ECommandLineError,
      "Option \"" + name() + "\" cannot be set to \"" + value + "\"");
  }

  std::string string() override
  {
    for (auto i : mPossibleValues)
      if (this->value() == i.second)
        return i.first;

    AURORA_UNREACHABLE;
  }

  std::vector<std::string> valueList() const
  {
    std::vector<std::string> result;
    for (auto i : mPossibleValues)
      result.push_back(i.first);

    return result;
  }

private:
  /// function to add recusively the possible values
  template<class... Values>
  void addValue(const std::string& s, DataType value, Values... values)
  {
    std::size_t pos = s.find_last_of(':');
    assert(std::string::npos != pos);
    mPossibleValues.push_back(std::make_pair(s.substr(pos + 1), value));
    addValue(values...);
  }

  // End of the recursion
  void addValue() {}

  std::vector<std::pair<std::string, DataType>> mPossibleValues;
};

template<class DataType>
class Option : public OptionHelper<DataType, std::is_enum<DataType>::value>
{
public:
  using OptionHelper<DataType, std::is_enum<DataType>::value>::OptionHelper;
};

AURORA_API void parseConfigFile(const std::string& filename);
AURORA_API bool parse(int argc, char** argv,
                      const std::string& description = std::string());
AURORA_API std::string configString();

} // namespace CommandLine

/// The current verbosity level.
AURORA_API extern CL::Option<int64_t> verbosity;

} // namespace Aurora

#endif
