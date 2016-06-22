//--- Aurora/AuroraLib/Factory.hpp --------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_FACTORY_HPP
#define AURORA_AURORA_LIB_FACTORY_HPP

#include "AuroraLib/Exception.hpp"

#include <functional>
#include <map>
#include <string>
#include <vector>

namespace Aurora
{

/// Template-class used to instantiate preregistered classes with a common
/// base class.
/// @remarks
///  The factory takes at least one template parameter specifying the base
///  class of the classes that can be instantiated. The other template
///  parameters are the parameters needed by the constructors of the child
/// classes.
/// @remarks
///  To instantiate a class one has to pass the name under which it has been
///  registered to the function @ref create(). The additional parameters of
///  create() are identical with the parameters needed by the constructors of
///  the child classes.
/// @remarks
///  To register a class one has to call @ref registerClass().
template<class Base, class... Parameters>
class Factory
{
public:
  /// The pointer type of the Base.
  typedef std::unique_ptr<Base> BasePtr;

  /// Registers a new class with the base class Base.
  /// @tparam Derived
  ///  The type of the class to be registered.
  /// @param name
  ///  The name under which the new class should be registered. If a class with
  ///  this name is already registered an EInvalidParameter exception will be
  ///  thrown.
  template<class Derived>
  static void registerClass(const std::string& name)
  {
    Creators& creators = Base::creators();

    // check if a class with this name has already been registered
    if (!creators.insert(std::make_pair(name,
                  static_cast<CreateFunction>(createFunction<Derived>))).second)
    {
      AURORA_THROW(EInvalidParameter,
                   "\"" + name + "\" has already been registered.");
    }
  }

  /// Create a new instance to the class specified by the name @p name.
  /// @param name
  ///  Specifies the class to be instantiated. A class with this name has to be
  ///  registered with a call to @ref registerClass(). If no matching class can
  ///  be found an EInvalidParameter exception will be thrown.
  /// @param parameters
  ///  Additional parameters passed through to the constructor of the class.
  /// @retval
  ///  A pointer to the newly instantiated class.
  static BasePtr create(const std::string& name, Parameters... parameters)
  {
    Creators& creators = Base::creators();

    // try to find the map entry corresponding to name
    auto it = creators.find(name);

    // If no map entry with this name can be found, raise an exception
    if (creators.end() == it)
      AURORA_THROW(EInvalidParameter, "\"" + name + "\" is unknown");

    // call the registered function to instantiate the class
    return it->second(std::forward<Parameters>(parameters)...);
  }

  /// Returns the names of all registered device implementations.
  static std::vector<std::string> namesList()
  {
    Creators& creators(Base::creators());

    std::vector<std::string> result;
    for (auto i : creators)
      result.push_back(i.first);

    return result;
  }

  static std::string namesString()
  {
    Creators& creators(Base::creators());

    std::string result = "{";
    for (auto i : creators)
      result += i.first + ", ";
    // delete the last comma
    result.erase(result.length() - 2, 2);
    result += "}";

    return result;
  }

protected:
  typedef std::function<BasePtr (Parameters...)> CreateFunction;
  typedef std::map<std::string, CreateFunction> Creators;

private:
  // Template function with the purpose of creating an instance of @p Derived.
  // All the parameters are passed through to the constructor of @p Derived.
  // Functors pointing to this function are stored in mCreators.
  template<class Derived>
  static BasePtr createFunction(Parameters... parameters)
  {
    return std::make_unique<Derived>(std::forward<Parameters>(parameters)...);
  }
};

} // namespace Aurora

#endif
