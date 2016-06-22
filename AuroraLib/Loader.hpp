//--- Aurora/AuroraLib/Loader.hpp --------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_LOADER_HPP
#define AURORA_AURORA_LIB_LOADER_HPP

#include "AuroraLib/Exception.hpp"

#include <map>

namespace Aurora
{

/// Class to load an object of type T from a file
template<class T, class... Params>
class Loader
{
public:
  virtual ~Loader() { }

  virtual void load(const std::string& filename, T& t, Params... params) = 0;
};

template<class T, class... Params>
class Loadable
{
public:
  static void registerLoader(const std::string& extension,
                            const std::shared_ptr<Loader<T, Params...>>& loader)
  {
    Loaders& loaders = T::loaders();

    if (!loaders.insert(std::make_pair(extension, loader)).second)
    {
      AURORA_THROW(EInvalidParameter,
                   "\"" + extension + "\" has already been registered.");
    }
  }

  void load(const std::string& filename, Params... params)
  {
    size_t pos = filename.find_last_of(".");
    if ((std::string::npos == pos) || (filename.length() - 1 == pos))
      AURORA_THROW(EInvalidParameter, "Invalid filename.");

    std::string extension = filename.substr(pos + 1);

    Loaders& loaders = T::loaders();
    auto loader = loaders.find(extension);
    if (loaders.end() == loader)
    {
      AURORA_THROW(EInvalidParameter,
                   "Unsupported filetype \"" + extension + "\".");
    }

    loader->second->load(filename, *static_cast<T*>(this),
                         std::forward<Params>(params)...);
  }

protected:
  typedef std::map<std::string, std::shared_ptr<Loader<T, Params...>>> Loaders;
};

}

#endif
