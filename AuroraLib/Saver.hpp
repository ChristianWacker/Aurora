//--- Aurora/AuroraLib/Saver.hpp -----------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_SAVER_HPP
#define AURORA_AURORA_LIB_SAVER_HPP

#include "AuroraLib/Exception.hpp"

#include <map>

namespace Aurora
{

/// Class to save an object of type T to a file
template<class T, class... Params>
class Saver
{
public:
  virtual ~Saver() { }

  virtual void save(const T& object, Params... params,
                    const std::string& filename) = 0;
};

template<class T, class... Params>
class Savable
{
public:
  static void registerSaver(const std::string& extension,
                            const std::shared_ptr<Saver<T, Params...>>& saver)
  {
    Savers& savers = T::savers();

    if (!savers.insert(std::make_pair(extension, saver)).second)
    {
      AURORA_THROW(EInvalidParameter,
                   "\"" + extension + "\" has already been registered.");
    }
  }

  void save(const std::string& filename, Params... params) const
  {
    size_t pos = filename.find_last_of(".");
    if ((std::string::npos == pos) || (filename.length() - 1 == pos))
      AURORA_THROW(EInvalidParameter, "Invalid filename.");

    std::string extension = filename.substr(pos + 1);

    Savers& savers = T::savers();
    auto saver = savers.find(extension);
    if (savers.end() == saver)
    {
      AURORA_THROW(EInvalidParameter,
                   "Unsupported filetype: \"" + extension + "\"");
    }

    saver->second->save(*static_cast<const T*>(this),
                        std::forward<Params>(params)..., filename);
  }

protected:
  typedef std::map<std::string, std::shared_ptr<Saver<T, Params...>>> Savers;
};

}

#endif
