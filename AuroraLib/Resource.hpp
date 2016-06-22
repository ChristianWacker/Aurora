//--- Aurora/AuroraLib/Resource.hpp --------------------------------------------
//
//   This file is distributed under the University of Illinois Open Source
//   License. See LICENSE.TXT for details.
//
//------------------------------------------------------------------------------
///
/// Inspired by http://accu.org/index.php/journals/2086
///
//------------------------------------------------------------------------------

#ifndef AURORA_AURORA_LIB_RESOURCE_HPP
#define AURORA_AURORA_LIB_RESOURCE_HPP

#include "AuroraLib/Prerequisites.hpp"

namespace Aurora
{

/// Class for the management of resources, handles, files and so on. This class
/// is a transparent layer that ensure the release of the resource as soon as it
/// goes out of scope.
template<class ResourceTag, class ResourceType>
class UniqueResource
{
public:
  UniqueResource() noexcept = default;
  explicit UniqueResource(const ResourceType& resource) noexcept :
    mResource(resource)
  { }

  UniqueResource(const UniqueResource&) = delete;
  UniqueResource& operator=(const ResourceType) = delete;

  UniqueResource(UniqueResource&& other) noexcept :
    mResource(std::move(other.mResource))
  {
    other.mResource = ResourceType();
  }

  UniqueResource& operator=(UniqueResource&& other) noexcept
  {
    assert(this != std::addressof(other) && "Self-assignment");
    release();
    mResource = std::move(other.mResource);
    other.mResource = ResourceType();
    return *this;
  }

  ~UniqueResource()
  {
    release();
  }

  operator const ResourceType&() const noexcept { return mResource; }

private:
  static constexpr bool allwaysFalse() noexcept { return false; }

  void release() noexcept
  {
    static_assert(allwaysFalse(),
                  "This function must be explicitly specialized");
  }

  ResourceType mResource = ResourceType();
};

} // namespace Aurora

#endif
