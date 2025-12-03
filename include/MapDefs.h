#ifndef TRACEON_MAPDEFS_H
#define TRACEON_MAPDEFS_H

// Check if the robin_hood header is available
#if __has_include("third_party/robin_hood.h")
    #include "third_party/robin_hood.h"
    #define TRACEON_USE_ROBIN_HOOD
#endif

#include <unordered_map>

namespace TracEon {

#ifdef TRACEON_USE_ROBIN_HOOD
    // OPTIMIZATION: Use flat_map for better cache locality (open addressing)
    template<typename Key, typename Value>
    using HashMap = robin_hood::unordered_flat_map<Key, Value>;
#else
    template<typename Key, typename Value>
    using HashMap = std::unordered_map<Key, Value>;
#endif

} // namespace TracEon

#endif // TRACEON_MAPDEFS_H