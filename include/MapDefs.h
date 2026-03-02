#ifndef TRACEON_MAPDEFS_H
#define TRACEON_MAPDEFS_H

// Check if the robin_hood header is available
#if __has_include("third_party/robin_hood.h")
    #include "third_party/robin_hood.h"
    #define TRACEON_USE_ROBIN_HOOD
#endif

#include <unordered_map>
#include <string>
#include <string_view>
#include <functional>

namespace TracEon {

#ifdef TRACEON_USE_ROBIN_HOOD
    // OPTIMIZATION: Use flat_map for better cache locality (open addressing)
    template<typename Key, typename Value>
    using HashMap = robin_hood::unordered_flat_map<Key, Value>;

    // Transparent hasher: enables heterogeneous lookup (find by string_view,
    // no std::string allocation on the hot read path).
    struct TransparentStringHash {
        using is_transparent = void;
        size_t operator()(std::string_view sv) const noexcept {
            return robin_hood::hash<std::string_view>{}(sv);
        }
        size_t operator()(const std::string& s) const noexcept {
            return robin_hood::hash<std::string_view>{}(s);
        }
    };

    // StringHashMap: std::string keys, but find/count/contains accept string_view.
    template<typename Value>
    using StringHashMap = robin_hood::unordered_flat_map<
        std::string, Value, TransparentStringHash, std::equal_to<>>;

#else
    template<typename Key, typename Value>
    using HashMap = std::unordered_map<Key, Value>;

    template<typename Value>
    using StringHashMap = std::unordered_map<std::string, Value>;
#endif

} // namespace TracEon

#endif // TRACEON_MAPDEFS_H