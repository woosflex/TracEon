#ifndef TRACEON_MAPDEFS_H
#define TRACEON_MAPDEFS_H

#include <ankerl/unordered_dense.h>
#include <string>
#include <string_view>
#include <functional>

namespace TracEon {

// OPTIMIZATION: Use ankerl::unordered_dense for better cache locality and insert performance
template<typename Key, typename Value>
using HashMap = ankerl::unordered_dense::map<Key, Value>;

// Transparent hasher: enables heterogeneous lookup (find by string_view,
// no std::string allocation on the hot read path).
struct TransparentStringHash {
    using is_transparent = void;
    size_t operator()(std::string_view sv) const noexcept {
        return std::hash<std::string_view>{}(sv);
    }
    size_t operator()(const std::string& s) const noexcept {
        return std::hash<std::string_view>{}(s);
    }
};

// StringHashMap: std::string keys, but find/count/contains accept string_view.
template<typename Value>
using StringHashMap = ankerl::unordered_dense::map<
    std::string, Value, TransparentStringHash, std::equal_to<>>;

} // namespace TracEon

#endif // TRACEON_MAPDEFS_H