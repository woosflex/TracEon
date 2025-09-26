#include "Cache.h"

// We need to define the namespace we declared in the header file.
namespace TracEon {

    // The constructor. For now, it doesn't need to do anything.
    Cache::Cache() {
        // We'll initialize our strategies here later.
    }

    // The destructor. Also doesn't need to do anything yet.
    Cache::~Cache() = default;

    // The implementation for the size() method.
    // It just returns the current number of elements in our map.
    size_t Cache::size() const {
        return m_store.size();
    }

} // namespace TracEon
//
// Created by Adnan Raza on 9/26/2025.
//