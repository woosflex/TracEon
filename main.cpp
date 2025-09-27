#include <iostream>
#include "Cache.h"

// This is a simple executable for manual testing or demonstration.
// It links against our traceon_core library.
int main() {
    std::cout << ">> TracEon Sanity Check Executable <<" << std::endl;

    TracEon::Cache cache;
    std::cout << "Cache object created successfully." << std::endl;
    std::cout << "Initial size: " << cache.size() << std::endl;

    cache.set("hello", "world");
    std::cout << "After set('hello', 'world'), size is: " << cache.size() << std::endl;
    std::cout << "get('hello') returns: " << cache.get("hello") << std::endl;

    return 0;
}