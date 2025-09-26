#include <iostream>
#include "Cache.h"

int main() {
    std::cout << ">> TracEon Test Executable <<" << std::endl;
    TracEon::Cache cache;
    std::cout << "Cache object created successfully." << std::endl;
    std::cout << "Initial number of items in cache: " << cache.size() << std::endl;
    return 0;
}

//
// Created by Adnan Raza on 9/26/2025.
//