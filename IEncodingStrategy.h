//
// Created by Adnan Raza on 9/26/2025.
//

#ifndef TRACEON_IENCODINGSTRATEGY_H
#define TRACEON_IENCODINGSTRATEGY_H

#include <string>
#include <vector>

class IEncodingStrategy {
public:
    virtual ~IEncodingStrategy() = default;
    virtual std::vector<unsigned char> encode(const std::string& data) const = 0;
    virtual std::string decode(const std::vector<unsigned char>& data) const = 0;
};

#endif //TRACEON_IENCODINGSTRATEGY_H