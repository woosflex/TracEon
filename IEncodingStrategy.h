#ifndef TRACEON_IENCODINGSTRATEGY_H
#define TRACEON_IENCODINGSTRATEGY_H

#include <string>
#include <vector>

// An interface (abstract base class) that all encoding strategies must implement.
class IEncodingStrategy {
public:
    virtual ~IEncodingStrategy() = default;

    // Takes a raw string and returns a compressed/encoded vector of bytes.
    virtual std::vector<unsigned char> encode(const std::string& data) const = 0;
    // Takes a compressed/encoded vector of bytes and returns the original raw string.
    virtual std::string decode(const std::vector<unsigned char>& data) const = 0;
};

#endif //TRACEON_IENCODINGSTRATEGY_H