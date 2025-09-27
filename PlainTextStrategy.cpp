#ifndef TRACEON_PLAINTEXTSTRATEGY_H
#define TRACEON_PLAINTEXTSTRATEGY_H

#include "IEncodingStrategy.h"

// A simple strategy that performs no compression.
// It just converts between string and a vector of its bytes.
class PlainTextStrategy : public IEncodingStrategy {
public:
    std::vector<unsigned char> encode(const std::string& data) const override;
    std::string decode(const std::vector<unsigned char>& data) const override;
};

#endif //TRACEON_PLAINTEXTSTRATEGY_H