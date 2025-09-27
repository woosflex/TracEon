#ifndef TRACEON_PLAINTEXTSTRATEGY_H
#define TRACEON_PLAINTEXTSTRATEGY_H

#include "IEncodingStrategy.h"

class PlainTextStrategy : public IEncodingStrategy {
public:
    ~PlainTextStrategy() override;
    std::vector<unsigned char> encode(const std::string& data) const override;
    std::string decode(const std::vector<unsigned char>& data) const override;
};

#endif //TRACEON_PLAINTEXTSTRATEGY_H