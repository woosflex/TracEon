#ifndef TRACEON_TWOBITDNASTRATEGY_H
#define TRACEON_TWOBITDNASTRATEGY_H

#include "IEncodingStrategy.h"

// A strategy that compresses DNA sequences using 2 bits per base.
class TwoBitDnaStrategy : public IEncodingStrategy {
public:
    ~TwoBitDnaStrategy() override = default;
    std::vector<unsigned char> encode(const std::string& data) const override;
    std::string decode(const std::vector<unsigned char>& data) const override;
};

#endif //TRACEON_TWOBITDNASTRATEGY_H