#ifndef TRACEON_TWOBITDNASTRATEGY_H
#define TRACEON_TWOBITDNASTRATEGY_H

#include "IEncodingStrategy.h"

class TwoBitDnaStrategy : public IEncodingStrategy {
public:
    ~TwoBitDnaStrategy() override;
    std::vector<unsigned char> encode(const std::string& data) const override;
    std::string decode(const std::vector<unsigned char>& data) const override;
};

#endif //TRACEON_TWOBITDNASTRATEGY_H