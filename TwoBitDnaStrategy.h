//
// Created by Adnan Raza on 9/26/2025.
//

#ifndef TRACEON_TWOBITDNASTRATEGY_H
#define TRACEON_TWOBITDNASTRATEGY_H

#include "IEncodingStrategy.h"

class TwoBitDnaStrategy : public IEncodingStrategy {
public:
    std::vector<unsigned char> encode(const std::string& data) const override;
    std::string decode(const std::vector<unsigned char>& data) const override;
};

#endif //TRACEON_TWOBITDNASTRATEGY_H