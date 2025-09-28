#ifndef TRACEON_QUALITYSCORERLESTRATEGY_H
#define TRACEON_QUALITYSCORERLESTRATEGY_H

#include "IEncodingStrategy.h"

// A strategy that compresses quality scores using Run-Length Encoding (RLE).
class QualityScoreRleStrategy : public IEncodingStrategy {
public:
    ~QualityScoreRleStrategy() override;
    std::vector<unsigned char> encode(const std::string& data) const override;
    std::string decode(const std::vector<unsigned char>& data) const override;
};

#endif //TRACEON_QUALITYSCORERLESTRATEGY_H