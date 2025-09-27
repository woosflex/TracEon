#ifndef TRACEON_RECORDTYPES_H
#define TRACEON_RECORDTYPES_H

#include <string>
#include <vector>

namespace TracEon {

    // Defines the available encoding strategies for the Cache.
    enum class StrategyType {
        PlainText,
        TwoBit
    };

    // A structure to hold the components of a FASTQ record.
    // This will be used when we implement FASTQ parsing.
    struct FastqRecord {
        std::string identifier;
        std::vector<unsigned char> compressed_sequence;
        std::vector<unsigned char> compressed_quality;
    };

} // namespace TracEon

#endif //TRACEON_RECORDTYPES_H