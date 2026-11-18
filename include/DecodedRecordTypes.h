#ifndef TRACEON_DECODEDRECORDTYPES_H
#define TRACEON_DECODEDRECORDTYPES_H

#include <string>

namespace TracEon {

    // A simple struct to hold a fully decoded FASTQ record,
    // containing both the sequence and its quality string.
    struct DecodedFastqRecord {
        std::string sequence;
        std::string quality;
    };

    // Non-owning view of a record (sequence and quality views).
    struct FastqRecordView {
        std::string_view sequence;
        std::string_view quality;
    };

} // namespace TracEon

#endif //TRACEON_DECODEDRECORDTYPES_H