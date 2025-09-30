#include "FileReader.h"
#include <zlib.h>
#include <fstream>
#include <vector>

namespace TracEon {

// The actual implementation is hidden inside this struct.
// The header file only knows that this struct exists, not what's inside it.
struct FileReader::FileReaderImpl {
    gzFile gz_file = nullptr;
    std::ifstream regular_file;
    bool is_gzipped = false;

    // Destructor to clean up resources
    ~FileReaderImpl() {
        if (gz_file) {
            gzclose(gz_file);
        }
    }
};

// --- Public Method Implementations ---

FileReader::FileReader(const std::string& filepath)
    : pimpl_(std::make_unique<FileReaderImpl>()) {

    // Check if the filepath ends with .gz
    if (filepath.size() > 3 && filepath.substr(filepath.size() - 3) == ".gz") {
        pimpl_->is_gzipped = true;
        pimpl_->gz_file = gzopen(filepath.c_str(), "rb");
        if (!pimpl_->gz_file) {
            throw std::runtime_error("FileReader Error: Cannot open gzipped file: " + filepath);
        }
    } else {
        pimpl_->is_gzipped = false;
        pimpl_->regular_file.open(filepath);
        if (!pimpl_->regular_file.is_open()) {
            throw std::runtime_error("FileReader Error: Cannot open file: " + filepath);
        }
    }
}

FileReader::~FileReader() = default; // The unique_ptr will handle cleanup automatically.

bool FileReader::is_open() const {
    if (!pimpl_) return false;
    return pimpl_->is_gzipped ? (pimpl_->gz_file != nullptr) : pimpl_->regular_file.is_open();
}

bool FileReader::getline(std::string& line) {
    line.clear();
    if (!is_open()) {
        return false;
    }

    if (pimpl_->is_gzipped) {
        // For gzipped files, we need a buffer to read into.
        const int BUFFER_SIZE = 8192;
        char buffer[BUFFER_SIZE];

        if (gzgets(pimpl_->gz_file, buffer, BUFFER_SIZE) == nullptr) {
            return false; // End of file or error
        }
        line = buffer;
    } else {
        // For regular files, std::getline is sufficient.
        if (!std::getline(pimpl_->regular_file, line)) {
            return false;
        }
    }

    // Clean up newline characters from the end of the line for both cases.
    while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) {
        line.pop_back();
    }

    return true;
}

} // namespace TracEon

