#ifndef TRACEON_FILEREADER_H
#define TRACEON_FILEREADER_H

#include <string>
#include <memory> // Required for std::unique_ptr

namespace TracEon {

    /**
     * @class FileReader
     * @brief A unified reader for both plain text and gzipped files.
     *
     * This class abstracts away the complexity of handling different file types.
     * It provides a single, simple interface to read files line-by-line,
     * automatically detecting whether the input file is compressed with gzip.
     */
    class FileReader {
    public:
        // Constructor: Opens the specified file.
        // Throws std::runtime_error if the file cannot be opened.
        explicit FileReader(const std::string& filepath);

        // Destructor: Automatically closes the file.
        ~FileReader();

        // Prevent copying and assignment to avoid resource management issues.
        FileReader(const FileReader&) = delete;
        FileReader& operator=(const FileReader&) = delete;

        /**
         * @brief Reads the next line from the file.
         * @param line A string to store the read line.
         * @return true if a line was successfully read, false otherwise (e.g., end of file).
         */
        bool getline(std::string& line);

        /**
         * @brief Checks if the file is successfully opened.
         * @return true if the file is open and ready for reading.
         */
        bool is_open() const;

    private:
        // Forward-declaration of the implementation struct (PImpl idiom).
        // This hides all private members from the header file.
        struct FileReaderImpl;
        std::unique_ptr<FileReaderImpl> pimpl_;
    };

} // namespace TracEon

#endif // TRACEON_FILEREADER_H

