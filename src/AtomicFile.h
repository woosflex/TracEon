#pragma once

#include <cerrno>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <stdexcept>
#include <string>
#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#include <share.h>
#include <sys/stat.h>
#include <random>
#else
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

namespace TracEon {

// A complete, synchronized file replaces the destination at a single commit point.
// No destination truncation occurs on serialization, flush, or pre-rename failure.
class AtomicFile {
public:
    explicit AtomicFile(const std::string& destination) : target_(destination) {
        try {
        if (destination.empty() || target_.filename().empty())
            throw std::invalid_argument("cache destination must name a file");
        auto parent = target_.parent_path();
        if (parent.empty()) parent = ".";
        std::error_code error;
        const auto status = std::filesystem::symlink_status(target_, error);
        if (error && error != std::errc::no_such_file_or_directory)
            throw std::filesystem::filesystem_error("inspect cache destination", target_, error);
        if (std::filesystem::exists(status) && !std::filesystem::is_regular_file(status))
            throw std::runtime_error("cache destination must be a regular file (not a symlink/device)");
        int& fd = descriptor_;
#ifdef _WIN32
        std::random_device random;
        for (int i = 0; i < 32 && fd < 0; ++i) {
            temporary_ = parent / (".traceon-tmp-" + std::to_string(random()) + std::to_string(random()));
            const auto result = _wsopen_s(&fd, temporary_.c_str(), _O_RDWR | _O_CREAT | _O_EXCL | _O_BINARY,
                                         _SH_DENYRW, _S_IREAD | _S_IWRITE);
            if (result != 0 && result != EEXIST) break;
        }
        if (fd < 0) fail("create temporary cache");
        temporary_created_ = true;
        file_ = _fdopen(fd, "w+b");
        if (!file_) fail("open temporary cache stream");
        descriptor_ = -1; // The stdio stream now owns the descriptor.
#else
        temporary_name_ = (parent / ".traceon-tmp-XXXXXX").string();
        directory_ = ::open(parent.c_str(), O_RDONLY | O_DIRECTORY | O_CLOEXEC);
        if (directory_ < 0) fail("open cache directory");
        fd = ::mkstemp(temporary_name_.data());
        if (fd < 0) fail("create temporary cache");
        temporary_created_ = true;
        if (::fcntl(fd, F_SETFD, FD_CLOEXEC) < 0) fail("set temporary cache close-on-exec");
        if (std::filesystem::exists(status) &&
            ::fchmod(fd, static_cast<mode_t>(status.permissions()) & 0777) != 0)
            fail("preserve cache permissions");
        file_ = ::fdopen(fd, "w+b");
        if (!file_) fail("open temporary cache stream");
        descriptor_ = -1;
#endif
        } catch (...) { cleanup(); throw; }
    }
    AtomicFile(const AtomicFile&) = delete;
    AtomicFile& operator=(const AtomicFile&) = delete;
    ~AtomicFile() { cleanup(); }

    void write(const char* bytes, std::size_t length) {
        if (length && std::fwrite(bytes, 1, length, file_) != length) fail("write cache");
    }
    std::uint64_t tell() {
#ifdef _WIN32
        const auto position = _ftelli64(file_);
#else
        const auto position = ::ftello(file_);
#endif
        if (position < 0) fail("tell cache position");
        return static_cast<std::uint64_t>(position);
    }
    void seek(std::uint64_t position) {
#ifdef _WIN32
        if (_fseeki64(file_, static_cast<__int64>(position), SEEK_SET) != 0) fail("seek cache");
#else
        if (::fseeko(file_, static_cast<off_t>(position), SEEK_SET) != 0) fail("seek cache");
#endif
    }
    void commit() {
        if (std::fflush(file_) != 0) fail("flush cache");
#ifdef _WIN32
        if (_commit(_fileno(file_)) != 0) fail("synchronize cache data");
#else
        if (::fsync(::fileno(file_)) != 0) fail("synchronize cache data");
#endif
        const int closed = std::fclose(file_);
        file_ = nullptr;
        if (closed != 0) fail("close cache");
#ifdef _WIN32
        std::filesystem::rename(temporary_, target_);
        temporary_created_ = false;
#else
        if (::rename(temporary_name_.c_str(), target_.c_str()) != 0) fail("replace cache");
        temporary_created_ = false; // Commit point: destination is now complete.
        // Filesystems like NFS, CIFS, and FAT/exFAT do not implement directory fsync.
        if (::fsync(directory_) != 0 && errno != EINVAL && errno != ENOTSUP &&
            errno != EOPNOTSUPP && errno != EROFS && errno != EBADF)
            fail("cache replacement completed, but directory durability is uncertain");
#endif
    }

private:
    [[noreturn]] void fail(const char* operation) const {
        throw std::runtime_error(std::string(operation) + ": " + target_.string() + ": " + std::strerror(errno));
    }
    void cleanup() noexcept {
        if (file_) { std::fclose(file_); file_ = nullptr; }
#ifdef _WIN32
        if (descriptor_ >= 0) { _close(descriptor_); descriptor_ = -1; }
        if (temporary_created_) {
            std::error_code error;
            std::filesystem::remove(temporary_, error);
            temporary_created_ = false;
        }
#else
        if (descriptor_ >= 0) { ::close(descriptor_); descriptor_ = -1; }
        if (temporary_created_) { ::unlink(temporary_name_.c_str()); temporary_created_ = false; }
        if (directory_ >= 0) { ::close(directory_); directory_ = -1; }
#endif
    }
    std::filesystem::path target_;
    std::FILE* file_ = nullptr;
    int descriptor_ = -1;
    bool temporary_created_ = false;
#ifdef _WIN32
    std::filesystem::path temporary_;
#else
    std::string temporary_name_;
    int directory_ = -1;
#endif
};
} // namespace TracEon
