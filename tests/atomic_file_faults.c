#define _GNU_SOURCE
#include <dlfcn.h>
#include <errno.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

static int mode(const char *expected) {
    const char *value = getenv("TRACEON_TEST_ATOMIC_FAULT");
    return value && strcmp(value, expected) == 0;
}

int fsync(int fd) {
    int (*real_fsync)(int) = (int (*)(int))dlsym(RTLD_NEXT, "fsync");
    struct stat status;
    if (fstat(fd, &status) == 0) {
        if ((S_ISREG(status.st_mode) && mode("file-sync")) ||
            (S_ISDIR(status.st_mode) && mode("directory-sync"))) {
            errno = EIO;
            return -1;
        }
        if (S_ISDIR(status.st_mode) && mode("directory-sync-einval")) {
            errno = EINVAL;
            return -1;
        }
        if (S_ISREG(status.st_mode) && mode("stop-before-rename")) {
            int result = real_fsync(fd);
            raise(SIGSTOP);
            return result;
        }
    }
    return real_fsync(fd);
}

int rename(const char *old, const char *new_path) {
    int (*real_rename)(const char *, const char *) =
        (int (*)(const char *, const char *))dlsym(RTLD_NEXT, "rename");
    if (mode("rename")) { errno = EIO; return -1; }
    return real_rename(old, new_path);
}

size_t fwrite(const void *ptr, size_t size, size_t count, FILE *stream) {
    size_t (*real_fwrite)(const void *, size_t, size_t, FILE *) =
        (size_t (*)(const void *, size_t, size_t, FILE *))dlsym(RTLD_NEXT, "fwrite");
    struct stat status;
    if (mode("write") && fstat(fileno(stream), &status) == 0 && S_ISREG(status.st_mode)) {
        errno = EIO;
        return 0;
    }
    return real_fwrite(ptr, size, count, stream);
}
