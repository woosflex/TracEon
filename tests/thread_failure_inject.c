// Linux integration probe: a real loader must catch failure of its second launch.
#define _GNU_SOURCE
#include <pthread.h>
#include <dlfcn.h>
#include <errno.h>
#include <stdatomic.h>

int pthread_create(pthread_t* t, const pthread_attr_t* a,
                   void* (*fn)(void*), void* arg) {
    static atomic_int count;
    if (atomic_fetch_add(&count, 1) == 1) return EAGAIN;
    int (*real_create)(pthread_t*, const pthread_attr_t*, void* (*)(void*), void*) =
        dlsym(RTLD_NEXT, "pthread_create");
    return real_create(t, a, fn, arg);
}
