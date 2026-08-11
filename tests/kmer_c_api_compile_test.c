/*
 * C99 compilation + linkage proof for the kmer C API boundary.
 *
 * This file is deliberately compiled as C (C_STANDARD 99) and linked
 * against the C++ traceon_kmer static library. If kmerindex_c_api.h ever
 * leaked a C++-only construct (references, templates, namespaces, member
 * access on the opaque handles, ...), this translation unit would fail to
 * compile; if the ABI (symbol linkage, opaque struct layout, exceptions
 * crossing the boundary) broke, it would fail to link or would crash.
 *
 * The checks exercised here mirror the Catch2 suite (happy path, minimap2
 * ignore-bit0 semantics, freeze enforcement, error model) so the C-facing
 * behavior is verified from an actual C caller's perspective.
 */
#include "kmerindex_c_api.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>

static int fail(const char* what, int code) {
    fprintf(stderr, "kmer_c_compile_test: FAIL: %s (last_error='%s')\n",
            what, kmerindex_last_error(NULL));
    return code;
}

int main(void) {
    kmerindex_t* h = kmerindex_create();
    if (!h) return fail("create returned NULL", 1);

    /* reserve: happy path + null-handle */
    if (kmerindex_reserve(h, 16) != 1) return fail("reserve failed", 2);
    if (kmerindex_reserve(NULL, 16) != 0) return fail("reserve(NULL) must return 0", 3);
    if (strlen(kmerindex_last_error(h)) == 0)
        return fail("reserve(NULL) must set last_error", 4);

    /* insert with minimap2's final-key semantics: bit0 = singleton flag */
    if (kmerindex_insert(h, (uint64_t)2 | 1u, 10) != 1)
        return fail("insert(final_key=3) failed", 5);
    /* same base key (>>1 == 1) -> collision, returns 0 */
    if (kmerindex_insert(h, 2, 99) != 0)
        return fail("collision insert must return 0", 6);
    if (kmerindex_size(h) != 1) return fail("size != 1 after inserts", 7);

    /* get: query_key with bit0=0 must find the entry stored with bit0=1 */
    uint64_t matched = 0;
    const uint64_t* v = kmerindex_get(h, 2, &matched);
    if (!v || *v != 10) return fail("get(2) must return value 10", 8);
    if (matched != 3) return fail("matched_key must be the stored final_key (3)", 9);

    /* get froze the index -> insert must now fail with 'frozen' error */
    if (kmerindex_insert(h, 4, 20) != 0)
        return fail("insert after get must return 0 (frozen)", 10);
    if (strstr(kmerindex_last_error(h), "frozen") == NULL)
        return fail("frozen error message expected", 11);

    /* caller-owned iteration over the frozen index */
    kmerindex_iter_t it;
    if (kmerindex_iter_begin(h, &it) != 1) return fail("iter_begin failed", 12);
    uint64_t k = 0, val = 0;
    int n = 0;
    while (kmerindex_iter_next(&it, &k, &val) == 1) ++n;
    if (n != 1) return fail("iteration must visit exactly 1 entry", 13);
    if (k != 3 || val != 10) return fail("iteration must yield (3, 10)", 14);

    /* iter_next on an uninitialized iterator must fail cleanly, not crash */
    kmerindex_iter_t uninit;
    memset(&uninit, 0, sizeof(uninit));
    if (kmerindex_iter_next(&uninit, &k, &val) != 0)
        return fail("iter_next on uninitialized iterator must return 0", 15);

    /* explicit freeze + null-handle sweep */
    kmerindex_t* h2 = kmerindex_create();
    if (!h2) return fail("second create returned NULL", 16);
    if (kmerindex_freeze(h2) != 1) return fail("freeze failed", 17);
    if (kmerindex_insert(h2, 8, 30) != 0) return fail("insert after freeze must return 0", 18);
    kmerindex_destroy(h2);

    if (kmerindex_size(NULL) != 0) return fail("size(NULL) must return 0", 19);
    if (kmerindex_get(NULL, 1, &matched) != NULL) return fail("get(NULL) must return NULL", 20);
    if (kmerindex_iter_begin(h, NULL) != 0) return fail("iter_begin(NULL out) must return 0", 21);
    if (kmerindex_iter_next(NULL, &k, &val) != 0) return fail("iter_next(NULL) must return 0", 22);
    if (kmerindex_get(h, 1, NULL) != NULL) return fail("get(NULL out) must return NULL", 23);

    kmerindex_destroy(h);
    puts("kmer_c_compile_test: OK (C99 compile + link + runtime contract)");
    return 0;
}
