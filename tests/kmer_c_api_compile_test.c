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

    /* reads still work on the frozen index; miss freezes too */
    v = kmerindex_get(h, 2, &matched);
    if (!v || *v != 10) return fail("get after freeze must still hit", 12);

    if (kmerindex_size(NULL) != 0) return fail("size(NULL) must return 0", 13);
    if (kmerindex_get(NULL, 1, &matched) != NULL) return fail("get(NULL) must return NULL", 14);
    if (kmerindex_get(h, 1, NULL) != NULL) return fail("get(NULL out) must return NULL", 15);

    kmerindex_destroy(h);
    puts("kmer_c_compile_test: OK (C99 compile + link + runtime contract)");
    return 0;
}
