#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdio.h>
#include <stdlib.h>
/* Process-local Pango setting: retain fractional glyph positions in vector PDF.
 * No changes to installed R, Pango, fonts, or the system environment. */
void *pango_cairo_create_layout(void *cr) {
    static void *(*create)(void *) = NULL;
    static void *(*context)(void *) = NULL;
    static void (*round_positions)(void *, int) = NULL;
    if (!create) {
        void *pc = dlopen("/home/xd/y/Benchmark/.conda-env/lib/libpangocairo-1.0.so.0", RTLD_LAZY | RTLD_LOCAL);
        void *pg = dlopen("/home/xd/y/Benchmark/.conda-env/lib/libpango-1.0.so.0", RTLD_LAZY | RTLD_LOCAL);
        create = pc ? dlsym(pc, "pango_cairo_create_layout") : NULL;
        context = pg ? dlsym(pg, "pango_layout_get_context") : NULL;
        round_positions = pg ? dlsym(pg, "pango_context_set_round_glyph_positions") : NULL;
        if (!create || !context || !round_positions) {
            fputs("Required Pango fractional positioning API unavailable\n", stderr);
            abort();
        }
    }
    void *layout = create(cr);
    round_positions(context(layout), 0);
    return layout;
}
