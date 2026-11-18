/* Linux-only experiment: one fully initialized parent, forked mapping workers.
 * No exec in workers: immutable index pages are inherited through fork/COW.
 * ponytail: controlled local stdin protocol, at most four jobs per cohort;
 * this is a benchmark driver, not a general-purpose mapping service. */
#include <sys/wait.h>
#include <unistd.h>
#include <fcntl.h>
#include <signal.h>

static int resident_pilot(const mm_idx_t *mi, const mm_mapopt_t *opt,
                          int threads, const char *rg, int argc, char **argv)
{
    char line[16384];
    int n, i, status, failed;
    pid_t children[4];
    char *queries[4], *outputs[4];
    fprintf(stdout, "READY %ld\n", (long)getpid());
    if (fflush(stdout)) return 1;
    while (fgets(line, sizeof(line), stdin)) {
        char extra;
        if (sscanf(line, "%d %c", &n, &extra) != 1 || n < 1 || n > 4) return 1;
        memset(queries, 0, sizeof(queries));
        memset(outputs, 0, sizeof(outputs));
        memset(children, 0, sizeof(children));
        failed = 0;
        for (i = 0; i < n; ++i) {
            char *tab, *end;
            if (!fgets(line, sizeof(line), stdin) || !(end = strchr(line, '\n')) ||
                !(tab = strchr(line, '\t')) || tab > end) { failed = 1; break; }
            *end = 0; *tab++ = 0;
            if (!*line || !*tab || strchr(tab, '\t')) { failed = 1; break; }
            queries[i] = strdup(line); outputs[i] = strdup(tab);
            if (!queries[i] || !outputs[i]) { failed = 1; break; }
        }
        if (!failed && fflush(NULL)) failed = 1;
        for (i = 0; !failed && i < n; ++i) {
            children[i] = fork();
            if (children[i] < 0) { failed = 1; break; }
            if (children[i] == 0) {
                int result, fd = open(outputs[i], O_WRONLY|O_CREAT|O_EXCL, 0600);
                if (fd < 0 || dup2(fd, STDOUT_FILENO) < 0) _exit(2);
                close(fd);
                if (mm_write_sam_hdr(mi, rg, MM_VERSION, argc, argv)) _exit(3);
                result = mm_map_file(mi, queries[i], opt, threads);
                if (fflush(stdout)) result = -1;
                /* Never destroy the inherited index: that would dirty COW pages. */
                _exit(result < 0 ? 4 : 0);
            }
        }
        if (failed)
            for (i = 0; i < n; ++i) if (children[i] > 0) kill(children[i], SIGTERM);
        for (i = 0; i < n; ++i) if (children[i] > 0) {
            struct rusage usage;
            pid_t waited;
            do { waited = wait4(children[i], &status, 0, &usage); } while (waited < 0 && errno == EINTR);
            if (waited < 0 || !WIFEXITED(status) || WEXITSTATUS(status)) failed = 1;
            if (waited > 0)
                fprintf(stdout, "CHILD %ld %d %.6f %.6f %ld\n", (long)children[i],
                    WIFEXITED(status) ? WEXITSTATUS(status) : 128 + WTERMSIG(status),
                    usage.ru_utime.tv_sec + usage.ru_utime.tv_usec / 1e6,
                    usage.ru_stime.tv_sec + usage.ru_stime.tv_usec / 1e6, usage.ru_maxrss);
        }
        for (i = 0; i < n; ++i) { free(queries[i]); free(outputs[i]); }
        fprintf(stdout, "DONE %d\n", failed);
        if (fflush(stdout) || failed) return 1;
    }
    return ferror(stdin) ? 1 : 0;
}
