#define _GNU_SOURCE

#include <dlfcn.h>
#include <errno.h>
#include <limits.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

typedef size_t (*fwrite_fn)(const void *, size_t, size_t, FILE *);
typedef int (*fseek_fn)(FILE *, long, int);
typedef FILE *(*fopen_fn)(const char *, const char *);

static int injected = 0;

static int IsRestartPartial(FILE *stream) {
  char fd_path[64];
  char path[PATH_MAX + 1];
  const int fd = fileno(stream);
  if (fd < 0) return 0;
  snprintf(fd_path, sizeof(fd_path), "/proc/self/fd/%d", fd);
  const ssize_t length = readlink(fd_path, path, PATH_MAX);
  if (length < 0) return 0;
  path[length] = '\0';
  return strstr(path, ".rst.partial") != NULL;
}

static fwrite_fn RealFwrite(void) {
  static fwrite_fn real_fwrite = NULL;
  if (real_fwrite == NULL) {
    real_fwrite = (fwrite_fn)dlsym(RTLD_NEXT, "fwrite");
  }
  return real_fwrite;
}

static fseek_fn RealFseek(void) {
  static fseek_fn real_fseek = NULL;
  if (real_fseek == NULL) {
    real_fseek = (fseek_fn)dlsym(RTLD_NEXT, "fseek");
  }
  return real_fseek;
}

static fopen_fn RealFopen(void) {
  static fopen_fn real_fopen = NULL;
  if (real_fopen == NULL) {
    real_fopen = (fopen_fn)dlsym(RTLD_NEXT, "fopen");
  }
  return real_fopen;
}

FILE *fopen(const char *path, const char *mode) {
  const char *fault = getenv("ATHENAK_RESTART_FAULT");
  if (fault != NULL && strcmp(fault, "rank_one_fopen_failure") == 0 &&
      strstr(path, "rank_00000001") != NULL &&
      strstr(path, ".rst.partial") != NULL) {
    errno = EACCES;
    return NULL;
  }
  return RealFopen()(path, mode);
}

size_t fwrite(const void *ptr, size_t size, size_t count, FILE *stream) {
  fwrite_fn real_fwrite = RealFwrite();
  const char *fault = getenv("ATHENAK_RESTART_FAULT");
  if (!injected && fault != NULL && IsRestartPartial(stream)) {
    if (strcmp(fault, "short_header_write") == 0 && count > 1) {
      injected = 1;
      return real_fwrite(ptr, size, count - 1, stream);
    }
    if (strcmp(fault, "kill_writer") == 0) {
      injected = 1;
      const size_t written = real_fwrite(ptr, size, count, stream);
      fflush(stream);
      fsync(fileno(stream));
      kill(getpid(), SIGKILL);
      _exit(137);
      return written;
    }
  }
  return real_fwrite(ptr, size, count, stream);
}

int fseek(FILE *stream, long offset, int whence) {
  const char *fault = getenv("ATHENAK_RESTART_FAULT");
  if (!injected && fault != NULL && strcmp(fault, "fseek_failure") == 0 &&
      IsRestartPartial(stream)) {
    injected = 1;
    errno = EIO;
    return -1;
  }
  return RealFseek()(stream, offset, whence);
}
