#ifndef E_LIB_H_
#define E_LIB_H_

#include <stddef.h>


void eprintf(const char *fmt, ...);
void error(const char *msg);

char *estrdup(const char *s);
void *emalloc(size_t n);
void *ecalloc(size_t nmemb, size_t size);
void *efopen(const char *fileName, const char *mode);

void setProgName(const char *str);
const char *getProgName(void);
void printSpace(void);

#endif /* end of elib.h */
