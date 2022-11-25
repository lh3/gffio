#ifndef GIO_PRIV_H
#define GIO_PRIV_H

#include <stddef.h>
#include "gffio.h"

#define GIO_MALLOC(ptr, len) ((ptr) = (__typeof__(ptr))malloc((len) * sizeof(*(ptr))))
#define GIO_CALLOC(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))
#define GIO_REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))
#define GIO_BZERO(ptr, len) memset((ptr), 0, (len) * sizeof(*(ptr)))
#define GIO_EXPAND(a, m) do { \
		(m) += ((m)>>1) + 16; \
		GIO_REALLOC((a), (m)); \
	} while (0)

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif

gio_dict_t *gio_dict_init(void);
void gio_dict_destroy(gio_dict_t *d);
const char *gio_dict_put(gio_dict_t *d, const char *s);
const char *gio_dict_get(const gio_dict_t *d, const char *s);

void gio_attr_append(gio_gff_t *g, gio_feat_t *f, const char *key, const char *val);
const char *gio_attr_find(const gio_gff_t *g, const gio_feat_t *f, const char *key);
void gio_label(gio_gff_t *gff);

#ifdef __cplusplus
}
#endif

#endif
