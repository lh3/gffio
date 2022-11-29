#ifndef GF_PRIV_H
#define GF_PRIV_H

#include <stdarg.h>
#include <stdlib.h>
#include "gffio.h"

#define GF_MALLOC(ptr, len) ((ptr) = (__typeof__(ptr))malloc((len) * sizeof(*(ptr))))
#define GF_CALLOC(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))
#define GF_REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))
#define GF_BZERO(ptr, len) memset((ptr), 0, (len) * sizeof(*(ptr)))
#define GF_EXPAND(a, m) do { \
		(m) += ((m)>>1) + 16; \
		GF_REALLOC((a), (m)); \
	} while (0)
#define GF_PUSH_BACK(n_, m_, a_, v_) do { \
		if ((n_) == (m_)) { \
			(m_) += ((m_)>>1) + 16; \
			GF_REALLOC((a_), (m_)); \
		} \
		(a_)[(n_)++] = (v_); \
	} while (0)

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

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

gf_dict_t *gf_dict_init(void);
void gf_dict_destroy(gf_dict_t *d);
const char *gf_dict_put(gf_dict_t *d, const char *s);
const char *gf_dict_get(const gf_dict_t *d, const char *s);
char *gf_strndup(const char *src, size_t n);

void *gf_id_init(void);
void gf_id_destroy(void *d);
int32_t gf_id_put(void *d, const char *id, int32_t idx);
int32_t gf_id_get(void *d, const char *id);
int32_t gf_id_size(void *d);

void gf_attr_append(gf_gff_t *g, gf_feat_t *f, const char *key, const char *val);
const char *gf_attr_find(const gf_gff_t *g, const gf_feat_t *f, const char *key);

void gf_label(gf_gff_t *gff);
void gf_connect(gf_gff_t *gff);

const gf_feat_t **gf_toposort(const gf_gff_t *gff);;

void gf_mrna_init(gf_mrna_t *t);
int32_t gf_mrna_gen(gf_qbuf_t *b, const gf_gff_t *gff, const gf_feat_t *f, gf_mrna_t *t);
void gf_mrna_free(gf_mrna_t *t);

void gf_sprintf_lite(kstring_t *s, const char *fmt, ...);

#ifdef __cplusplus
}
#endif

#endif
