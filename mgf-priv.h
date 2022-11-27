#ifndef MGF_PRIV_H
#define MGF_PRIV_H

#include <stdarg.h>
#include <stdlib.h>
#include "minigff.h"

#define MGF_MALLOC(ptr, len) ((ptr) = (__typeof__(ptr))malloc((len) * sizeof(*(ptr))))
#define MGF_CALLOC(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))
#define MGF_REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))
#define MGF_BZERO(ptr, len) memset((ptr), 0, (len) * sizeof(*(ptr)))
#define MGF_EXPAND(a, m) do { \
		(m) += ((m)>>1) + 16; \
		MGF_REALLOC((a), (m)); \
	} while (0)
#define MGF_PUSH_BACK(n_, m_, a_, v_) do { \
		if ((n_) == (m_)) { \
			(m_) += ((m_)>>1) + 16; \
			MGF_REALLOC((a_), (m_)); \
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

mgf_dict_t *mgf_dict_init(void);
void mgf_dict_destroy(mgf_dict_t *d);
const char *mgf_dict_put(mgf_dict_t *d, const char *s);
const char *mgf_dict_get(const mgf_dict_t *d, const char *s);
char *mgf_strndup(const char *src, size_t n);

void *mgf_id_init(void);
void mgf_id_destroy(void *d);
int32_t mgf_id_put(void *d, const char *id, int32_t idx);
int32_t mgf_id_get(void *d, const char *id);
int32_t mgf_id_size(void *d);

void mgf_attr_append(mgf_gff_t *g, mgf_feat_t *f, const char *key, const char *val);
const char *mgf_attr_find(const mgf_gff_t *g, const mgf_feat_t *f, const char *key);

void mgf_label(mgf_gff_t *gff);
void mgf_connect(mgf_gff_t *gff);

const mgf_feat_t **mgf_toposort(const mgf_gff_t *gff);;

void mgf_mrna_init(mgf_mrna_t *t);
int32_t mgf_mrna_gen(mgf_qbuf_t *b, const mgf_gff_t *gff, const mgf_feat_t *f, mgf_mrna_t *t);

void mgf_sprintf_lite(kstring_t *s, const char *fmt, ...);

#ifdef __cplusplus
}
#endif

#endif
