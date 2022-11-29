#include "gfpriv.h"
#include "khashl.h"
KHASHL_CMAP_INIT(KH_LOCAL, gf_strhash_t, gf_sh, const char*, int32_t, kh_hash_str, kh_eq_str)

struct gf_dict_s {
	int32_t n_str, m_str;
	char **str;
	gf_strhash_t *h;
};

gf_dict_t *gf_dict_init(void)
{
	gf_dict_t *d;
	GF_CALLOC(d, 1);
	d->h = gf_sh_init();
	return d;
}

char *gf_strdup(const char *src)
{
	int32_t len;
	char *dst;
	len = strlen(src);
	dst = (char*)malloc(len + 1);
	memcpy(dst, src, len + 1);
	return dst;
}

char *gf_strndup(const char *src, size_t n)
{
	char *dst;
	dst = (char*)malloc(n + 1);
	strncpy(dst, src, n);
	dst[n] = 0;
	return dst;
}

const char *gf_dict_put(gf_dict_t *d, const char *s)
{
	int absent;
	khint_t k;
	k = gf_sh_put(d->h, s, &absent);
	if (absent) {
		if (d->n_str == d->m_str)
			GF_EXPAND(d->str, d->m_str);
		kh_key(d->h, k) = d->str[d->n_str] = gf_strdup(s);
		kh_val(d->h, k) = d->n_str++;
	}
	return kh_key(d->h, k);
}

const char *gf_dict_get(const gf_dict_t *d, const char *s)
{
	khint_t k;
	k = gf_sh_get(d->h, s);
	return k == kh_end(d->h)? 0 : kh_key(d->h, k);
}

void gf_dict_destroy(gf_dict_t *d)
{
	int32_t i;
	if (d == 0) return;
	gf_sh_destroy(d->h);
	for (i = 0; i < d->n_str; ++i)
		free(d->str[i]);
	free(d->str);
	free(d);
}

void *gf_id_init(void)
{
	return (void*)gf_sh_init();
}

void gf_id_destroy(void *d)
{
	gf_sh_destroy((gf_strhash_t*)d);
}

int32_t gf_id_put(void *d, const char *id, int32_t idx)
{
	int absent;
	khint_t k;
	gf_strhash_t *h = (gf_strhash_t*)d;
	k = gf_sh_put(h, id, &absent);
	if (absent) kh_val(h, k) = idx;
	return absent;
}

int32_t gf_id_get(void *d, const char *id)
{
	khint_t k;
	gf_strhash_t *h = (gf_strhash_t*)d;
	k = gf_sh_get(h, id);
	return k == kh_end(h)? -1 : kh_val(h, k);
}

int32_t gf_id_size(void *d)
{
	gf_strhash_t *h = (gf_strhash_t*)d;
	return kh_size(h);
}
