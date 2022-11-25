#include "mgf-priv.h"
#include "khashl.h"

KHASHL_CMAP_INIT(KH_LOCAL, mgf_strhash_t, mgf_sh, const char*, int32_t, kh_hash_str, kh_eq_str)

struct mgf_dict_s {
	int32_t n_str, m_str;
	char **str;
	mgf_strhash_t *h;
};

mgf_dict_t *mgf_dict_init(void)
{
	mgf_dict_t *d;
	MGF_CALLOC(d, 1);
	d->h = mgf_sh_init();
	return d;
}

char *mgf_strdup(const char *src)
{
	int32_t len;
	char *dst;
	len = strlen(src);
	dst = (char*)malloc(len + 1);
	memcpy(dst, src, len + 1);
	return dst;
}

const char *mgf_dict_put(mgf_dict_t *d, const char *s)
{
	int absent;
	khint_t k;
	k = mgf_sh_put(d->h, s, &absent);
	if (absent) {
		if (d->n_str == d->m_str)
			MGF_EXPAND(d->str, d->m_str);
		kh_key(d->h, k) = d->str[d->n_str] = mgf_strdup(s);
		kh_val(d->h, k) = d->n_str++;
	}
	return kh_key(d->h, k);
}

const char *mgf_dict_get(const mgf_dict_t *d, const char *s)
{
	khint_t k;
	k = mgf_sh_get(d->h, s);
	return k == kh_end(d->h)? 0 : kh_key(d->h, k);
}

void mgf_dict_destroy(mgf_dict_t *d)
{
	int32_t i;
	if (d == 0) return;
	mgf_sh_destroy(d->h);
	for (i = 0; i < d->n_str; ++i)
		free(d->str[i]);
	free(d->str);
	free(d);
}
