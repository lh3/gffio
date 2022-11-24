#include "gio-priv.h"
#include "khashl.h"

KHASHL_CMAP_INIT(KH_LOCAL, gio_strhash_t, gio_sh, const char*, int32_t, kh_hash_str, kh_eq_str)

struct gio_dict_s {
	int32_t n_str, m_str;
	char **str;
	gio_strhash_t *h;
};

gio_dict_t *gio_dict_init(void)
{
	gio_dict_t *d;
	GIO_CALLOC(d, 1);
	d->h = gio_sh_init();
	return d;
}

char *gio_strdup(const char *src)
{
	int32_t len;
	char *dst;
	len = strlen(src);
	dst = (char*)malloc(len + 1);
	memcpy(dst, src, len + 1);
	return dst;
}

int32_t gio_dict_put(gio_dict_t *d, const char *s)
{
	int absent;
	khint_t k;
	k = gio_sh_put(d->h, s, &absent);
	if (absent) {
		if (d->n_str == d->m_str)
			GIO_EXPAND(d->str, d->m_str);
		kh_key(d->h, k) = d->str[d->n_str] = gio_strdup(s);
		kh_val(d->h, k) = d->n_str++;
	}
	return kh_val(d->h, k);
}

int32_t gio_dict_get(const gio_dict_t *d, const char *s)
{
	khint_t k;
	k = gio_sh_get(d->h, s);
	return k == kh_end(d->h)? -1 : kh_val(d->h, k);
}

void gio_dict_destroy(gio_dict_t *d)
{
	int32_t i;
	if (d == 0) return;
	gio_sh_destroy(d->h);
	for (i = 0; i < d->n_str; ++i)
		free(d->str[i]);
	free(d->str);
	free(d);
}
