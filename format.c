#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "mgf-priv.h"

#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

static inline void str_enlarge(kstring_t *s, int l)
{
	if (s->l + l + 1 > s->m) {
		s->m = s->l + l + 1;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
}

static inline void str_copy(kstring_t *s, const char *st, const char *en)
{
	str_enlarge(s, en - st);
	memcpy(&s->s[s->l], st, en - st);
	s->l += en - st;
}

void mgf_sprintf_lite(kstring_t *s, const char *fmt, ...)
{
	char buf[32]; // for integer to string conversion
	const char *p, *q;
	va_list ap;
	va_start(ap, fmt);
	for (q = p = fmt; *p; ++p) {
		if (*p == '%') {
			if (p > q) str_copy(s, q, p);
			++p;
			if (*p == 'd') {
				int c, i, l = 0;
				unsigned int x;
				c = va_arg(ap, int);
				x = c >= 0? c : -c;
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				if (c < 0) buf[l++] = '-';
				str_enlarge(s, l);
				for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
			} else if (*p == 'l' && *(p+1) == 'd') {
				int c, i, l = 0;
				unsigned long x;
				c = va_arg(ap, long);
				x = c >= 0? c : -c;
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				if (c < 0) buf[l++] = '-';
				str_enlarge(s, l);
				for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
				++p;
			} else if (*p == 'u') {
				int i, l = 0;
				uint32_t x;
				x = va_arg(ap, uint32_t);
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				str_enlarge(s, l);
				for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
			} else if (*p == 's') {
				char *r = va_arg(ap, char*);
				str_copy(s, r, r + strlen(r));
			} else if (*p == 'c') {
				str_enlarge(s, 1);
				s->s[s->l++] = va_arg(ap, int);
			} else {
				fprintf(stderr, "ERROR: unrecognized type '%%%c'\n", *p);
				abort();
			}
			q = p + 1;
		}
	}
	if (p > q) str_copy(s, q, p);
	va_end(ap);
	s->s[s->l] = 0;
}

static void write_comment(FILE *fp, const mgf_gff_t *gff, kstring_t *str)
{
	int32_t i;
	for (i = 0; i < gff->n_comm; ++i) {
		str->l = 0;
		mgf_sprintf_lite(str, "%s\n", gff->comm[i].line);
		fwrite(str->s, 1, str->l, fp);
	}
}

static void write_feat(kstring_t *str, const mgf_gff_t *gff, const mgf_feat_t *f, int32_t fmt)
{
	int32_t i;
	str->l = 0;
	mgf_sprintf_lite(str, "%s\t%s\t%s\t%ld\t%ld\t", f->ctg, f->src, f->feat_ori, f->st + 1, f->en);
	if (!isnan(f->score)) {
		char buf[32];
		snprintf(buf, 32, "%g", f->score);
		mgf_sprintf_lite(str, "%s", buf);
	} else mgf_sprintf_lite(str, ".");
	mgf_sprintf_lite(str, "\t%c", f->strand < 0? '-' : f->strand > 0? '+' : '.');
	if (f->frame < 0) mgf_sprintf_lite(str, "\t.\t");
	else mgf_sprintf_lite(str, "\t%d\t", f->frame);
	if (fmt == MGF_FMT_GFF3) {
		if (f->n_attr == 0) mgf_sprintf_lite(str, ".");
		for (i = 0; i < f->n_attr; ++i) {
			mgf_sprintf_lite(str, "%s=%s", f->attr[i].key, f->attr[i].val);
			if (i != f->n_attr - 1) mgf_sprintf_lite(str, ";");
		}
	} else if (fmt == MGF_FMT_GTF) {
	}
	mgf_sprintf_lite(str, "\n");
}

void mgf_write_gff_stream(FILE *fp, const mgf_gff_t *gff, int32_t fmt)
{
	int32_t i;
	kstring_t str = {0,0,0};
	write_comment(fp, gff, &str);
	for (i = 0; i < gff->n_feat; ++i) {
		write_feat(&str, gff, gff->feat_view? gff->feat_view[i] : &gff->feat[i], fmt);
		fwrite(str.s, 1, str.l, fp);
	}
	free(str.s);
}

void mgf_write(const char *fn, const mgf_gff_t *gff, int32_t fmt)
{
	FILE *fp;
	fp = fn && strcmp(fn, "-")? fopen(fn, "w") : fdopen(1, "w");
	mgf_write_gff_stream(fp, gff, fmt);
	fclose(fp);
}
