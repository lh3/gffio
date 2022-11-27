#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include "mgf-priv.h"

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

void mgf_write_gff_comment(FILE *fp, const mgf_gff_t *gff, kstring_t *str)
{
	int32_t i;
	for (i = 0; i < gff->n_comm; ++i) {
		str->l = 0;
		mgf_sprintf_lite(str, "%s\n", gff->comm[i].line);
		fwrite(str->s, 1, str->l, fp);
	}
}

void write_feat(kstring_t *str, const mgf_gff_t *gff, const mgf_feat_t *f, int32_t fmt)
{
	int32_t i;
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

void mgf_write_feat(char **str, int32_t *len, int32_t *cap, const mgf_gff_t *gff, const mgf_feat_t *f, int32_t fmt)
{
	kstring_t s;
	s.s = *str, s.l = *len, s.m = *cap;
	write_feat(&s, gff, f, fmt);
	*str = s.s, *len = s.l, *cap = s.m;
}

void mgf_write_gff_stream(FILE *fp, const mgf_gff_t *gff, int32_t fmt)
{
	int32_t i;
	kstring_t str = {0,0,0};
	mgf_write_gff_comment(fp, gff, &str);
	for (i = 0; i < gff->n_feat; ++i) {
		str.l = 0;
		write_feat(&str, gff, gff->feat_view? gff->feat_view[i] : &gff->feat[i], fmt);
		fwrite(str.s, 1, str.l, fp);
	}
	free(str.s);
}

void mgf_write_bed12_stream(FILE *fp, const mgf_gff_t *gff, int32_t fmt)
{
	int32_t i, j;
	kstring_t str = {0,0,0};
	mgf_qbuf_t *b;
	mgf_mrna_t t;
	b = mgf_qbuf_init(gff);
	mgf_mrna_init(&t);
	for (i = 0; i < gff->n_feat; ++i) {
		int32_t ret;
		const mgf_feat_t *f = &gff->feat[i];
		const char *name = 0;
		if (f->feat != MGF_FEAT_MRNA) continue;
		ret = mgf_mrna_gen(b, gff, f, &t);
		if (ret < 0) continue; // error
		str.l = 0;
		name = fmt == MGF_FMT_BED12L? t.name : f->id? f->id : "NA";
		mgf_sprintf_lite(&str, "%s\t%ld\t%ld\t%s", f->ctg, t.st, t.en, name);
		mgf_sprintf_lite(&str, "\t0\t%c\t", t.strand < 0? '-' : t.strand > 0? '+' : '.');
		mgf_sprintf_lite(&str, "%ld\t%ld\t0,0,0\t%d\t", t.st_cds, t.en_cds, t.n_exon);
		for (j = 0; j < t.n_exon; ++j)
			mgf_sprintf_lite(&str, "%ld,", t.exon[j].st);
		mgf_sprintf_lite(&str, "\t");
		for (j = 0; j < t.n_exon; ++j)
			mgf_sprintf_lite(&str, "%ld,", t.exon[j].en - t.exon[j].st);
		mgf_sprintf_lite(&str, "\n");
		fwrite(str.s, 1, str.l, fp);
	}
	free(str.s);
	mgf_qbuf_destroy(b);
}

void mgf_write(const char *fn, const mgf_gff_t *gff, int32_t fmt)
{
	FILE *fp;
	fp = fn && strcmp(fn, "-")? fopen(fn, "w") : fdopen(1, "w");
	if (fmt == MGF_FMT_GFF3 || fmt == MGF_FMT_GTF)
		mgf_write_gff_stream(fp, gff, fmt);
	else if (fmt == MGF_FMT_BED12 || fmt == MGF_FMT_BED12L)
		mgf_write_bed12_stream(fp, gff, fmt);
	fclose(fp);
}

void mgf_write_list(const char *fn, const mgf_gff_t *gff, int32_t fmt, int32_t n, char **list)
{
	int32_t i, len = 0, cap = 0;
	char *str = 0;
	mgf_qbuf_t *b;
	FILE *fp;
	if (n <= 0) return;
	fp = fn && strcmp(fn, "-")? fopen(fn, "w") : fdopen(1, "w");
	assert(fp);
	b = mgf_qbuf_init(gff);
	for (i = 0; i < n; ++i) {
		const mgf_feat_t *f;
		f = mgf_get_by_id(gff, list[i]);
		if (f == 0) {
			if (mgf_verbose >= 2)
				fprintf(stderr, "WARNING: failed to find ID '%s'\n", list[i]);
		} else {
			const mgf_feat_t **fs;
			int32_t j, n_fs;
			fs = mgf_descend(b, f, &n_fs);
			for (j = 0; j < n_fs; ++j) {
				len = 0;
				mgf_write_feat(&str, &len, &cap, gff, fs[j], MGF_FMT_GFF3);
				fwrite(str, 1, len, fp);
			}
		}
	}
	mgf_qbuf_destroy(b);
	fclose(fp);
}
