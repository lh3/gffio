#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include "gfpriv.h"

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

void gf_sprintf_lite(kstring_t *s, const char *fmt, ...)
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

void gf_write_gff_comment(FILE *fp, const gf_gff_t *gff, kstring_t *str)
{
	int32_t i;
	for (i = 0; i < gff->n_comm; ++i) {
		str->l = 0;
		gf_sprintf_lite(str, "%s\n", gff->comm[i].line);
		fwrite(str->s, 1, str->l, fp);
	}
}

void write_feat(kstring_t *str, const gf_gff_t *gff, const gf_feat_t *f, int32_t fmt)
{
	int32_t i;
	gf_sprintf_lite(str, "%s\t%s\t", f->ctg, f->src);
	if (f->feat == GF_FEAT_MRNA) {
		if (fmt == GF_FMT_GFF3) gf_sprintf_lite(str, "mRNA");
		else gf_sprintf_lite(str, "transcript");
	} else gf_sprintf_lite(str, "%s", f->feat_ori);
	gf_sprintf_lite(str, "\t%ld\t%ld\t", f->st + 1, f->en);
	if (!isnan(f->score)) {
		char buf[32];
		snprintf(buf, 32, "%g", f->score);
		gf_sprintf_lite(str, "%s", buf);
	} else gf_sprintf_lite(str, ".");
	gf_sprintf_lite(str, "\t%c", f->strand < 0? '-' : f->strand > 0? '+' : '.');
	if (f->frame < 0) gf_sprintf_lite(str, "\t.\t");
	else gf_sprintf_lite(str, "\t%d\t", f->frame);
	if (fmt == GF_FMT_GFF3) {
		int32_t l0 = str->l;
		if (!f->has_id && f->id) gf_sprintf_lite(str, "ID=%s;", f->id);
		if (!f->has_parent && f->parent1) gf_sprintf_lite(str, "Parent=%s;", f->parent1);
		if (!f->has_name && f->name) gf_sprintf_lite(str, "Name=%s;", f->name);
		if (str->l == l0 && f->n_attr == 0) gf_sprintf_lite(str, ".");
		for (i = 0; i < f->n_attr; ++i) {
			gf_sprintf_lite(str, "%s=%s", f->attr[i].key, f->attr[i].val);
			if (i != f->n_attr - 1) gf_sprintf_lite(str, ";");
		}
	} else if (fmt == GF_FMT_GTF) {
		if (f->n_attr == 0) gf_sprintf_lite(str, ".");
		for (i = 0; i < f->n_attr; ++i) {
			gf_sprintf_lite(str, "%s \"%s\";", f->attr[i].key, f->attr[i].val);
			if (i != f->n_attr - 1) gf_sprintf_lite(str, " ");
		}
	}
	gf_sprintf_lite(str, "\n");
}

void gf_write_feat(char **str, int32_t *len, int32_t *cap, const gf_gff_t *gff, const gf_feat_t *f, int32_t fmt)
{
	kstring_t s;
	s.s = *str, s.l = *len, s.m = *cap;
	write_feat(&s, gff, f, fmt);
	*str = s.s, *len = s.l, *cap = s.m;
}

void gf_write_gff_stream(FILE *fp, const gf_gff_t *gff, int32_t fmt)
{
	int64_t n_view = gff->n_feat_view > 0? gff->n_feat_view : gff->n_feat;
	int32_t i;
	kstring_t str = {0,0,0};
	gf_write_gff_comment(fp, gff, &str);
	for (i = 0; i < n_view; ++i) {
		str.l = 0;
		write_feat(&str, gff, gff->feat_view? gff->feat_view[i] : &gff->feat[i], fmt);
		fwrite(str.s, 1, str.l, fp);
	}
	free(str.s);
}

void gf_write_bed12_stream(FILE *fp, const gf_gff_t *gff, int32_t fmt)
{
	int32_t i, j;
	kstring_t str = {0,0,0};
	gf_qbuf_t *b;
	gf_mrna_t t;
	b = gf_qbuf_init(gff);
	gf_mrna_init(&t);
	for (i = 0; i < gff->n_feat; ++i) {
		int32_t ret;
		const gf_feat_t *f = &gff->feat[i];
		const char *name = 0;
		if (f->feat != GF_FEAT_MRNA) continue;
		ret = gf_mrna_gen(b, gff, f, &t);
		if (ret < 0) continue; // error
		str.l = 0;
		name = fmt == GF_FMT_BED12L? t.name : f->id? f->id : "NA";
		gf_sprintf_lite(&str, "%s\t%ld\t%ld\t%s", f->ctg, t.st, t.en, name);
		gf_sprintf_lite(&str, "\t0\t%c\t", t.strand < 0? '-' : t.strand > 0? '+' : '.');
		gf_sprintf_lite(&str, "%ld\t%ld\t0,0,0\t%d\t", t.st_cds, t.en_cds, t.n_exon);
		for (j = 0; j < t.n_exon; ++j)
			gf_sprintf_lite(&str, "%ld,", t.exon[j].st);
		gf_sprintf_lite(&str, "\t");
		for (j = 0; j < t.n_exon; ++j)
			gf_sprintf_lite(&str, "%ld,", t.exon[j].en - t.exon[j].st);
		gf_sprintf_lite(&str, "\n");
		fwrite(str.s, 1, str.l, fp);
	}
	gf_mrna_free(&t);
	free(str.s);
	gf_qbuf_destroy(b);
}

void gf_write_bed6_stream(FILE *fp, const gf_gff_t *gff, int32_t fmt)
{
	int32_t i, j;
	kstring_t str = {0,0,0};
	gf_qbuf_t *b;
	gf_mrna_t t;
	b = gf_qbuf_init(gff);
	gf_mrna_init(&t);
	for (i = 0; i < gff->n_feat; ++i) {
		int32_t ret, strand;
		const gf_feat_t *f = &gff->feat[i];
		if (f->feat != GF_FEAT_MRNA) continue;
		ret = gf_mrna_gen(b, gff, f, &t);
		if (ret < 0) continue; // error
		str.l = 0;
		strand = t.strand < 0? '-' : t.strand > 0? '+' : '.';
		if (fmt == GF_FMT_BED_EXON) {
			for (j = 0; j < t.n_exon; ++j)
				gf_sprintf_lite(&str, "%s\t%ld\t%ld\t%s\t0\t%c\n", f->ctg, t.exon[j].st, t.exon[j].en, t.name, strand);
		} else if (fmt == GF_FMT_BED_INTRON) {
			for (j = 0; j < t.n_exon - 1; ++j)
				gf_sprintf_lite(&str, "%s\t%ld\t%ld\t%s\t0\t%c\n", f->ctg, t.exon[j].en, t.exon[j+1].st, t.name, strand);
		} else if (fmt == GF_FMT_BED_CDS && t.has_cds) {
			for (j = 0; j < t.n_exon; ++j) {
				int64_t st, en;
				st = t.exon[j].st > t.st_cds? t.exon[j].st : t.st_cds;
				en = t.exon[j].en < t.en_cds? t.exon[j].en : t.en_cds;
				if (st >= en) continue;
				gf_sprintf_lite(&str, "%s\t%ld\t%ld\t%s\t0\t%c\n", f->ctg, st, en, t.name, strand);
			}
		}
		if (str.l > 0)
			fwrite(str.s, 1, str.l, fp);
	}
	gf_mrna_free(&t);
	free(str.s);
	gf_qbuf_destroy(b);
}

void gf_write(const char *fn, const gf_gff_t *gff, int32_t fmt)
{
	FILE *fp;
	fp = fn && strcmp(fn, "-")? fopen(fn, "w") : stdout;
	if (fmt == GF_FMT_GFF3 || fmt == GF_FMT_GTF)
		gf_write_gff_stream(fp, gff, fmt);
	else if (fmt == GF_FMT_BED12S || fmt == GF_FMT_BED12L)
		gf_write_bed12_stream(fp, gff, fmt);
	else if (fmt == GF_FMT_BED_EXON || fmt == GF_FMT_BED_INTRON || fmt == GF_FMT_BED_CDS)
		gf_write_bed6_stream(fp, gff, fmt);
	fclose(fp);
}

void gf_write_list(const char *fn, const gf_gff_t *gff, int32_t fmt, int32_t n, char **list)
{
	int32_t i, len = 0, cap = 0;
	char *str = 0;
	gf_qbuf_t *b;
	FILE *fp;
	if (n <= 0) return;
	fp = fn && strcmp(fn, "-")? fopen(fn, "w") : stdout;
	assert(fp);
	b = gf_qbuf_init(gff);
	for (i = 0; i < n; ++i) {
		const gf_feat_t *f;
		f = gf_get_by_id(gff, list[i]);
		if (f == 0) {
			if (gf_verbose >= 2)
				fprintf(stderr, "WARNING: failed to find ID '%s'\n", list[i]);
		} else {
			const gf_feat_t **fs;
			int32_t j, n_fs;
			fs = gf_descend(b, f, &n_fs);
			for (j = 0; j < n_fs; ++j) {
				len = 0;
				gf_write_feat(&str, &len, &cap, gff, fs[j], GF_FMT_GFF3);
				fwrite(str, 1, len, fp);
			}
		}
	}
	gf_qbuf_destroy(b);
	fclose(fp);
}

void gf_write_fasta_stream(FILE *fp, const gf_gff_t *gff, const gf_seqs_t *seq, int32_t fmt)
{
	int32_t i, cap = 0;
	char *s = 0;
	kstring_t str = {0,0,0};
	gf_qbuf_t *b;
	gf_mrna_t t;
	b = gf_qbuf_init(gff);
	gf_mrna_init(&t);
	for (i = 0; i < gff->n_feat; ++i) {
		int32_t ret;
		const gf_feat_t *f = &gff->feat[i];
		if (f->feat != GF_FEAT_MRNA) continue;
		ret = gf_mrna_gen(b, gff, f, &t);
		if (ret < 0) continue; // error
		ret = gf_extract_seq(gff, seq, &t, fmt, &s, &cap);
		if (ret < 0) continue; // error
		str.l = 0;
		gf_sprintf_lite(&str, ">%s\n%s\n", t.name, s);
		fwrite(str.s, 1, str.l, fp);
	}
	free(s); free(str.s);
	gf_mrna_free(&t);
	gf_qbuf_destroy(b);
}

void gf_write_fasta(const char *fn, const gf_gff_t *gff, const gf_seqs_t *seq, int32_t fmt)
{
	FILE *fp;
	fp = fn && strcmp(fn, "-")? fopen(fn, "w") : stdout;
	if (fmt == GF_FMT_FA_MRNA || fmt == GF_FMT_FA_CDS || fmt == GF_FMT_FA_PROTEIN)
		gf_write_fasta_stream(fp, gff, seq, fmt);
	fclose(fp);
}
