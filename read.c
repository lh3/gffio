#include <string.h>
#include <stdio.h>
#include <math.h>
#include <zlib.h>
#include "gfpriv.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

static void gf_parse_attr(gf_gff_t *g, gf_feat_t *f, char *str)
{
	char *q = str;
	int32_t len;
	len = strlen(str);
	while (*q && q - str < len) {
		char *p = q, *key = 0, *val = 0;
		while (*p && (*p == ';' || *p == ' ')) ++p; // skip leading ; or SPACE
		key = p;
		while (*p && *p != '=' && *p != ' ') ++p; // get the key string
		if (*p == 0) break; // no value
		*p++ = 0;
		while (*p && (*p == ';' || *p == ' ')) ++p; // skip ; or SPACE
		if (*p == '\'' || *p == '\"') {
			int32_t c = *p++;
			val = p;
			while (*p && *p != c) ++p; // TODO: support escape
			if (*p == 0) break; // missing the ending quotation mark
		} else {
			val = p;
			while (*p && *p != ';') ++p; // TODO: support escape
		}
		*p = 0;
		q = p + 1;
		gf_attr_append(g, f, key, val);
	}
}

static void gf_parse_feat(gf_gff_t *gff, char *str, int64_t lineoff)
{
	int32_t i, tmp;
	char *p, *q;
	gf_feat_t *f;
	if (gff->n_feat == gff->m_feat) {
		int32_t oldm = gff->m_feat;
		GF_EXPAND(gff->feat, gff->m_feat);
		memset(&gff->feat[oldm], 0, sizeof(gf_feat_t) * (gff->m_feat - oldm));
	}
	f = &gff->feat[gff->n_feat++];
	f->lineoff = lineoff;
	for (p = q = str, i = 0;; ++p) {
		if (*p == '\t' || *p == 0) {
			int32_t c = *p;
			*p = 0;
			if (i == 0) { // contig name
				f->ctg = gf_dict_put(gff->dict, q);
			} else if (i == 1) {
				f->src = gf_dict_put(gff->dict, q);
			} else if (i == 2) { // source
				f->feat_ori = gf_dict_put(gff->dict, q);
			} else if (i == 3) { // start
				f->st = atol(q) - 1;
			} else if (i == 4) { // end
				f->en = atol(q);
			} else if (i == 5) { // score
				char *r;
				f->score = strcmp(q, ".") == 0? nan("1") : strtod(q, &r);
			} else if (i == 6) { // strand
				f->strand = *q == '+'? 1 : *q == '-'? -1 : 0;
			} else if (i == 7) { // frame
				tmp = *q >= '0' && *q <= '9'? atoi(q) : -1;
				f->frame = tmp < 0? -1 : tmp % 3;
			} else if (i == 8) { // attributes
				gf_parse_attr(gff, f, q);
			}
			q = p + 1, ++i;
			if (c == 0 || i == 9) break;
		}
	}
}

gf_gff_t *gf_read_ks(kstream_t *ks)
{
	gf_gff_t *gff = 0;
	kstring_t str = {0,0,0};
	int dret;
	int64_t lineoff = 0;
	GF_CALLOC(gff, 1);
	gff->dict = gf_dict_init();
	while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
		if (str.s[0] == '#') {
			gf_comm_t *p;
			if (gff->n_comm == gff->m_comm)
				GF_EXPAND(gff->comm, gff->m_comm);
			p = &gff->comm[gff->n_comm++];
			p->lineoff = lineoff;
			GF_CALLOC(p->line, str.l + 1);
			memcpy(p->line, str.s, str.l + 1);
		} else {
			gf_parse_feat(gff, str.s, lineoff);
		}
		++lineoff;
	}
	free(str.s);
	gf_label(gff);
	gf_connect(gff);
	return gff;
}

gf_gff_t *gf_read(const char *fn)
{
	gzFile fp;
	kstream_t *ks;
	gf_gff_t *gff;
	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(0, "r");
	if (fp == 0) return 0;
	ks = ks_init(fp);
	gff = gf_read_ks(ks);
	ks_destroy(ks);
	gzclose(fp);
	return gff;
}

void gf_destroy(gf_gff_t *gff)
{
	int32_t i;
	if (gff == 0) return;
	for (i = 0; i < gff->n_feat; ++i) {
		gf_feat_t *f = &gff->feat[i];
		free(f->attr); free(f->parent); // f->child is allocated together with parent
	}
	free(gff->feat); free(gff->feat_view);
	for (i = 0; i < gff->n_comm; ++i)
		free(gff->comm[i].line);
	free(gff->comm);
	gf_id_destroy(gff->dict_id);
	gf_dict_destroy(gff->dict);
	free(gff);
}

gf_seqs_t *gf_seqs_read(const char *fn)
{
	gzFile fp;
	kseq_t *ks;
	gf_seqs_t *s;
	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(0, "r");
	if (fp == 0) return 0;
	ks = kseq_init(fp);
	GF_CALLOC(s, 1);
	s->h = gf_id_init();
	while (kseq_read(ks) >= 0) {
		int32_t absent;
		char *name;
		name = gf_strndup(ks->name.s, ks->name.l);
		absent = gf_id_put(s->h, name, s->n_seq);
		if (!absent) {
			if (gf_verbose >= 2)
				fprintf(stderr, "WARNING: skipped duplicated sequence '%s'\n", name);
			free(name);
			continue;
		}
		if (s->n_seq == s->m_seq) {
			s->m_seq += (s->m_seq>>1) + 16;
			GF_REALLOC(s->name, s->m_seq);
			GF_REALLOC(s->seq, s->m_seq);
			GF_REALLOC(s->len, s->m_seq);
		}
		s->name[s->n_seq] = name;
		s->seq[s->n_seq] = gf_strndup(ks->seq.s, ks->seq.l);
		s->len[s->n_seq++] = ks->seq.l;
	}
	kseq_destroy(ks);
	gzclose(fp);
	return s;
}

void gf_seqs_destroy(gf_seqs_t *s)
{
	int32_t i;
	if (s == 0) return;
	gf_id_destroy(s->h);
	for (i = 0; i < s->n_seq; ++i) {
		free(s->name[i]); free(s->seq[i]);
	}
	free(s->name); free(s->seq); free(s->len);
	free(s);
}

char **gf_read_list(const char *o, int *n_) // from gfatools
{
	int n = 0, m = 0;
	char **s = 0;
	*n_ = 0;
	if (*o != '@') {
		const char *q = o, *p;
		for (p = q;; ++p) {
			if (*p == ',' || *p == 0) {
				if (n == m) {
					m = m? m<<1 : 16;
					s = (char**)realloc(s, m * sizeof(char*));
				}
				s[n++] = gf_strndup(q, p - q);
				if (*p == 0) break;
				q = p + 1;
			}
		}
	} else {
		gzFile fp;
		kstream_t *ks;
		kstring_t str = {0,0,0};
		int dret;

		fp = gzopen(o + 1, "r");
		if (fp == 0) return 0;
		ks = ks_init(fp);
		while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
			char *p;
			for (p = str.s; *p && !isspace(*p); ++p);
			if (n == m) {
				m = m? m<<1 : 16;
				s = (char**)realloc(s, m * sizeof(char*));
			}
			s[n++] = gf_strndup(str.s, p - str.s);
		}
		ks_destroy(ks);
		gzclose(fp);
	}
	if (s) s = (char**)realloc(s, n * sizeof(char*));
	*n_ = n;
	return s;
}
