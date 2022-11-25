#include <string.h>
#include <stdio.h>
#include <math.h>
#include <zlib.h>
#include "mgf-priv.h"
#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 0x10000)

static void mgf_parse_attr(mgf_gff_t *g, mgf_feat_t *f, char *str)
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
		mgf_attr_append(g, f, key, val);
	}
}

static void mgf_parse_feat(mgf_gff_t *gff, char *str)
{
	int32_t i;
	char *p, *q;
	mgf_feat_t *f;
	if (gff->n_feat == gff->m_feat) {
		int32_t oldm = gff->m_feat;
		MGF_EXPAND(gff->feat, gff->m_feat);
		memset(&gff->feat[oldm], 0, sizeof(mgf_feat_t) * (gff->m_feat - oldm));
	}
	f = &gff->feat[gff->n_feat++];
	for (p = q = str, i = 0;; ++p) {
		if (*p == '\t' || *p == 0) {
			int32_t c = *p;
			*p = 0;
			if (i == 0) { // contig name
				f->ctg = mgf_dict_put(gff->dict, q);
			} else if (i == 1) {
				f->src = mgf_dict_put(gff->dict, q);
			} else if (i == 2) { // source
				f->feat_ori = mgf_dict_put(gff->dict, q);
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
				f->frame = *q >= '0' && *q <= '9'? atoi(q) : -1;
			} else if (i == 8) { // attributes
				mgf_parse_attr(gff, f, q);
			}
			q = p + 1, ++i;
			if (c == 0 || i == 9) break;
		}
	}
}

mgf_gff_t *mgf_read_ks(kstream_t *ks)
{
	mgf_gff_t *gff = 0;
	kstring_t str = {0,0,0};
	int dret;
	int64_t lineoff = 0;
	MGF_CALLOC(gff, 1);
	gff->dict = mgf_dict_init();
	while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
		if (str.s[0] == '#') {
			mgf_comm_t *p;
			if (gff->n_comm == gff->m_comm)
				MGF_EXPAND(gff->comm, gff->m_comm);
			p = &gff->comm[gff->n_comm++];
			p->lineoff = lineoff;
			MGF_CALLOC(p->line, str.l + 1);
			memcpy(p->line, str.s, str.l + 1);
		} else {
			mgf_parse_feat(gff, str.s);
		}
		++lineoff;
	}
	mgf_label(gff);
	return gff;
}

mgf_gff_t *mgf_read(const char *fn)
{
	gzFile fp;
	kstream_t *ks;
	mgf_gff_t *gff;
	fp = fn || strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(0, "r");
	if (fp == 0) return 0;
	ks = ks_init(fp);
	gff = mgf_read_ks(ks);
	ks_destroy(ks);
	gzclose(fp);
	return gff;
}

void mgf_destroy(mgf_gff_t *gff)
{
}
