#include "mgf-priv.h"

const mgf_feat_t **mgf_get_by_id(const mgf_gff_t *gff, const char *id, const char *feat, int32_t *n_feat_)
{
	int32_t i, k, n_feat;
	const char *s_feat = 0;
	const mgf_feat_t *f;
	const mgf_feat_t **r;
	*n_feat_ = 0;
	k = mgf_id_get(gff->dict_id, id);
	if (k < 0) return 0; // the ID not found
	if (feat) {
		s_feat = mgf_dict_get(gff->dict, feat);
		if (s_feat == 0) return 0; // invalid feat
	}
	f = &gff->feat[k];
	if (f->n_child == 0) return 0; // no child features
	if (s_feat == 0) {
		n_feat = f->n_child;
		MGF_CALLOC(r, n_feat);
		for (i = 0; i < f->n_child; ++i)
			r[i] = f->child[i];
	} else {
		for (i = 0, n_feat = 0; i < f->n_child; ++i)
			if (f->child[i]->feat_ori == s_feat)
				++n_feat;
		MGF_CALLOC(r, n_feat);
		for (i = 0, n_feat = 0; i < f->n_child; ++i)
			if (f->child[i]->feat_ori == s_feat)
				r[n_feat++] = f->child[i];
	}
	*n_feat_ = n_feat;
	return r;
}
