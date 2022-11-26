#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "minigff.h"
#include "ketopt.h"

int main_view(int argc, char *argv[])
{
	mgf_gff_t *gff;
	ketopt_t o = KETOPT_INIT;
	int32_t c, to_group = 0;
	char *id = 0;
	while ((c = ketopt(&o, argc, argv, 1, "gv:d:", 0)) >= 0) {
		if (c == 'v') mgf_verbose = atoi(o.arg);
		else if (c == 'g') to_group = 1;
		else if (c == 'd') id = o.arg;
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: minigff view [options] <in.gff>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -g        group by hierarchy\n");
		fprintf(stderr, "  -v INT    verbose level [%d]\n", mgf_verbose);
		fprintf(stderr, "  -d STR    extract records descended from ID 'STR' []\n");
		return 1;
	}
	gff = mgf_read(argv[o.ind]);

	if (id) {
		const mgf_feat_t *f;
		f = mgf_get_by_id(gff, id);
		if (f) {
			char *str = 0;
			int32_t i, n_fs, len = 0, cap = 0;
			const mgf_feat_t **fs;
			mgf_qbuf_t *b;
			b = mgf_qbuf_init(gff);
			fs = mgf_descend(b, f, &n_fs);
			for (i = 0; i < n_fs; ++i) {
				len = 0;
				mgf_write_feat(&str, &len, &cap, gff, fs[i], MGF_FMT_GFF3);
				fputs(str, stdout);
			}
			mgf_qbuf_destroy(b);
		}
	} else {
		if (to_group) mgf_group(gff);
		mgf_write(0, gff, MGF_FMT_GFF3);
	}

	mgf_destroy(gff);
	return 0;
}

static int usage(FILE *fp)
{
	fprintf(fp, "Usage: minigff <command> <arguments>\n");
	fprintf(fp, "Commands:\n");
	fprintf(fp, "  view      read GTF/GFF3\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc < 2) return usage(stderr);
	if (strcmp(argv[1], "view") == 0) return main_view(argc-1, argv+1);
	else {
		fprintf(stderr, "ERROR: unrecognized command '%s'\n", argv[1]);
		return 1;
	}
	return 0;
}
