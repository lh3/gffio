#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "minigff.h"
#include "ketopt.h"

int main_view(int argc, char *argv[])
{
	mgf_gff_t *gff;
	ketopt_t o = KETOPT_INIT;
	int32_t c;
	while ((c = ketopt(&o, argc, argv, 1, "v:", 0)) >= 0) {
		if (c == 'v') mgf_verbose = atoi(o.arg);
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: minigff view [options] <in.gff>\n");
		return 1;
	}
	gff = mgf_read(argv[o.ind]);

	/*
	char *str = 0;
	int32_t i, n_fs, len = 0, cap = 0;
	const mgf_feat_t **fs = mgf_get_by_id(gff, "ENST00000525758.1", 0, &n_fs);
	for (i = 0; i < n_fs; ++i)
		mgf_write_feat(&str, &len, &cap, gff, fs[i], MGF_FMT_GFF3);
	fputs(str, stdout);
	*/

	char *str = 0;
	int32_t i, n_fs, len = 0, cap = 0;
	const mgf_feat_t **fs;
	fs = mgf_toposort(gff);
	for (i = 0; i < gff->n_feat; ++i) {
		len = 0;
		mgf_write_feat(&str, &len, &cap, gff, fs[i], MGF_FMT_GFF3);
		fputs(str, stdout);
	}

	mgf_write(0, gff, MGF_FMT_GFF3);

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
