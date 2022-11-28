#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "minigff.h"
#include "ketopt.h"

int main_view(int argc, char *argv[])
{
	mgf_gff_t *gff;
	ketopt_t o = KETOPT_INIT;
	int32_t c, to_group = 0, fmt = MGF_FMT_GFF3;
	char *id_list = 0;
	while ((c = ketopt(&o, argc, argv, 1, "gv:l:t", 0)) >= 0) {
		if (c == 'v') mgf_verbose = atoi(o.arg);
		else if (c == 'g') to_group = 1;
		else if (c == 'l') id_list = o.arg;
		else if (c == 't') fmt = MGF_FMT_GTF;
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: minigff view [options] <in.gff>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -t        output GTF\n");
		fprintf(stderr, "  -g        group by hierarchy\n");
		fprintf(stderr, "  -l STR    extract records descended from ID list []\n");
		fprintf(stderr, "  -v INT    verbose level [%d]\n", mgf_verbose);
		return 1;
	}
	gff = mgf_read(argv[o.ind]);

	if (id_list) {
		int32_t n;
		char **list;
		list = mgf_read_list(id_list, &n);
		mgf_write_list(0, gff, fmt, n, list);
	} else {
		if (to_group) mgf_group(gff);
		mgf_write(0, gff, fmt);
	}

	mgf_destroy(gff);
	return 0;
}

int main_gff2bed(int argc, char *argv[])
{
	mgf_gff_t *gff;
	ketopt_t o = KETOPT_INIT;
	int32_t c, fmt = MGF_FMT_BED12L;
	while ((c = ketopt(&o, argc, argv, 1, "sv:eci", 0)) >= 0) {
		if (c == 'v') mgf_verbose = atoi(o.arg);
		else if (c == 's') fmt = MGF_FMT_BED12S;
		else if (c == 'e') fmt = MGF_FMT_BED_EXON;
		else if (c == 'i') fmt = MGF_FMT_BED_INTRON;
		else if (c == 'c') fmt = MGF_FMT_BED_CDS;
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: minigff gff2bed [options] <in.gff>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -s        output in BED12 with transcript ID only\n");
		fprintf(stderr, "  -e        output exons in 6-column BED (aka BED6)\n");
		fprintf(stderr, "  -i        output introns in BED6\n");
		fprintf(stderr, "  -c        output CDS in BED6\n");
		fprintf(stderr, "  -v INT    verbose level [%d]\n", mgf_verbose);
		return 1;
	}
	gff = mgf_read(argv[o.ind]);
	mgf_write(0, gff, fmt);
	mgf_destroy(gff);
	return 0;
}

int main_gff2fa(int argc, char *argv[])
{
	mgf_gff_t *gff;
	mgf_seqs_t *seq;
	ketopt_t o = KETOPT_INIT;
	int32_t c, fmt = MGF_FMT_FA_MRNA;
	while ((c = ketopt(&o, argc, argv, 1, "v:tcp", 0)) >= 0) {
		if (c == 'v') mgf_verbose = atoi(o.arg);
		else if (c == 't') fmt = MGF_FMT_FA_MRNA;
		else if (c == 'c') fmt = MGF_FMT_FA_CDS;
		else if (c == 'p') fmt = MGF_FMT_FA_PROTEIN;
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: minigff gff2fa [options] <in.gff> <ref.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -c        extract CDS (mRNA/transcript by default)\n");
		fprintf(stderr, "  -p        extract protein sequences\n");
		fprintf(stderr, "  -v INT    verbose level [%d]\n", mgf_verbose);
		return 1;
	}
	gff = mgf_read(argv[o.ind]);
	seq = mgf_seqs_read(argv[o.ind+1]);
	mgf_write_fasta(0, gff, seq, fmt);
	mgf_seqs_destroy(seq);
	mgf_destroy(gff);
	return 0;
}

static int usage(FILE *fp)
{
	fprintf(fp, "Usage: minigff <command> <arguments>\n");
	fprintf(fp, "Commands:\n");
	fprintf(fp, "  view      read GFF3/GTF\n");
	fprintf(fp, "  gff2bed   convert GFF3/GTF to BED\n");
	fprintf(fp, "  gff2fa    extract sequences\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc < 2) return usage(stderr);
	if (strcmp(argv[1], "view") == 0) return main_view(argc-1, argv+1);
	else if (strcmp(argv[1], "gff2bed") == 0) return main_gff2bed(argc-1, argv+1);
	else if (strcmp(argv[1], "gff2fa") == 0) return main_gff2fa(argc-1, argv+1);
	else {
		fprintf(stderr, "ERROR: unrecognized command '%s'\n", argv[1]);
		return 1;
	}
	return 0;
}
