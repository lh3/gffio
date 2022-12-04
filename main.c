#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "gffio.h"
#include "ketopt.h"

int main_view(int argc, char *argv[])
{
	gf_gff_t *gff;
	ketopt_t o = KETOPT_INIT;
	int32_t c, to_group = 0, fmt = GF_FMT_GFF3, sel_long = 0;
	char *id_list = 0;
	while ((c = ketopt(&o, argc, argv, 1, "gv:l:tL", 0)) >= 0) {
		if (c == 'v') gf_verbose = atoi(o.arg);
		else if (c == 'g') to_group = 1;
		else if (c == 'l') id_list = o.arg;
		else if (c == 't') fmt = GF_FMT_GTF;
		else if (c == 'L') sel_long = 1;
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: gffio view [options] <in.gff>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -t        output GTF\n");
		fprintf(stderr, "  -g        group by hierarchy\n");
		fprintf(stderr, "  -L        choose the longest CDS/tanscript\n");
		fprintf(stderr, "  -l STR    extract records descended from ID list []\n");
		fprintf(stderr, "  -v INT    verbose level [%d]\n", gf_verbose);
		return 1;
	}
	gff = gf_read(argv[o.ind]);
	if (to_group) gf_group(gff);
	if (sel_long) gf_mrna_choose_long(gff);

	if (id_list) {
		int32_t n;
		char **list;
		list = gf_read_list(id_list, &n);
		gf_write_list(0, gff, fmt, n, list);
	} else {
		gf_write(0, gff, fmt);
	}

	gf_destroy(gff);
	return 0;
}

int main_gff2bed(int argc, char *argv[])
{
	gf_gff_t *gff;
	ketopt_t o = KETOPT_INIT;
	int32_t c, fmt = GF_FMT_BED12L;
	while ((c = ketopt(&o, argc, argv, 1, "sv:eci", 0)) >= 0) {
		if (c == 'v') gf_verbose = atoi(o.arg);
		else if (c == 's') fmt = GF_FMT_BED12S;
		else if (c == 'e') fmt = GF_FMT_BED_EXON;
		else if (c == 'i') fmt = GF_FMT_BED_INTRON;
		else if (c == 'c') fmt = GF_FMT_BED_CDS;
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: gffio gff2bed [options] <in.gff>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -s        output in BED12 with transcript ID only\n");
		fprintf(stderr, "  -e        output exons in 6-column BED (aka BED6)\n");
		fprintf(stderr, "  -i        output introns in BED6\n");
		fprintf(stderr, "  -c        output CDS in BED6\n");
		fprintf(stderr, "  -v INT    verbose level [%d]\n", gf_verbose);
		return 1;
	}
	gff = gf_read(argv[o.ind]);
	gf_write(0, gff, fmt);
	gf_destroy(gff);
	return 0;
}

int main_gff2fa(int argc, char *argv[])
{
	gf_gff_t *gff;
	gf_seqs_t *seq;
	ketopt_t o = KETOPT_INIT;
	int32_t c, fmt = GF_FMT_FA_MRNA;
	while ((c = ketopt(&o, argc, argv, 1, "v:tcp", 0)) >= 0) {
		if (c == 'v') gf_verbose = atoi(o.arg);
		else if (c == 't') fmt = GF_FMT_FA_MRNA;
		else if (c == 'c') fmt = GF_FMT_FA_CDS;
		else if (c == 'p') fmt = GF_FMT_FA_PROTEIN;
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: gffio gff2fa [options] <in.gff> <ref.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -c        extract CDS (mRNA/transcript by default)\n");
		fprintf(stderr, "  -p        extract protein sequences\n");
		fprintf(stderr, "  -v INT    verbose level [%d]\n", gf_verbose);
		return 1;
	}
	gff = gf_read(argv[o.ind]);
	seq = gf_seqs_read(argv[o.ind+1]);
	gf_write_fasta(0, gff, seq, fmt);
	gf_seqs_destroy(seq);
	gf_destroy(gff);
	return 0;
}

static int usage(FILE *fp)
{
	fprintf(fp, "Usage: gffio <command> <arguments>\n");
	fprintf(fp, "Commands:\n");
	fprintf(fp, "  view      read GFF3/GTF\n");
	fprintf(fp, "  gff2bed   convert GFF3/GTF to BED\n");
	fprintf(fp, "  gff2fa    extract sequences\n");
	fprintf(fp, "  version   print version number\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc < 2) return usage(stderr);
	if (strcmp(argv[1], "view") == 0) return main_view(argc-1, argv+1);
	else if (strcmp(argv[1], "gff2bed") == 0) return main_gff2bed(argc-1, argv+1);
	else if (strcmp(argv[1], "gff2fa") == 0) return main_gff2fa(argc-1, argv+1);
	else if (strcmp(argv[1], "version") == 0) {
		puts(GF_VERSION);
		return 0;
	} else {
		fprintf(stderr, "ERROR: unrecognized command '%s'\n", argv[1]);
		return 1;
	}
	return 0;
}
