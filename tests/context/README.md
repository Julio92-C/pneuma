# Real-data genomic-context test

`genomic_context_real.nf.test` (tag `context_real`, not run in CI) runs the GENOMIC_CONTEXT subworkflow on public chicken reads (SRR37217670) with the real geNomad and MOB-suite databases and the Medaka model.

Build the fixture with `bash tests/make_context_fixture.sh`, then follow the `Database` and `Run` lines in the test header.

What it checks:

- every read hit of CARD, VFDB and PlasmidFinder gets one of the four context categories; replicons are never `chromosome`;
- `context/hit_context.tsv` has one row per read hit across all databases, with the final call, `context_conflict`, the Kraken2 taxon fields and the contig host vote (`weak_vote`);
- `context/summary/` holds `arg_context_by_sample.tsv`, `arg_context_by_class.tsv`, `gene_context.tsv` and `taxon_resolution.tsv`.

The database build uses `--skip_pathogen_download`, which skips only the large pathogen databases (CheckM2, Bakta, AMRFinderPlus); the geNomad database, MOB-suite and the Medaka model are still built with the mobilome.
