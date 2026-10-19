# Julio92-C/pneuma: Citations

If you use PNEUMA, please also cite the tools that produced your results. The list below covers the tools and databases the pipeline actually runs (Phase 1). Tools of the planned Phase 2 and Phase 3 will be added when they are implemented.

## Framework

- [nf-core](https://pubmed.ncbi.nlm.nih.gov/32055031/)

  > Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. doi: [10.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x).

- [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/)

  > Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: [10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820).

## Preprocessing

- [Porechop_ABI](https://github.com/bonsai-team/Porechop_ABI) (optional adapter trimming, `--run_porechop`)

  > Bonenfant Q, Noé L, Touzet H. Porechop_ABI: discovering unknown adapters in Oxford Nanopore Technology sequencing reads for downstream trimming. Bioinformatics Advances. 2022;3(1):vbac085. doi: [10.1093/bioadv/vbac085](https://doi.org/10.1093/bioadv/vbac085).

- [nanoq](https://github.com/esteinig/nanoq)

  > Steinig E, Coin L. Nanoq: ultra-fast quality control for nanopore reads. Journal of Open Source Software. 2022;7(69):2991. doi: [10.21105/joss.02991](https://doi.org/10.21105/joss.02991).

- [NanoPlot](https://github.com/wdecoster/NanoPlot) (NanoPack)

  > De Coster W, D'Hert S, Schultz DT, Cruts M, Van Broeckhoven C. NanoPack: visualizing and processing long-read sequencing data. Bioinformatics. 2018;34(15):2666-2669. doi: [10.1093/bioinformatics/bty149](https://doi.org/10.1093/bioinformatics/bty149).

  > De Coster W, Rademakers R. NanoPack2: population-scale evaluation of long-read sequencing data. Bioinformatics. 2023;39(5):btad311. doi: [10.1093/bioinformatics/btad311](https://doi.org/10.1093/bioinformatics/btad311).

- [SeqKit](https://bioinf.shenwei.me/seqkit/)

  > Shen W, Le S, Li Y, Hu F. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PLoS ONE. 2016;11(10):e0163962. doi: [10.1371/journal.pone.0163962](https://doi.org/10.1371/journal.pone.0163962).

- [Hostile](https://github.com/bede/hostile) (host read removal; `human-t2t-hla` index)

  > Constantinides B, Hunt M, Crook DW. Hostile: accurate decontamination of microbial host sequences. Bioinformatics. 2023;39(12):btad728. doi: [10.1093/bioinformatics/btad728](https://doi.org/10.1093/bioinformatics/btad728).

  > Nurk S, Koren S, Rhie A, et al. The complete sequence of a human genome. Science. 2022;376(6588):44-53. doi: [10.1126/science.abj6987](https://doi.org/10.1126/science.abj6987). (T2T-CHM13 reference underlying the Hostile index.)

- [MultiQC](https://pubmed.ncbi.nlm.nih.gov/27312411/)

  > Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8. doi: [10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354).

## Taxonomy

- [Kraken 2](https://github.com/DerrickWood/kraken2) (PlusPF index, [Langmead lab index collection](https://benlangmead.github.io/aws-indexes/k2))

  > Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2. Genome Biology. 2019;20:257. doi: [10.1186/s13059-019-1891-0](https://doi.org/10.1186/s13059-019-1891-0).

- [Bracken](https://github.com/jenniferlu717/Bracken)

  > Lu J, Breitwieser FP, Thielen P, Salzberg SL. Bracken: estimating species abundance in metagenomics data. PeerJ Computer Science. 2017;3:e104. doi: [10.7717/peerj-cs.104](https://doi.org/10.7717/peerj-cs.104).

- [KrakenTools](https://github.com/jenniferlu717/KrakenTools) (report combination)

  > Lu J, Rincon N, Wood DE, Breitwieser FP, Pockrandt C, Langmead B, Salzberg SL, Steinegger M. Metagenome analysis using the Kraken software suite. Nature Protocols. 2022;17:2815-2839. doi: [10.1038/s41596-022-00738-y](https://doi.org/10.1038/s41596-022-00738-y).

- [Krona](https://github.com/marbl/Krona)

  > Ondov BD, Bergman NH, Phillippy AM. Interactive metagenomic visualization in a Web browser. BMC Bioinformatics. 2011;12:385. doi: [10.1186/1471-2105-12-385](https://doi.org/10.1186/1471-2105-12-385).

- [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) (taxdump)

  > Schoch CL, Ciufo S, Domrachev M, et al. NCBI Taxonomy: a comprehensive update on curation, resources and tools. Database. 2020;2020:baaa062. doi: [10.1093/database/baaa062](https://doi.org/10.1093/database/baaa062).

## Contamination

- [Recentrifuge](https://github.com/khyox/recentrifuge)

  > Martí JM. Recentrifuge: robust comparative analysis and contamination removal for metagenomics. PLoS Computational Biology. 2019;15(4):e1006967. doi: [10.1371/journal.pcbi.1006967](https://doi.org/10.1371/journal.pcbi.1006967).

## Functional annotation (resistome, virulome, plasmids)

- [ABRicate](https://github.com/tseemann/abricate)

  > Seemann T. ABRicate: mass screening of contigs for antimicrobial resistance or virulence genes. GitHub: <https://github.com/tseemann/abricate> (no publication).

- [CARD](https://card.mcmaster.ca/)

  > Alcock BP, Huynh W, Chalil R, et al. CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. Nucleic Acids Research. 2023;51(D1):D690-D699. doi: [10.1093/nar/gkac920](https://doi.org/10.1093/nar/gkac920).

- [VFDB](http://www.mgc.ac.cn/VFs/)

  > Liu B, Zheng D, Zhou S, Chen L, Yang J. VFDB 2022: a general classification scheme for bacterial virulence factors. Nucleic Acids Research. 2022;50(D1):D912-D917. doi: [10.1093/nar/gkab1107](https://doi.org/10.1093/nar/gkab1107).

- [PlasmidFinder](https://cge.food.dtu.dk/services/PlasmidFinder/)

  > Carattoli A, Zankari E, García-Fernández A, Voldby Larsen M, Lund O, Villa L, Møller Aarestrup F, Hasman H. In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing. Antimicrobial Agents and Chemotherapy. 2014;58(7):3895-3903. doi: [10.1128/AAC.02412-14](https://doi.org/10.1128/AAC.02412-14).

## Downstream analysis

- [MERIDIAN](https://github.com/Julio92-C/meridian) (container `ghcr.io/julio92-c/meridian:v1.0.1`)

  > Ortega Cambara JC. MERIDIAN: Metagenomic Evaluation of Resistance, Identity & Diversity through Integrated Analysis of Nanopore sequencing. GitHub: <https://github.com/Julio92-C/meridian>. Manuscript in preparation.

## Software packaging/containerisation tools

- [Bioconda](https://pubmed.ncbi.nlm.nih.gov/29967506/)

  > Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018 Jul;15(7):475-476. doi: [10.1038/s41592-018-0046-7](https://doi.org/10.1038/s41592-018-0046-7).

- [BioContainers](https://pubmed.ncbi.nlm.nih.gov/28379341/)

  > da Veiga Leprevost F, Grüning B, Aflitos SA, Röst HL, Uszkoreit J, Barsnes H, Vaudel M, Moreno P, Gatto L, Weber J, Bai M, Jimenez RC, Sachsenberg T, Pfeuffer J, Alvarez RV, Griss J, Nesvizhskii AI, Perez-Riverol Y. BioContainers: an open-source and community-driven framework for software standardization. Bioinformatics. 2017 Aug 15;33(16):2580-2582. doi: [10.1093/bioinformatics/btx192](https://doi.org/10.1093/bioinformatics/btx192).

- [Docker](https://dl.acm.org/doi/10.5555/2600239.2600241)

  > Merkel D. Docker: lightweight Linux containers for consistent development and deployment. Linux Journal. 2014;2014(239):2. doi: 10.5555/2600239.2600241.

- [Singularity](https://pubmed.ncbi.nlm.nih.gov/28494014/)

  > Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. doi: [10.1371/journal.pone.0177459](https://doi.org/10.1371/journal.pone.0177459).
