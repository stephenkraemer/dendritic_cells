#!/usr/bin/env bash

cd /icgc/dkfzlsdf/analysis/hs_ontogeny/databases/region_set_profiler_databases
# manually downloaded genesets from msigdb webpage
# names are quite cryptic, rename
mv h.all.v6.2.symbols.gmt hallmarks.gmt
mv c2.all.v6.2.symbols.gmt curated-genesets.gmt
mv c2.cp.v6.2.symbols.gmt canonical-pathways.gmt
mv c5.bp.v6.2.symbols.gmt go_biological-process.gmt
mv c5.cc.v6.2.symbols.gmt go_cellular-component.gmt
mv c5.mf.v6.2.symbols.gmt go_molecular-function.gmt
mv c6.all.v6.2.symbols.gmt oncogenic-signatures.gmt
mv c7.all.v6.2.symbols.gmt immunologic-signatures.gmt

