#!/usr/bin/env bash

cd /icgc/dkfzlsdf/analysis/hs_ontogeny/databases/gene_annotations
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M19/gencode.vM19.annotation.gtf.gz
zcat gencode.vM19.annotation.gtf.gz > gencode.vM19.annotation.gtf
grep 'tag "appris_principal_' gencode.vM19.annotation.gtf > gencode.vM19.annotation.appris-principal.gtf
wc -l gencode.vM19.annotation.gtf
wc -l gencode.vM19.annotation.appris-principal.gtf
sed 's:^chr::' gencode.vM19.annotation.gtf > gencode.vM19.annotation.no-prefix.gtf
sed 's:^chr::' gencode.vM19.annotation.appris-principal.gtf > gencode.vM19.annotation.appris-principal.no-prefix.gtf
