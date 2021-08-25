#!/bin/sh

rm covizu/data/lineages.csv covizu/data/problematic_sites_sarsCov2.vcf 
wget -O covizu/data/problematic_sites_sarsCov2.vcf https://github.com/W-L/ProblematicSites_SARS-CoV2/blob/master/problematic_sites_sarsCov2.vcf?raw=true
wget -O covizu/data/lineages.csv https://github.com/cov-lineages/pango-designation/blob/master/lineages.csv?raw=true
