import rpy2.robjects as robjects
from rpy2.robjects import r, ListVector
from rpy2.robjects.packages import importr
import rpy2.robjects.vectors as rvectors
from rpy2.rinterface_lib.sexp import NALogicalType
import glob
import os
import re
import csv

nwk_files = glob.glob('ctree/*.nwk')
nwk_files_abs = [os.path.abspath(file) for file in nwk_files]
lineages = [os.path.basename(file).split('.nwk')[0] for file in nwk_files]

parallel = importr("parallel")

robjects.r('set.seed(123456)')

r_lineages = rvectors.StrVector(lineages)
r_data = rvectors.StrVector(nwk_files_abs)

r('''
    find_ne <- function(file) {
        library(ape)
        library(LambdaSkyline)

        tryCatch({
            tree <- read.tree(file)
            tree <- multi2di(tree)
            tree <- collapse.singles(tree)
            alpha <- betacoal.maxlik(tree)
            sky <- skyline.multi.phylo(tree, alpha$p1)
            pop_sizes <- head(sky$population.size, n = 5)
            return(mean(pop_sizes, na.rm = TRUE))
        }, error = function(e) {
            return(NA)
        })
    }
''')

r_results = parallel.mclapply(r_data, r('find_ne'), mc_cores=20)
r_named_results = robjects.r['setNames'](r_results, r_lineages)

results_dict = dict(zip(lineages, r_named_results))
res = {k: 'NA' if isinstance(v[0], NALogicalType) else float(v[0]) for k, v in results_dict.items()}

with open("hunepi_ne.csv", "w", newline='') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(['Lineage', 'Ne'])
    for lineage, ne in res.items():
        writer.writerow([lineage, ne])