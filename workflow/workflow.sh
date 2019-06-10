##mkdir ../results/diff_tests

Rscript ../scripts/diff_test.r --response ../data/clinical_data.csv --expression ../data/metab_conc_data.csv --logic_variable adhd_16_case --id_variable id --o ../results/diff_tests/metab_conc_diff.csv

Rscript ../scripts/diff_test.r --response ../data/clinical_data.csv --expression ../data/metab_pct_data.csv --logic_variable adhd_16_case --id_variable id --o ../results/diff_tests/metab_pct_diff.csv

Rscript diff_test.r --response ../data/clinical_data.csv --expression ../data/methyl_data.csv --logic_variable adhd_16_case --id_variable id --o ../results/diff_tests/methyl_diff.csv

#conc
mkdir ../results/conc

# CONC WGCNA
mkdir ../results/conc/wgcna

Rscript ../scripts/WGCNA_script.r --matrix ../data/metab_conc_data.csv --power 8 --module_size 5 --o ../results/conc/wgcna/metab_conc_module

Rscript ../scripts/enrichment_PCA.r --matrix ../data/metab_conc_data.csv --module ../results/conc/wgcna/metab_conc_module.csv --metadata ../data/metab_meta.csv --variable Final_Group --id id --o ../results/conc/wgcna/metab_conc

echo "\nWGCNA_O2PLS\n"
Rscript ../scripts/O2_PLS.r --response ../results/conc/wgcna/metab_conc_module_res.csv --predictor ../data/methyl_data_sig.csv --o ../results/conc/wgcna/metab_conc_module

echo "\nWGCNA_SPLS\n"
Rscript ../scripts/sparse_PLS.r --response ../results/conc/wgcna/metab_conc_module_res.csv --predictor ../data/methyl_data_sig.csv --o ../results/conc/wgcna/metab_conc_module

Rscript ../scripts/CCA.r --response ../results/conc/wgcna/metab_conc_module_res.csv --predictor ../data/methyl_data_sig.csv --type single_single --reps 5 --p 4 --o ../results/conc/wgcna/metab_conc --it 10

Rscript ../scripts/CCA.r --response ../results/conc/wgcna/metab_conc_module_res.csv --predictor ../data/methyl_data_sig.csv --type many_single --reps 3 --p 4 --o ../results/conc/wgcna/metab_conc --it 5

echo "\nWGCNA_CCA_MM\n"
Rscript ../scripts/CCA.r --response ../results/conc/wgcna/metab_conc_module_res.csv --predictor ../data/methyl_data_sig.csv --type many_many --reps 2 --p 7 --o ../results/conc/wgcna/metab_conc --it 20

# CONC Biological group
mkdir ../results/conc/biological

Rscript ../scripts/biological_group.r

Rscript ../scripts/enrichment_PCA.r --matrix ../data/metab_conc_data.csv --module ../results/conc/biological/metab_bio_group.csv --metadata ../data/metab_meta.csv --variable Final_Group --id id --o ../results/conc/biological/metab_bio_group

echo "\nBIO_O2PLS\n"
Rscript ../scripts/O2_PLS.r --response ../results/conc/biological/metab_bio_group_module_res.csv --predictor ../data/methyl_data_sig.csv --o ../results/conc/biological/metab_bio

echo "\nBIO_SPLS\n"
Rscript ../scripts/sparse_PLS.r --response ../results/conc/biological/metab_bio_group_module_res.csv --predictor ../data/methyl_data_sig.csv --o ../results/conc/biological/metab_bio

Rscript ../scripts/CCA.r --response ../results/conc/biological/metab_bio_group_module_res.csv --predictor ../data/methyl_data_sig.csv --type single_single --reps 5 --p 4 --o ../results/conc/biological/metab_bio --it 10

Rscript ../scripts/CCA.r --response ../results/conc/biological/metab_bio_group_module_res.csv --predictor ../data/methyl_data_sig.csv --type many_single --reps 3 --p 4 --o ../results/conc/biological/metab_bio --it 5

echo "\nBIO_CCA_MM\n"
Rscript ../scripts/CCA.r --response ../results/conc/biological/metab_bio_group_module_res.csv --predictor ../data/methyl_data_sig.csv --type many_many --reps 2 --p 4 --o ../results/conc/biological/metab_bio --it 20

Rscript ../scripts/combine_groups.r --biological ../results/conc/biological/metab_bio_group_module_res.csv --WGCNA ../results/conc/wgcna/metab_conc_module_res.csv --o ../results/conc/wgcna_biological_comb.csv

# CONC COLINEAR
mkdir ../results/conc/colinear

Rscript ../scripts/colinearity_script.r --matrix ../data/metab_conc_data.csv --o ../results/conc/colinear/metab_conc_res.csv

echo "\nCOL_O2PLS\n"
Rscript ../scripts/O2_PLS.r --response ../results/conc/colinear/metab_conc_res.csv --predictor ../data/methyl_data_sig.csv --o ../results/conc/colinear/metab_conc

echo "\nCOL_SPLS\n"
Rscript ../scripts/sparse_PLS.r --response ../results/conc/colinear/metab_conc_res.csv --predictor ../data/methyl_data_sig.csv --o ../results/conc/colinear/metab_conc

Rscript ../scripts/CCA.r --response ../results/conc/colinear/metab_conc_res.csv --predictor ../data/methyl_data_sig.csv --type single_single --reps 5 --p 4 --o ../results/conc/colinear/metab_conc_ --it 10

Rscript ../scripts/CCA.r --response ../results/conc/colinear/metab_conc_res.csv --predictor ../data/methyl_data_sig.csv --type many_single --reps 3 --p 4 --o ../results/conc/colinear/metab_conc --it 5

echo "\nCOL_CCA_MM\n"
Rscript ../scripts/CCA.r --response ../results/conc/colinear/metab_conc_res.csv --predictor ../data/methyl_data_sig.csv --type many_many --reps 2 --p 4 --o ../results/conc/colinear/metab_conc --it 20

#conc
mkdir ../results/pct

# PCT WGCNA
mkdir ../results/pct/wgcna

Rscript ../scripts/WGCNA_script.r --matrix ../data/metab_pct_data.csv --power 8 --module_size 5 --o ../results/pct/wgcna/metab_pct_module

Rscript ../scripts/enrichment_PCA.r --matrix ../data/metab_pct_data.csv --module ../results/pct/wgcna/metab_pct_module.csv --metadata ../data/metab_meta.csv --variable Final_Group --id id --o ../results/pct/wgcna/metab_pct

Rscript ../scripts/O2_PLS.r --response ../results/pct/wgcna/metab_pct_module_res.csv --predictor ../data/methyl_data_sig.csv --o ../results/pct/wgcna/metab_pct_module

Rscript ../scripts/sparse_PLS.r --response ../results/pct/wgcna/metab_pct_module_res.csv --predictor ../data/methyl_data_sig.csv --o ../results/pct/wgcna/metab_pct_module

Rscript ../scripts/CCA.r --response ../results/pct/wgcna/metab_pct_module_res.csv --predictor ../data/methyl_data_sig.csv --type single_single --reps 5 --p 4 --o ../results/pct/wgcna/metab_pct_ --it 10

Rscript ../scripts/CCA.r --response ../results/pct/wgcna/metab_pct_module_PC1_res--predictor ../data/methyl_data_sig.csv --type many_single --reps 5 --p 4 --o ../results/pct/wgcna/metab_pct_ --it 10

Rscript ../scripts/CCA.r --response ../results/pct/wgcna/metab_pct_module_res.csv --predictor ../data/methyl_data_sig.csv --type many_many --reps 5 --p 4 --o ../results/pct/wgcna/metab_pct_ --it 10

# # PCT COLINEAR
mkdir ../results/pct/colinear

Rscript ../scripts/colinearity_script.r --matrix ../data/metab_pct_data.csv --o ../results/pct/colinear/metab_pct_res.csv

Rscript ../scripts/O2_PLS.r --response ../results/pct/colinear/metab_pct_res.csv --predictor ../data/methyl_data_sig.csv --o ../results/pct/colinear/metab_pct

Rscript ../scripts/sparse_PLS.r --response ../results/pct/colinear/metab_pct_res.csv --predictor ../data/methyl_data_sig.csv --o ../results/pct/colinear/metab_pct

Rscript ../scripts/CCA.r --response ../results/pct/colinear/metab_pct_res.csv --predictor ../data/methyl_data_sig.csv --type single_single --reps 5 --p 4 --o ../results/pct/colinear/metab_pct_ --it 10

Rscript ../scripts/CCA.r --response ../results/pct/colinear/metab_pct_res.csv --predictor ../data/methyl_data_sig.csv --type many_many --reps 5 --p 4 --o ../results/pct/colinear/metab_pct_ --it 3
