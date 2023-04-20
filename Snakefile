rule biomart:      
    output:
        "data/tables/biomart.csv"
    shell:
        "Rscript R/start.R"

rule rna_download:
    output:
        "data/RDS/rna-raw.rds"
    shell:
        "Rscript R/rna-download.R"

rule rna_norm:
    input:
        "data/RDS/rna-raw.rds",
        "data/tables/biomart.csv"
    output:
        "data/RDS/rna-clean.rds",
        "data/RDS/rna-norm.rds"
    shell:
        "Rscript R/rna-norm.R"

rule cnv_download:
    output:
        "data/RDS/cnvs-raw.rds"
    shell:
        "Rscript R/cnv-download.R"

rule cnv_clean:
    input:
        "data/RDS/cnvs-raw.rds",
    output:
        "data/RDS/cnvs-clean.rds",
    shell:
        "Rscript R/cnv-clean.R"

rule pairing:
    input:
        "data/RDS/cnvs-clean.rds",
        "data/RDS/rna-norm.rds",
    output:
        "data/RDS/cnvs-paired.rds",
        "data/RDS/rna-paired.rds",
    shell:
        "Rscript R/pairing.R"

rule rna_deg:
    input:
        "data/RDS/rna-paired.rds",
    output:
        "data/plots/rna-volc-deg.png",
	      "data/plots/rna-log2fc-whole.png",
        "data/RDS/rna-deg-NtvsPT.rds",
    shell:
        "Rscript R/deg.R"


