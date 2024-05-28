
# admixture results
'admixture_files' folder is downloaded from 
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/admixture_files/

# plink file used in admixture results
ped/map files are converted to plink 2 format using
plink2 --pedmap ALL.wgs.phase3_shapeit2_filtered.20141217.maf0.05 --make-pgen --out ALL.wgs.phase3_shapeit2_filtered.20141217.maf0.05
plink2 --pedmap admixture_files/ALL.wgs.phase3_shapeit2_filtered.20141217.maf0.05 --make-bed --out ALL.wgs.phase3_shapeit2_filtered.20141217.maf0.05


# all_phase3_ns files
downloaded from https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg
(2016-05-05 primary release (build 37, 2504 samples)); no singleton option

## decompress 
plink2 --zst-decompress all_phase3_ns.pvar.zst > all_phase3_ns.pvar
plink2 --zst-decompress all_phase3_ns.pgen.zst > all_phase3_ns.pgen

## data with only SNPS used in admixture results
plink2 --pfile 'vzs' all_phase3_ns --snps-only --maf 0.05 --make-pgen --out all_phase3_ns_maf0.05
run 'code/build_AA.r' to create 'data/AA.rds'
