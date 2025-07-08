tabix -p vcf gnomad.exomes.v4.1.sites.chr20.vcf.bgz
bcftools view -r chr20:5114953-5114953 gnomad.exomes.v4.1.sites.chr20.vcf.bgz -Oz -o region_subset.vcf.gz
tabix -p vcf region_subset.vcf.gz
bcftools view -H region_subset.vcf.gz | awk '{print $1"\t"($2-1)"\t"$2"\t"$3}' > region_subset.bed
