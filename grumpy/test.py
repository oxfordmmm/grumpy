import grumpy

ref = grumpy.Genome("reference/NC_000962.3.gbk")

vcf = grumpy.VCFFile("/home/jeremy/Documents/work/gnomonicus/e4b54393-e037-435e-b942-527f0e2b1616.merged.vcf", False, 5)

sample = grumpy.mutate(ref, vcf)

diff = grumpy.GenomeDifference(ref, sample, grumpy.MinorType.COV)

for variant in diff.variants:
    print(variant.variant, variant.evidence.position, variant.evidence.reference, variant.evidence.alternative, variant.evidence.filter, variant.evidence.fields)
    for key, val in variant.evidence.fields.items():
        print(key, val)

# for variant in diff.minor_variants:
#     print(variant.variant)
