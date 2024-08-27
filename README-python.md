# grumpy - Re-implementation of [gumpy](https://github.com/oxfordmmm/gumpy) in Rust for speed
Does not implement the same interface, but rather the same premise. 

Up to 100x faster than `gumpy`

## Install
```
pip install bio-grumpy
```

## Usage
```
import grumpy

# Parse a genbank file
ref = grumpy.Genome("some/path/to/a/genbank/file.gbk")

# Parse a VCF file, respecting filter fails with a MIN_DP of 3 reads to make a call
vcf = grumpy.VCFFile("some/path/to/a/vcf/file.vcf", False, 3)

# Apply the VCF's mutations to the genome
sample = grumpy.mutate(ref, vcf)

# Get the genome level differences
genome_diff = grumpy.GenomeDifference(ref, sample)
for variant in genome_diff.variants:
    print(variant.variant)
# And minor alleles
for variant in genome_diff.minor_variants:
    print(variant.variant)

# Get gene level differences for all genes with mutations
for gene_name in sample.genes_with_mutations:
    print(gene_name)
    gene_diff = grumpy.GeneDifference(
            ref.get_gene(gene_name),
            sample.get_gene(gene_name),
            grumpy.MinorType.COV,
        )
    for mutation in gene_diff.mutations:
        print(mutation.mutation)
    # And minor alleles
    for mutation in gene_diff.minor_mutations:
        print(mutation.mutation)
```

## Threads
By default, VCF parsing will run with the same number of threads as CPU cores. This can give significant speed improvements for complex VCF files

To force single threading (or other number of threads), run `grumpy.thread_setup(<thread count>)`. Note that should be run before **any** other function in this library, and should not be called more than once! Not following this will cause errors!