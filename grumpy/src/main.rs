
pub mod common;
pub mod genome;
pub mod gene;
pub mod vcf;
pub mod difference;

use difference::GeneDifference;
use genome::mutate;
use genome::Genome;
use vcf::VCFFile;
use difference::GenomeDifference;
use common::MinorType;


fn main() {
    let mut vcf = VCFFile::new("/home/jeremy/Documents/work/gnomonicus/e4b54393-e037-435e-b942-527f0e2b1616.merged.vcf".to_string(), false);
    let mut reference = Genome::new("reference/NC_000962.3.gbk");
    // reference.build_all_genes();
    let mut sample = mutate(&reference, vcf);
    // sample.build_all_genes();

    let difference = GenomeDifference::new(reference.clone(), sample.clone(), MinorType::COV);
    println!("{:?}", difference.variants.iter().map(|variant| variant.variant.clone()).collect::<Vec<String>>());

    for gene_name in sample.genes_with_mutations.clone().iter(){
        println!("{}", gene_name);
        let gene_diff = GeneDifference::new(reference.get_gene(gene_name.clone()), sample.get_gene(gene_name.clone()), MinorType::COV);
        println!("{:?}", gene_diff.mutations.iter().map(|mutation| mutation.mutation.clone()).collect::<Vec<String>>());
        println!("{:?}\n", gene_diff.minor_mutations.iter().map(|mutation| mutation.mutation.clone()).collect::<Vec<String>>());
    }
    // let gene_name = "fabG1".to_string();
    // let rrs_diff = GeneDifference::new(reference.get_gene(gene_name.clone()), sample.get_gene(gene_name.clone()), MinorType::COV);
    // println!("{:?}", rrs_diff.mutations.iter().map(|mutation| mutation.mutation.clone()).collect::<Vec<String>>());
    // println!("{:?}", rrs_diff.minor_mutations.iter().map(|mutation| mutation.mutation.clone()).collect::<Vec<String>>());



}