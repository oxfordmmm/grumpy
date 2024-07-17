use std::time::SystemTime;

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
    let vcf_start = SystemTime::now();
    let vcf = VCFFile::new("/home/jeremy/Documents/work/gnomonicus/e4b54393-e037-435e-b942-527f0e2b1616.merged.vcf".to_string(), false, 5);
    // let vcf = VCFFile::new("/home/jeremy/Downloads/e4b54393-e037-435e-b942-527f0e2b1616_final.full.vcf".to_string(), false, 5);
    let vcf_end = SystemTime::now();

    let reference_start = SystemTime::now();
    let mut reference = Genome::new("reference/NC_000962.3.gbk");
    let reference_end = SystemTime::now();

    let sample_start = SystemTime::now();
    let mut sample = mutate(&reference, vcf);
    let sample_end = SystemTime::now();
    // sample.build_all_genes();

    // let gene_name = "ponA1".to_string();

    // let gene_diff = GeneDifference::new(reference.get_gene(gene_name.clone()), sample.get_gene(gene_name.clone()), MinorType::COV);
    // println!("{:?}", gene_diff.mutations.iter().map(|mutation| mutation.mutation.clone()).collect::<Vec<String>>());
    // println!("{:?}\n", gene_diff.minor_mutations.iter().map(|mutation| mutation.mutation.clone()).collect::<Vec<String>>());
    let genome_start = SystemTime::now();
    let difference = GenomeDifference::new(reference.clone(), sample.clone(), MinorType::COV);
    println!("{:?}", difference.variants.iter().map(|variant| variant.variant.clone()).collect::<Vec<String>>());
    let genome_end = SystemTime::now();

    let gene_start = SystemTime::now();
    for gene_name in sample.genes_with_mutations.clone().iter(){
        println!("{}", gene_name);
        let gene_diff = GeneDifference::new(reference.get_gene(gene_name.clone()), sample.get_gene(gene_name.clone()), MinorType::COV);
        println!("{:?}", gene_diff.mutations.iter().map(|mutation| mutation.mutation.clone()).collect::<Vec<String>>());
        println!("{:?}\n", gene_diff.minor_mutations.iter().map(|mutation| mutation.mutation.clone()).collect::<Vec<String>>());
    }
    let gene_end = SystemTime::now();

    println!("VCF took {:?}", vcf_end.duration_since(vcf_start).unwrap());
    println!("Reference took {:?}", reference_end.duration_since(reference_start).unwrap());
    println!("Sample took {:?}", sample_end.duration_since(sample_start).unwrap());
    println!("Genome diff took {:?}", genome_end.duration_since(genome_start).unwrap());
    println!("Gene diff took {:?}", gene_end.duration_since(gene_start).unwrap());


    // let gene_name = "fabG1".to_string();
    // let rrs_diff = GeneDifference::new(reference.get_gene(gene_name.clone()), sample.get_gene(gene_name.clone()), MinorType::COV);
    // println!("{:?}", rrs_diff.mutations.iter().map(|mutation| mutation.mutation.clone()).collect::<Vec<String>>());
    // println!("{:?}", rrs_diff.minor_mutations.iter().map(|mutation| mutation.mutation.clone()).collect::<Vec<String>>());



}