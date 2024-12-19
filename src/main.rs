use std::env;
use std::time::SystemTime;

pub mod common;
pub mod difference;
pub mod gene;
pub mod genome;
pub mod vcf;

use common::MinorType;
use difference::GeneDifference;
use difference::GenomeDifference;
use genome::mutate;
use genome::Genome;
use vcf::VCFFile;

#[cfg(not(tarpaulin_include))]
fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        println!("Usage: grumpy <reference_genome> <vcf_file>");
        return;
    }
    let reference_path = &args[1];
    let vcf_path = &args[2];

    let vcf_start = SystemTime::now();
    let vcf = VCFFile::new(vcf_path.to_string(), false, 3);
    let vcf_end = SystemTime::now();

    let reference_start = SystemTime::now();
    let mut reference = Genome::new(reference_path);
    let reference_end = SystemTime::now();

    let sample_start = SystemTime::now();
    let mut sample = mutate(&reference, vcf);
    let sample_end = SystemTime::now();

    let target_gene = "PE1";

    let genome_start = SystemTime::now();
    let mut difference = GenomeDifference::new(reference.clone(), sample.clone(), MinorType::COV);
    let genome_end = SystemTime::now();
    // println!("{:?}", difference.variants.iter().map(|variant| variant.variant.clone()).collect::<Vec<String>>());
    for variant in difference.variants.iter_mut() {
        if variant.gene_name.clone().is_none()
            || (variant.gene_name.clone().is_some()
                && variant.gene_name.clone().unwrap() != target_gene)
        {
            continue;
        }
        println!(
            "{:?}@{:?} --> {:?}",
            variant.gene_name,
            variant.gene_position,
            variant.variant.clone()
        );
    }
    for variant in difference.minor_variants.iter_mut() {
        if variant.gene_name.clone().is_none()
            || (variant.gene_name.clone().is_some()
                && variant.gene_name.clone().unwrap() != target_gene)
        {
            continue;
        }
        println!(
            "{:?}@{:?} --> {:?}",
            variant.gene_name,
            variant.gene_position,
            variant.variant.clone()
        );
    }

    let gene_start = SystemTime::now();
    for gene_name in sample.genes_with_mutations.clone().iter() {
        if gene_name != target_gene {
            continue;
        }
        println!("{}", gene_name);
        let gene_diff = GeneDifference::new(
            reference.get_gene(gene_name.clone()),
            sample.get_gene(gene_name.clone()),
            MinorType::COV,
        );
        println!(
            "{:?}",
            gene_diff
                .mutations
                .iter()
                .map(|mutation| mutation.mutation.clone())
                .collect::<Vec<String>>()
        );
        println!(
            "{:?}\n",
            gene_diff
                .minor_mutations
                .iter()
                .map(|mutation| mutation.mutation.clone())
                .collect::<Vec<String>>()
        );
    }
    let gene_end = SystemTime::now();

    println!("\n-----------------------------------\n");
    println!("VCF took {:?}", vcf_end.duration_since(vcf_start).unwrap());
    println!(
        "Reference took {:?}",
        reference_end.duration_since(reference_start).unwrap()
    );
    println!(
        "Sample took {:?}",
        sample_end.duration_since(sample_start).unwrap()
    );
    println!(
        "Genome diff took {:?}",
        genome_end.duration_since(genome_start).unwrap()
    );
    println!(
        "Gene diff took {:?}",
        gene_end.duration_since(gene_start).unwrap()
    );
}
