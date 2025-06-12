use std::time::SystemTime;

use clap::Parser;

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

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Path to the reference genbank file
    #[arg(long)]
    reference: String,

    /// Path to the VCF file
    #[arg(long)]
    vcf: String,

    /// Target gene name. If given, pull out all variants and mutations within this gene
    #[arg(long)]
    gene: Option<String>,

    /// Path to the output variable length FASTA file. If given, the sample genome will be written to this file
    #[arg(long)]
    fasta: Option<String>,
}

#[cfg(not(tarpaulin_include))]
fn main() {
    let args = Args::parse();

    let reference_path = args.reference;
    let vcf_path = args.vcf;

    let vcf_start = SystemTime::now();
    let vcf = VCFFile::new(vcf_path.to_string(), false, 3);
    let vcf_end = SystemTime::now();

    let reference_start = SystemTime::now();
    let mut reference = Genome::new(&reference_path);
    let reference_end = SystemTime::now();

    let sample_start = SystemTime::now();
    let mut sample = mutate(&reference, vcf);
    let sample_end = SystemTime::now();

    // If given a gene name, pull out the genome and gene level differences
    if let Some(target_gene) = args.gene {
        let genome_start = SystemTime::now();
        let mut difference =
            GenomeDifference::new(reference.clone(), sample.clone(), MinorType::COV);
        let genome_end = SystemTime::now();
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
            if gene_name != &target_gene {
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

    // If a FASTA path is provided, write the sample genome to it
    if let Some(fasta_path) = args.fasta {
        sample.write_fasta(&fasta_path);
    }
}
