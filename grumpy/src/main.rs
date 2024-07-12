
pub mod common;
pub mod genome;
pub mod gene;
pub mod vcf;

use genome::mutate;
use genome::Genome;
use vcf::VCFFile;


fn main() {
    let mut vcf = VCFFile::new("/home/jeremy/Documents/work/gnomonicus/e4b54393-e037-435e-b942-527f0e2b1616.merged.vcf".to_string(), false);
    let mut reference = Genome::new("reference/NC_000962.3.gbk");
    let mut sample = mutate(&reference, vcf);

    let ref_nc = reference.clone().nucleotide_sequence.clone();
    let mut sample_chars = sample.nucleotide_sequence.chars();
    let mut idx = 0;
    println!("{} {}", ref_nc.len(), sample.nucleotide_sequence.len());
    for (r, a) in ref_nc.chars().zip(sample_chars){
        // let a = sample_chars.nth(idx).unwrap().clone();
        if r != a{
            println!("{} {} {}", idx+1, r, a);
        }
        idx += 1;
    }

    // reference.assign_promoters();
    // for gene in reference.gene_definitions.iter(){
    //     // if gene.promoter_start == -1{
    //     if gene.name == "rpoB" || gene.name == "katG" {
    //         println!("Name {}", gene.name);
    //         println!("Is reverse complement {}", gene.reverse_complement);
    //         println!("Is coding {}", gene.coding);
    //         println!("Start pos {}", gene.start);
    //         println!("End pos {}", gene.end);
    //         println!("Promoter start pos {}", gene.promoter_start);
    //         println!("Promoter size {}", gene.promoter_size);
    //         println!("Ribosomal shifts {:?}", gene.ribosomal_shifts);
    //         println!("");
    //     }
    // }
    // let orf1ab = reference.build_gene("orf1ab".to_string());
    // println!("{:?}", orf1ab);
    // println!("{:?}", orf1ab.nucleotide_sequence.len());
    // println!("{:?}", orf1ab.nucleotide_index.len());
    // println!("{:?}", orf1ab.nucleotide_number.len());

    // let katG = reference.build_gene("katG".to_string());
    // println!("{:?}", katG);
    // println!("{:?}", katG.nucleotide_sequence.len());
    // println!("{:?}", katG.nucleotide_index.len());
    // println!("{:?}", katG.nucleotide_number.len());

    // let rpoB = reference.build_gene("rpoB".to_string());
    // println!("{:?}", rpoB);
    // println!("{:?}", rpoB.nucleotide_sequence.len());
    // println!("{:?}", rpoB.nucleotide_index.len());
    // println!("{:?}", rpoB.nucleotide_number.len());

    // let rrs = reference.build_gene("rrs".to_string());
    // println!("{:?}", rrs);
    // println!("{:?}", rrs.nucleotide_sequence.len());
    // println!("{:?}", rrs.nucleotide_index.len());
    // println!("{:?}", rrs.nucleotide_number.len());

    // let embC = reference.build_gene("embC".to_string());
    // println!("{:?}", embC);
    // println!("{:?}", embC.nucleotide_sequence.len());
    // println!("{:?}", embC.nucleotide_index.len());
    // println!("{:?}", embC.nucleotide_number.len());


}