//! Grumpy, genetic analysis in Rust.
//! 
//! This library provides a set of tools for genetic analysis, including:
//! - Genome representation
//! - Gene representation
//! - VCF file parsing
//! - Finding effects of a given VCF file at both genome and gene levels
//! 
//! # Example
//! ```
//! use grumpy::genome::{Genome, mutate};
//! use grumpy::vcf::VCFFile;
//! use grumpy::difference::{GenomeDifference, GeneDifference};
//! 
//! let reference = Genome::new("some/path/to/a/genbank/file.gb");
//! let vcf = VCFFile::new("some/path/to/a/vcf/file.vcf", false, 3);
//! let sample = mutate(&reference, vcf);
//! 
//! let genome_diff = GenomeDifference::new(reference, sample);
//! for variant in genome_diff.variants.iter(){
//!    println!("{}", variant.variant);
//! }
//! 
//! for gene_name in sample.genes_with_mutations.clone().iter(){
//!   let gene_diff = GeneDifference::new(reference.get_gene(gene_name.clone()), sample.get_gene(gene_name.clone()));
//!   for mutation in gene_diff.mutations.iter(){
//!     println!("{}", mutation.mutation);
//!   }
//! }
//! 
//! ```
//! 
//! Also provides an interface to this library as a Python module using PyO3.
//! `pip install bio-grumpy`
use pyo3::prelude::*;

pub mod common;
pub mod genome;
pub mod gene;
pub mod vcf;
pub mod difference;

/// A Python module implemented in Rust.
#[pymodule]
fn grumpy(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<common::AltType>()?;
    m.add_class::<common::MinorType>()?;
    m.add_class::<common::VCFRow>()?;
    m.add_class::<common::Evidence>()?;

    m.add_class::<genome::Genome>()?;

    m.add_class::<gene::Gene>()?;

    m.add_class::<vcf::VCFFile>()?;

    m.add_class::<difference::GenomeDifference>()?;
    m.add_class::<difference::GeneDifference>()?;
    
    m.add_function(wrap_pyfunction!(genome::mutate, m)?)?;

    Ok(())
}


