use pyo3::prelude::*;

pub mod common;
pub mod genome;
pub mod gene;
pub mod vcf;
pub mod difference;

// use difference::GeneDifference;
// use genome::mutate;
// use genome::Genome;
// use vcf::VCFFile;
// use difference::GenomeDifference;
// use common::MinorType;

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


