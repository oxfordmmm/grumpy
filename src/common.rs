//! Module of common structs and enums used throughout the program
use pyo3::prelude::*;
use std::collections::HashMap;

use ordered_float::OrderedFloat;

#[pyclass(eq, eq_int)]
#[derive(Debug, Clone, PartialEq, Eq)]
/// Enum for the types alts can take
pub enum AltType{
    SNP,
    REF,
    HET,
    NULL,
    INS,
    DEL
}

#[pyclass(eq, eq_int)]
#[derive(Debug, Clone, PartialEq, Eq)]
/// Enum for the types of minor evidence
pub enum MinorType{
    COV,
    FRS
}

#[pyclass(eq)]
#[derive(Clone, Debug, Eq, PartialEq)]
/// Struct to hold the information from a VCF row
pub struct VCFRow{
    #[pyo3(get, set)]
    /// Genome position specified in the VCF row
    pub position: i64,

    #[pyo3(get, set)]
    /// Reference base
    pub reference: String,

    #[pyo3(get, set)]
    /// Alt bases
    pub alternative: Vec<String>,

    #[pyo3(get, set)]
    /// Items in the filter column
    pub filter: Vec<String>,

    #[pyo3(get, set)]
    /// Mapping of FORMAT -> SAMPLE fields
    pub fields: HashMap<String, Vec<String>>,

    #[pyo3(get, set)]
    /// True if the filter column passes
    pub is_filter_pass: bool,
}

#[pyclass(eq)]
#[derive(Clone, Debug, Eq, PartialEq)]
/// Struct to hold the information parsed for a call
pub struct Evidence{
    #[pyo3(get, set)]
    /// Coverage at this position. i.e how many reads cover this call
    pub cov: Option<i32>,

    /// FRS at this position. i.e proportion of reads at this position have this alt
    /// 
    /// Annoyingly f32 doesn't implement Eq so use OrderedFloat
    pub frs: Option<OrderedFloat<f32>>,

    #[pyo3(get, set)]
    /// Genotype string from VCF row; 1/1 or 0/0 or 0/1 etc
    pub genotype: String,

    #[pyo3(get, set)]
    /// What type of call this is
    pub call_type: AltType,

    #[pyo3(get, set)]
    /// Reference base
    pub reference: String,

    #[pyo3(get, set)]
    /// Alt call. If SNP/HET/NULL this will be a single base. If INS/DEL this will be the sequence inserted/deleted
    pub alt: String,

    #[pyo3(get, set)]
    /// 1-based genome index this refers to
    pub genome_index: i64,

    #[pyo3(get, set)]
    /// Whether this is a minor call
    pub is_minor: bool,

    #[pyo3(get, set)]
    /// VCF row which this call originated from
    pub vcf_row: VCFRow,

    #[pyo3(get, set)]
    /// Index of the COV field in the VCF row which this call originated from
    pub vcf_idx: i64,
}

#[pymethods]
impl Evidence{
    #[getter]
    fn frs(&self) -> PyResult<i32> {
        match self.frs{
            Some(frs) => Ok(frs.into_inner() as i32),
            None => Ok(0)
        }
    }
}

#[pyclass]
#[derive(Clone)]
/// Struct to hold the information to construct a gene
pub struct GeneDef{
    #[pyo3(get, set)]
    /// Gene name
    pub name: String,

    #[pyo3(get, set)]
    /// Whether this gene codes protein
    pub coding: bool,

    #[pyo3(get, set)]
    /// Whether this gene is reverse complement
    pub reverse_complement: bool,

    #[pyo3(get, set)]
    /// Genome index of the gene start
    pub start: i64,

    #[pyo3(get, set)]
    /// Genome index of the gene end
    pub end: i64,

    #[pyo3(get, set)]
    /// Genome index of the gene promoter start
    pub promoter_start: i64,

    #[pyo3(get, set)]
    /// Number of bases in the promoter
    pub promoter_size: i64,

    #[pyo3(get, set)]
    /// Vec of duplicated positions due to ribosomal shifts
    pub ribosomal_shifts: Vec<i64>,
}

#[pyclass(eq)]
#[derive(Clone, Debug, Eq, PartialEq)]
/// Struct to hold the information of a call
pub struct Alt{
    #[pyo3(get, set)]
    /// Type of the call
    pub alt_type: AltType,

    #[pyo3(get, set)]
    /// Alt call. If SNP/HET/NULL this will be a single base. If INS/DEL this will be the sequence inserted/deleted
    pub base: String,

    #[pyo3(get, set)]
    /// Evidence associated with this call
    pub evidence: Evidence
}
