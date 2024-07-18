use pyo3::prelude::*;
use std::collections::HashMap;

use ordered_float::OrderedFloat;

#[pyclass(eq, eq_int)]
#[derive(Debug, Clone, PartialEq, Eq)]
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
pub enum MinorType{
    COV,
    FRS
}

#[pyclass(eq)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct VCFRow{
    #[pyo3(get, set)]
    pub position: i64,

    #[pyo3(get, set)]
    pub reference: String,

    #[pyo3(get, set)]
    pub alternative: Vec<String>,

    #[pyo3(get, set)]
    pub filter: Vec<String>,

    #[pyo3(get, set)]
    pub fields: HashMap<String, Vec<String>>,

    #[pyo3(get, set)]
    pub is_filter_pass: bool,
}

#[pyclass(eq)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Evidence{
    #[pyo3(get, set)]
    pub cov: Option<i32>,

    pub frs: Option<OrderedFloat<f32>>, // Annoyingly f32 doesn't implement Eq so use OrderedFloat

    #[pyo3(get, set)]
    pub genotype: String, // 1/1 or 0/0 or 0/1 etc

    #[pyo3(get, set)]
    pub call_type: AltType,

    #[pyo3(get, set)]
    pub reference: String,

    #[pyo3(get, set)]
    pub alt: String,

    #[pyo3(get, set)]
    pub genome_index: i64, // 1-based genome index this refers to

    #[pyo3(get, set)]
    pub is_minor: bool,

    #[pyo3(get, set)]
    pub vcf_row: VCFRow,

    #[pyo3(get, set)]
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
pub struct GeneDef{
    #[pyo3(get, set)]
    pub name: String,

    #[pyo3(get, set)]
    pub coding: bool,

    #[pyo3(get, set)]
    pub reverse_complement: bool,

    #[pyo3(get, set)]
    pub start: i64,

    #[pyo3(get, set)]
    pub end: i64,

    #[pyo3(get, set)]
    pub promoter_start: i64,

    #[pyo3(get, set)]
    pub promoter_size: i64,

    #[pyo3(get, set)]
    pub ribosomal_shifts: Vec<i64>,
}

#[pyclass(eq)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Alt{
    #[pyo3(get, set)]
    pub alt_type: AltType,

    #[pyo3(get, set)]
    pub base: String,

    #[pyo3(get, set)]
    pub evidence: Evidence
}
