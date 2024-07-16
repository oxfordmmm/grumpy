use std::collections::HashMap;

use ordered_float::OrderedFloat;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AltType{
    SNP,
    REF,
    HET,
    NULL,
    INS,
    DEL
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MinorType{
    COV,
    FRS
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct VCFRow{
    pub position: i64,
    pub reference: String,
    pub alternative: Vec<String>,
    pub filter: Vec<String>,
    pub fields: HashMap<String, Vec<String>>,
    pub is_filter_pass: bool,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Evidence{
    pub cov: Option<i32>,
    pub frs: Option<OrderedFloat<f32>>, // Annoyingly f32 doesn't implement Eq so use OrderedFloat
    pub genotype: String, // 1/1 or 0/0 or 0/1 etc
    pub call_type: AltType,
    pub reference: String,
    pub alt: String,
    pub genome_index: i64, // 1-based genome index this refers to
    pub is_minor: bool,
    pub vcf_row: VCFRow,
    pub vcf_idx: i64,
}

#[derive(Clone)]
pub struct GeneDef{
    pub name: String,
    pub coding: bool,
    pub reverse_complement: bool,
    pub start: i64,
    pub end: i64,
    pub promoter_start: i64,
    pub promoter_size: i64,
    pub ribosomal_shifts: Vec<i64>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Alt{
    pub alt_type: AltType,
    pub base: String,
    pub evidence: Evidence
}
