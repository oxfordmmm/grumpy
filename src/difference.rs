//! Module for handling differences between genomes and genes
use std::collections::{HashMap, HashSet};

use pyo3::prelude::*;

use ordered_float::{Float, OrderedFloat};

use crate::common::{Alt, AltType, Evidence, MinorType};
use crate::gene::{codon_to_aa, Gene, GenePos};
use crate::genome::Genome;

#[pyclass]
#[derive(Clone, Debug, Eq, PartialEq)]
/// Genome level variant
pub struct Variant {
    #[pyo3(get, set)]
    /// GARC for genome level variant
    pub variant: String,

    #[pyo3(get, set)]
    /// Nucleotide index of this variant
    pub nucleotide_index: i64,

    #[pyo3(get, set)]
    /// VCF row index for this variant's evidence
    pub evidence: usize,

    #[pyo3(get, set)]
    /// Index of the VCF row. i.e COV at vcf_idx == coverage for this variant
    pub vcf_idx: Option<i64>,

    #[pyo3(get, set)]
    /// Length of the indel
    pub indel_length: i64,

    #[pyo3(get, set)]
    /// Bases for the indel (if this is an indel)
    pub indel_nucleotides: Option<String>,

    #[pyo3(get, set)]
    /// Gene name this variant lies in. None if not in a gene or in multiple genes
    pub gene_name: Option<String>,

    #[pyo3(get, set)]
    /// Gene position of this variant. None if not in a gene or in multiple genes
    pub gene_position: Option<i64>,

    #[pyo3(get, set)]
    /// Codon index of this variant. None if not in exactly 1 gene and not within a codon
    pub codon_idx: Option<i64>,
}

#[pyclass]
#[derive(Clone, Debug, Eq, PartialEq)]
/// Gene level mutation
pub struct Mutation {
    #[pyo3(get, set)]
    /// GARC for gene level mutation
    pub mutation: String,

    #[pyo3(get, set)]
    /// Gene name
    pub gene: String,

    #[pyo3(get, set)]
    /// Evidence to support this mutation
    pub evidence: Vec<Evidence>,

    #[pyo3(get, set)]
    /// Gene position of this mutation. None if large deletion
    pub gene_position: Option<i64>,

    #[pyo3(get, set)]
    /// Whether this mutation lies within a region which codes protein. None if large deletion
    pub codes_protein: Option<bool>,

    #[pyo3(get, set)]
    /// Reference nucleotides for this mutation. None if not a SNP
    pub ref_nucleotides: Option<String>,

    #[pyo3(get, set)]
    /// Alternate nucleotides for this mutation. None if not a SNP
    pub alt_nucleotides: Option<String>,

    #[pyo3(get, set)]
    /// Nucleotide number for this mutation. None if not referring to a single nucleotide
    pub nucleotide_number: Option<i64>,

    #[pyo3(get, set)]
    /// Nucleotide index for this mutation. None if not referring to a single nucleotide
    pub nucleotide_index: Option<i64>,

    #[pyo3(get, set)]
    /// Length of the indel. None if not an indel
    pub indel_length: Option<i64>,

    #[pyo3(get, set)]
    /// Bases for the indel (if this is an indel). None if not an indel
    pub indel_nucleotides: Option<String>,

    #[pyo3(get, set)]
    /// Amino acid number this mutation refers to. None if not an amino acid SNP
    pub amino_acid_number: Option<i64>,

    #[pyo3(get, set)]
    /// Amino acid alt. None if not an amino acid SNP
    pub amino_acid_sequence: Option<char>,
}

#[pyclass]
#[derive(Debug)]
/// Struct to hold the difference between two genomes
pub struct GenomeDifference {
    #[pyo3(get, set)]
    /// Variants in the genome
    pub variants: Vec<Variant>,

    #[pyo3(get, set)]
    /// Minor variants in the genome
    pub minor_variants: Vec<Variant>,
}

#[pyclass]
#[derive(Debug)]
/// Struct to hold the difference between two genes
pub struct GeneDifference {
    #[pyo3(get, set)]
    /// Mutations in the gene
    pub mutations: Vec<Mutation>,

    #[pyo3(get, set)]
    /// Minor mutations in the gene
    pub minor_mutations: Vec<Mutation>,
}

// Type definitions to appease clippy
type MinorSnp = (
    char,
    Option<i32>,
    Option<OrderedFloat<f32>>,
    Option<Vec<Evidence>>,
);
type NcSnpArgs = (String, i64, bool, char, char, i64, i64, Vec<Alt>);
type MixedIndelArgs = (String, i64, bool, i64, i64, Vec<Alt>, MinorType, String);

#[pymethods]
impl GenomeDifference {
    #[new]
    /// Create a new GenomeDifference object
    ///
    /// # Arguments
    /// - `ref_genome` - Reference genome
    /// - `alt_genome` - Alternate genome
    /// - `minor_type` - Type of minor allele evidence to use in minor variants
    pub fn new(ref_genome: Genome, mut alt_genome: Genome, minor_type: MinorType) -> Self {
        let mut variants: Vec<Variant> = Vec::new();
        let mut minor_variants: Vec<Variant> = Vec::new();

        for idx in 0..ref_genome.genome_positions.len() {
            let ref_pos = &ref_genome.genome_positions[idx];
            let alt_pos = alt_genome.genome_positions[idx].clone();

            if !alt_pos.alts.is_empty() {
                // Alt has a variant at this position, so figure out what it is
                for alt in alt_pos.alts.iter() {
                    let mut garc = "".to_string();
                    let mut indel_bases = None;
                    let mut indel_length = 0;

                    if alt.alt_type == AltType::REF {
                        // Skip ref calls
                        continue;
                    }

                    let mut gene_names = Vec::new();
                    let mut gene_positions = Vec::new();
                    let mut codon_idxs = Vec::new();
                    let mut genes_to_check = HashSet::new();
                    let mut gene_to_genome = HashMap::new();
                    if alt.alt_type == AltType::DEL {
                        // This is a deletion, so check for all genes within each base of this deletion
                        for (i, _) in alt.base.chars().enumerate() {
                            for g in alt_genome.genome_positions[idx + i].genes.clone() {
                                let gene_def = alt_genome.gene_name_to_def.get(&g).unwrap();

                                if !genes_to_check.contains(&g) {
                                    genes_to_check.insert(g.clone());
                                    // Get the genome index of the first base of the deletion
                                    gene_to_genome.insert(g.clone(), (idx + i + 1) as i64);
                                }

                                if gene_def.reverse_complement {
                                    // If this gene is reverse complemented, we want the last base rather than first
                                    gene_to_genome.insert(g, (idx + i + 1) as i64);
                                }
                            }
                        }
                    }
                    for gene in alt_pos.genes.iter() {
                        genes_to_check.insert(gene.clone());
                        gene_to_genome.insert(gene.clone(), ref_pos.genome_idx);
                    }
                    let mut genes = genes_to_check.iter().cloned().collect::<Vec<String>>();
                    genes.sort();
                    for gene in genes.iter() {
                        let gene_name = Some(gene.clone());
                        let g = alt_genome.get_gene(gene.clone());
                        let mut __gene_position = None;
                        let (_gene_position, _codon_idx) = g
                            .genome_idx_map
                            .get(gene_to_genome.get(gene).unwrap())
                            .unwrap();
                        codon_idxs.push(*_codon_idx);
                        if alt.alt_type == AltType::SNP
                            || alt.alt_type == AltType::HET
                            || alt.alt_type == AltType::NULL
                            || _gene_position < &0
                        {
                            // Use gene position for these as it should cover cases of being the codon index
                            __gene_position = Some(*_gene_position);
                        } else {
                            // Use nucleotide number for indels as they shouldn't be construed as codon indices
                            __gene_position = GenomeDifference::get_nucleotide_number(&g, alt);
                        }
                        gene_names.push(gene_name.clone());
                        gene_positions.push(__gene_position);
                    }

                    if alt.alt_type == AltType::SNP
                        || alt.alt_type == AltType::HET
                        || alt.alt_type == AltType::NULL
                    {
                        garc = ref_pos.genome_idx.to_string()
                            + &ref_pos.reference.to_string()
                            + ">"
                            + &alt.base;
                    }
                    if alt.alt_type == AltType::INS {
                        garc = ref_pos.genome_idx.to_string() + "_ins_" + &alt.base;
                        indel_bases = Some(alt.base.clone());
                        indel_length = alt.base.len() as i64;
                    }
                    if alt.alt_type == AltType::DEL {
                        garc = ref_pos.genome_idx.to_string() + "_del_" + &alt.base;
                        indel_bases = Some(alt.base.clone());
                        indel_length = -(alt.base.len() as i64);
                    }

                    if alt.evidence.is_minor {
                        // Append coverage to the variant
                        if minor_type == MinorType::COV && alt.evidence.cov.is_some() {
                            garc = garc + ":" + &alt.evidence.cov.unwrap().to_string();
                        } else if minor_type == MinorType::FRS && alt.evidence.frs.is_some() {
                            garc = garc
                                + ":"
                                + &trim_float_string(format!("{:.3}", alt.evidence.frs.unwrap()));
                        } else {
                            println!("Missing evidence for minor allele: {}, skipping!", garc);
                            continue;
                        }
                    }

                    if gene_names.is_empty() {
                        gene_names.push(None);
                        gene_positions.push(None);
                        codon_idxs.push(None);
                    }

                    for ((gene_name, gene_position), codon_idx) in gene_names
                        .iter()
                        .zip(gene_positions.iter())
                        .zip(codon_idxs.iter())
                    {
                        let variant = Variant {
                            variant: garc.clone(),
                            nucleotide_index: ref_pos.genome_idx,
                            evidence: alt.evidence.vcf_row,
                            vcf_idx: alt.evidence.vcf_idx,
                            indel_length,
                            indel_nucleotides: indel_bases.clone(),
                            gene_position: *gene_position,
                            codon_idx: *codon_idx,
                            gene_name: gene_name.clone(),
                        };

                        if alt.evidence.is_minor {
                            minor_variants.push(variant);
                        } else {
                            variants.push(variant);
                        }
                    }
                }
            }
        }

        GenomeDifference {
            variants,
            minor_variants,
        }
    }

    #[staticmethod]
    /// Find the nucleotide number (at the gene level) of a given alt (at the genome level)
    /// Easy enough to pull out a gene position, but for indels we need to find the nucleotide number
    ///
    /// # Arguments
    /// - `gene` - Gene to search
    /// - `genome_alt` - Alt to search for
    ///
    /// # Returns
    /// Nucleotide number of the alt in the gene. None if not found
    pub fn get_nucleotide_number(gene: &Gene, genome_alt: &Alt) -> Option<i64> {
        for gene_pos in gene.gene_positions.iter() {
            match &gene_pos.gene_position_data {
                GenePos::Codon(c) => {
                    for codon in c.codon.iter() {
                        for alt in codon.alts.iter() {
                            // Match on evidence's VCF row and index
                            // This is to catch cases where a deletion starts outside of this gene
                            if genome_alt.evidence.vcf_row == alt.evidence.vcf_row
                                && genome_alt.evidence.vcf_idx == alt.evidence.vcf_idx
                                && genome_alt.alt_type == alt.alt_type
                                && genome_alt.evidence.is_minor == alt.evidence.is_minor
                            {
                                return Some(codon.nucleotide_number);
                            }
                        }
                    }
                }
                GenePos::Nucleotide(n) => {
                    for alt in n.alts.iter() {
                        // Match on evidence's VCF row and index
                        // This is to catch cases where a deletion starts outside of this gene
                        if genome_alt.evidence.vcf_row == alt.evidence.vcf_row
                            && genome_alt.evidence.vcf_idx == alt.evidence.vcf_idx
                            && genome_alt.alt_type == alt.alt_type
                            && genome_alt.evidence.is_minor == alt.evidence.is_minor
                        {
                            return Some(n.nucleotide_number);
                        }
                    }
                }
            }
        }
        // I hate implicit returns, but appease clippy
        None
    }
}

#[pymethods]
impl GeneDifference {
    #[new]
    /// Create a new GeneDifference object
    ///
    /// # Arguments
    /// - `ref_gene` - Reference gene
    /// - `alt_gene` - Alternate gene
    /// - `minor_type` - Type of minor allele evidence to use in minor variants
    pub fn new(ref_gene: Gene, alt_gene: Gene, minor_type: MinorType) -> Self {
        if ref_gene.name != alt_gene.name {
            panic!("Gene names do not match!");
        }
        let mut mutations = Vec::new();
        let mut minor_mutations = Vec::new();
        let mut deleted_bases = 0;
        let mut minor_deleted_bases = 0;
        let mut minor_deleted_cov = 0;
        let mut minor_deleted_frs = OrderedFloat(0.0);
        let gene_name = ref_gene.name.clone();

        for idx in 0..ref_gene.gene_positions.len() {
            let ref_pos = &ref_gene.gene_positions[idx];
            let alt_pos = &alt_gene.gene_positions[idx];

            if ref_pos != alt_pos {
                // We have a difference so unpack it
                let gene_position = ref_pos.gene_position;
                let codes_protein = ref_gene.coding && ref_pos.gene_position > 0;
                match &alt_pos.gene_position_data {
                    GenePos::Codon(alt_codon) => {
                        // Unpack the ref codon too
                        let mut _maybe_ref_codon = None;
                        match &ref_pos.gene_position_data {
                            GenePos::Codon(x) => {
                                _maybe_ref_codon = Some(x);
                            }
                            _ => panic!("Reference gene position is not a codon"),
                        }
                        let ref_codon = _maybe_ref_codon.unwrap();

                        let mut _ref_nucleotides = None;
                        let mut _alt_nucleotides = None;
                        let mut nucleotide_number = None;
                        let mut nucleotide_index = None;
                        let mut indel_length = None;
                        let mut indel_nucleotides = None;
                        let mut _amino_acid_number = None;
                        let mut _amino_acid_sequence = None;

                        let mut _mutation = "".to_string();
                        let mut evidence: Vec<Evidence> = Vec::new();

                        // Codons are annoying as they can be split across multiple nucleotides
                        // So let's first check for a SNP as that's simple
                        let ref_codon_str = ref_codon
                            .codon
                            .iter()
                            .map(|x| x.reference.to_string())
                            .collect::<String>();
                        _ref_nucleotides = Some(ref_codon_str.clone());
                        _alt_nucleotides = Some(
                            alt_codon
                                .codon
                                .iter()
                                .map(|x| x.reference.to_string())
                                .collect::<String>(),
                        );
                        if alt_codon.amino_acid != ref_codon.amino_acid {
                            // SNP at this amino acid so mutation is simple
                            _mutation = ref_codon.amino_acid.to_string()
                                + &ref_pos.gene_position.to_string()
                                + &alt_codon.amino_acid.to_string();
                            _amino_acid_number = Some(gene_position);
                            _amino_acid_sequence = Some(alt_codon.amino_acid);
                            // Filter evidence only for SNP
                            for ev in alt_codon.codon.clone() {
                                if !ev.alts.is_empty() {
                                    for e in ev.alts.iter() {
                                        if !e.evidence.is_minor
                                            && (e.alt_type == AltType::SNP
                                                || e.alt_type == AltType::HET
                                                || e.alt_type == AltType::NULL)
                                        {
                                            evidence.push(e.evidence.clone());
                                        }
                                    }
                                }
                            }
                            mutations.push(Mutation {
                                mutation: _mutation,
                                gene: gene_name.clone(),
                                evidence,
                                gene_position: Some(gene_position),
                                codes_protein: Some(codes_protein),
                                ref_nucleotides: _ref_nucleotides,
                                alt_nucleotides: _alt_nucleotides,
                                nucleotide_number,
                                nucleotide_index,
                                indel_length,
                                indel_nucleotides: indel_nucleotides.clone(),
                                amino_acid_number: _amino_acid_number,
                                amino_acid_sequence: _amino_acid_sequence,
                            });
                        } else {
                            // Pull out nucleotide variants in cases of synonymous mutations
                            let mut synon_snp = false;
                            for (ref_cd, alt_cd) in
                                ref_codon.codon.iter().zip(alt_codon.codon.clone())
                            {
                                if ref_cd.reference != alt_cd.reference {
                                    synon_snp = true;
                                    mutations.push(GeneDifference::nc_snp((
                                        gene_name.clone(),
                                        alt_cd.nucleotide_number,
                                        codes_protein,
                                        ref_cd.reference,
                                        alt_cd.reference,
                                        ref_cd.nucleotide_number,
                                        ref_cd.nucleotide_index,
                                        alt_cd
                                            .alts
                                            .iter()
                                            .filter(|x| !x.evidence.is_minor)
                                            .map(|x| (*x).clone())
                                            .collect::<Vec<Alt>>(),
                                    )))
                                }
                            }
                            if synon_snp {
                                _mutation = ref_codon.amino_acid.to_string()
                                    + &ref_pos.gene_position.to_string()
                                    + &alt_codon.amino_acid.to_string();
                                _amino_acid_number = Some(gene_position);
                                _amino_acid_sequence = Some(alt_codon.amino_acid);
                                // Filter evidence only for SNP
                                for ev in alt_codon.codon.clone() {
                                    if !ev.alts.is_empty() {
                                        for e in ev.alts.iter() {
                                            if !e.evidence.is_minor && e.alt_type == AltType::SNP {
                                                evidence.push(e.evidence.clone());
                                            }
                                        }
                                    }
                                }
                                mutations.push(Mutation {
                                    mutation: _mutation,
                                    gene: gene_name.clone(),
                                    evidence,
                                    gene_position: Some(gene_position),
                                    codes_protein: Some(codes_protein),
                                    ref_nucleotides: _ref_nucleotides,
                                    alt_nucleotides: _alt_nucleotides,
                                    nucleotide_number,
                                    nucleotide_index,
                                    indel_length,
                                    indel_nucleotides: indel_nucleotides.clone(),
                                    amino_acid_number: _amino_acid_number,
                                    amino_acid_sequence: _amino_acid_sequence,
                                });
                            }
                        }

                        // Then, check for indels and minor mutations in a codon

                        // Codons being fun mean we need to iter the codon and check for minor mutations
                        // An overall minor amino acid can be found, but the individual minor mutations need to be checked
                        // In cases of >1 minor SNP at a nucleotide, treat it as a het call with the minimum coverage
                        // Vec of (alt, cov, frs, evidence).
                        // If no minor mutation at nucleotide, give ref nucleotide and None for other values
                        let mut minor_snps: Vec<MinorSnp> = Vec::new();
                        let mut minor_snp_exists = false;
                        for (ref_cd, alt_cd) in ref_codon.codon.iter().zip(alt_codon.codon.clone())
                        {
                            if alt_cd.is_deleted {
                                deleted_bases += 1;
                            }
                            if alt_cd.is_deleted_minor {
                                minor_deleted_bases += 1;
                            }
                            let mut these_minor_snps = Vec::new();
                            let mut these_minor_indels = Vec::new();

                            // Make sure these are empty for each nc in the codon
                            _ref_nucleotides = None;
                            _alt_nucleotides = None;
                            nucleotide_number = None;
                            nucleotide_index = None;
                            indel_length = None;
                            indel_nucleotides = None;
                            _amino_acid_number = None;
                            _amino_acid_sequence = None;
                            _mutation = "".to_string();
                            evidence = Vec::new();
                            for e in alt_cd.alts.iter() {
                                if e.evidence.is_minor {
                                    if e.alt_type == AltType::REF {
                                        //Not sure why we get minor ref calls but skip
                                        continue;
                                    }
                                    if e.alt_type == AltType::SNP
                                        || e.alt_type == AltType::HET
                                        || e.alt_type == AltType::NULL
                                    {
                                        these_minor_snps.push(e);
                                        minor_snp_exists = true;
                                    } else {
                                        if e.alt_type == AltType::DEL {
                                            if e.evidence.cov.unwrap() > minor_deleted_cov {
                                                minor_deleted_cov = e.evidence.cov.unwrap();
                                            }
                                            if e.evidence.frs.unwrap() > minor_deleted_frs {
                                                minor_deleted_frs = e.evidence.frs.unwrap();
                                            }
                                        }

                                        these_minor_indels.push(e);
                                    }
                                } else {
                                    if e.alt_type == AltType::INS {
                                        _mutation = alt_cd.nucleotide_number.to_string()
                                            + "_ins_"
                                            + &e.base;
                                        indel_length = Some(e.base.len() as i64);
                                        indel_nucleotides = Some(e.base.clone());
                                        evidence = vec![e.evidence.clone()];
                                    }
                                    if e.alt_type == AltType::DEL {
                                        _mutation = alt_cd.nucleotide_number.to_string()
                                            + "_del_"
                                            + &e.base;
                                        indel_length = Some(-(e.base.len() as i64));
                                        indel_nucleotides = Some(e.base.clone());
                                        evidence = vec![e.evidence.clone()];
                                    }
                                }
                            }
                            if _mutation != *"" {
                                // We picked up a mutation so lets append it
                                mutations.push(Mutation {
                                    mutation: _mutation.clone(),
                                    gene: gene_name.clone(),
                                    evidence: evidence.clone(),
                                    gene_position: Some(alt_cd.nucleotide_number),
                                    codes_protein: Some(codes_protein),
                                    ref_nucleotides: _ref_nucleotides.clone(),
                                    alt_nucleotides: _alt_nucleotides.clone(),
                                    nucleotide_number,
                                    nucleotide_index,
                                    indel_length,
                                    indel_nucleotides: indel_nucleotides.clone(),
                                    amino_acid_number: _amino_acid_number,
                                    amino_acid_sequence: _amino_acid_sequence,
                                });
                            }
                            if !these_minor_indels.is_empty() && !these_minor_snps.is_empty() {
                                // Mix of indel and SNP at this position
                                let mut these_minors = these_minor_indels.clone();
                                for snp in these_minor_snps.iter() {
                                    these_minors.push(snp);
                                }
                                minor_mutations.push(GeneDifference::mixed_indel((
                                    gene_name.clone(),
                                    alt_cd.nucleotide_number,
                                    codes_protein,
                                    alt_cd.nucleotide_number,
                                    alt_cd.nucleotide_index,
                                    these_minors
                                        .iter()
                                        .map(|x| (*x).clone())
                                        .collect::<Vec<Alt>>(),
                                    minor_type.clone(),
                                    "mixed".to_string(),
                                )));
                            } else {
                                if !these_minor_indels.is_empty() {
                                    // Minor indel
                                    if these_minor_indels.len() > 1 {
                                        // We have a mixed minor indel, so treat it as such
                                        minor_mutations.push(GeneDifference::mixed_indel((
                                            gene_name.clone(),
                                            alt_cd.nucleotide_number,
                                            codes_protein,
                                            alt_cd.nucleotide_number,
                                            alt_cd.nucleotide_index,
                                            these_minor_indels
                                                .iter()
                                                .map(|x| (*x).clone())
                                                .collect::<Vec<Alt>>(),
                                            minor_type.clone(),
                                            "indel".to_string(),
                                        )));
                                    } else {
                                        let e = these_minor_indels[0];
                                        // We have a single minor indel which is much easier
                                        if e.alt_type == AltType::INS {
                                            _mutation = alt_cd.nucleotide_number.to_string()
                                                + "_ins_"
                                                + &e.base;
                                            indel_length = Some(e.base.len() as i64);
                                            indel_nucleotides = Some(e.base.clone());
                                            evidence = vec![e.evidence.clone()];
                                        }
                                        if e.alt_type == AltType::DEL {
                                            _mutation = alt_cd.nucleotide_number.to_string()
                                                + "_del_"
                                                + &e.base;
                                            indel_length = Some(-(e.base.len() as i64));
                                            indel_nucleotides = Some(e.base.clone());
                                            evidence = vec![e.evidence.clone()];
                                        }
                                        if minor_type == MinorType::COV && e.evidence.cov.is_some()
                                        {
                                            _mutation = _mutation
                                                + ":"
                                                + &e.evidence.cov.unwrap().to_string();
                                        } else if minor_type == MinorType::FRS
                                            && e.evidence.frs.is_some()
                                        {
                                            _mutation = _mutation
                                                + ":"
                                                + &trim_float_string(format!(
                                                    "{:.3}",
                                                    e.evidence.frs.unwrap()
                                                ));
                                        } else {
                                            println!(
                                                "Missing evidence for minor allele: {}, skipping!",
                                                _mutation
                                            );
                                            continue;
                                        }
                                        minor_mutations.push(Mutation {
                                            mutation: _mutation.clone(),
                                            gene: gene_name.clone(),
                                            evidence: evidence.clone(),
                                            gene_position: Some(alt_cd.nucleotide_number),
                                            codes_protein: Some(codes_protein),
                                            ref_nucleotides: None,
                                            alt_nucleotides: None,
                                            nucleotide_number,
                                            nucleotide_index,
                                            indel_length,
                                            indel_nucleotides: indel_nucleotides.clone(),
                                            amino_acid_number: None,
                                            amino_acid_sequence: None,
                                        });
                                    }
                                }
                                if !these_minor_snps.is_empty() {
                                    // Minor SNP
                                    if these_minor_snps.len() > 1 {
                                        // Mixed minor SNP
                                        let min_cov = these_minor_snps
                                            .iter()
                                            .filter_map(|x| x.evidence.cov)
                                            .max()
                                            .unwrap();
                                        let min_frs = these_minor_snps
                                            .iter()
                                            .filter_map(|x| x.evidence.frs)
                                            .max()
                                            .unwrap();
                                        minor_snps.push((
                                            'z',
                                            Some(min_cov),
                                            Some(min_frs),
                                            Some(
                                                these_minor_snps
                                                    .iter()
                                                    .map(|x| x.evidence.clone())
                                                    .collect(),
                                            ),
                                        ));
                                    } else {
                                        // Single minor SNP
                                        minor_snps.push((
                                            these_minor_snps[0].base.chars().nth(0).unwrap(),
                                            these_minor_snps[0].evidence.cov,
                                            these_minor_snps[0].evidence.frs,
                                            Some(vec![these_minor_snps[0].evidence.clone()]),
                                        ));
                                    }
                                } else {
                                    // No minor SNP at this position
                                    minor_snps.push((ref_cd.reference, None, None, None));
                                }
                            }
                        }
                        if minor_snp_exists && minor_snps.len() == 3 {
                            // Construct the minor amino acid change
                            let mut codon = "".to_string();
                            let mut minor_cov = 0;
                            let mut minor_frs = OrderedFloat::min_value();
                            let mut minor_evidence = Vec::new();
                            for (nc, cov, frs, ev) in minor_snps.iter() {
                                codon += &nc.to_string();
                                // Evidence is the only field which is None in the case of no minor SNP at this nc
                                if let Some(e) = ev {
                                    for e in e.iter() {
                                        minor_evidence.push(e.clone());
                                    }
                                    if cov.is_some() && cov.unwrap() > minor_cov {
                                        minor_cov = cov.unwrap();
                                    }
                                    if frs.is_some() && frs.unwrap() > minor_frs {
                                        minor_frs = frs.unwrap();
                                    }
                                }
                            }
                            let aa = codon_to_aa(codon.clone());
                            let mut mutation = "".to_string();
                            if minor_type == MinorType::COV {
                                mutation = ref_codon.amino_acid.to_string()
                                    + &ref_pos.gene_position.to_string()
                                    + &aa.to_string()
                                    + ":"
                                    + &minor_cov.to_string();
                            }
                            if minor_type == MinorType::FRS {
                                mutation = ref_codon.amino_acid.to_string()
                                    + &ref_pos.gene_position.to_string()
                                    + &aa.to_string()
                                    + ":"
                                    + &trim_float_string(format!("{:.3}", minor_frs));
                            }
                            minor_mutations.push(Mutation {
                                mutation,
                                gene: gene_name.clone(),
                                evidence: minor_evidence,
                                gene_position: Some(gene_position),
                                codes_protein: Some(codes_protein),
                                ref_nucleotides: Some(ref_codon_str.clone()),
                                alt_nucleotides: Some(codon.clone()),
                                nucleotide_number: None,
                                nucleotide_index: None,
                                indel_length: None,
                                indel_nucleotides: None,
                                amino_acid_number: Some(gene_position),
                                amino_acid_sequence: Some(aa),
                            });
                        }
                    }
                    GenePos::Nucleotide(alt_nc) => {
                        // Unpack the ref nc too
                        let mut _maybe_ref_nc = None;
                        match &ref_pos.gene_position_data {
                            GenePos::Nucleotide(x) => {
                                _maybe_ref_nc = Some(x);
                            }
                            _ => panic!("Reference gene position is not a nucleotide"),
                        }
                        let ref_nc = _maybe_ref_nc.unwrap();

                        // Nucleotide mutations are much more simple
                        for alt in alt_nc.alts.iter() {
                            let mut mutation = "".to_string();
                            let evidence: Vec<Evidence> = vec![alt.evidence.clone()];
                            let mut ref_nucleotides = None;
                            let mut alt_nucleotides = None;
                            let nucleotide_number = Some(ref_nc.nucleotide_number);
                            let nucleotide_index = Some(ref_nc.nucleotide_index);
                            let mut indel_length: Option<i64> = None;
                            let mut indel_nucleotides = None;
                            let amino_acid_number = None;
                            let amino_acid_sequence = None;

                            if alt.alt_type == AltType::REF {
                                // Skip ref calls
                                continue;
                            }

                            if alt.alt_type == AltType::SNP
                                || alt.alt_type == AltType::HET
                                || alt.alt_type == AltType::NULL
                            {
                                mutation = ref_nc.reference.to_string()
                                    + &ref_nc.nucleotide_number.to_string()
                                    + &alt.base;
                                ref_nucleotides = Some(ref_nc.reference.to_string());
                                alt_nucleotides = Some(alt.base.to_string());
                            }
                            if alt.alt_type == AltType::INS {
                                mutation =
                                    ref_nc.nucleotide_number.to_string() + "_ins_" + &alt.base;
                                indel_length = Some(alt.base.len() as i64);
                                indel_nucleotides = Some(alt.base.clone());
                            }
                            if alt.alt_type == AltType::DEL {
                                mutation =
                                    ref_nc.nucleotide_number.to_string() + "_del_" + &alt.base;
                                indel_length = Some(-(alt.base.len() as i64));
                                indel_nucleotides = Some(alt.base.clone());
                            }

                            if alt.evidence.is_minor {
                                // Append coverage to the variant
                                if minor_type == MinorType::COV && alt.evidence.cov.is_some() {
                                    mutation =
                                        mutation + ":" + &alt.evidence.cov.unwrap().to_string();
                                } else if minor_type == MinorType::FRS && alt.evidence.frs.is_some()
                                {
                                    mutation = mutation
                                        + ":"
                                        + &trim_float_string(format!(
                                            "{:.3}",
                                            alt.evidence.frs.unwrap()
                                        ));
                                }
                                if alt.alt_type == AltType::DEL {
                                    if alt.evidence.cov.is_some()
                                        && alt.evidence.cov.unwrap() > minor_deleted_cov
                                    {
                                        minor_deleted_cov = alt.evidence.cov.unwrap();
                                    }
                                    if alt.evidence.frs.is_some()
                                        && alt.evidence.frs.unwrap() > minor_deleted_frs
                                    {
                                        minor_deleted_frs = alt.evidence.frs.unwrap();
                                    }
                                }
                                if alt.alt_type == AltType::DEL || alt.alt_type == AltType::INS {
                                    minor_deleted_bases += 1;
                                }
                            } else if alt.alt_type == AltType::DEL || alt.alt_type == AltType::INS {
                                deleted_bases += 1;
                            }
                            let m = Mutation {
                                mutation,
                                gene: gene_name.clone(),
                                evidence,
                                gene_position: Some(gene_position),
                                codes_protein: Some(codes_protein),
                                ref_nucleotides,
                                alt_nucleotides,
                                nucleotide_number,
                                nucleotide_index,
                                indel_length,
                                indel_nucleotides,
                                amino_acid_number,
                                amino_acid_sequence,
                            };

                            if alt.evidence.is_minor {
                                minor_mutations.push(m);
                            } else {
                                mutations.push(m);
                            }
                        }
                    }
                }
            }
        }

        let deleted_percent = deleted_bases as f64 / ref_gene.nucleotide_number.len() as f64;
        let minor_deleted_percent =
            minor_deleted_bases as f64 / ref_gene.nucleotide_number.len() as f64;
        if deleted_percent >= 0.5 {
            mutations.push(Mutation {
                mutation: "del_".to_string()
                    + trim_float_string(format!("{:.2}", deleted_percent)).as_str(),
                gene: gene_name.clone(),
                evidence: Vec::new(),
                gene_position: None,
                codes_protein: None,
                ref_nucleotides: None,
                alt_nucleotides: None,
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: None,
                amino_acid_sequence: None,
            });
        }
        if minor_deleted_percent >= 0.5 {
            let mut mutation = "del_".to_string()
                + trim_float_string(format!("{:.2}", minor_deleted_percent)).as_str();
            if minor_type == MinorType::COV {
                mutation = mutation + ":" + &minor_deleted_cov.to_string();
            }
            if minor_type == MinorType::FRS {
                mutation = mutation + ":" + &trim_float_string(format!("{:.3}", minor_deleted_frs));
            }
            minor_mutations.push(Mutation {
                mutation,
                gene: gene_name.clone(),
                evidence: Vec::new(),
                gene_position: None,
                codes_protein: None,
                ref_nucleotides: None,
                alt_nucleotides: None,
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: None,
                amino_acid_sequence: None,
            });
        }

        // I hate implicit returns, but appease clippy
        GeneDifference {
            mutations,
            minor_mutations,
        }
    }

    #[staticmethod]
    /// Create a new Mutation object for a nucleotide SNP
    ///
    /// # Arguments (as NcSnpArgs tuple)
    /// - `gene_name` - Name of the gene
    /// - `gene_position` - Position of the SNP in the gene
    /// - `codes_protein` - Whether this SNP codes protein
    /// - `ref_nc` - Reference nucleotide
    /// - `alt_nc` - Alternate nucleotide
    /// - `nc_num` - Nucleotide number
    /// - `nc_idx` - Nucleotide index
    /// - `evidence` - Evidence for this SNP
    ///
    /// # Returns
    /// - Mutation of the SNP
    fn nc_snp(args: NcSnpArgs) -> Mutation {
        let (gene_name, gene_position, codes_protein, ref_nc, alt_nc, nc_num, nc_idx, evidence) =
            args;
        let mutation = ref_nc.to_string() + &nc_num.to_string() + &alt_nc.to_string();
        let ref_nucleotides = Some(ref_nc.to_string());
        let alt_nucleotides = Some(alt_nc.to_string());
        let nucleotide_number = Some(nc_num);
        let nucleotide_index = Some(nc_idx);
        let evidence = evidence
            .iter()
            .filter(|x| x.alt_type == AltType::SNP)
            .map(|x| x.evidence.clone())
            .collect();

        // I hate implicit returns, but appease clippy
        Mutation {
            mutation,
            gene: gene_name.clone(),
            evidence,
            gene_position: Some(gene_position),
            codes_protein: Some(codes_protein),
            ref_nucleotides,
            alt_nucleotides,
            nucleotide_number,
            nucleotide_index,
            indel_length: None,
            indel_nucleotides: None,
            amino_acid_number: None,
            amino_acid_sequence: None,
        }
    }

    #[staticmethod]
    /// Create a new Mutation object for a mixed indel. This could either be >1 indel at this position, or >=1 indel and >=1 SNP
    ///
    /// # Arguments (as MixedIndelArgs tuple)
    /// - `gene_name` - Name of the gene
    /// - `gene_position` - Position of the indel in the gene
    /// - `codes_protein` - Whether this indel codes protein
    /// - `nc_num` - Nucleotide number
    /// - `nc_idx` - Nucleotide index
    /// - `these_minors` - Minor alleles at this position
    /// - `minor_type` - Type of minor allele evidence to use
    /// - `mutation_name` - 'indel' or 'mixed' depending on the type of mixed indel
    ///
    /// # Returns
    /// - Mutation object for the mixed indel
    fn mixed_indel(args: MixedIndelArgs) -> Mutation {
        let (
            gene_name,
            gene_position,
            codes_protein,
            nc_num,
            nc_idx,
            these_minors,
            minor_type,
            mutation_name,
        ) = args;
        let mut min_coverage = "".to_string();
        if minor_type == MinorType::COV {
            min_coverage = these_minors
                .iter()
                .filter_map(|x| x.evidence.cov)
                .max()
                .unwrap()
                .to_string();
        }
        if minor_type == MinorType::FRS {
            min_coverage = trim_float_string(format!(
                "{:.3}",
                these_minors
                    .iter()
                    .filter_map(|x| x.evidence.frs)
                    .max()
                    .unwrap()
            ));
        }
        let mutation = nc_num.to_string() + "_" + &mutation_name + ":" + &min_coverage;
        let evidence = these_minors.iter().map(|x| x.evidence.clone()).collect();

        let ref_nucleotides = None;
        let alt_nucleotides = None;
        let nucleotide_number = Some(nc_num);
        let nucleotide_index = Some(nc_idx);

        // I hate implicit returns, but appease clippy
        Mutation {
            mutation,
            gene: gene_name.clone(),
            evidence,
            gene_position: Some(gene_position),
            codes_protein: Some(codes_protein),
            ref_nucleotides,
            alt_nucleotides,
            nucleotide_number,
            nucleotide_index,
            indel_length: None,
            indel_nucleotides: None,
            amino_acid_number: None,
            amino_acid_sequence: None,
        }
    }
}

#[pyfunction]
/// Given a string of a float, trim trailling `0` chars
///
/// # Arguments
/// - `float_string` - String representation of a float
///
/// # Returns
/// - Trimmed string
fn trim_float_string(mut float_string: String) -> String {
    while float_string.ends_with('0') {
        float_string.pop();
        if float_string.ends_with("1.0") {
            // Keep trailing 0 for 1.0
            return float_string;
        }
    }
    // I hate implicit returns, but appease clippy
    float_string
}
