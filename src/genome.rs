//! Module for handling genome data
use pyo3::prelude::*;

use std::collections::{HashMap, HashSet};
use std::fs::File;

use gb_io::reader::SeqReader;
use gb_io::seq::Location::Complement;
use gb_io::seq::Location::Join;
use gb_io::seq::Location::Range;

use crate::common::{Alt, AltType, Evidence, GeneDef, VCFRow};
use crate::gene::Gene;
use crate::vcf::VCFFile;

#[pyclass]
#[derive(Clone, Debug, Eq, PartialEq)]
/// Struct to hold the information of a genome position
pub struct GenomePosition {
    // Updated during mutation
    #[pyo3(get, set)]
    /// Nucleotide at this position
    pub reference: char,

    #[pyo3(get, set)]
    /// Whether this position has been deleted
    pub is_deleted: bool,

    #[pyo3(get, set)]
    /// Whether this position has been deleted in a minor allele
    pub is_deleted_minor: bool,

    #[pyo3(get, set)]
    /// Added for evidence of a deletion which didn't start at this position but affects it
    /// Includes evidence of minor deletions for simplicity
    pub deleted_evidence: Vec<Evidence>,

    #[pyo3(get, set)]
    /// Used to store calls from VCF
    pub alts: Vec<Alt>,

    #[pyo3(get, set)]
    /// 1-indexed genome index
    pub genome_idx: i64,

    #[pyo3(get, set)]
    /// Names of genes present at this position
    pub genes: Vec<String>,
}

#[pyclass]
#[derive(Clone)]
/// Struct to hold the information about a genome
pub struct Genome {
    #[pyo3(get, set)]
    /// Name of the genome
    pub name: String,

    #[pyo3(get, set)]
    /// Nucleotide sequence which comprises this genome
    pub nucleotide_sequence: String,

    #[pyo3(get, set)]
    /// Definitions for each gene in the genome
    pub gene_definitions: Vec<GeneDef>,

    #[pyo3(get, set)]
    /// Positions in the genome
    pub genome_positions: Vec<GenomePosition>,

    #[pyo3(get, set)]
    /// Names of all of the genes with definitions
    pub gene_names: Vec<String>,

    #[pyo3(get, set)]
    /// HashMap of gene names to definitions
    pub gene_name_to_def: HashMap<String, GeneDef>,

    #[pyo3(get, set)]
    /// HashMap of genes which have been built
    pub genes: HashMap<String, Gene>,

    #[pyo3(get, set)]
    /// Set of genes with mutations
    pub genes_with_mutations: HashSet<String>,

    #[pyo3(get, set)]
    /// List of the VCF records associated with this genome (if any)
    pub vcf_records: Option<Vec<VCFRow>>,
}

#[pymethods]
impl Genome {
    #[new]
    /// Create a new genome from a GenBank file
    ///
    /// # Arguments
    /// - `filename` - A string reference to the GenBank file
    pub fn new(filename: &str) -> Self {
        let file = File::open(filename).unwrap();
        let mut _gene_definitions = Vec::new();
        let mut _nucleotide_sequence: String = "".to_string();
        let mut genome_name: String = "".to_string();
        let mut gene_names: Vec<String> = Vec::new();
        for seq in SeqReader::new(file) {
            let seq = seq.unwrap();
            _nucleotide_sequence = match String::from_utf8(seq.seq) {
                Ok(s) => s.to_lowercase(),
                Err(e) => panic!("Problem reading sequence data: {:?}", e),
            };
            genome_name = seq.name.unwrap();
            for feature in seq.features {
                let mut name: String = "".to_string();
                let start: i64;
                let end: i64;
                let mut coding: bool = false;
                let mut reverse_complement: bool = false;
                let mut ribosomal_shifts: Vec<i64> = Vec::new();
                if feature.kind == *"CDS" || feature.kind == *"rRNA" {
                    if feature.kind != *"rRNA" {
                        coding = true;
                    }
                    match feature.location {
                        Complement(x) => {
                            reverse_complement = true;
                            match *x {
                                Range(s, e) => {
                                    start = e.0;
                                    end = s.0;
                                }
                                _ => panic!("Complement location not range"),
                            }
                        }
                        Range(s, e) => {
                            start = s.0;
                            end = e.0
                        }
                        Join(ranges) => {
                            // Checking for PRFS
                            let mut start_pos: i64 = 0;
                            let mut end_pos = 0;
                            let mut first = true;
                            for range in ranges {
                                match range {
                                    Range(s, e) => {
                                        if first {
                                            start_pos = s.0;
                                            first = false;
                                        } else {
                                            // Adjust the position to match the 1-indexed genome
                                            ribosomal_shifts.push(s.0 + 1);
                                        }
                                        end_pos = e.0;
                                    }
                                    _ => panic!("Join location not range"),
                                }
                            }
                            start = start_pos;
                            end = end_pos;
                        }
                        _ => panic!("Location not range, complement or join"),
                    }
                    for qual in feature.qualifiers {
                        match qual.1 {
                            Some(val) => {
                                if qual.0 == *"gene" {
                                    // Use gene name if it exists
                                    name = val.clone();
                                }
                                if name.is_empty() && qual.0 == *"locus_tag" {
                                    // Else default to locus tag
                                    name = val.clone();
                                }
                            }
                            None => continue,
                        }
                    }
                    while gene_names.contains(&name) {
                        // Duplicate gene names can exist :(
                        // Repeatedly add _2 to the end of the name until it is unique
                        // not ideal but exact mirror of gumpy
                        name += "_2";
                    }
                    gene_names.push(name.clone());
                    _gene_definitions.push(GeneDef {
                        name,
                        reverse_complement,
                        coding,
                        start,
                        end,
                        promoter_start: -1,
                        promoter_size: 0,
                        ribosomal_shifts,
                    });
                }
            }
        }
        let mut genome_positions = Vec::new();
        for (genome_idx, c) in _nucleotide_sequence.chars().enumerate() {
            genome_positions.push(GenomePosition {
                reference: c,
                genome_idx: genome_idx as i64 + 1,
                alts: Vec::new(),
                genes: Vec::new(),
                is_deleted: false,
                is_deleted_minor: false,
                deleted_evidence: Vec::new(),
            });
        }
        let mut genome = Genome {
            name: genome_name,
            nucleotide_sequence: _nucleotide_sequence,
            gene_definitions: _gene_definitions,
            genome_positions,
            genes: HashMap::new(),
            gene_names,
            genes_with_mutations: HashSet::new(),
            gene_name_to_def: HashMap::new(),
            vcf_records: None,
        };
        genome.assign_promoters();

        // I hate implicit returns, but appease clippy
        genome
    }

    /// Assign promoters to genes
    fn assign_promoters(&mut self) {
        /*
        Assigns promoters to genes iteratively.
        Expand out to a given distance from the start of each gene without overlapping with other genes.
        */
        let max_promoter_length = 99;
        for gene in self.gene_definitions.iter_mut() {
            // First pass to add gene names to all positions they exist in
            let mut start_idx = gene.start;
            let mut end_idx = gene.end;
            if gene.reverse_complement {
                start_idx = gene.end;
                end_idx = gene.start;
            }
            for i in start_idx..end_idx {
                self.genome_positions[i as usize]
                    .genes
                    .push(gene.name.clone());
            }
        }
        for gene in self.gene_definitions.iter_mut() {
            // Check for overlapping genes, assigning promoters to non-overlapping genes
            if self.genome_positions[gene.start as usize].genes.len() > 1 {
                continue;
            } else if gene.start == 0 {
                // Catch edge case of gene starting at genome index 0
                // this couldn't have a promoter so mark as such
                gene.promoter_start = -1;
            } else if gene.reverse_complement {
                gene.promoter_start = gene.start - 1;
            } else {
                gene.promoter_start = gene.start;
            }
        }

        let mut complete = false;
        while !complete {
            let mut this_complete = true;
            for gene in self.gene_definitions.iter_mut() {
                let mut expanding = -1;
                if gene.reverse_complement {
                    // This pushes `gene.promoter_start += expanding` the right way for reverse_complement
                    expanding = 1;
                }

                if gene.promoter_start == -1
                    || gene.promoter_start == 0
                    || gene.promoter_size == max_promoter_length
                    || (gene.promoter_start + expanding) >= self.genome_positions.len() as i64
                {
                    continue;
                }
                if !self.genome_positions[(gene.promoter_start + expanding) as usize]
                    .genes
                    .is_empty()
                {
                    // Pre-existing gene at this position so ignore expanding!
                    continue;
                } else {
                    // No gene in the new position so expand
                    gene.promoter_start += expanding;
                    gene.promoter_size += 1;
                    self.genome_positions[gene.promoter_start as usize]
                        .genes
                        .push(gene.name.clone());
                    this_complete = false;
                }
            }
            complete = this_complete;
        }

        // Add gene definitions to hashmap for easy lookup
        for gene in self.gene_definitions.iter() {
            self.gene_name_to_def
                .insert(gene.name.clone(), gene.clone());
        }
    }

    /// Build a gene from the genome
    ///
    /// # Arguments
    /// - `gene_name` - Name of the gene to build
    ///
    /// # Returns
    /// Corresponding Gene object
    pub fn build_gene(&self, gene_name: String) -> Gene {
        let mut valid = false;
        let mut maybe_gene_def: Option<GeneDef> = None;
        for gene in self.gene_definitions.iter() {
            if gene.name == gene_name {
                valid = true;
                maybe_gene_def = Some(gene.clone());
                break;
            }
        }
        if !valid {
            panic!("Gene {} not found in genome {}", gene_name, self.name);
        }

        let gene_def = maybe_gene_def.unwrap();
        let mut nucleotide_sequence = "".to_string();
        let mut nucleotide_index = Vec::new();
        let mut genome_positions: Vec<GenomePosition> = Vec::new();
        if gene_def.reverse_complement {
            let mut last_idx = gene_def.promoter_start;
            if gene_def.promoter_start == -1 {
                last_idx = gene_def.start;
            }
            for i in gene_def.end..last_idx + 1 {
                nucleotide_sequence.push(self.genome_positions[i as usize].reference);
                nucleotide_index.push(self.genome_positions[i as usize].genome_idx);
                genome_positions.push(self.genome_positions[i as usize].clone());
            }
        } else {
            let mut first_idx = gene_def.promoter_start - 1;
            if gene_def.promoter_start == -1 {
                first_idx = gene_def.start;
            }
            if first_idx == -1 {
                first_idx = 0;
            }
            for i in first_idx..gene_def.end {
                nucleotide_sequence.push(self.genome_positions[i as usize].reference);
                nucleotide_index.push(self.genome_positions[i as usize].genome_idx);
                genome_positions.push(self.genome_positions[i as usize].clone());
            }
        }

        // I hate implicit returns, but appease clippy
        Gene::new(
            gene_def,
            nucleotide_sequence,
            nucleotide_index,
            genome_positions,
        )
    }

    /// Build all genes in the genome, storing them in the genes hashmap
    pub fn build_all_genes(&mut self) {
        for gene_name in self.gene_names.iter() {
            let gene = self.build_gene(gene_name.clone());
            self.genes.insert(gene_name.clone(), gene);
        }
    }

    /// Get a gene from the genome, building as required
    pub fn get_gene(&mut self, gene_name: String) -> Gene {
        // Return a gene from the hashmap if it exists, else build and cache it
        if !self.genes.contains_key(&gene_name) {
            let gene = self.build_gene(gene_name.clone());
            self.genes.insert(gene_name.clone(), gene);
        }

        // I hate implicit returns, but appease clippy
        self.genes.get(&gene_name).unwrap().clone()
    }

    /// Get the data at a given genome index
    ///
    /// # Arguments
    /// - `index` - 1-indexed genome index
    ///
    /// # Returns
    /// GenomePosition at the given index
    pub fn at_genome_index(&self, index: i64) -> GenomePosition {
        // 1-indexed genome index
        if index < 1 || index > self.genome_positions.len() as i64 {
            panic!("Genome index {} out of range", index);
        }
        self.genome_positions[(index - 1) as usize].clone()
    }

    /// Get the VCFRow associated with a given VCF row's index for this sample
    ///
    /// # Arguments
    /// - `index` - Index of the VCF row to retrieve (i.e row number of the records in a VCF)
    ///
    /// # Returns
    /// VCFRow associated with the given index
    pub fn get_vcf_row(&self, index: usize) -> VCFRow {
        if self.vcf_records.is_none() {
            panic!("No VCF records associated with this genome");
        }
        self.vcf_records.as_ref().unwrap()[index].clone()
    }
}

#[pyfunction]
/// Mutate a genome using a VCF file
///
/// # Arguments
/// - `reference` - Reference genome to mutate
/// - `vcf` - VCF file to use for mutation
///
/// # Returns
/// Mutated genome
pub fn mutate(reference: &Genome, vcf: VCFFile) -> Genome {
    let mut new_genome = reference.clone();
    let mut deleted_positions = HashSet::new();
    let mut deleted_evidence: HashMap<usize, Evidence> = HashMap::new();
    for (idx, position) in new_genome.genome_positions.iter_mut().enumerate() {
        if vcf.calls.contains_key(&position.genome_idx) {
            for call in vcf.calls.get(&position.genome_idx).unwrap() {
                // Set the new base etc from the call
                let c = call.clone();
                position.alts.push(Alt {
                    alt_type: c.call_type,
                    base: c.alt,
                    evidence: call.clone(),
                });

                // Mark containing gene as containing a mutation
                // Should make finding genes level mutations easier
                for gene_name in position.genes.iter() {
                    new_genome.genes_with_mutations.insert(gene_name.clone());
                }

                if call.call_type == AltType::DEL {
                    // Mark all associated bases as deleted
                    for del_idx in 0..call.alt.len() {
                        deleted_positions.insert(idx + del_idx);
                        deleted_evidence.insert(idx + del_idx, call.clone());
                    }
                }

                let mut base = call.clone().alt.chars().next().unwrap();
                if call.call_type == AltType::HET {
                    base = 'z';
                }
                if call.call_type == AltType::NULL {
                    base = 'x';
                }

                if call.call_type == AltType::SNP
                    || call.call_type == AltType::HET
                    || call.call_type == AltType::NULL
                {
                    // Update nucleotide for SNP/het/null
                    position.reference = base;
                    new_genome
                        .nucleotide_sequence
                        .replace_range(idx..idx + 1, &base.to_string());
                }
            }
        }

        if vcf.minor_calls.contains_key(&position.genome_idx) {
            for call in vcf.minor_calls.get(&position.genome_idx).unwrap() {
                // Set the new base etc from the call
                let c = call.clone();
                position.alts.push(Alt {
                    alt_type: c.call_type,
                    base: c.alt,
                    evidence: call.clone(),
                });

                // Mark containing gene as containing a mutation
                // Should make finding genes level mutations easier
                for gene_name in position.genes.iter() {
                    new_genome.genes_with_mutations.insert(gene_name.clone());
                }

                if call.call_type == AltType::DEL {
                    position.is_deleted_minor = true;
                    for del_idx in 0..call.alt.len() {
                        deleted_positions.insert(idx + del_idx);
                        deleted_evidence.insert(idx + del_idx, call.clone());
                    }
                }
            }
        }

        if deleted_positions.contains(&idx) {
            for ev in deleted_evidence.get(&idx).iter() {
                position.deleted_evidence.push((*ev).clone());
                // Set is deleted as appropriate from this evidence
                position.is_deleted = position.is_deleted || !ev.is_minor;
                position.is_deleted_minor = position.is_deleted_minor || ev.is_minor;
            }
            // Ensure that we mark all deleted positions as having mutations
            // this ensures we pick up cases where a deletion starts in an upstream gene
            for gene_name in position.genes.iter() {
                new_genome.genes_with_mutations.insert(gene_name.clone());
            }
        }
    }

    // Reset the gene hashmap as the nucleotide sequence has changed
    new_genome.genes = HashMap::new();

    // Keep track of the VCF records for ease of pulling out VCF rows later
    new_genome.vcf_records = Some(vcf.records.clone());

    // I hate implicit returns, but appease clippy
    new_genome
}

#[cfg(test)]
mod tests {
    use ordered_float::OrderedFloat;

    use crate::{
        common::{MinorType, VCFRow},
        difference::{GeneDifference, GenomeDifference, Mutation, Variant},
        gene::{codon_to_aa, complement_base},
    };
    use pretty_assertions::assert_eq;

    use super::*;

    macro_rules! assert_panics {
        ($expression:expr) => {
            let result = std::panic::catch_unwind(|| $expression);
            assert!(result.is_err());
        };
    }

    #[test]
    fn test_tb_genome() {
        let mut genome = Genome::new("reference/NC_000962.3.gbk");
        assert_eq!(genome.name, "NC_000962");
        assert_eq!(genome.nucleotide_sequence.len(), 4411532);
        assert_eq!(genome.gene_definitions.len(), 3909);
        assert_eq!(genome.gene_names.len(), 3909);

        // Check overlapping genes are fine
        assert_eq!(
            genome.genome_positions[2288680].genes,
            vec!["Rv2042c".to_string(), "pncA".to_string()]
        );

        let kat_g = genome.get_gene("katG".to_string());
        assert_eq!(kat_g.name, "katG");
        assert_eq!(kat_g.amino_acid_sequence,"VPEQHPPITETTTGAASNGCPVVGHMKYPVEGGGNQDWWPNRLNLKVLHQNPAVADPMGAAFDYAAEVATIDVDALTRDIEEVMTTSQPWWPADYGHYGPLFIRMAWHAAGTYRIHDGRGGAGGGMQRFAPLNSWPDNASLDKARRLLWPVKKKYGKKLSWADLIVFAGNCALESMGFKTFGFGFGRVDQWEPDEVYWGKEATWLGDERYSGKRDLENPLAAVQMGLIYVNPEGPNGNPDPMAAAVDIRETFRRMAMNDVETAALIVGGHTFGKTHGAGPADLVGPEPEAAPLEQMGLGWKSSYGTGTGKDAITSGIEVVWTNTPTKWDNSFLEILYGYEWELTKSPAGAWQYTAKDGAGAGTIPDPFGGPGRSPTMLATDLSLRVDPIYERITRRWLEHPEELADEFAKAWYKLIHRDMGPVARYLGPLVPKQTLLWQDPVPAVSHDLVGEAEIASLKSQIRASGLTVSQLVSTAWAAASSFRGSDKRGGANGGRIRLQPQVGWEVNDPDGDLRKVIRTLEEIQESFNSAAPGNIKVSFADLVVLGGCAAIEKAAKAAGHNITVPFTPGRTDASQEQTDVESFAVLEPKADGFRNYLGKGNPLPAEYMLLDKANLLTLSAPEMTVLVGGLRVLGANYKRLPLGVFTEASESLTNDFFVNLLDMGITWEPSPADDGTYQGKDGSGKVKWTGSRVDLVFGSNSELRALVEVYGADDAQPKFVQDFVAAWDKVMNLDRFDVR!".to_string());
        assert_eq!(
            kat_g.at_promoter(&kat_g.nucleotide_sequence.chars().collect::<Vec<char>>()),
            [
                't', 'c', 'a', 'c', 'a', 'g', 'c', 'c', 'c', 'g', 'a', 't', 'a', 'a', 'c', 'a',
                'c', 'c', 'a', 'a', 'c', 't', 'c', 'c', 't', 'g', 'g', 'a', 'a', 'g', 'g', 'a',
                'a', 't', 'g', 'c', 't'
            ]
        );
        assert_eq!(
            kat_g.not_promoter(&kat_g.nucleotide_sequence.chars().collect::<Vec<char>>()),
            [
                'g', 't', 'g', 'c', 'c', 'c', 'g', 'a', 'g', 'c', 'a', 'a', 'c', 'a', 'c', 'c',
                'c', 'a', 'c', 'c', 'c', 'a', 't', 't', 'a', 'c', 'a', 'g', 'a', 'a', 'a', 'c',
                'c', 'a', 'c', 'c', 'a', 'c', 'c', 'g', 'g', 'a', 'g', 'c', 'c', 'g', 'c', 't',
                'a', 'g', 'c', 'a', 'a', 'c', 'g', 'g', 'c', 't', 'g', 't', 'c', 'c', 'c', 'g',
                't', 'c', 'g', 't', 'g', 'g', 'g', 't', 'c', 'a', 't', 'a', 't', 'g', 'a', 'a',
                'a', 't', 'a', 'c', 'c', 'c', 'c', 'g', 't', 'c', 'g', 'a', 'g', 'g', 'g', 'c',
                'g', 'g', 'c', 'g', 'g', 'a', 'a', 'a', 'c', 'c', 'a', 'g', 'g', 'a', 'c', 't',
                'g', 'g', 't', 'g', 'g', 'c', 'c', 'c', 'a', 'a', 'c', 'c', 'g', 'g', 'c', 't',
                'c', 'a', 'a', 't', 'c', 't', 'g', 'a', 'a', 'g', 'g', 't', 'a', 'c', 't', 'g',
                'c', 'a', 'c', 'c', 'a', 'a', 'a', 'a', 'c', 'c', 'c', 'g', 'g', 'c', 'c', 'g',
                't', 'c', 'g', 'c', 't', 'g', 'a', 'c', 'c', 'c', 'g', 'a', 't', 'g', 'g', 'g',
                't', 'g', 'c', 'g', 'g', 'c', 'g', 't', 't', 'c', 'g', 'a', 'c', 't', 'a', 't',
                'g', 'c', 'c', 'g', 'c', 'g', 'g', 'a', 'g', 'g', 't', 'c', 'g', 'c', 'g', 'a',
                'c', 'c', 'a', 't', 'c', 'g', 'a', 'c', 'g', 't', 't', 'g', 'a', 'c', 'g', 'c',
                'c', 'c', 't', 'g', 'a', 'c', 'g', 'c', 'g', 'g', 'g', 'a', 'c', 'a', 't', 'c',
                'g', 'a', 'g', 'g', 'a', 'a', 'g', 't', 'g', 'a', 't', 'g', 'a', 'c', 'c', 'a',
                'c', 'c', 't', 'c', 'g', 'c', 'a', 'g', 'c', 'c', 'g', 't', 'g', 'g', 't', 'g',
                'g', 'c', 'c', 'c', 'g', 'c', 'c', 'g', 'a', 'c', 't', 'a', 'c', 'g', 'g', 'c',
                'c', 'a', 'c', 't', 'a', 'c', 'g', 'g', 'g', 'c', 'c', 'g', 'c', 't', 'g', 't',
                't', 't', 'a', 't', 'c', 'c', 'g', 'g', 'a', 't', 'g', 'g', 'c', 'g', 't', 'g',
                'g', 'c', 'a', 'c', 'g', 'c', 't', 'g', 'c', 'c', 'g', 'g', 'c', 'a', 'c', 'c',
                't', 'a', 'c', 'c', 'g', 'c', 'a', 't', 'c', 'c', 'a', 'c', 'g', 'a', 'c', 'g',
                'g', 'c', 'c', 'g', 'c', 'g', 'g', 'c', 'g', 'g', 'c', 'g', 'c', 'c', 'g', 'g',
                'g', 'g', 'g', 'c', 'g', 'g', 'c', 'a', 't', 'g', 'c', 'a', 'g', 'c', 'g', 'g',
                't', 't', 'c', 'g', 'c', 'g', 'c', 'c', 'g', 'c', 't', 't', 'a', 'a', 'c', 'a',
                'g', 'c', 't', 'g', 'g', 'c', 'c', 'c', 'g', 'a', 'c', 'a', 'a', 'c', 'g', 'c',
                'c', 'a', 'g', 'c', 't', 't', 'g', 'g', 'a', 'c', 'a', 'a', 'g', 'g', 'c', 'g',
                'c', 'g', 'c', 'c', 'g', 'g', 'c', 't', 'g', 'c', 't', 'g', 't', 'g', 'g', 'c',
                'c', 'g', 'g', 't', 'c', 'a', 'a', 'g', 'a', 'a', 'g', 'a', 'a', 'g', 't', 'a',
                'c', 'g', 'g', 'c', 'a', 'a', 'g', 'a', 'a', 'g', 'c', 't', 'c', 't', 'c', 'a',
                't', 'g', 'g', 'g', 'c', 'g', 'g', 'a', 'c', 'c', 't', 'g', 'a', 't', 't', 'g',
                't', 't', 't', 't', 'c', 'g', 'c', 'c', 'g', 'g', 'c', 'a', 'a', 'c', 't', 'g',
                'c', 'g', 'c', 'g', 'c', 't', 'g', 'g', 'a', 'a', 't', 'c', 'g', 'a', 't', 'g',
                'g', 'g', 'c', 't', 't', 'c', 'a', 'a', 'g', 'a', 'c', 'g', 't', 't', 'c', 'g',
                'g', 'g', 't', 't', 'c', 'g', 'g', 'c', 't', 't', 'c', 'g', 'g', 'c', 'c', 'g',
                'g', 'g', 't', 'c', 'g', 'a', 'c', 'c', 'a', 'g', 't', 'g', 'g', 'g', 'a', 'g',
                'c', 'c', 'c', 'g', 'a', 't', 'g', 'a', 'g', 'g', 't', 'c', 't', 'a', 't', 't',
                'g', 'g', 'g', 'g', 'c', 'a', 'a', 'g', 'g', 'a', 'a', 'g', 'c', 'c', 'a', 'c',
                'c', 't', 'g', 'g', 'c', 't', 'c', 'g', 'g', 'c', 'g', 'a', 't', 'g', 'a', 'g',
                'c', 'g', 't', 't', 'a', 'c', 'a', 'g', 'c', 'g', 'g', 't', 'a', 'a', 'g', 'c',
                'g', 'g', 'g', 'a', 't', 'c', 't', 'g', 'g', 'a', 'g', 'a', 'a', 'c', 'c', 'c',
                'g', 'c', 't', 'g', 'g', 'c', 'c', 'g', 'c', 'g', 'g', 't', 'g', 'c', 'a', 'g',
                'a', 't', 'g', 'g', 'g', 'g', 'c', 't', 'g', 'a', 't', 'c', 't', 'a', 'c', 'g',
                't', 'g', 'a', 'a', 'c', 'c', 'c', 'g', 'g', 'a', 'g', 'g', 'g', 'g', 'c', 'c',
                'g', 'a', 'a', 'c', 'g', 'g', 'c', 'a', 'a', 'c', 'c', 'c', 'g', 'g', 'a', 'c',
                'c', 'c', 'c', 'a', 't', 'g', 'g', 'c', 'c', 'g', 'c', 'g', 'g', 'c', 'g', 'g',
                't', 'c', 'g', 'a', 'c', 'a', 't', 't', 'c', 'g', 'c', 'g', 'a', 'g', 'a', 'c',
                'g', 't', 't', 't', 'c', 'g', 'g', 'c', 'g', 'c', 'a', 't', 'g', 'g', 'c', 'c',
                'a', 't', 'g', 'a', 'a', 'c', 'g', 'a', 'c', 'g', 't', 'c', 'g', 'a', 'a', 'a',
                'c', 'a', 'g', 'c', 'g', 'g', 'c', 'g', 'c', 't', 'g', 'a', 't', 'c', 'g', 't',
                'c', 'g', 'g', 'c', 'g', 'g', 't', 'c', 'a', 'c', 'a', 'c', 't', 't', 't', 'c',
                'g', 'g', 't', 'a', 'a', 'g', 'a', 'c', 'c', 'c', 'a', 't', 'g', 'g', 'c', 'g',
                'c', 'c', 'g', 'g', 'c', 'c', 'c', 'g', 'g', 'c', 'c', 'g', 'a', 't', 'c', 't',
                'g', 'g', 't', 'c', 'g', 'g', 'c', 'c', 'c', 'c', 'g', 'a', 'a', 'c', 'c', 'c',
                'g', 'a', 'g', 'g', 'c', 't', 'g', 'c', 't', 'c', 'c', 'g', 'c', 't', 'g', 'g',
                'a', 'g', 'c', 'a', 'g', 'a', 't', 'g', 'g', 'g', 'c', 't', 't', 'g', 'g', 'g',
                'c', 't', 'g', 'g', 'a', 'a', 'g', 'a', 'g', 'c', 't', 'c', 'g', 't', 'a', 't',
                'g', 'g', 'c', 'a', 'c', 'c', 'g', 'g', 'a', 'a', 'c', 'c', 'g', 'g', 't', 'a',
                'a', 'g', 'g', 'a', 'c', 'g', 'c', 'g', 'a', 't', 'c', 'a', 'c', 'c', 'a', 'g',
                'c', 'g', 'g', 'c', 'a', 't', 'c', 'g', 'a', 'g', 'g', 't', 'c', 'g', 't', 'a',
                't', 'g', 'g', 'a', 'c', 'g', 'a', 'a', 'c', 'a', 'c', 'c', 'c', 'c', 'g', 'a',
                'c', 'g', 'a', 'a', 'a', 't', 'g', 'g', 'g', 'a', 'c', 'a', 'a', 'c', 'a', 'g',
                't', 't', 't', 'c', 'c', 't', 'c', 'g', 'a', 'g', 'a', 't', 'c', 'c', 't', 'g',
                't', 'a', 'c', 'g', 'g', 'c', 't', 'a', 'c', 'g', 'a', 'g', 't', 'g', 'g', 'g',
                'a', 'g', 'c', 't', 'g', 'a', 'c', 'g', 'a', 'a', 'g', 'a', 'g', 'c', 'c', 'c',
                't', 'g', 'c', 't', 'g', 'g', 'c', 'g', 'c', 't', 't', 'g', 'g', 'c', 'a', 'a',
                't', 'a', 'c', 'a', 'c', 'c', 'g', 'c', 'c', 'a', 'a', 'g', 'g', 'a', 'c', 'g',
                'g', 'c', 'g', 'c', 'c', 'g', 'g', 't', 'g', 'c', 'c', 'g', 'g', 'c', 'a', 'c',
                'c', 'a', 't', 'c', 'c', 'c', 'g', 'g', 'a', 'c', 'c', 'c', 'g', 't', 't', 'c',
                'g', 'g', 'c', 'g', 'g', 'g', 'c', 'c', 'a', 'g', 'g', 'g', 'c', 'g', 'c', 't',
                'c', 'c', 'c', 'c', 'g', 'a', 'c', 'g', 'a', 't', 'g', 'c', 't', 'g', 'g', 'c',
                'c', 'a', 'c', 't', 'g', 'a', 'c', 'c', 't', 'c', 't', 'c', 'g', 'c', 't', 'g',
                'c', 'g', 'g', 'g', 't', 'g', 'g', 'a', 't', 'c', 'c', 'g', 'a', 't', 'c', 't',
                'a', 't', 'g', 'a', 'g', 'c', 'g', 'g', 'a', 't', 'c', 'a', 'c', 'g', 'c', 'g',
                't', 'c', 'g', 'c', 't', 'g', 'g', 'c', 't', 'g', 'g', 'a', 'a', 'c', 'a', 'c',
                'c', 'c', 'c', 'g', 'a', 'g', 'g', 'a', 'a', 't', 't', 'g', 'g', 'c', 'c', 'g',
                'a', 'c', 'g', 'a', 'g', 't', 't', 'c', 'g', 'c', 'c', 'a', 'a', 'g', 'g', 'c',
                'c', 't', 'g', 'g', 't', 'a', 'c', 'a', 'a', 'g', 'c', 't', 'g', 'a', 't', 'c',
                'c', 'a', 'c', 'c', 'g', 'a', 'g', 'a', 'c', 'a', 't', 'g', 'g', 'g', 't', 'c',
                'c', 'c', 'g', 't', 't', 'g', 'c', 'g', 'a', 'g', 'a', 't', 'a', 'c', 'c', 't',
                't', 'g', 'g', 'g', 'c', 'c', 'g', 'c', 't', 'g', 'g', 't', 'c', 'c', 'c', 'c',
                'a', 'a', 'g', 'c', 'a', 'g', 'a', 'c', 'c', 'c', 't', 'g', 'c', 't', 'g', 't',
                'g', 'g', 'c', 'a', 'g', 'g', 'a', 't', 'c', 'c', 'g', 'g', 't', 'c', 'c', 'c',
                't', 'g', 'c', 'g', 'g', 't', 'c', 'a', 'g', 'c', 'c', 'a', 'c', 'g', 'a', 'c',
                'c', 't', 'c', 'g', 't', 'c', 'g', 'g', 'c', 'g', 'a', 'a', 'g', 'c', 'c', 'g',
                'a', 'g', 'a', 't', 't', 'g', 'c', 'c', 'a', 'g', 'c', 'c', 't', 't', 'a', 'a',
                'g', 'a', 'g', 'c', 'c', 'a', 'g', 'a', 't', 'c', 'c', 'g', 'g', 'g', 'c', 'a',
                't', 'c', 'g', 'g', 'g', 'a', 't', 't', 'g', 'a', 'c', 't', 'g', 't', 'c', 't',
                'c', 'a', 'c', 'a', 'g', 'c', 't', 'a', 'g', 't', 't', 't', 'c', 'g', 'a', 'c',
                'c', 'g', 'c', 'a', 't', 'g', 'g', 'g', 'c', 'g', 'g', 'c', 'g', 'g', 'c', 'g',
                't', 'c', 'g', 't', 'c', 'g', 't', 't', 'c', 'c', 'g', 't', 'g', 'g', 't', 'a',
                'g', 'c', 'g', 'a', 'c', 'a', 'a', 'g', 'c', 'g', 'c', 'g', 'g', 'c', 'g', 'g',
                'c', 'g', 'c', 'c', 'a', 'a', 'c', 'g', 'g', 't', 'g', 'g', 't', 'c', 'g', 'c',
                'a', 't', 'c', 'c', 'g', 'c', 'c', 't', 'g', 'c', 'a', 'g', 'c', 'c', 'a', 'c',
                'a', 'a', 'g', 't', 'c', 'g', 'g', 'g', 't', 'g', 'g', 'g', 'a', 'g', 'g', 't',
                'c', 'a', 'a', 'c', 'g', 'a', 'c', 'c', 'c', 'c', 'g', 'a', 'c', 'g', 'g', 'g',
                'g', 'a', 't', 'c', 't', 'g', 'c', 'g', 'c', 'a', 'a', 'g', 'g', 't', 'c', 'a',
                't', 't', 'c', 'g', 'c', 'a', 'c', 'c', 'c', 't', 'g', 'g', 'a', 'a', 'g', 'a',
                'g', 'a', 't', 'c', 'c', 'a', 'g', 'g', 'a', 'g', 't', 'c', 'a', 't', 't', 'c',
                'a', 'a', 'c', 't', 'c', 'c', 'g', 'c', 'g', 'g', 'c', 'g', 'c', 'c', 'g', 'g',
                'g', 'g', 'a', 'a', 'c', 'a', 't', 'c', 'a', 'a', 'a', 'g', 't', 'g', 't', 'c',
                'c', 't', 't', 'c', 'g', 'c', 'c', 'g', 'a', 'c', 'c', 't', 'c', 'g', 't', 'c',
                'g', 't', 'g', 'c', 't', 'c', 'g', 'g', 't', 'g', 'g', 'c', 't', 'g', 't', 'g',
                'c', 'c', 'g', 'c', 'c', 'a', 't', 'a', 'g', 'a', 'g', 'a', 'a', 'a', 'g', 'c',
                'a', 'g', 'c', 'a', 'a', 'a', 'g', 'g', 'c', 'g', 'g', 'c', 't', 'g', 'g', 'c',
                'c', 'a', 'c', 'a', 'a', 'c', 'a', 't', 'c', 'a', 'c', 'g', 'g', 't', 'g', 'c',
                'c', 'c', 't', 't', 'c', 'a', 'c', 'c', 'c', 'c', 'g', 'g', 'g', 'c', 'c', 'g',
                'c', 'a', 'c', 'g', 'g', 'a', 't', 'g', 'c', 'g', 't', 'c', 'g', 'c', 'a', 'g',
                'g', 'a', 'a', 'c', 'a', 'a', 'a', 'c', 'c', 'g', 'a', 'c', 'g', 't', 'g', 'g',
                'a', 'a', 't', 'c', 'c', 't', 't', 't', 'g', 'c', 'c', 'g', 't', 'g', 'c', 't',
                'g', 'g', 'a', 'g', 'c', 'c', 'c', 'a', 'a', 'g', 'g', 'c', 'a', 'g', 'a', 't',
                'g', 'g', 'c', 't', 't', 'c', 'c', 'g', 'a', 'a', 'a', 'c', 't', 'a', 'c', 'c',
                't', 'c', 'g', 'g', 'a', 'a', 'a', 'g', 'g', 'g', 'c', 'a', 'a', 'c', 'c', 'c',
                'g', 't', 't', 'g', 'c', 'c', 'g', 'g', 'c', 'c', 'g', 'a', 'g', 't', 'a', 'c',
                'a', 't', 'g', 'c', 't', 'g', 'c', 't', 'c', 'g', 'a', 'c', 'a', 'a', 'g', 'g',
                'c', 'g', 'a', 'a', 'c', 'c', 't', 'g', 'c', 't', 't', 'a', 'c', 'g', 'c', 't',
                'c', 'a', 'g', 't', 'g', 'c', 'c', 'c', 'c', 't', 'g', 'a', 'g', 'a', 't', 'g',
                'a', 'c', 'g', 'g', 't', 'g', 'c', 't', 'g', 'g', 't', 'a', 'g', 'g', 't', 'g',
                'g', 'c', 'c', 't', 'g', 'c', 'g', 'c', 'g', 't', 'c', 'c', 't', 'c', 'g', 'g',
                'c', 'g', 'c', 'a', 'a', 'a', 'c', 't', 'a', 'c', 'a', 'a', 'g', 'c', 'g', 'c',
                't', 't', 'a', 'c', 'c', 'g', 'c', 't', 'g', 'g', 'g', 'c', 'g', 't', 'g', 't',
                't', 'c', 'a', 'c', 'c', 'g', 'a', 'g', 'g', 'c', 'c', 't', 'c', 'c', 'g', 'a',
                'g', 't', 'c', 'a', 'c', 't', 'g', 'a', 'c', 'c', 'a', 'a', 'c', 'g', 'a', 'c',
                't', 't', 'c', 't', 't', 'c', 'g', 't', 'g', 'a', 'a', 'c', 'c', 't', 'g', 'c',
                't', 'c', 'g', 'a', 'c', 'a', 't', 'g', 'g', 'g', 't', 'a', 't', 'c', 'a', 'c',
                'c', 't', 'g', 'g', 'g', 'a', 'g', 'c', 'c', 'c', 't', 'c', 'g', 'c', 'c', 'a',
                'g', 'c', 'a', 'g', 'a', 't', 'g', 'a', 'c', 'g', 'g', 'g', 'a', 'c', 'c', 't',
                'a', 'c', 'c', 'a', 'g', 'g', 'g', 'c', 'a', 'a', 'g', 'g', 'a', 't', 'g', 'g',
                'c', 'a', 'g', 't', 'g', 'g', 'c', 'a', 'a', 'g', 'g', 't', 'g', 'a', 'a', 'g',
                't', 'g', 'g', 'a', 'c', 'c', 'g', 'g', 'c', 'a', 'g', 'c', 'c', 'g', 'c', 'g',
                't', 'g', 'g', 'a', 'c', 'c', 't', 'g', 'g', 't', 'c', 't', 't', 'c', 'g', 'g',
                'g', 't', 'c', 'c', 'a', 'a', 'c', 't', 'c', 'g', 'g', 'a', 'g', 't', 't', 'g',
                'c', 'g', 'g', 'g', 'c', 'g', 'c', 't', 't', 'g', 't', 'c', 'g', 'a', 'g', 'g',
                't', 'c', 't', 'a', 't', 'g', 'g', 'c', 'g', 'c', 'c', 'g', 'a', 't', 'g', 'a',
                'c', 'g', 'c', 'g', 'c', 'a', 'g', 'c', 'c', 'g', 'a', 'a', 'g', 't', 't', 'c',
                'g', 't', 'g', 'c', 'a', 'g', 'g', 'a', 'c', 't', 't', 'c', 'g', 't', 'c', 'g',
                'c', 't', 'g', 'c', 'c', 't', 'g', 'g', 'g', 'a', 'c', 'a', 'a', 'g', 'g', 't',
                'g', 'a', 't', 'g', 'a', 'a', 'c', 'c', 't', 'c', 'g', 'a', 'c', 'a', 'g', 'g',
                't', 't', 'c', 'g', 'a', 'c', 'g', 't', 'g', 'c', 'g', 'c', 't', 'g', 'a'
            ]
        );

        let rpo_b = genome.get_gene("rpoB".to_string());
        assert_eq!(rpo_b.name, "rpoB");
        assert_eq!(rpo_b.amino_acid_sequence,"LADSRQSKTAASPSPSRPQSSSNNSVPGAPNRVSFAKLREPLEVPGLLDVQTDSFEWLIGSPRWRESAAERGDVNPVGGLEEVLYELSPIEDFSGSMSLSFSDPRFDDVKAPVDECKDKDMTYAAPLFVTAEFINNNTGEIKSQTVFMGDFPMMTEKGTFIINGTERVVVSQLVRSPGVYFDETIDKSTDKTLHSVKVIPSRGAWLEFDVDKRDTVGVRIDRKRRQPVTVLLKALGWTSEQIVERFGFSEIMRSTLEKDNTVGTDEALLDIYRKLRPGEPPTKESAQTLLENLFFKEKRYDLARVGRYKVNKKLGLHVGEPITSSTLTEEDVVATIEYLVRLHEGQTTMTVPGGVEVPVETDDIDHFGNRRLRTVGELIQNQIRVGMSRMERVVRERMTTQDVEAITPQTLINIRPVVAAIKEFFGTSQLSQFMDQNNPLSGLTHKRRLSALGPGGLSRERAGLEVRDVHPSHYGRMCPIETPEGPNIGLIGSLSVYARVNPFGFIETPYRKVVDGVVSDEIVYLTADEEDRHVVAQANSPIDADGRFVEPRVLVRRKAGEVEYVPSSEVDYMDVSPRQMVSVATAMIPFLEHDDANRALMGANMQRQAVPLVRSEAPLVGTGMELRAAIDAGDVVVAEESGVIEEVSADYITVMHDNGTRRTYRMRKFARSNHGTCANQCPIVDAGDRVEAGQVIADGPCTDDGEMALGKNLLVAIMPWEGHNYEDAIILSNRLVEEDVLTSIHIEEHEIDARDTKLGAEEITRDIPNISDEVLADLDERGIVRIGAEVRDGDILVGKVTPKGETELTPEERLLRAIFGEKAREVRDTSLKVPHGESGKVIGIRVFSREDEDELPAGVNELVRVYVAQKRKISDGDKLAGRHGNKGVIGKILPVEDMPFLADGTPVDIILNTHGVPRRMNIGQILETHLGWCAHSGWKVDAAKGVPDWAARLPDELLEAQPNAIVSTPVFDGAQEAELQGLLSCTLPNRDGDVLVDADGKAMLFDGRSGEPFPYPVTVGYMYIMKLHHLVDDKIHARSTGPYSMITQQPLGGKAQFGGQRFGEMECWAMQAYGAAYTLQELLTIKSDDTVGRVKVYEAIVKGENIPEPGIPESFKVLLKELQSLCLNVEVLSSDGAAIELREGEDEDLERAAANLGINLSRNESASVEDLA!".to_string());

        let pnc_a = genome.get_gene("pncA".to_string());
        assert_eq!(pnc_a.name, "pncA");
        assert_eq!(pnc_a.amino_acid_sequence,"MRALIIVDVQNDFCEGGSLAVTGGAALARAISDYLAEAADYHHVVATKDFHIDPGDHFSGTPDYSSSWPPHCVSGTPGADFHPSLDTSAIEAVFYKGAYTGAYSGFEGVDENGTPLLNWLRQRGVDEVDVVGIATDHCVRQTAEDAVRNGLATRVLVDLTAGVSADTTVAALEEMRTASVELVCSS!".to_string());

        let rv2042_c = genome.get_gene("Rv2042c".to_string());
        assert_eq!(rv2042_c.name, "Rv2042c");
        assert_eq!(rv2042_c.amino_acid_sequence,"MAPPNRDELLAAVERSPQAAAAHDRAGWVGLFTGDARVEDPVGSQPQVGHEAIGRFYDTFIGPRDITFHRDLDIVSGTVVLRDLELEVAMDSAVTVFIPAFLRYDLRPVTGEWQIAALRAYWELPAMMLQFLRTGSGATRPALQLSRALLGNQGLGGTAGFLTGFRRAGRRHKKLVETFLNAASRADKSAAYHALSRTATMTLGEDELLDIVELFEQLRGASWTKVTGAGSTVAVSLASDHRRGIMFADVPWRGNRINRIRYFPA!".to_string());

        let rrs = genome.get_gene("rrs".to_string());
        assert_eq!(rrs.name, "rrs");
        assert_eq!(rrs.amino_acid_sequence, "".to_string());
        assert_eq!(
            rrs.at_promoter(&rrs.nucleotide_sequence.chars().collect::<Vec<char>>()),
            [
                'g', 'g', 'c', 'c', 'a', 't', 'g', 'c', 't', 'c', 't', 't', 'g', 'a', 't', 'g',
                'c', 'c', 'c', 'c', 'g', 't', 't', 'g', 't', 'c', 'g', 'g', 'g', 'g', 'g', 'c',
                'g', 't', 'g', 'g', 'c', 'c', 'g', 't', 't', 't', 'g', 't', 't', 't', 't', 'g',
                't', 'c', 'a', 'g', 'g', 'a', 't', 'a', 't', 't', 't', 'c', 't', 'a', 'a', 'a',
                't', 'a', 'c', 'c', 't', 't', 't', 'g', 'g', 'c', 't', 'c', 'c', 'c', 't', 't',
                't', 't', 'c', 'c', 'a', 'a', 'a', 'g', 'g', 'g', 'a', 'g', 't', 'g', 't', 't',
                't', 'g', 'g', 'g'
            ]
        );
        assert_eq!(
            rrs.not_promoter(&rrs.nucleotide_sequence.chars().collect::<Vec<char>>()).iter().collect::<String>(),
            "ttttgtttggagagtttgatcctggctcaggacgaacgctggcggcgtgcttaacacatgcaagtcgaacggaaaggtctcttcggagatactcgagtggcgaacgggtgagtaacacgtgggtgatctgccctgcacttcgggataagcctgggaaactgggtctaataccggataggaccacgggatgcatgtcttgtggtggaaagcgctttagcggtgtgggatgagcccgcggcctatcagcttgttggtggggtgacggcctaccaaggcgacgacgggtagccggcctgagagggtgtccggccacactgggactgagatacggcccagactcctacgggaggcagcagtggggaatattgcacaatgggcgcaagcctgatgcagcgacgccgcgtgggggatgacggccttcgggttgtaaacctctttcaccatcgacgaaggtccgggttctctcggattgacggtaggtggagaagaagcaccggccaactacgtgccagcagccgcggtaatacgtagggtgcgagcgttgtccggaattactgggcgtaaagagctcgtaggtggtttgtcgcgttgttcgtgaaatctcacggcttaactgtgagcgtgcgggcgatacgggcagactagagtactgcaggggagactggaattcctggtgtagcggtggaatgcgcagatatcaggaggaacaccggtggcgaaggcgggtctctgggcagtaactgacgctgaggagcgaaagcgtggggagcgaacaggattagataccctggtagtccacgccgtaaacggtgggtactaggtgtgggtttccttccttgggatccgtgccgtagctaacgcattaagtaccccgcctggggagtacggccgcaaggctaaaactcaaaggaattgacgggggcccgcacaagcggcggagcatgtggattaattcgatgcaacgcgaagaaccttacctgggtttgacatgcacaggacgcgtctagagataggcgttcccttgtggcctgtgtgcaggtggtgcatggctgtcgtcagctcgtgtcgtgagatgttgggttaagtcccgcaacgagcgcaacccttgtctcatgttgccagcacgtaatggtggggactcgtgagagactgccggggtcaactcggaggaaggtggggatgacgtcaagtcatcatgccccttatgtccagggcttcacacatgctacaatggccggtacaaagggctgcgatgccgcgaggttaagcgaatccttaaaagccggtctcagttcggatcggggtctgcaactcgaccccgtgaagtcggagtcgctagtaatcgcagatcagcaacgctgcggtgaatacgttcccgggccttgtacacaccgcccgtcacgtcatgaaagtcggtaacacccgaagccagtggcctaaccctcgggagggagctgtcgaaggtgggatcggcgattgggacgaagtcgtaacaaggtagccgtaccggaaggtgcggctggatcacctcctttct".to_string()
        );
    }

    #[test]
    fn test_covid_genome() {
        let mut genome = Genome::new("reference/MN908947.3.gb");
        assert_eq!(genome.name, "MN908947");
        assert_eq!(genome.nucleotide_sequence.len(), 29903);
        assert_eq!(genome.gene_definitions.len(), 10);
        assert_eq!(genome.gene_names.len(), 10);

        // orf1ab has ribosomal shifts, so ensure these work by comparing to a truth sequence of protein
        let orf1ab = genome.get_gene("orf1ab".to_string());
        assert_eq!(orf1ab.name, "orf1ab");
        assert_eq!(
            orf1ab.amino_acid_sequence,
            "MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGIQYGRSGETLGVLVPHVGEIPVAYRKVLLRKNGNKGAGGHSYGADLKSFDLGDELGTDPYEDFQENWNTKHSSGVTRELMRELNGGAYTRYVDNNFCGPDGYPLECIKDLLARAGKASCTLSEQLDFIDTKRGVYCCREHEHEIAWYTERSEKSYELQTPFEIKLAKKFDTFNGECPNFVFPLNSIIKTIQPRVEKKKLDGFMGRIRSVYPVASPNECNQMCLSTLMKCDHCGETSWQTGDFVKATCEFCGTENLTKEGATTCGYLPQNAVVKIYCPACHNSEVGPEHSLAEYHNESGLKTILRKGGRTIAFGGCVFSYVGCHNKCAYWVPRASANIGCNHTGVVGEGSEGLNDNLLEILQKEKVNINIVGDFKLNEEIAIILASFSASTSAFVETVKGLDYKAFKQIVESCGNFKVTKGKAKKGAWNIGEQKSILSPLYAFASEAARVVRSIFSRTLETAQNSVRVLQKAAITILDGISQYSLRLIDAMMFTSDLATNNLVVMAYITGGVVQLTSQWLTNIFGTVYEKLKPVLDWLEEKFKEGVEFLRDGWEIVKFISTCACEIVGGQIVTCAKEIKESVQTFFKLVNKFLALCADSIIIGGAKLKALNLGETFVTHSKGLYRKCVKSREETGLLMPLKAPKEIIFLEGETLPTEVLTEEVVLKTGDLQPLEQPTSEAVEAPLVGTPVCINGLMLLEIKDTEKYCALAPNMMVTNNTFTLKGGAPTKVTFGDDTVIEVQGYKSVNITFELDERIDKVLNEKCSAYTVELGTEVNEFACVVADAVIKTLQPVSELLTPLGIDLDEWSMATYYLFDESGEFKLASHMYCSFYPPDEDEEEGDCEEEEFEPSTQYEYGTEDDYQGKPLEFGATSAALQPEEEQEEDWLDDDSQQTVGQQDGSEDNQTTTIQTIVEVQPQLEMELTPVVQTIEVNSFSGYLKLTDNVYIKNADIVEEAKKVKPTVVVNAANVYLKHGGGVAGALNKATNNAMQVESDDYIATNGPLKVGGSCVLSGHNLAKHCLHVVGPNVNKGEDIQLLKSAYENFNQHEVLLAPLLSAGIFGADPIHSLRVCVDTVRTNVYLAVFDKNLYDKLVSSFLEMKSEKQVEQKIAEIPKEEVKPFITESKPSVEQRKQDDKKIKACVEEVTTTLEETKFLTENLLLYIDINGNLHPDSATLVSDIDITFLKKDAPYIVGDVVQEGVLTAVVIPTKKAGGTTEMLAKALRKVPTDNYITTYPGQGLNGYTVEEAKTVLKKCKSAFYILPSIISNEKQEILGTVSWNLREMLAHAEETRKLMPVCVETKAIVSTIQRKYKGIKIQEGVVDYGARFYFYTSKTTVASLINTLNDLNETLVTMPLGYVTHGLNLEEAARYMRSLKVPATVSVSSPDAVTAYNGYLTSSSKTPEEHFIETISLAGSYKDWSYSGQSTQLGIEFLKRGDKSVYYTSNPTTFHLDGEVITFDNLKTLLSLREVRTIKVFTTVDNINLHTQVVDMSMTYGQQFGPTYLDGADVTKIKPHNSHEGKTFYVLPNDDTLRVEAFEYYHTTDPSFLGRYMSALNHTKKWKYPQVNGLTSIKWADNNCYLATALLTLQQIELKFNPPALQDAYYRARAGEAANFCALILAYCNKTVGELGDVRETMSYLFQHANLDSCKRVLNVVCKTCGQQQTTLKGVEAVMYMGTLSYEQFKKGVQIPCTCGKQATKYLVQQESPFVMMSAPPAQYELKHGTFTCASEYTGNYQCGHYKHITSKETLYCIDGALLTKSSEYKGPITDVFYKENSYTTTIKPVTYKLDGVVCTEIDPKLDNYYKKDNSYFTEQPIDLVPNQPYPNASFDNFKFVCDNIKFADDLNQLTGYKKPASRELKVTFFPDLNGDVVAIDYKHYTPSFKKGAKLLHKPIVWHVNNATNKATYKPNTWCIRCLWSTKPVETSNSFDVLKSEDAQGMDNLACEDLKPVSEEVVENPTIQKDVLECNVKTTEVVGDIILKPANNSLKITEEVGHTDLMAAYVDNSSLTIKKPNELSRVLGLKTLATHGLAAVNSVPWDTIANYAKPFLNKVVSTTTNIVTRCLNRVCTNYMPYFFTLLLQLCTFTRSTNSRIKASMPTTIAKNTVKSVGKFCLEASFNYLKSPNFSKLINIIIWFLLLSVCLGSLIYSTAALGVLMSNLGMPSYCTGYREGYLNSTNVTIATYCTGSIPCSVCLSGLDSLDTYPSLETIQITISSFKWDLTAFGLVAEWFLAYILFTRFFYVLGLAAIMQLFFSYFAVHFISNSWLMWLIINLVQMAPISAMVRMYIFFASFYYVWKSYVHVVDGCNSSTCMMCYKRNRATRVECTTIVNGVRRSFYVYANGGKGFCKLHNWNCVNCDTFCAGSTFISDEVARDLSLQFKRPINPTDQSSYIVDSVTVKNGSIHLYFDKAGQKTYERHSLSHFVNLDNLRANNTKGSLPINVIVFDGKSKCEESSAKSASVYYSQLMCQPILLLDQALVSDVGDSAEVAVKMFDAYVNTFSSTFNVPMEKLKTLVATAEAELAKNVSLDNVLSTFISAARQGFVDSDVETKDVVECLKLSHQSDIEVTGDSCNNYMLTYNKVENMTPRDLGACIDCSARHINAQVAKSHNIALIWNVKDFMSLSEQLRKQIRSAAKKNNLPFKLTCATTRQVVNVVTTKIALKGGKIVNNWLKQLIKVTLVFLFVAAIFYLITPVHVMSKHTDFSSEIIGYKAIDGGVTRDIASTDTCFANKHADFDTWFSQRGGSYTNDKACPLIAAVITREVGFVVPGLPGTILRTTNGDFLHFLPRVFSAVGNICYTPSKLIEYTDFATSACVLAAECTIFKDASGKPVPYCYDTNVLEGSVAYESLRPDTRYVLMDGSIIQFPNTYLEGSVRVVTTFDSEYCRHGTCERSEAGVCVSTSGRWVLNNDYYRSLPGVFCGVDAVNLLTNMFTPLIQPIGALDISASIVAGGIVAIVVTCLAYYFMRFRRAFGEYSHVVAFNTLLFLMSFTVLCLTPVYSFLPGVYSVIYLYLTFYLTNDVSFLAHIQWMVMFTPLVPFWITIAYIICISTKHFYWFFSNYLKRRVVFNGVSFSTFEEAALCTFLLNKEMYLKLRSDVLLPLTQYNRYLALYNKYKYFSGAMDTTSYREAACCHLAKALNDFSNSGSDVLYQPPQTSITSAVLQSGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDVVYCPRHVICTSEDMLNPNYEDLLIRKSNHNFLVQAGNVQLRVIGHSMQNCVLKLKVDTANPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNFTIKGSFLNGSCGSVGFNIDYDCVSFCYMHHMELPTGVHAGTDLEGNFYGPFVDRQTAQAAGTDTTITVNVLAWLYAAVINGDRWFLNRFTTTLNDFNLVAMKYNYEPLTQDHVDILGPLSAQTGIAVLDMCASLKELLQNGMNGRTILGSALLEDEFTPFDVVRQCSGVTFQSAVKRTIKGTHHWLLLTILTSLLVLVQSTQWSLFFFLYENAFLPFAMGIIAMSAFAMMFVKHKHAFLCLFLLPSLATVAYFNMVYMPASWVMRIMTWLDMVDTSLSGFKLKDCVMYASAVVLLILMTARTVYDDGARRVWTLMNVLTLVYKVYYGNALDQAISMWALIISVTSNYSGVVTTVMFLARGIVFMCVEYCPIFFITGNTLQCIMLVYCFLGYFCTCYFGLFCLLNRYFRLTLGVYDYLVSTQEFRYMNSQGLLPPKNSIDAFKLNIKLLGVGGKPCIKVATVQSKMSDVKCTSVVLLSVLQQLRVESSSKLWAQCVQLHNDILLAKDTTEAFEKMVSLLSVLLSMQGAVDINKLCEEMLDNRATLQAIASEFSSLPSYAAFATAQEAYEQAVANGDSEVVLKKLKKSLNVAKSEFDRDAAMQRKLEKMADQAMTQMYKQARSEDKRAKVTSAMQTMLFTMLRKLDNDALNNIINNARDGCVPLNIIPLTTAAKLMVVIPDYNTYKNTCDGTTFTYASALWEIQQVVDADSKIVQLSEISMDNSPNLAWPLIVTALRANSAVKLQNNELSPVALRQMSCAAGTTQTACTDDNALAYYNTTKGGRFVLALLSDLQDLKWARFPKSDGTGTIYTELEPPCRFVTDTPKGPKVKYLYFIKGLNNLNRGMVLGSLAATVRLQAGNATEVPANSTVLSFCAFAVDAAKAYKDYLASGGQPITNCVKMLCTHTGTGQAITVTPEANMDQESFGGASCCLYCRCHIDHPNPKGFCDLKGKYVQIPTTCANDPVGFTLKNTVCTVCGMWKGYGCSCDQLREPMLQSADAQSFLNRVCGVSAARLTPCGTGTSTDVVYRAFDIYNDKVAGFAKFLKTNCCRFQEKDEDDNLIDSYFVVKRHTFSNYQHEETIYNLLKDCPAVAKHDFFKFRIDGDMVPHISRQRLTKYTMADLVYALRHFDEGNCDTLKEILVTYNCCDDDYFNKKDWYDFVENPDILRVYANLGERVRQALLKTVQFCDAMRNAGIVGVLTLDNQDLNGNWYDFGDFIQTTPGSGVPVVDSYYSLLMPILTLTRALTAESHVDTDLTKPYIKWDLLKYDFTEERLKLFDRYFKYWDQTYHPNCVNCLDDRCILHCANFNVLFSTVFPPTSFGPLVRKIFVDGVPFVVSTGYHFRELGVVHNQDVNLHSSRLSFKELLVYAADPAMHAASGNLLLDKRTTCFSVAALTNNVAFQTVKPGNFNKDFYDFAVSKGFFKEGSSVELKHFFFAQDGNAAISDYDYYRYNLPTMCDIRQLLFVVEVVDKYFDCYDGGCINANQVIVNNLDKSAGFPFNKWGKARLYYDSMSYEDQDALFAYTKRNVIPTITQMNLKYAISAKNRARTVAGVSICSTMTNRQFHQKLLKSIAATRGATVVIGTSKFYGGWHNMLKTVYSDVENPHLMGWDYPKCDRAMPNMLRIMASLVLARKHTTCCSLSHRFYRLANECAQVLSEMVMCGGSLYVKPGGTSSGDATTAYANSVFNICQAVTANVNALLSTDGNKIADKYVRNLQHRLYECLYRNRDVDTDFVNEFYAYLRKHFSMMILSDDAVVCFNSTYASQGLVASIKNFKSVLYYQNNVFMSEAKCWTETDLTKGPHEFCSQHTMLVKQGDDYVYLPYPDPSRILGAGCFVDDIVKTDGTLMIERFVSLAIDAYPLTKHPNQEYADVFHLYLQYIRKLHDELTGHMLDMYSVMLTNDNTSRYWEPEFYEAMYTPHTVLQAVGACVLCNSQTSLRCGACIRRPFLCCKCCYDHVISTSHKLVLSVNPYVCNAPGCDVTDVTQLYLGGMSYYCKSHKPPISFPLCANGQVFGLYKNTCVGSDNVTDFNAIATCDWTNAGDYILANTCTERLKLFAAETLKATEETFKLSYGIATVREVLSDRELHLSWEVGKPRPPLNRNYVFTGYRVTKNSKVQIGEYTFEKGDYGDAVVYRGTTTYKLNVGDYFVLTSHTVMPLSAPTLVPQEHYVRITGLYPTLNISDEFSSNVANYQKVGMQKYSTLQGPPGTGKSHFAIGLALYYPSARIVYTACSHAAVDALCEKALKYLPIDKCSRIIPARARVECFDKFKVNSTLEQYVFCTVNALPETTADIVVFDEISMATNYDLSVVNARLRAKHYVYIGDPAQLPAPRTLLTKGTLEPEYFNSVCRLMKTIGPDMFLGTCRRCPAEIVDTVSALVYDNKLKAHKDKSAQCFKMFYKGVITHDVSSAINRPQIGVVREFLTRNPAWRKAVFISPYNSQNAVASKILGLPTQTVDSSQGSEYDYVIFTQTTETAHSCNVNRFNVAITRAKVGILCIMSDRDLYDKLQFTSLEIPRRNVATLQAENVTGLFKDCSKVITGLHPTQAPTHLSVDTKFKTEGLCVDIPGIPKDMTYRRLISMMGFKMNYQVNGYPNMFITREEAIRHVRAWIGFDVEGCHATREAVGTNLPLQLGFSTGVNLVAVPTGYVDTPNNTDFSRVSAKPPPGDQFKHLIPLMYKGLPWNVVRIKIVQMLSDTLKNLSDRVVFVLWAHGFELTSMKYFVKIGPERTCCLCDRRATCFSTASDTYACWHHSIGFDYVYNPFMIDVQQWGFTGNLQSNHDLYCQVHGNAHVASCDAIMTRCLAVHECFVKRVDWTIEYPIIGDELKINAACRKVQHMVVKAALLADKFPVLHDIGNPKAIKCVPQADVEWKFYDAQPCSDKAYKIEELFYSYATHSDKFTDGVCLFWNCNVDRYPANSIVCRFDTRVLSNLNLPGCDGGSLYVNKHAFHTPAFDKSAFVNLKQLPFFYYSDSPCESHGKQVVSDIDYVPLKSATCITRCNLGGAVCRHHANEYRLYLDAYNMMISAGFSLWVYKQFDTYNLWNTFTRLQSLENVAFNVVNKGHFDGQQGEVPVSIINNTVYTKVDGVDVELFENKTTLPVNVAFELWAKRNIKPVPEVKILNNLGVDIAANTVIWDYKRDAPAHISTIGVCSMTDIAKKPTETICAPLTVFFDGRVDGQVDLFRNARNGVLITEGSVKGLQPSVGPKQASLNGVTLIGEAVKTQFNYYKKVDGVVQQLPETYFTQSRNLQEFKPRSQMEIDFLELAMDEFIERYKLEGYAFEHIVYGDFSHSQLGGLHLLIGLAKRFKESPFELEDFIPMDSTVKNYFITDAQTGSSKCVCSVIDLLLDDFVEIIKSQDLSVVSKVVKVTIDYTEISFMLWCKDGHVETFYPKLQSSQAWQPGVAMPNLYKMQRMLLEKCDLQNYGDSATLPKGIMMNVAKYTQLCQYLNTLTLAVPYNMRVIHFGAGSDKGVAPGTAVLRQWLPTGTLLVDSDLNDFVSDADSTLIGDCATVHTANKWDLIISDMYDPKTKNVTKENDSKEGFFTYICGFIQQKLALGGSVAIKITEHSWNADLYKLMGHFAWWTAFVTNVNASSSEAFLIGCNYLGKPREQIDGYVMHANYIFWRNTNPIQLSSYSLFDMSKFPLKLRGTAVMSLKEGQINDMILSLLSKGRLIIRENNRVVISSDVLVNN!".to_string()
        );
    }

    #[test]
    fn test_invalid_args() {
        assert_panics!(Genome::new("not/a/path"));
    }

    #[test]
    fn test_dna_genome() {
        // Dummy DNA genbank file, used as it's short enough to know the expected values
        let genome = Genome::new("reference/TEST-DNA.gbk");
        assert_eq!(genome.name, "TEST_DNA");
        assert_eq!(genome.nucleotide_sequence, "aaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccc");
        let mut offset = 0;
        for (idx, pos) in genome.genome_positions.iter().enumerate() {
            if idx - offset < 10 {
                assert_eq!(pos.reference, 'a')
            } else if idx - offset < 20 {
                assert_eq!(pos.reference, 'c');
            } else if idx - offset < 30 {
                assert_eq!(pos.reference, 'g');
            } else if idx - offset < 40 {
                assert_eq!(pos.reference, 't');
            } else {
                assert_eq!(pos.reference, 'a');
                offset = idx;
            }

            if (0..30).filter(|x| x == &idx).count() > 0 {
                assert!(pos.genes.contains(&"A".to_string()));
            }
            if (28..60).filter(|x| x == &idx).count() > 0 {
                assert!(pos.genes.contains(&"B".to_string()));
            }
            if (91..99).filter(|x| x == &idx).count() > 0 {
                assert!(pos.genes.contains(&"C".to_string()));
            }
        }
    }

    #[test]
    fn test_dna_genes() {
        let mut genome = Genome::new("reference/TEST-DNA.gbk");
        let gene_a = genome.get_gene("A".to_string());
        assert_eq!(gene_a.name, "A");
        assert_eq!(gene_a.nucleotide_sequence, "aaaaaaaaaaccccccccccgggggggggg");
        assert_eq!(gene_a.amino_acid_sequence, "KKTPPPGGG".to_string());
        assert_eq!(
            gene_a.codons,
            vec!["aaa", "aaa", "acc", "ccc", "ccc", "ccg", "ggg", "ggg", "ggg"]
        );
        assert!(gene_a.coding);
        assert!(!gene_a.reverse_complement);

        let gene_b = genome.get_gene("B".to_string());
        assert_eq!(gene_b.name, "B");
        assert_eq!(
            gene_b.nucleotide_sequence,
            "gggttttttttttaaaaaaaaaacccccccccc"
        );
        assert_eq!(gene_b.amino_acid_sequence, "GFFF!KKNPPP".to_string());
        assert_eq!(
            gene_b.codons,
            vec!["ggg", "ttt", "ttt", "ttt", "taa", "aaa", "aaa", "aac", "ccc", "ccc", "ccc"]
        );
        assert!(gene_b.coding);
        assert!(!gene_b.reverse_complement);

        let gene_c = genome.get_gene("C".to_string());
        assert_eq!(gene_c.name, "C");
        assert_eq!(gene_c.nucleotide_sequence, "ggggggggg");
        assert_eq!(gene_c.amino_acid_sequence, "GG".to_string());
        assert_eq!(gene_c.codons, vec!["ggg", "ggg"]);
        assert!(gene_c.coding);
        assert!(gene_c.reverse_complement);
    }

    #[test]
    fn test_mutate_genome() {
        let mut reference = Genome::new("reference/MN908947.3.gb");
        // Note that this VCF is completely dummy data
        // so variants aren't necessarily matching the reference in the VCF row
        let vcf = VCFFile::new("test/dummy.vcf".to_string(), false, 3);
        let mut sample = mutate(&reference, vcf);

        let genome_diff = GenomeDifference::new(reference.clone(), sample.clone(), MinorType::COV);

        let expected_vcf_rows = vec![
            VCFRow {
                // 0
                position: 4687,
                reference: "t".to_string(),
                alternative: vec!["c".to_string()],
                filter: vec!["PASS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["1/1".to_string()]),
                    ("DP".to_string(), vec!["68".to_string()]),
                    ("COV".to_string(), vec!["0".to_string(), "68".to_string()]),
                    ("GT_CONF".to_string(), vec!["613.77".to_string()]),
                ]),
                is_filter_pass: true,
            },
            VCFRow {
                // 1
                position: 4725,
                reference: "t".to_string(),
                alternative: vec!["c".to_string()],
                filter: vec![],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["1/1".to_string()]),
                    ("DP".to_string(), vec!["68".to_string()]),
                    ("COV".to_string(), vec!["0".to_string(), "68".to_string()]),
                    ("GT_CONF".to_string(), vec!["613.77".to_string()]),
                ]),
                is_filter_pass: false,
            },
            VCFRow {
                // 2
                position: 4730,
                reference: "c".to_string(),
                alternative: vec!["t".to_string(), "g".to_string()],
                filter: vec!["PASS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["1/2".to_string()]),
                    ("DP".to_string(), vec!["200".to_string()]),
                    (
                        "COV".to_string(),
                        vec!["1".to_string(), "99".to_string(), "100".to_string()],
                    ),
                    ("GT_CONF".to_string(), vec!["613.77".to_string()]),
                ]),
                is_filter_pass: true,
            },
            VCFRow {
                // 3
                position: 4735,
                reference: "g".to_string(),
                alternative: vec!["gcc".to_string()],
                filter: vec!["PASS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["1/1".to_string()]),
                    ("DP".to_string(), vec!["68".to_string()]),
                    ("COV".to_string(), vec!["0".to_string(), "68".to_string()]),
                    ("GT_CONF".to_string(), vec!["63.77".to_string()]),
                ]),
                is_filter_pass: true,
            },
            VCFRow {
                // 4
                position: 4740,
                reference: "c".to_string(),
                alternative: vec!["gtt".to_string()],
                filter: vec![],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["1/1".to_string()]),
                    ("DP".to_string(), vec!["68".to_string()]),
                    ("COV".to_string(), vec!["0".to_string(), "68".to_string()]),
                    ("GT_CONF".to_string(), vec!["63.77".to_string()]),
                ]),
                is_filter_pass: false,
            },
            VCFRow {
                // 5
                position: 13148,
                reference: "t".to_string(),
                alternative: vec!["g".to_string()],
                filter: vec![],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["./.".to_string()]),
                    ("DP".to_string(), vec!["68".to_string()]),
                    ("COV".to_string(), vec!["0".to_string(), "68".to_string()]),
                    ("GT_CONF".to_string(), vec!["613.77".to_string()]),
                ]),
                is_filter_pass: false,
            },
            VCFRow {
                // 6
                position: 13149,
                reference: "g".to_string(),
                alternative: vec!["t".to_string()],
                filter: vec!["PASS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["./.".to_string()]),
                    ("DP".to_string(), vec!["68".to_string()]),
                    ("COV".to_string(), vec!["0".to_string(), "68".to_string()]),
                    ("GT_CONF".to_string(), vec!["613.77".to_string()]),
                ]),
                is_filter_pass: true,
            },
            VCFRow {
                // 7
                position: 13150,
                reference: "a".to_string(),
                alternative: vec!["tcg".to_string()],
                filter: vec!["PASS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["./.".to_string()]),
                    ("DP".to_string(), vec!["68".to_string()]),
                    ("COV".to_string(), vec!["0".to_string(), "68".to_string()]),
                    ("GT_CONF".to_string(), vec!["613.77".to_string()]),
                ]),
                is_filter_pass: true,
            }, // Latter rows omitted as are all filter fails
        ];

        for (idx, row) in expected_vcf_rows.iter().enumerate() {
            assert_eq!(*row, sample.get_vcf_row(idx));
        }

        let expected_genome_variants = vec![
            Variant {
                variant: "4687a>c".to_string(),
                nucleotide_index: 4687,
                evidence: 0,
                vcf_idx: Some(1),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(1474),
                codon_idx: Some(2),
                gene_name: Some("orf1ab".to_string()),
            },
            Variant {
                variant: "4730t>z".to_string(),
                nucleotide_index: 4730,
                evidence: 2,
                vcf_idx: None,
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(1489),
                codon_idx: Some(0),
                gene_name: Some("orf1ab".to_string()),
            },
            Variant {
                variant: "4735_ins_cc".to_string(),
                nucleotide_index: 4735,
                evidence: 3,
                vcf_idx: Some(1),
                indel_length: 2,
                indel_nucleotides: Some("cc".to_string()),
                gene_position: Some(4470),
                codon_idx: Some(2),
                gene_name: Some("orf1ab".to_string()),
            },
            Variant {
                variant: "13148g>x".to_string(),
                nucleotide_index: 13148,
                evidence: 5,
                vcf_idx: None,
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(4295),
                codon_idx: Some(0),
                gene_name: Some("orf1ab".to_string()),
            },
            Variant {
                variant: "13149t>x".to_string(),
                nucleotide_index: 13149,
                evidence: 6,
                vcf_idx: None,
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(4295),
                codon_idx: Some(1),
                gene_name: Some("orf1ab".to_string()),
            },
            Variant {
                variant: "13150t>x".to_string(),
                nucleotide_index: 13150,
                evidence: 7,
                vcf_idx: None,
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(4295),
                codon_idx: Some(2),
                gene_name: Some("orf1ab".to_string()),
            },
        ];

        for (idx, variant) in genome_diff.variants.iter().enumerate() {
            assert_eq!(*variant, expected_genome_variants[idx]);
        }

        let expected_genome_minor_variants = vec![
            Variant {
                // Looks weird, but the reference in the genbank = 't' and in vcf is 'c'
                variant: "4730t>t:99".to_string(),
                nucleotide_index: 4730,
                evidence: 2,
                vcf_idx: Some(1),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(1489),
                codon_idx: Some(0),
                gene_name: Some("orf1ab".to_string()),
            },
            Variant {
                variant: "4730t>g:100".to_string(),
                nucleotide_index: 4730,
                evidence: 2,
                vcf_idx: Some(2),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(1489),
                codon_idx: Some(0),
                gene_name: Some("orf1ab".to_string()),
            },
        ];

        for (idx, variant) in genome_diff.minor_variants.iter().enumerate() {
            assert_eq!(*variant, expected_genome_minor_variants[idx]);
        }

        assert_eq!(
            sample.genes_with_mutations,
            HashSet::from(["orf1ab".to_string()])
        );

        let expected_orf1ab_mutations = vec![
            Mutation {
                mutation: "a4422c".to_string(),
                gene: "orf1ab".to_string(),
                evidence: vec![Evidence {
                    cov: Some(68),
                    frs: Some(ordered_float::OrderedFloat(1.0)),
                    genotype: "1/1".to_string(),
                    call_type: AltType::SNP,
                    vcf_row: 0,
                    reference: "t".to_string(),
                    alt: "c".to_string(),
                    genome_index: 4687,
                    is_minor: false,
                    vcf_idx: Some(1),
                }],
                gene_position: Some(4422),
                codes_protein: Some(true),
                ref_nucleotides: Some("a".to_string()),
                alt_nucleotides: Some("c".to_string()),
                nucleotide_number: Some(4422),
                nucleotide_index: Some(4687),
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: None,
                amino_acid_sequence: None,
            },
            Mutation {
                mutation: "T1474T".to_string(),
                gene: "orf1ab".to_string(),
                evidence: vec![Evidence {
                    cov: Some(68),
                    frs: Some(ordered_float::OrderedFloat(1.0)),
                    genotype: "1/1".to_string(),
                    call_type: AltType::SNP,
                    vcf_row: 0,
                    reference: "t".to_string(),
                    alt: "c".to_string(),
                    genome_index: 4687,
                    is_minor: false,
                    vcf_idx: Some(1),
                }],
                gene_position: Some(1474),
                codes_protein: Some(true),
                ref_nucleotides: Some("aca".to_string()),
                alt_nucleotides: Some("acc".to_string()),
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: Some(1474),
                amino_acid_sequence: Some('T'),
            },
            Mutation {
                mutation: "Y1489Z".to_string(),
                gene: "orf1ab".to_string(),
                evidence: vec![Evidence {
                    cov: None,
                    frs: None,
                    genotype: "1/2".to_string(),
                    call_type: AltType::HET,
                    vcf_row: 2,
                    reference: "c".to_string(),
                    alt: "z".to_string(),
                    genome_index: 4730,
                    is_minor: false,
                    vcf_idx: None,
                }],
                gene_position: Some(1489),
                codes_protein: Some(true),
                ref_nucleotides: Some("tat".to_string()),
                alt_nucleotides: Some("zat".to_string()),
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: Some(1489),
                amino_acid_sequence: Some('Z'),
            },
            Mutation {
                mutation: "4470_ins_cc".to_string(),
                gene: "orf1ab".to_string(),
                evidence: vec![Evidence {
                    cov: Some(68),
                    frs: Some(ordered_float::OrderedFloat(1.0)),
                    genotype: "1/1".to_string(),
                    call_type: AltType::INS,
                    vcf_row: 3,
                    reference: "g".to_string(),
                    alt: "cc".to_string(),
                    genome_index: 4735,
                    is_minor: false,
                    vcf_idx: Some(1),
                }],
                gene_position: Some(4470),
                codes_protein: Some(true),
                ref_nucleotides: None,
                alt_nucleotides: None,
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: Some(2),
                indel_nucleotides: Some("cc".to_string()),
                amino_acid_number: None,
                amino_acid_sequence: None,
            },
            Mutation {
                mutation: "V4295X".to_string(),
                gene: "orf1ab".to_string(),
                evidence: vec![
                    Evidence {
                        cov: None,
                        frs: None,
                        genotype: "./.".to_string(),
                        call_type: AltType::NULL,
                        reference: "t".to_string(),
                        alt: "x".to_string(),
                        genome_index: 13148,
                        is_minor: false,
                        vcf_row: 5,
                        vcf_idx: None,
                    },
                    Evidence {
                        cov: None,
                        frs: None,
                        genotype: "./.".to_string(),
                        call_type: AltType::NULL,
                        reference: "g".to_string(),
                        alt: "x".to_string(),
                        genome_index: 13149,
                        is_minor: false,
                        vcf_row: 6,
                        vcf_idx: None,
                    },
                    Evidence {
                        cov: None,
                        frs: None,
                        genotype: "./.".to_string(),
                        call_type: AltType::NULL,
                        reference: "a".to_string(),
                        alt: "x".to_string(),
                        genome_index: 13150,
                        is_minor: false,
                        vcf_row: 7,
                        vcf_idx: None,
                    },
                ],
                gene_position: Some(4295),
                codes_protein: Some(true),
                ref_nucleotides: Some("gtt".to_string()),
                alt_nucleotides: Some("xxx".to_string()),
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: Some(4295),
                amino_acid_sequence: Some('X'),
            },
        ];

        let gene_diff = GeneDifference::new(
            reference.get_gene("orf1ab".to_string()),
            sample.get_gene("orf1ab".to_string()),
            MinorType::COV,
        );
        for (idx, mutation) in gene_diff.mutations.iter().enumerate() {
            assert_eq!(*mutation, expected_orf1ab_mutations[idx]);
        }

        let expected_orf1ab_minor_mutations = vec![Mutation {
            mutation: "Y1489Z:100".to_string(),
            gene: "orf1ab".to_string(),
            evidence: vec![
                Evidence {
                    cov: Some(99),
                    frs: Some(ordered_float::OrderedFloat(0.495)),
                    genotype: "1/2".to_string(),
                    call_type: AltType::SNP,
                    vcf_row: 2,
                    reference: "c".to_string(),
                    alt: "t".to_string(),
                    genome_index: 4730,
                    is_minor: true,
                    vcf_idx: Some(1),
                },
                Evidence {
                    cov: Some(100),
                    frs: Some(ordered_float::OrderedFloat(0.5)),
                    genotype: "1/2".to_string(),
                    call_type: AltType::SNP,
                    vcf_row: 2,
                    reference: "c".to_string(),
                    alt: "g".to_string(),
                    genome_index: 4730,
                    is_minor: true,
                    vcf_idx: Some(2),
                },
            ],
            gene_position: Some(1489),
            codes_protein: Some(true),
            ref_nucleotides: Some("tat".to_string()),
            alt_nucleotides: Some("zat".to_string()),
            nucleotide_number: None,
            nucleotide_index: None,
            indel_length: None,
            indel_nucleotides: None,
            amino_acid_number: Some(1489),
            amino_acid_sequence: Some('Z'),
        }];

        for (idx, mutation) in gene_diff.minor_mutations.iter().enumerate() {
            assert_eq!(*mutation, expected_orf1ab_minor_mutations[idx]);
        }
    }

    #[test]
    fn test_mutate_test_dna() {
        let mut reference = Genome::new("reference/TEST-DNA.gbk");
        // Note that this VCF is completely dummy data
        // so variants aren't necessarily matching the reference in the VCF row
        let vcf = VCFFile::new("test/TEST-DNA.vcf".to_string(), false, 3);
        let mut sample = mutate(&reference, vcf);

        let genome_diff = GenomeDifference::new(reference.clone(), sample.clone(), MinorType::COV);

        let expected_vcf_rows = vec![
            VCFRow {
                // 0
                position: 2,
                reference: "a".to_string(),
                alternative: vec!["g".to_string()],
                filter: vec!["PASS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["./.".to_string()]),
                    ("DP".to_string(), vec!["2".to_string()]),
                    ("COV".to_string(), vec!["1".to_string(), "1".to_string()]),
                    ("GT_CONF".to_string(), vec!["2.05".to_string()]),
                ]),
                is_filter_pass: true,
            },
            VCFRow {
                // 1
                position: 4,
                reference: "a".to_string(),
                alternative: vec!["g".to_string(), "t".to_string()],
                filter: vec![
                    "MIN_GCP".to_string(),
                    "MIN_DP".to_string(),
                    "MIN_FRS".to_string(),
                ],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["./.".to_string()]),
                    ("DP".to_string(), vec!["4".to_string()]),
                    (
                        "COV".to_string(),
                        vec!["1".to_string(), "2".to_string(), "1".to_string()],
                    ),
                    ("GT_CONF".to_string(), vec!["3.77".to_string()]),
                ]),
                is_filter_pass: false,
            },
            VCFRow {
                // 2
                position: 6,
                reference: "aaa".to_string(),
                alternative: vec!["ggt".to_string(), "gta".to_string(), "ata".to_string()],
                filter: vec!["PASS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["./.".to_string()]),
                    ("DP".to_string(), vec!["4".to_string()]),
                    (
                        "COV".to_string(),
                        vec![
                            "1".to_string(),
                            "1".to_string(),
                            "1".to_string(),
                            "1".to_string(),
                        ],
                    ),
                    ("GT_CONF".to_string(), vec!["2.76".to_string()]),
                ]),
                is_filter_pass: true,
            },
            VCFRow {
                // 3
                position: 12,
                reference: "c".to_string(),
                alternative: vec!["t".to_string()],
                filter: vec!["PASS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["1/1".to_string()]),
                    ("DP".to_string(), vec!["50".to_string()]),
                    ("COV".to_string(), vec!["0".to_string(), "50".to_string()]),
                    ("GT_CONF".to_string(), vec!["200.58".to_string()]),
                ]),
                is_filter_pass: true,
            },
            VCFRow {
                // 4
                position: 14,
                reference: "c".to_string(),
                alternative: vec!["t".to_string(), "g".to_string()],
                filter: vec!["PASS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["2/2".to_string()]),
                    ("DP".to_string(), vec!["45".to_string()]),
                    (
                        "COV".to_string(),
                        vec!["0".to_string(), "2".to_string(), "43".to_string()],
                    ),
                    ("GT_CONF".to_string(), vec!["155.58".to_string()]),
                ]),
                is_filter_pass: true,
            },
            VCFRow {
                // 5
                position: 16,
                reference: "ccc".to_string(),
                alternative: vec!["tgc".to_string(), "gtg".to_string()],
                filter: vec!["PASS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["1/1".to_string()]),
                    ("DP".to_string(), vec!["70".to_string()]),
                    (
                        "COV".to_string(),
                        vec!["0".to_string(), "68".to_string(), "8".to_string()],
                    ),
                    ("GT_CONF".to_string(), vec!["300.25".to_string()]),
                ]),
                is_filter_pass: true,
            },
            VCFRow {
                // 6
                position: 22,
                reference: "g".to_string(),
                alternative: vec!["t".to_string(), "c".to_string(), "a".to_string()],
                filter: vec!["PASS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["1/2".to_string()]),
                    ("DP".to_string(), vec!["202".to_string()]),
                    (
                        "COV".to_string(),
                        vec![
                            "1".to_string(),
                            "99".to_string(),
                            "100".to_string(),
                            "2".to_string(),
                        ],
                    ),
                    ("GT_CONF".to_string(), vec!["613.77".to_string()]),
                ]),
                is_filter_pass: true,
            },
            VCFRow {
                // 7
                position: 24,
                reference: "g".to_string(),
                alternative: vec!["t".to_string(), "c".to_string(), "a".to_string()],
                filter: vec!["PASS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["0/2".to_string()]),
                    ("DP".to_string(), vec!["202".to_string()]),
                    (
                        "COV".to_string(),
                        vec![
                            "99".to_string(),
                            "1".to_string(),
                            "100".to_string(),
                            "2".to_string(),
                        ],
                    ),
                    ("GT_CONF".to_string(), vec!["613.77".to_string()]),
                ]),
                is_filter_pass: true,
            },
            VCFRow {
                // 8
                position: 26,
                reference: "gg".to_string(),
                alternative: vec!["aa".to_string(), "ct".to_string(), "at".to_string()],
                filter: vec!["PASS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["1/2".to_string()]),
                    ("DP".to_string(), vec!["100".to_string()]),
                    (
                        "COV".to_string(),
                        vec![
                            "0".to_string(),
                            "48".to_string(),
                            "50".to_string(),
                            "2".to_string(),
                        ],
                    ),
                    ("GT_CONF".to_string(), vec!["475.54".to_string()]),
                ]),
                is_filter_pass: true,
            },
            VCFRow {
                // 9
                position: 28,
                reference: "gg".to_string(),
                alternative: vec!["aa".to_string(), "t".to_string(), "a".to_string()],
                filter: vec!["PASS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["1/3".to_string()]),
                    ("DP".to_string(), vec!["100".to_string()]),
                    (
                        "COV".to_string(),
                        vec![
                            "0".to_string(),
                            "48".to_string(),
                            "2".to_string(),
                            "50".to_string(),
                        ],
                    ),
                    ("GT_CONF".to_string(), vec!["315.11".to_string()]),
                ]),
                is_filter_pass: true,
            },
            VCFRow {
                // 10
                position: 33,
                reference: "t".to_string(),
                alternative: vec!["ttt".to_string()],
                filter: vec!["PASS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["1/1".to_string()]),
                    ("DP".to_string(), vec!["200".to_string()]),
                    ("COV".to_string(), vec!["1".to_string(), "199".to_string()]),
                    ("GT_CONF".to_string(), vec!["145.21".to_string()]),
                ]),
                is_filter_pass: true,
            },
            VCFRow {
                // 11
                position: 36,
                reference: "tt".to_string(),
                alternative: vec!["t".to_string()],
                filter: vec!["PASS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["1/1".to_string()]),
                    ("DP".to_string(), vec!["200".to_string()]),
                    ("COV".to_string(), vec!["1".to_string(), "199".to_string()]),
                    ("GT_CONF".to_string(), vec!["145.21".to_string()]),
                ]),
                is_filter_pass: true,
            },
            VCFRow {
                // 12
                position: 39,
                reference: "tt".to_string(),
                alternative: vec!["agt".to_string()],
                filter: vec!["PASS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["1/1".to_string()]),
                    ("DP".to_string(), vec!["200".to_string()]),
                    ("COV".to_string(), vec!["1".to_string(), "199".to_string()]),
                    ("GT_CONF".to_string(), vec!["145.21".to_string()]),
                ]),
                is_filter_pass: true,
            },
            VCFRow {
                // 13
                position: 65,
                reference: "gg".to_string(),
                alternative: vec!["cagg".to_string()],
                filter: vec!["PASS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["1/1".to_string()]),
                    ("DP".to_string(), vec!["200".to_string()]),
                    ("COV".to_string(), vec!["1".to_string(), "199".to_string()]),
                    ("GT_CONF".to_string(), vec!["145.21".to_string()]),
                ]),
                is_filter_pass: true,
            },
            VCFRow {
                // 14
                position: 69,
                reference: "gg".to_string(),
                alternative: vec!["gg".to_string()],
                filter: vec!["PASS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["0/0".to_string()]),
                    ("DP".to_string(), vec!["200".to_string()]),
                    ("COV".to_string(), vec!["1".to_string(), "199".to_string()]),
                    ("GT_CONF".to_string(), vec!["145.21".to_string()]),
                ]),
                is_filter_pass: true,
            },
            VCFRow {
                // 15
                position: 73,
                reference: "t".to_string(),
                alternative: vec!["ta".to_string(), "at".to_string()],
                filter: vec!["PASS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["1/1".to_string()]),
                    ("DP".to_string(), vec!["200".to_string()]),
                    (
                        "COV".to_string(),
                        vec!["1".to_string(), "198".to_string(), "1".to_string()],
                    ),
                    ("GT_CONF".to_string(), vec!["145.21".to_string()]),
                ]),
                is_filter_pass: true,
            },
        ];

        for (idx, row) in expected_vcf_rows.iter().enumerate() {
            assert_eq!(*row, sample.get_vcf_row(idx));
        }

        let expected_variants = vec![
            Variant {
                variant: "2a>x".to_string(),
                nucleotide_index: 2,
                evidence: 0,
                vcf_idx: None,
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(-2),
                codon_idx: None,
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "4a>x".to_string(),
                nucleotide_index: 4,
                evidence: 1,
                vcf_idx: None,
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(1),
                codon_idx: Some(0),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "6a>x".to_string(),
                nucleotide_index: 6,
                evidence: 2,
                vcf_idx: None,
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(1),
                codon_idx: Some(2),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "7a>x".to_string(),
                nucleotide_index: 7,
                evidence: 2,
                vcf_idx: None,
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(2),
                codon_idx: Some(0),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "8a>x".to_string(),
                nucleotide_index: 8,
                evidence: 2,
                vcf_idx: None,
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(2),
                codon_idx: Some(1),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "12c>t".to_string(),
                nucleotide_index: 12,
                evidence: 3,
                vcf_idx: Some(1),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(3),
                codon_idx: Some(2),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "14c>g".to_string(),
                nucleotide_index: 14,
                evidence: 4,
                vcf_idx: Some(2),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(4),
                codon_idx: Some(1),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "16c>t".to_string(),
                nucleotide_index: 16,
                evidence: 5,
                vcf_idx: Some(1),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(5),
                codon_idx: Some(0),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "17c>g".to_string(),
                nucleotide_index: 17,
                evidence: 5,
                vcf_idx: Some(1),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(5),
                codon_idx: Some(1),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "22g>z".to_string(),
                nucleotide_index: 22,
                evidence: 6,
                vcf_idx: None,
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(7),
                codon_idx: Some(0),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "24g>z".to_string(),
                nucleotide_index: 24,
                evidence: 7,
                vcf_idx: None,
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(7),
                codon_idx: Some(2),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "26g>z".to_string(),
                nucleotide_index: 26,
                evidence: 8,
                vcf_idx: None,
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(8),
                codon_idx: Some(1),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "27g>z".to_string(),
                nucleotide_index: 27,
                evidence: 8,
                vcf_idx: None,
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(8),
                codon_idx: Some(2),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "28g>z".to_string(),
                nucleotide_index: 28,
                evidence: 9,
                vcf_idx: None,
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(9),
                codon_idx: Some(0),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "28g>z".to_string(),
                nucleotide_index: 28,
                evidence: 9,
                vcf_idx: None,
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(1),
                codon_idx: Some(0),
                gene_name: Some("B".to_string()),
            },
            Variant {
                variant: "29g>z".to_string(),
                nucleotide_index: 29,
                evidence: 9,
                vcf_idx: None,
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(9),
                codon_idx: Some(1),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "29g>z".to_string(),
                nucleotide_index: 29,
                evidence: 9,
                vcf_idx: None,
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(1),
                codon_idx: Some(1),
                gene_name: Some("B".to_string()),
            },
            Variant {
                variant: "33_ins_tt".to_string(),
                nucleotide_index: 33,
                evidence: 10,
                vcf_idx: Some(1),
                indel_length: 2,
                indel_nucleotides: Some("tt".to_string()),
                gene_position: Some(6),
                codon_idx: Some(2),
                gene_name: Some("B".to_string()),
            },
            Variant {
                variant: "37_del_t".to_string(),
                nucleotide_index: 37,
                evidence: 11,
                vcf_idx: Some(1),
                indel_length: -1,
                indel_nucleotides: Some("t".to_string()),
                gene_position: Some(10),
                codon_idx: Some(0),
                gene_name: Some("B".to_string()),
            },
            Variant {
                variant: "39_ins_g".to_string(),
                nucleotide_index: 39,
                evidence: 12,
                vcf_idx: Some(1),
                indel_length: 1,
                indel_nucleotides: Some("g".to_string()),
                gene_position: Some(12),
                codon_idx: Some(2),
                gene_name: Some("B".to_string()),
            },
            Variant {
                variant: "39t>a".to_string(),
                nucleotide_index: 39,
                evidence: 12,
                vcf_idx: Some(1),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(4),
                codon_idx: Some(2),
                gene_name: Some("B".to_string()),
            },
            Variant {
                variant: "64_ins_ca".to_string(),
                nucleotide_index: 64,
                evidence: 13,
                vcf_idx: Some(1),
                indel_length: 2,
                indel_nucleotides: Some("ca".to_string()),
                gene_position: None,
                codon_idx: None,
                gene_name: None,
            },
            Variant {
                variant: "69g>x".to_string(),
                nucleotide_index: 69,
                evidence: 14,
                vcf_idx: Some(0),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: None,
                codon_idx: None,
                gene_name: None,
            },
            Variant {
                variant: "70g>x".to_string(),
                nucleotide_index: 70,
                evidence: 14,
                vcf_idx: Some(0),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: None,
                codon_idx: None,
                gene_name: None,
            },
            Variant {
                variant: "73_ins_a".to_string(),
                nucleotide_index: 73,
                evidence: 15,
                vcf_idx: Some(1),
                indel_length: 1,
                indel_nucleotides: Some("a".to_string()),
                gene_position: None,
                codon_idx: None,
                gene_name: None,
            },
        ];

        for (idx, variant) in genome_diff.variants.iter().enumerate() {
            assert_eq!(variant, &expected_variants[idx])
        }

        let expected_minor_variants = vec![
            Variant {
                variant: "16c>g:8".to_string(),
                nucleotide_index: 16,
                evidence: 5,
                vcf_idx: Some(2),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(5),
                codon_idx: Some(0),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "17c>t:8".to_string(),
                nucleotide_index: 17,
                evidence: 5,
                vcf_idx: Some(2),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(5),
                codon_idx: Some(1),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "18c>g:8".to_string(),
                nucleotide_index: 18,
                evidence: 5,
                vcf_idx: Some(2),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(5),
                codon_idx: Some(2),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "22g>t:99".to_string(),
                nucleotide_index: 22,
                evidence: 6,
                vcf_idx: Some(1),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(7),
                codon_idx: Some(0),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "22g>c:100".to_string(),
                nucleotide_index: 22,
                evidence: 6,
                vcf_idx: Some(2),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(7),
                codon_idx: Some(0),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "24g>c:100".to_string(),
                nucleotide_index: 24,
                evidence: 7,
                vcf_idx: Some(2),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(7),
                codon_idx: Some(2),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "26g>a:48".to_string(),
                nucleotide_index: 26,
                evidence: 8,
                vcf_idx: Some(1),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(8),
                codon_idx: Some(1),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "26g>c:50".to_string(),
                nucleotide_index: 26,
                evidence: 8,
                vcf_idx: Some(2),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(8),
                codon_idx: Some(1),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "27g>a:48".to_string(),
                nucleotide_index: 27,
                evidence: 8,
                vcf_idx: Some(1),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(8),
                codon_idx: Some(2),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "27g>t:50".to_string(),
                nucleotide_index: 27,
                evidence: 8,
                vcf_idx: Some(2),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(8),
                codon_idx: Some(2),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "28g>a:48".to_string(),
                nucleotide_index: 28,
                evidence: 9,
                vcf_idx: Some(1),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(9),
                codon_idx: Some(0),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "28g>a:48".to_string(),
                nucleotide_index: 28,
                evidence: 9,
                vcf_idx: Some(1),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(1),
                codon_idx: Some(0),
                gene_name: Some("B".to_string()),
            },
            Variant {
                variant: "28g>a:50".to_string(),
                nucleotide_index: 28,
                evidence: 9,
                vcf_idx: Some(3),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(9),
                codon_idx: Some(0),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "28g>a:50".to_string(),
                nucleotide_index: 28,
                evidence: 9,
                vcf_idx: Some(3),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(1),
                codon_idx: Some(0),
                gene_name: Some("B".to_string()),
            },
            Variant {
                variant: "29g>a:48".to_string(),
                nucleotide_index: 29,
                evidence: 9,
                vcf_idx: Some(1),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(9),
                codon_idx: Some(1),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "29g>a:48".to_string(),
                nucleotide_index: 29,
                evidence: 9,
                vcf_idx: Some(1),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(1),
                codon_idx: Some(1),
                gene_name: Some("B".to_string()),
            },
            Variant {
                variant: "29_del_g:50".to_string(),
                nucleotide_index: 29,
                evidence: 9,
                vcf_idx: Some(3),
                indel_length: -1,
                indel_nucleotides: Some("g".to_string()),
                gene_position: Some(26),
                codon_idx: Some(1),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "29_del_g:50".to_string(),
                nucleotide_index: 29,
                evidence: 9,
                vcf_idx: Some(3),
                indel_length: -1,
                indel_nucleotides: Some("g".to_string()),
                gene_position: Some(2),
                codon_idx: Some(1),
                gene_name: Some("B".to_string()),
            },
        ];
        for (idx, variant) in genome_diff.minor_variants.iter().enumerate() {
            assert_eq!(variant, &expected_minor_variants[idx])
        }

        let genome_diff = GenomeDifference::new(reference.clone(), sample.clone(), MinorType::FRS);
        let expected_minor_variants = vec![
            Variant {
                variant: "16c>g:0.105".to_string(),
                nucleotide_index: 16,
                evidence: 5,
                vcf_idx: Some(2),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(5),
                codon_idx: Some(0),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "17c>t:0.105".to_string(),
                nucleotide_index: 17,
                evidence: 5,
                vcf_idx: Some(2),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(5),
                codon_idx: Some(1),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "18c>g:0.105".to_string(),
                nucleotide_index: 18,
                evidence: 5,
                vcf_idx: Some(2),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(5),
                codon_idx: Some(2),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "22g>t:0.49".to_string(),
                nucleotide_index: 22,
                evidence: 6,
                vcf_idx: Some(1),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(7),
                codon_idx: Some(0),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "22g>c:0.495".to_string(),
                nucleotide_index: 22,
                evidence: 6,
                vcf_idx: Some(2),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(7),
                codon_idx: Some(0),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "24g>c:0.495".to_string(),
                nucleotide_index: 24,
                evidence: 7,
                vcf_idx: Some(2),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(7),
                codon_idx: Some(2),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "26g>a:0.48".to_string(),
                nucleotide_index: 26,
                evidence: 8,
                vcf_idx: Some(1),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(8),
                codon_idx: Some(1),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "26g>c:0.5".to_string(),
                nucleotide_index: 26,
                evidence: 8,
                vcf_idx: Some(2),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(8),
                codon_idx: Some(1),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "27g>a:0.48".to_string(),
                nucleotide_index: 27,
                evidence: 8,
                vcf_idx: Some(1),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(8),
                codon_idx: Some(2),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "27g>t:0.5".to_string(),
                nucleotide_index: 27,
                evidence: 8,
                vcf_idx: Some(2),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(8),
                codon_idx: Some(2),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "28g>a:0.48".to_string(),
                nucleotide_index: 28,
                evidence: 9,
                vcf_idx: Some(1),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(9),
                codon_idx: Some(0),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "28g>a:0.48".to_string(),
                nucleotide_index: 28,
                evidence: 9,
                vcf_idx: Some(1),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(1),
                codon_idx: Some(0),
                gene_name: Some("B".to_string()),
            },
            Variant {
                variant: "28g>a:0.5".to_string(),
                nucleotide_index: 28,
                evidence: 9,
                vcf_idx: Some(3),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(9),
                codon_idx: Some(0),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "28g>a:0.5".to_string(),
                nucleotide_index: 28,
                evidence: 9,
                vcf_idx: Some(3),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(1),
                codon_idx: Some(0),
                gene_name: Some("B".to_string()),
            },
            Variant {
                variant: "29g>a:0.48".to_string(),
                nucleotide_index: 29,
                evidence: 9,
                vcf_idx: Some(1),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(9),
                codon_idx: Some(1),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "29g>a:0.48".to_string(),
                nucleotide_index: 29,
                evidence: 9,
                vcf_idx: Some(1),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(1),
                codon_idx: Some(1),
                gene_name: Some("B".to_string()),
            },
            Variant {
                variant: "29_del_g:0.5".to_string(),
                nucleotide_index: 29,
                evidence: 9,
                vcf_idx: Some(3),
                indel_length: -1,
                indel_nucleotides: Some("g".to_string()),
                gene_position: Some(26),
                codon_idx: Some(1),
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "29_del_g:0.5".to_string(),
                nucleotide_index: 29,
                evidence: 9,
                vcf_idx: Some(3),
                indel_length: -1,
                indel_nucleotides: Some("g".to_string()),
                gene_position: Some(2),
                codon_idx: Some(1),
                gene_name: Some("B".to_string()),
            },
        ];
        for (idx, variant) in genome_diff.minor_variants.iter().enumerate() {
            assert_eq!(variant, &expected_minor_variants[idx])
        }

        let expected_a_mutations = vec![
            Mutation {
                mutation: "a-2x".to_string(),
                gene: "A".to_string(),
                evidence: vec![Evidence {
                    cov: None,
                    frs: None,
                    genotype: "./.".to_string(),
                    call_type: AltType::NULL,
                    vcf_row: 0,
                    reference: "a".to_string(),
                    alt: "x".to_string(),
                    genome_index: 2,
                    is_minor: false,
                    vcf_idx: None,
                }],
                gene_position: Some(-2),
                codes_protein: Some(false),
                ref_nucleotides: Some("a".to_string()),
                alt_nucleotides: Some("x".to_string()),
                nucleotide_number: Some(-2),
                nucleotide_index: Some(2),
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: None,
                amino_acid_sequence: None,
            },
            Mutation {
                mutation: "K1X".to_string(),
                gene: "A".to_string(),
                evidence: vec![
                    Evidence {
                        cov: None,
                        frs: None,
                        genotype: "./.".to_string(),
                        call_type: AltType::NULL,
                        vcf_row: 1,
                        reference: "a".to_string(),
                        alt: "x".to_string(),
                        genome_index: 4,
                        is_minor: false,
                        vcf_idx: None,
                    },
                    Evidence {
                        cov: None,
                        frs: None,
                        genotype: "./.".to_string(),
                        call_type: AltType::NULL,
                        vcf_row: 2,
                        reference: "a".to_string(),
                        alt: "x".to_string(),
                        genome_index: 6,
                        is_minor: false,
                        vcf_idx: None,
                    },
                ],
                gene_position: Some(1),
                codes_protein: Some(true),
                ref_nucleotides: Some("aaa".to_string()),
                alt_nucleotides: Some("xax".to_string()),
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: Some(1),
                amino_acid_sequence: Some('X'),
            },
            Mutation {
                mutation: "K2X".to_string(),
                gene: "A".to_string(),
                evidence: vec![
                    Evidence {
                        cov: None,
                        frs: None,
                        genotype: "./.".to_string(),
                        call_type: AltType::NULL,
                        vcf_row: 2,
                        reference: "a".to_string(),
                        alt: "x".to_string(),
                        genome_index: 7,
                        is_minor: false,
                        vcf_idx: None,
                    },
                    Evidence {
                        cov: None,
                        frs: None,
                        genotype: "./.".to_string(),
                        call_type: AltType::NULL,
                        vcf_row: 2,
                        reference: "a".to_string(),
                        alt: "x".to_string(),
                        genome_index: 8,
                        is_minor: false,
                        vcf_idx: None,
                    },
                ],
                gene_position: Some(2),
                codes_protein: Some(true),
                ref_nucleotides: Some("aaa".to_string()),
                alt_nucleotides: Some("xxa".to_string()),
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: Some(2),
                amino_acid_sequence: Some('X'),
            },
            Mutation {
                mutation: "c9t".to_string(),
                gene: "A".to_string(),
                evidence: vec![Evidence {
                    cov: Some(50),
                    frs: Some(OrderedFloat(1.0)),
                    genotype: "1/1".to_string(),
                    call_type: AltType::SNP,
                    vcf_row: 3,
                    reference: "c".to_string(),
                    alt: "t".to_string(),
                    genome_index: 12,
                    is_minor: false,
                    vcf_idx: Some(1),
                }],
                gene_position: Some(9),
                codes_protein: Some(true),
                ref_nucleotides: Some("c".to_string()),
                alt_nucleotides: Some("t".to_string()),
                nucleotide_number: Some(9),
                nucleotide_index: Some(12),
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: None,
                amino_acid_sequence: None,
            },
            Mutation {
                mutation: "T3T".to_string(),
                gene: "A".to_string(),
                evidence: vec![Evidence {
                    cov: Some(50),
                    frs: Some(OrderedFloat(1.0)),
                    genotype: "1/1".to_string(),
                    call_type: AltType::SNP,
                    vcf_row: 3,
                    reference: "c".to_string(),
                    alt: "t".to_string(),
                    genome_index: 12,
                    is_minor: false,
                    vcf_idx: Some(1),
                }],
                gene_position: Some(3),
                codes_protein: Some(true),
                ref_nucleotides: Some("acc".to_string()),
                alt_nucleotides: Some("act".to_string()),
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: Some(3),
                amino_acid_sequence: Some('T'),
            },
            Mutation {
                mutation: "P4R".to_string(),
                gene: "A".to_string(),
                evidence: vec![Evidence {
                    cov: Some(43),
                    frs: Some(OrderedFloat(43.0 / 45.0)),
                    genotype: "2/2".to_string(),
                    call_type: AltType::SNP,
                    vcf_row: 4,
                    reference: "c".to_string(),
                    alt: "g".to_string(),
                    genome_index: 14,
                    is_minor: false,
                    vcf_idx: Some(2),
                }],
                gene_position: Some(4),
                codes_protein: Some(true),
                ref_nucleotides: Some("ccc".to_string()),
                alt_nucleotides: Some("cgc".to_string()),
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: Some(4),
                amino_acid_sequence: Some('R'),
            },
            Mutation {
                mutation: "P5C".to_string(),
                gene: "A".to_string(),
                evidence: vec![
                    Evidence {
                        cov: Some(68),
                        frs: Some(OrderedFloat(68.0 / 76.0)),
                        genotype: "1/1".to_string(),
                        call_type: AltType::SNP,
                        vcf_row: 5,
                        reference: "c".to_string(),
                        alt: "t".to_string(),
                        genome_index: 16,
                        is_minor: false,
                        vcf_idx: Some(1),
                    },
                    Evidence {
                        cov: Some(68),
                        frs: Some(OrderedFloat(68.0 / 76.0)),
                        genotype: "1/1".to_string(),
                        call_type: AltType::SNP,
                        vcf_row: 5,
                        reference: "c".to_string(),
                        alt: "g".to_string(),
                        genome_index: 17,
                        is_minor: false,
                        vcf_idx: Some(1),
                    },
                ],
                gene_position: Some(5),
                codes_protein: Some(true),
                ref_nucleotides: Some("ccc".to_string()),
                alt_nucleotides: Some("tgc".to_string()),
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: Some(5),
                amino_acid_sequence: Some('C'),
            },
            Mutation {
                mutation: "G7Z".to_string(),
                gene: "A".to_string(),
                evidence: vec![
                    Evidence {
                        cov: None,
                        frs: None,
                        genotype: "1/2".to_string(),
                        call_type: AltType::HET,
                        vcf_row: 6,
                        reference: "g".to_string(),
                        alt: "z".to_string(),
                        genome_index: 22,
                        is_minor: false,
                        vcf_idx: None,
                    },
                    Evidence {
                        cov: None,
                        frs: None,
                        genotype: "0/2".to_string(),
                        call_type: AltType::HET,
                        vcf_row: 7,
                        reference: "g".to_string(),
                        alt: "z".to_string(),
                        genome_index: 24,
                        is_minor: false,
                        vcf_idx: None,
                    },
                ],
                gene_position: Some(7),
                codes_protein: Some(true),
                ref_nucleotides: Some("ggg".to_string()),
                alt_nucleotides: Some("zgz".to_string()),
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: Some(7),
                amino_acid_sequence: Some('Z'),
            },
            Mutation {
                mutation: "G8Z".to_string(),
                gene: "A".to_string(),
                evidence: vec![
                    Evidence {
                        cov: None,
                        frs: None,
                        genotype: "1/2".to_string(),
                        call_type: AltType::HET,
                        vcf_row: 8,
                        reference: "g".to_string(),
                        alt: "z".to_string(),
                        genome_index: 26,
                        is_minor: false,
                        vcf_idx: None,
                    },
                    Evidence {
                        cov: None,
                        frs: None,
                        genotype: "1/2".to_string(),
                        call_type: AltType::HET,
                        vcf_row: 8,
                        reference: "g".to_string(),
                        alt: "z".to_string(),
                        genome_index: 27,
                        is_minor: false,
                        vcf_idx: None,
                    },
                ],
                gene_position: Some(8),
                codes_protein: Some(true),
                ref_nucleotides: Some("ggg".to_string()),
                alt_nucleotides: Some("gzz".to_string()),
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: Some(8),
                amino_acid_sequence: Some('Z'),
            },
            Mutation {
                mutation: "G9Z".to_string(),
                gene: "A".to_string(),
                evidence: vec![
                    Evidence {
                        cov: None,
                        frs: None,
                        genotype: "1/3".to_string(),
                        call_type: AltType::HET,
                        vcf_row: 9,
                        reference: "g".to_string(),
                        alt: "z".to_string(),
                        genome_index: 28,
                        is_minor: false,
                        vcf_idx: None,
                    },
                    Evidence {
                        cov: None,
                        frs: None,
                        genotype: "1/3".to_string(),
                        call_type: AltType::HET,
                        vcf_row: 9,
                        reference: "g".to_string(),
                        alt: "z".to_string(),
                        genome_index: 29,
                        is_minor: false,
                        vcf_idx: None,
                    },
                ],
                gene_position: Some(9),
                codes_protein: Some(true),
                ref_nucleotides: Some("ggg".to_string()),
                alt_nucleotides: Some("zzg".to_string()),
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: Some(9),
                amino_acid_sequence: Some('Z'),
            },
        ];
        let a_diff = GeneDifference::new(
            reference.get_gene("A".to_string()),
            sample.get_gene("A".to_string()),
            MinorType::COV,
        );

        for (idx, mutation) in a_diff.mutations.iter().enumerate() {
            assert_eq!(mutation, &expected_a_mutations[idx])
        }

        let expected_a_minor_mutations = vec![
            Mutation {
                mutation: "P5V:8".to_string(),
                gene: "A".to_string(),
                evidence: vec![
                    Evidence {
                        cov: Some(8),
                        frs: Some(OrderedFloat(8.0 / 76.0)),
                        genotype: "1/1".to_string(),
                        call_type: AltType::SNP,
                        vcf_row: 5,
                        reference: "c".to_string(),
                        alt: "g".to_string(),
                        genome_index: 16,
                        is_minor: true,
                        vcf_idx: Some(2),
                    },
                    Evidence {
                        cov: Some(8),
                        frs: Some(OrderedFloat(8.0 / 76.0)),
                        genotype: "1/1".to_string(),
                        call_type: AltType::SNP,
                        vcf_row: 5,
                        reference: "c".to_string(),
                        alt: "t".to_string(),
                        genome_index: 17,
                        is_minor: true,
                        vcf_idx: Some(2),
                    },
                    Evidence {
                        cov: Some(8),
                        frs: Some(OrderedFloat(8.0 / 76.0)),
                        genotype: "1/1".to_string(),
                        call_type: AltType::SNP,
                        vcf_row: 5,
                        reference: "c".to_string(),
                        alt: "g".to_string(),
                        genome_index: 18,
                        is_minor: true,
                        vcf_idx: Some(2),
                    },
                ],
                gene_position: Some(5),
                codes_protein: Some(true),
                ref_nucleotides: Some("ccc".to_string()),
                alt_nucleotides: Some("gtg".to_string()),
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: Some(5),
                amino_acid_sequence: Some('V'),
            },
            Mutation {
                mutation: "G7Z:100".to_string(),
                gene: "A".to_string(),
                evidence: vec![
                    Evidence {
                        cov: Some(99),
                        frs: Some(OrderedFloat(99.0 / 202.0)),
                        genotype: "1/2".to_string(),
                        call_type: AltType::SNP,
                        vcf_row: 6,
                        reference: "g".to_string(),
                        alt: "t".to_string(),
                        genome_index: 22,
                        is_minor: true,
                        vcf_idx: Some(1),
                    },
                    Evidence {
                        cov: Some(100),
                        frs: Some(OrderedFloat(100.0 / 202.0)),
                        genotype: "1/2".to_string(),
                        call_type: AltType::SNP,
                        vcf_row: 6,
                        reference: "g".to_string(),
                        alt: "c".to_string(),
                        genome_index: 22,
                        is_minor: true,
                        vcf_idx: Some(2),
                    },
                    Evidence {
                        cov: Some(100),
                        frs: Some(OrderedFloat(100.0 / 202.0)),
                        genotype: "0/2".to_string(),
                        call_type: AltType::SNP,
                        vcf_row: 7,
                        reference: "g".to_string(),
                        alt: "c".to_string(),
                        genome_index: 24,
                        is_minor: true,
                        vcf_idx: Some(2),
                    },
                ],
                gene_position: Some(7),
                codes_protein: Some(true),
                ref_nucleotides: Some("ggg".to_string()),
                alt_nucleotides: Some("zgc".to_string()),
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: Some(7),
                amino_acid_sequence: Some('Z'),
            },
            Mutation {
                mutation: "G8Z:50".to_string(),
                gene: "A".to_string(),
                evidence: vec![
                    Evidence {
                        cov: Some(48),
                        frs: Some(OrderedFloat(48.0 / 100.0)),
                        genotype: "1/2".to_string(),
                        call_type: AltType::SNP,
                        vcf_row: 8,
                        reference: "g".to_string(),
                        alt: "a".to_string(),
                        genome_index: 26,
                        is_minor: true,
                        vcf_idx: Some(1),
                    },
                    Evidence {
                        cov: Some(50),
                        frs: Some(OrderedFloat(50.0 / 100.0)),
                        genotype: "1/2".to_string(),
                        call_type: AltType::SNP,
                        vcf_row: 8,
                        reference: "g".to_string(),
                        alt: "c".to_string(),
                        genome_index: 26,
                        is_minor: true,
                        vcf_idx: Some(2),
                    },
                    Evidence {
                        cov: Some(48),
                        frs: Some(OrderedFloat(48.0 / 100.0)),
                        genotype: "1/2".to_string(),
                        call_type: AltType::SNP,
                        vcf_row: 8,
                        reference: "g".to_string(),
                        alt: "a".to_string(),
                        genome_index: 27,
                        is_minor: true,
                        vcf_idx: Some(1),
                    },
                    Evidence {
                        cov: Some(50),
                        frs: Some(OrderedFloat(50.0 / 100.0)),
                        genotype: "1/2".to_string(),
                        call_type: AltType::SNP,
                        vcf_row: 8,
                        reference: "g".to_string(),
                        alt: "t".to_string(),
                        genome_index: 27,
                        is_minor: true,
                        vcf_idx: Some(2),
                    },
                ],
                gene_position: Some(8),
                codes_protein: Some(true),
                ref_nucleotides: Some("ggg".to_string()),
                alt_nucleotides: Some("gzz".to_string()),
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: Some(8),
                amino_acid_sequence: Some('Z'),
            },
            Mutation {
                mutation: "26_mixed:50".to_string(),
                gene: "A".to_string(),
                evidence: vec![
                    Evidence {
                        cov: Some(50),
                        frs: Some(OrderedFloat(50.0 / 100.0)),
                        genotype: "1/3".to_string(),
                        call_type: AltType::DEL,
                        vcf_row: 9,
                        reference: "g".to_string(),
                        alt: "g".to_string(),
                        genome_index: 29,
                        is_minor: true,
                        vcf_idx: Some(3),
                    },
                    Evidence {
                        cov: Some(48),
                        frs: Some(OrderedFloat(48.0 / 100.0)),
                        genotype: "1/3".to_string(),
                        call_type: AltType::SNP,
                        vcf_row: 9,
                        reference: "g".to_string(),
                        alt: "a".to_string(),
                        genome_index: 29,
                        is_minor: true,
                        vcf_idx: Some(1),
                    },
                ],
                gene_position: Some(26),
                codes_protein: Some(true),
                ref_nucleotides: None,
                alt_nucleotides: None,
                nucleotide_number: Some(26),
                nucleotide_index: Some(29),
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: None,
                amino_acid_sequence: None,
            },
        ];

        for (idx, mutation) in a_diff.minor_mutations.iter().enumerate() {
            assert_eq!(mutation, &expected_a_minor_mutations[idx])
        }

        let expected_a_minor_mutations_frs = vec![
            Mutation {
                mutation: "P5V:0.105".to_string(),
                gene: "A".to_string(),
                evidence: vec![
                    Evidence {
                        cov: Some(8),
                        frs: Some(OrderedFloat(8.0 / 76.0)),
                        genotype: "1/1".to_string(),
                        call_type: AltType::SNP,
                        vcf_row: 5,
                        reference: "c".to_string(),
                        alt: "g".to_string(),
                        genome_index: 16,
                        is_minor: true,
                        vcf_idx: Some(2),
                    },
                    Evidence {
                        cov: Some(8),
                        frs: Some(OrderedFloat(8.0 / 76.0)),
                        genotype: "1/1".to_string(),
                        call_type: AltType::SNP,
                        vcf_row: 5,
                        reference: "c".to_string(),
                        alt: "t".to_string(),
                        genome_index: 17,
                        is_minor: true,
                        vcf_idx: Some(2),
                    },
                    Evidence {
                        cov: Some(8),
                        frs: Some(OrderedFloat(8.0 / 76.0)),
                        genotype: "1/1".to_string(),
                        call_type: AltType::SNP,
                        vcf_row: 5,
                        reference: "c".to_string(),
                        alt: "g".to_string(),
                        genome_index: 18,
                        is_minor: true,
                        vcf_idx: Some(2),
                    },
                ],
                gene_position: Some(5),
                codes_protein: Some(true),
                ref_nucleotides: Some("ccc".to_string()),
                alt_nucleotides: Some("gtg".to_string()),
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: Some(5),
                amino_acid_sequence: Some('V'),
            },
            Mutation {
                mutation: "G7Z:0.495".to_string(),
                gene: "A".to_string(),
                evidence: vec![
                    Evidence {
                        cov: Some(99),
                        frs: Some(OrderedFloat(99.0 / 202.0)),
                        genotype: "1/2".to_string(),
                        call_type: AltType::SNP,
                        vcf_row: 6,
                        reference: "g".to_string(),
                        alt: "t".to_string(),
                        genome_index: 22,
                        is_minor: true,
                        vcf_idx: Some(1),
                    },
                    Evidence {
                        cov: Some(100),
                        frs: Some(OrderedFloat(100.0 / 202.0)),
                        genotype: "1/2".to_string(),
                        call_type: AltType::SNP,
                        vcf_row: 6,
                        reference: "g".to_string(),
                        alt: "c".to_string(),
                        genome_index: 22,
                        is_minor: true,
                        vcf_idx: Some(2),
                    },
                    Evidence {
                        cov: Some(100),
                        frs: Some(OrderedFloat(100.0 / 202.0)),
                        genotype: "0/2".to_string(),
                        call_type: AltType::SNP,
                        vcf_row: 7,
                        reference: "g".to_string(),
                        alt: "c".to_string(),
                        genome_index: 24,
                        is_minor: true,
                        vcf_idx: Some(2),
                    },
                ],
                gene_position: Some(7),
                codes_protein: Some(true),
                ref_nucleotides: Some("ggg".to_string()),
                alt_nucleotides: Some("zgc".to_string()),
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: Some(7),
                amino_acid_sequence: Some('Z'),
            },
            Mutation {
                mutation: "G8Z:0.5".to_string(),
                gene: "A".to_string(),
                evidence: vec![
                    Evidence {
                        cov: Some(48),
                        frs: Some(OrderedFloat(48.0 / 100.0)),
                        genotype: "1/2".to_string(),
                        call_type: AltType::SNP,
                        vcf_row: 8,
                        reference: "g".to_string(),
                        alt: "a".to_string(),
                        genome_index: 26,
                        is_minor: true,
                        vcf_idx: Some(1),
                    },
                    Evidence {
                        cov: Some(50),
                        frs: Some(OrderedFloat(50.0 / 100.0)),
                        genotype: "1/2".to_string(),
                        call_type: AltType::SNP,
                        vcf_row: 8,
                        reference: "g".to_string(),
                        alt: "c".to_string(),
                        genome_index: 26,
                        is_minor: true,
                        vcf_idx: Some(2),
                    },
                    Evidence {
                        cov: Some(48),
                        frs: Some(OrderedFloat(48.0 / 100.0)),
                        genotype: "1/2".to_string(),
                        call_type: AltType::SNP,
                        vcf_row: 8,
                        reference: "g".to_string(),
                        alt: "a".to_string(),
                        genome_index: 27,
                        is_minor: true,
                        vcf_idx: Some(1),
                    },
                    Evidence {
                        cov: Some(50),
                        frs: Some(OrderedFloat(50.0 / 100.0)),
                        genotype: "1/2".to_string(),
                        call_type: AltType::SNP,
                        vcf_row: 8,
                        reference: "g".to_string(),
                        alt: "t".to_string(),
                        genome_index: 27,
                        is_minor: true,
                        vcf_idx: Some(2),
                    },
                ],
                gene_position: Some(8),
                codes_protein: Some(true),
                ref_nucleotides: Some("ggg".to_string()),
                alt_nucleotides: Some("gzz".to_string()),
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: Some(8),
                amino_acid_sequence: Some('Z'),
            },
            Mutation {
                mutation: "26_mixed:0.5".to_string(),
                gene: "A".to_string(),
                evidence: vec![
                    Evidence {
                        cov: Some(50),
                        frs: Some(OrderedFloat(50.0 / 100.0)),
                        genotype: "1/3".to_string(),
                        call_type: AltType::DEL,
                        vcf_row: 9,
                        reference: "g".to_string(),
                        alt: "g".to_string(),
                        genome_index: 29,
                        is_minor: true,
                        vcf_idx: Some(3),
                    },
                    Evidence {
                        cov: Some(48),
                        frs: Some(OrderedFloat(48.0 / 100.0)),
                        genotype: "1/3".to_string(),
                        call_type: AltType::SNP,
                        vcf_row: 9,
                        reference: "g".to_string(),
                        alt: "a".to_string(),
                        genome_index: 29,
                        is_minor: true,
                        vcf_idx: Some(1),
                    },
                ],
                gene_position: Some(26),
                codes_protein: Some(true),
                ref_nucleotides: None,
                alt_nucleotides: None,
                nucleotide_number: Some(26),
                nucleotide_index: Some(29),
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: None,
                amino_acid_sequence: None,
            },
        ];

        let a_diff = GeneDifference::new(
            reference.get_gene("A".to_string()),
            sample.get_gene("A".to_string()),
            MinorType::FRS,
        );
        for (idx, mutation) in a_diff.minor_mutations.iter().enumerate() {
            assert_eq!(mutation, &expected_a_minor_mutations_frs[idx])
        }

        let b_diff = GeneDifference::new(
            reference.get_gene("B".to_string()),
            sample.get_gene("B".to_string()),
            MinorType::COV,
        );

        let expected_b_mutations = vec![
            Mutation {
                mutation: "G1Z".to_string(),
                gene: "B".to_string(),
                evidence: vec![
                    Evidence {
                        cov: None,
                        frs: None,
                        genotype: "1/3".to_string(),
                        call_type: AltType::HET,
                        vcf_row: 9,
                        reference: "g".to_string(),
                        alt: "z".to_string(),
                        genome_index: 28,
                        is_minor: false,
                        vcf_idx: None,
                    },
                    Evidence {
                        cov: None,
                        frs: None,
                        genotype: "1/3".to_string(),
                        call_type: AltType::HET,
                        vcf_row: 9,
                        reference: "g".to_string(),
                        alt: "z".to_string(),
                        genome_index: 29,
                        is_minor: false,
                        vcf_idx: None,
                    },
                ],
                gene_position: Some(1),
                codes_protein: Some(true),
                ref_nucleotides: Some("ggg".to_string()),
                alt_nucleotides: Some("zzg".to_string()),
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: Some(1),
                amino_acid_sequence: Some('Z'),
            },
            Mutation {
                mutation: "6_ins_tt".to_string(),
                gene: "B".to_string(),
                evidence: vec![Evidence {
                    cov: Some(199),
                    frs: Some(OrderedFloat(199.0 / 200.0)),
                    genotype: "1/1".to_string(),
                    call_type: AltType::INS,
                    vcf_row: 10,
                    reference: "t".to_string(),
                    alt: "tt".to_string(),
                    genome_index: 33,
                    is_minor: false,
                    vcf_idx: Some(1),
                }],
                gene_position: Some(6),
                codes_protein: Some(true),
                ref_nucleotides: None,
                alt_nucleotides: None,
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: Some(2),
                indel_nucleotides: Some("tt".to_string()),
                amino_acid_number: None,
                amino_acid_sequence: None,
            },
            Mutation {
                mutation: "F4L".to_string(),
                gene: "B".to_string(),
                evidence: vec![Evidence {
                    cov: Some(199),
                    frs: Some(OrderedFloat(199.0 / 200.0)),
                    genotype: "1/1".to_string(),
                    call_type: AltType::SNP,
                    vcf_row: 12,
                    reference: "t".to_string(),
                    alt: "a".to_string(),
                    genome_index: 39,
                    is_minor: false,
                    vcf_idx: Some(1),
                }],
                gene_position: Some(4),
                codes_protein: Some(true),
                ref_nucleotides: Some("ttt".to_string()),
                alt_nucleotides: Some("tta".to_string()),
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: Some(4),
                amino_acid_sequence: Some('L'),
            },
            Mutation {
                mutation: "10_del_t".to_string(),
                gene: "B".to_string(),
                evidence: vec![Evidence {
                    cov: Some(199),
                    frs: Some(OrderedFloat(199.0 / 200.0)),
                    genotype: "1/1".to_string(),
                    call_type: AltType::DEL,
                    vcf_row: 11,
                    reference: "t".to_string(),
                    alt: "t".to_string(),
                    genome_index: 37,
                    is_minor: false,
                    vcf_idx: Some(1),
                }],
                gene_position: Some(10),
                codes_protein: Some(true),
                ref_nucleotides: None,
                alt_nucleotides: None,
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: Some(-1),
                indel_nucleotides: Some("t".to_string()),
                amino_acid_number: None,
                amino_acid_sequence: None,
            },
            Mutation {
                mutation: "12_ins_g".to_string(),
                gene: "B".to_string(),
                evidence: vec![Evidence {
                    cov: Some(199),
                    frs: Some(OrderedFloat(199.0 / 200.0)),
                    genotype: "1/1".to_string(),
                    call_type: AltType::INS,
                    vcf_row: 12,
                    reference: "t".to_string(),
                    alt: "g".to_string(),
                    genome_index: 39,
                    is_minor: false,
                    vcf_idx: Some(1),
                }],
                gene_position: Some(12),
                codes_protein: Some(true),
                ref_nucleotides: None,
                alt_nucleotides: None,
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: Some(1),
                indel_nucleotides: Some("g".to_string()),
                amino_acid_number: None,
                amino_acid_sequence: None,
            },
        ];

        for (idx, mutation) in b_diff.mutations.iter().enumerate() {
            assert_eq!(mutation, &expected_b_mutations[idx])
        }

        let expected_b_minor_mutations = vec![Mutation {
            mutation: "2_mixed:50".to_string(),
            gene: "B".to_string(),
            evidence: vec![
                Evidence {
                    cov: Some(50),
                    frs: Some(OrderedFloat(50.0 / 100.0)),
                    genotype: "1/3".to_string(),
                    call_type: AltType::DEL,
                    vcf_row: 9,
                    reference: "g".to_string(),
                    alt: "g".to_string(),
                    genome_index: 29,
                    is_minor: true,
                    vcf_idx: Some(3),
                },
                Evidence {
                    cov: Some(48),
                    frs: Some(OrderedFloat(48.0 / 100.0)),
                    genotype: "1/3".to_string(),
                    call_type: AltType::SNP,
                    vcf_row: 9,
                    reference: "g".to_string(),
                    alt: "a".to_string(),
                    genome_index: 29,
                    is_minor: true,
                    vcf_idx: Some(1),
                },
            ],
            gene_position: Some(2),
            codes_protein: Some(true),
            ref_nucleotides: None,
            alt_nucleotides: None,
            nucleotide_number: Some(2),
            nucleotide_index: Some(29),
            indel_length: None,
            indel_nucleotides: None,
            amino_acid_number: None,
            amino_acid_sequence: None,
        }];

        for (idx, mutation) in b_diff.minor_mutations.iter().enumerate() {
            assert_eq!(mutation, &expected_b_minor_mutations[idx])
        }

        let c_diff = GeneDifference::new(
            reference.get_gene("C".to_string()),
            sample.get_gene("C".to_string()),
            MinorType::COV,
        );

        assert_eq!(c_diff.mutations.len(), 0);
        assert_eq!(c_diff.minor_mutations.len(), 0);
    }

    #[test]
    fn test_large_deletions() {
        let mut reference = Genome::new("reference/TEST-DNA.gbk");
        let vcf = VCFFile::new("test/TEST-DNA-4.vcf".to_string(), false, 3);

        let mut sample = mutate(&reference, vcf);

        let diff = GenomeDifference::new(reference.clone(), sample.clone(), MinorType::COV);

        let expected_vcf_rows = [
            VCFRow {
                position: 2,
                reference: "aaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc".to_string(),
                alternative: vec!["a".to_string()],
                filter: vec!["PASS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["1/1".to_string()]),
                    ("DP".to_string(), vec!["4".to_string()]),
                    ("COV".to_string(), vec!["1".to_string(), "3".to_string()]),
                    ("GT_CONF".to_string(), vec!["2.05".to_string()]),
                ]),
                is_filter_pass: true,
            }
        ];

        for (idx, row) in expected_vcf_rows.iter().enumerate() {
            assert_eq!(*row, sample.get_vcf_row(idx));
        }

        let expected_genome_variants = vec![
            Variant {
                variant: "3_del_aaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc".to_string(),
                nucleotide_index: 3,
                evidence: 0,
                vcf_idx: Some(1),
                indel_length: -91,
                indel_nucleotides: Some("aaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc".to_string()),
                gene_position: Some(-1),
                codon_idx: None,
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "3_del_aaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc".to_string(),
                nucleotide_index: 3,
                evidence: 0,
                vcf_idx: Some(1),
                indel_length: -91,
                indel_nucleotides: Some("aaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc".to_string()),
                gene_position: Some(1),
                codon_idx: Some(0),
                gene_name: Some("B".to_string()),
            },
            Variant {
                variant: "3_del_aaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc".to_string(),
                nucleotide_index: 3,
                evidence: 0,
                vcf_idx: Some(1),
                indel_length: -91,
                indel_nucleotides: Some("aaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc".to_string()),
                gene_position: Some(4),
                codon_idx: Some(0),
                gene_name: Some("C".to_string()),
            },
        ];

        for (idx, variant) in diff.variants.iter().enumerate() {
            assert_eq!(variant, &expected_genome_variants[idx]);
        }

        assert_eq!(diff.minor_variants.len(), 0);

        // This deletion should cover large deletions in A and B, as well as partial deletions in C
        let expected_a_mutations = vec![
            Mutation {
                mutation: "-1_del_aaaaaaaaccccccccccgggggggggg".to_string(),
                gene: "A".to_string(),
                evidence: vec![Evidence {
                    cov: Some(3),
                    frs: Some(OrderedFloat(0.75)),
                    genotype: "1/1".to_string(),
                    call_type: AltType::DEL,
                    vcf_row: 0,
                    reference: "a".to_string(),
                    alt: "aaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc".to_string(),
                    genome_index: 3,
                    is_minor: false,
                    vcf_idx: Some(1),
                }],
                gene_position: Some(-1),
                codes_protein: Some(false),
                ref_nucleotides: None,
                alt_nucleotides: None,
                nucleotide_number: Some(-1),
                nucleotide_index: Some(3),
                indel_length: Some(-28),
                indel_nucleotides: Some("aaaaaaaaccccccccccgggggggggg".to_string()),
                amino_acid_number: None,
                amino_acid_sequence: None,
            },
            Mutation {
                mutation: "del_0.93".to_string(),
                gene: "A".to_string(),
                evidence: vec![],
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
            },
        ];

        let a_diff = GeneDifference::new(
            reference.get_gene("A".to_string()),
            sample.get_gene("A".to_string()),
            MinorType::COV,
        );

        for (idx, mutation) in a_diff.mutations.iter().enumerate() {
            assert_eq!(mutation, &expected_a_mutations[idx]);
        }

        assert_eq!(a_diff.minor_mutations.len(), 0);

        let b_diff = GeneDifference::new(
            reference.get_gene("B".to_string()),
            sample.get_gene("B".to_string()),
            MinorType::COV,
        );

        let expected_b_mutations = vec![
            Mutation {
                mutation: "1_del_gggttttttttttaaaaaaaaaacccccccccc".to_string(),
                gene: "B".to_string(),
                evidence: vec![Evidence {
                    cov: Some(3),
                    frs: Some(OrderedFloat(0.75)),
                    genotype: "1/1".to_string(),
                    call_type: AltType::DEL,
                    vcf_row: 0,
                    reference: "a".to_string(),
                    alt: "gggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc"
                        .to_string(),
                    genome_index: 28,
                    is_minor: false,
                    vcf_idx: Some(1),
                }],
                gene_position: Some(1),
                codes_protein: Some(true),
                ref_nucleotides: None,
                alt_nucleotides: None,
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: Some(-33),
                indel_nucleotides: Some("gggttttttttttaaaaaaaaaacccccccccc".to_string()),
                amino_acid_number: None,
                amino_acid_sequence: None,
            },
            Mutation {
                mutation: "del_1.0".to_string(),
                gene: "B".to_string(),
                evidence: vec![],
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
            },
        ];

        for (idx, mutation) in b_diff.mutations.iter().enumerate() {
            assert_eq!(mutation, &expected_b_mutations[idx]);
        }

        assert_eq!(b_diff.minor_mutations.len(), 0);

        let c_diff = GeneDifference::new(
            reference.get_gene("C".to_string()),
            sample.get_gene("C".to_string()),
            MinorType::COV,
        );

        let expected_c_mutations = vec![Mutation {
            mutation: "4_del_ggg".to_string(),
            gene: "C".to_string(),
            evidence: vec![Evidence {
                cov: Some(3),
                frs: Some(OrderedFloat(0.75)),
                genotype: "1/1".to_string(),
                call_type: AltType::DEL,
                vcf_row: 0,
                reference: "g".to_string(),
                alt: "ggg".to_string(),
                genome_index: 94,
                is_minor: false,
                vcf_idx: Some(1),
            }],
            gene_position: Some(4),
            codes_protein: Some(true),
            ref_nucleotides: None,
            alt_nucleotides: None,
            nucleotide_number: None,
            nucleotide_index: None,
            indel_length: Some(-3),
            indel_nucleotides: Some("ggg".to_string()),
            amino_acid_number: None,
            amino_acid_sequence: None,
        }];

        for (idx, mutation) in c_diff.mutations.iter().enumerate() {
            assert_eq!(mutation, &expected_c_mutations[idx]);
        }

        assert_eq!(c_diff.minor_mutations.len(), 0);
    }

    #[test]
    fn test_large_deletions_minor() {
        // Mirror of test_large_deletions, but with a VCF where the deletion is minor instead
        let mut reference = Genome::new("reference/TEST-DNA.gbk");
        let vcf = VCFFile::new("test/TEST-DNA-4-minor.vcf".to_string(), false, 1);

        let mut sample = mutate(&reference, vcf);

        let diff = GenomeDifference::new(reference.clone(), sample.clone(), MinorType::COV);

        let expected_vcf_rows = [
            VCFRow {
                position: 2,
                reference: "aaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc".to_string(),
                alternative: vec!["a".to_string()],
                filter: vec!["PASS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["0/0".to_string()]),
                    ("DP".to_string(), vec!["4".to_string()]),
                    ("COV".to_string(), vec!["1".to_string(), "3".to_string()]),
                    ("GT_CONF".to_string(), vec!["2.05".to_string()]),
                ]),
                is_filter_pass: true,
            }

        ];

        for (idx, row) in expected_vcf_rows.iter().enumerate() {
            assert_eq!(*row, sample.get_vcf_row(idx));
        }

        let expected_genome_minor_variants = vec![
            Variant {
                variant: "3_del_aaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc:3".to_string(),
                nucleotide_index: 3,
                evidence: 0,
                vcf_idx: Some(1),
                indel_length: -91,
                indel_nucleotides: Some("aaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc".to_string()),
                gene_position: Some(-1),
                codon_idx: None,
                gene_name: Some("A".to_string()),
            },
            Variant {
                variant: "3_del_aaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc:3".to_string(),
                nucleotide_index: 3,
                evidence: 0,
                vcf_idx: Some(1),
                indel_length: -91,
                indel_nucleotides: Some("aaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc".to_string()),
                gene_position: Some(1),
                codon_idx: Some(0),
                gene_name: Some("B".to_string()),
            },
            Variant {
                variant: "3_del_aaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc:3".to_string(),
                nucleotide_index: 3,
                evidence: 0,
                vcf_idx: Some(1),
                indel_length: -91,
                indel_nucleotides: Some("aaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc".to_string()),
                gene_position: Some(4),
                codon_idx: Some(0),
                gene_name: Some("C".to_string()),
            },
        ];

        for (idx, variant) in diff.minor_variants.iter().enumerate() {
            assert_eq!(variant, &expected_genome_minor_variants[idx]);
        }

        assert_eq!(diff.variants.len(), 0);

        // This deletion should cover large deletions in A and B, as well as partial deletions in C
        let expected_a_minor_mutations = vec![
            Mutation {
                mutation: "-1_del_aaaaaaaaccccccccccgggggggggg:3".to_string(),
                gene: "A".to_string(),
                evidence: vec![Evidence {
                    cov: Some(3),
                    frs: Some(OrderedFloat(0.75)),
                    genotype: "0/0".to_string(),
                    call_type: AltType::DEL,
                    vcf_row: 0,
                    reference: "a".to_string(),
                    alt: "aaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc".to_string(),
                    genome_index: 3,
                    is_minor: true,
                    vcf_idx: Some(1),
                }],
                gene_position: Some(-1),
                codes_protein: Some(false),
                ref_nucleotides: None,
                alt_nucleotides: None,
                nucleotide_number: Some(-1),
                nucleotide_index: Some(3),
                indel_length: Some(-28),
                indel_nucleotides: Some("aaaaaaaaccccccccccgggggggggg".to_string()),
                amino_acid_number: None,
                amino_acid_sequence: None,
            },
            Mutation {
                mutation: "del_0.93:3".to_string(),
                gene: "A".to_string(),
                evidence: vec![],
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
            },
        ];

        let a_diff = GeneDifference::new(
            reference.get_gene("A".to_string()),
            sample.get_gene("A".to_string()),
            MinorType::COV,
        );

        for (idx, mutation) in a_diff.minor_mutations.iter().enumerate() {
            assert_eq!(mutation, &expected_a_minor_mutations[idx]);
        }

        assert_eq!(a_diff.mutations.len(), 0);

        let b_diff = GeneDifference::new(
            reference.get_gene("B".to_string()),
            sample.get_gene("B".to_string()),
            MinorType::COV,
        );

        let expected_b_minor_mutations = vec![
            Mutation {
                mutation: "1_del_gggttttttttttaaaaaaaaaacccccccccc:3".to_string(),
                gene: "B".to_string(),
                evidence: vec![Evidence {
                    cov: Some(3),
                    frs: Some(OrderedFloat(0.75)),
                    genotype: "0/0".to_string(),
                    call_type: AltType::DEL,
                    vcf_row: 0,
                    reference: "a".to_string(),
                    alt: "gggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc"
                        .to_string(),
                    genome_index: 28,
                    is_minor: true,
                    vcf_idx: Some(1),
                }],
                gene_position: Some(1),
                codes_protein: Some(true),
                ref_nucleotides: None,
                alt_nucleotides: None,
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: Some(-33),
                indel_nucleotides: Some("gggttttttttttaaaaaaaaaacccccccccc".to_string()),
                amino_acid_number: None,
                amino_acid_sequence: None,
            },
            Mutation {
                mutation: "del_1.0:3".to_string(),
                gene: "B".to_string(),
                evidence: vec![],
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
            },
        ];

        for (idx, mutation) in b_diff.minor_mutations.iter().enumerate() {
            assert_eq!(mutation, &expected_b_minor_mutations[idx]);
        }

        assert_eq!(b_diff.mutations.len(), 0);

        let b_diff = GeneDifference::new(
            reference.get_gene("B".to_string()),
            sample.get_gene("B".to_string()),
            MinorType::FRS,
        );

        let expected_b_minor_mutations = vec![
            Mutation {
                mutation: "1_del_gggttttttttttaaaaaaaaaacccccccccc:0.75".to_string(),
                gene: "B".to_string(),
                evidence: vec![Evidence {
                    cov: Some(3),
                    frs: Some(OrderedFloat(0.75)),
                    genotype: "0/0".to_string(),
                    call_type: AltType::DEL,
                    vcf_row: 0,
                    reference: "a".to_string(),
                    alt: "gggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc"
                        .to_string(),
                    genome_index: 28,
                    is_minor: true,
                    vcf_idx: Some(1),
                }],
                gene_position: Some(1),
                codes_protein: Some(true),
                ref_nucleotides: None,
                alt_nucleotides: None,
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: Some(-33),
                indel_nucleotides: Some("gggttttttttttaaaaaaaaaacccccccccc".to_string()),
                amino_acid_number: None,
                amino_acid_sequence: None,
            },
            Mutation {
                mutation: "del_1.0:0.75".to_string(),
                gene: "B".to_string(),
                evidence: vec![],
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
            },
        ];

        for (idx, mutation) in b_diff.minor_mutations.iter().enumerate() {
            assert_eq!(mutation, &expected_b_minor_mutations[idx]);
        }

        assert_eq!(b_diff.mutations.len(), 0);

        let c_diff = GeneDifference::new(
            reference.get_gene("C".to_string()),
            sample.get_gene("C".to_string()),
            MinorType::COV,
        );

        let expected_c_minor_mutations = vec![Mutation {
            mutation: "4_del_ggg:3".to_string(),
            gene: "C".to_string(),
            evidence: vec![Evidence {
                cov: Some(3),
                frs: Some(OrderedFloat(0.75)),
                genotype: "0/0".to_string(),
                call_type: AltType::DEL,
                vcf_row: 0,
                reference: "g".to_string(),
                alt: "ggg".to_string(),
                genome_index: 94,
                is_minor: true,
                vcf_idx: Some(1),
            }],
            gene_position: Some(4),
            codes_protein: Some(true),
            ref_nucleotides: None,
            alt_nucleotides: None,
            nucleotide_number: None,
            nucleotide_index: None,
            indel_length: Some(-3),
            indel_nucleotides: Some("ggg".to_string()),
            amino_acid_number: None,
            amino_acid_sequence: None,
        }];

        for (idx, mutation) in c_diff.minor_mutations.iter().enumerate() {
            assert_eq!(mutation, &expected_c_minor_mutations[idx]);
        }

        assert_eq!(c_diff.mutations.len(), 0);
    }

    #[test]
    fn test_revcomp_minor_population() {
        let mut reference = Genome::new("reference/NC_000962.3.gbk");
        let vcf = VCFFile::new("test/minor-populations-revcomp.vcf".to_string(), false, 3);
        let mut samples = mutate(&reference, vcf);

        let kat_g_diff = GeneDifference::new(
            reference.get_gene("katG".to_string()),
            samples.get_gene("katG".to_string()),
            MinorType::COV,
        );
        assert_eq!(kat_g_diff.mutations.len(), 0);

        let expected_vcf_rows = vec![
            VCFRow {
                // 0
                position: 2154397,
                reference: "g".to_string(),
                alternative: vec!["t".to_string(), "c".to_string()],
                filter: vec!["MIN_FRS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["0/0".to_string()]),
                    ("DP".to_string(), vec!["100".to_string()]),
                    (
                        "ALLELE_DP".to_string(),
                        vec!["75".to_string(), "15".to_string(), "10".to_string()],
                    ),
                    ("FRS".to_string(), vec!["0.75".to_string()]),
                    ("COV_TOTAL".to_string(), vec!["100".to_string()]),
                    (
                        "COV".to_string(),
                        vec!["75".to_string(), "15".to_string(), "10".to_string()],
                    ),
                    ("GT_CONF".to_string(), vec!["731.8".to_string()]),
                    ("GT_CONF_PERCENTILE".to_string(), vec!["77.89".to_string()]),
                ]),
                is_filter_pass: false,
            },
            VCFRow {
                // 1
                position: 2154400,
                reference: "cgg".to_string(),
                alternative: vec!["c".to_string()],
                filter: vec!["MIN_FRS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["0/0".to_string()]),
                    ("DP".to_string(), vec!["100".to_string()]),
                    (
                        "ALLELE_DP".to_string(),
                        vec!["75".to_string(), "25".to_string()],
                    ),
                    ("FRS".to_string(), vec!["0.75".to_string()]),
                    ("COV_TOTAL".to_string(), vec!["100".to_string()]),
                    ("COV".to_string(), vec!["75".to_string(), "25".to_string()]),
                    ("GT_CONF".to_string(), vec!["731.8".to_string()]),
                    ("GT_CONF_PERCENTILE".to_string(), vec!["77.89".to_string()]),
                ]),
                is_filter_pass: false,
            },
        ];

        for (idx, row) in expected_vcf_rows.iter().enumerate() {
            assert_eq!(*row, samples.get_vcf_row(idx));
        }

        let expected_katg_minor_mutations = vec![
            Mutation {
                mutation: "1710_del_cc:25".to_string(),
                gene: "katG".to_string(),
                evidence: vec![Evidence {
                    cov: Some(25),
                    frs: Some(OrderedFloat(25.0 / 100.0)),
                    genotype: "0/0".to_string(),
                    call_type: AltType::DEL,
                    vcf_row: 1,
                    reference: "g".to_string(),
                    alt: "gg".to_string(),
                    genome_index: 2154401,
                    is_minor: true,
                    vcf_idx: Some(1),
                }],
                gene_position: Some(1710),
                codes_protein: Some(true),
                ref_nucleotides: None,
                alt_nucleotides: None,
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: Some(-2),
                indel_nucleotides: Some("cc".to_string()),
                amino_acid_number: None,
                amino_acid_sequence: None,
            },
            Mutation {
                mutation: "T572Z:15".to_string(),
                gene: "katG".to_string(),
                evidence: vec![
                    Evidence {
                        cov: Some(15),
                        frs: Some(OrderedFloat(15.0 / 100.0)),
                        genotype: "0/0".to_string(),
                        call_type: AltType::SNP,
                        vcf_row: 0,
                        reference: "g".to_string(),
                        alt: "t".to_string(),
                        genome_index: 2154397,
                        is_minor: true,
                        vcf_idx: Some(1),
                    },
                    Evidence {
                        cov: Some(10),
                        frs: Some(OrderedFloat(10.0 / 100.0)),
                        genotype: "0/0".to_string(),
                        call_type: AltType::SNP,
                        vcf_row: 0,
                        reference: "g".to_string(),
                        alt: "c".to_string(),
                        genome_index: 2154397,
                        is_minor: true,
                        vcf_idx: Some(2),
                    },
                ],
                gene_position: Some(572),
                codes_protein: Some(true),
                ref_nucleotides: Some("acg".to_string()),
                alt_nucleotides: Some("azg".to_string()),
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: Some(572),
                amino_acid_sequence: Some('Z'),
            },
        ];

        for (idx, mutation) in kat_g_diff.minor_mutations.iter().enumerate() {
            assert_eq!(mutation, &expected_katg_minor_mutations[idx]);
        }
    }

    #[test]
    fn test_revcomp_minor_population_ad() {
        // Exact mirror of the previous test, but with AD tag instead of COV
        let mut reference = Genome::new("reference/NC_000962.3.gbk");
        let vcf = VCFFile::new(
            "test/minor-populations-revcomp-AD.vcf".to_string(),
            false,
            3,
        );
        let mut samples = mutate(&reference, vcf);

        let kat_g_diff = GeneDifference::new(
            reference.get_gene("katG".to_string()),
            samples.get_gene("katG".to_string()),
            MinorType::COV,
        );
        assert_eq!(kat_g_diff.mutations.len(), 0);

        let expected_vcf_rows = vec![
            VCFRow {
                // 0
                position: 2154397,
                reference: "g".to_string(),
                alternative: vec!["t".to_string(), "c".to_string()],
                filter: vec!["MIN_FRS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["0/0".to_string()]),
                    ("DP".to_string(), vec!["100".to_string()]),
                    (
                        "ALLELE_DP".to_string(),
                        vec!["75".to_string(), "15".to_string(), "10".to_string()],
                    ),
                    ("FRS".to_string(), vec!["0.75".to_string()]),
                    ("COV_TOTAL".to_string(), vec!["100".to_string()]),
                    (
                        "COV".to_string(),
                        vec!["75".to_string(), "15".to_string(), "10".to_string()],
                    ),
                    (
                        "AD".to_string(),
                        vec!["75".to_string(), "15".to_string(), "10".to_string()],
                    ),
                    ("GT_CONF".to_string(), vec!["731.8".to_string()]),
                    ("GT_CONF_PERCENTILE".to_string(), vec!["77.89".to_string()]),
                ]),
                is_filter_pass: false,
            },
            VCFRow {
                // 1
                position: 2154400,
                reference: "cgg".to_string(),
                alternative: vec!["c".to_string()],
                filter: vec!["MIN_FRS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["0/0".to_string()]),
                    ("DP".to_string(), vec!["100".to_string()]),
                    (
                        "ALLELE_DP".to_string(),
                        vec!["75".to_string(), "25".to_string()],
                    ),
                    ("FRS".to_string(), vec!["0.75".to_string()]),
                    ("COV_TOTAL".to_string(), vec!["100".to_string()]),
                    ("COV".to_string(), vec!["75".to_string(), "25".to_string()]),
                    ("AD".to_string(), vec!["75".to_string(), "25".to_string()]),
                    ("GT_CONF".to_string(), vec!["731.8".to_string()]),
                    ("GT_CONF_PERCENTILE".to_string(), vec!["77.89".to_string()]),
                ]),
                is_filter_pass: false,
            },
        ];

        for (idx, row) in expected_vcf_rows.iter().enumerate() {
            assert_eq!(*row, samples.get_vcf_row(idx));
        }

        let expected_katg_minor_mutations = vec![
            Mutation {
                mutation: "1710_del_cc:25".to_string(),
                gene: "katG".to_string(),
                evidence: vec![Evidence {
                    cov: Some(25),
                    frs: Some(OrderedFloat(25.0 / 100.0)),
                    genotype: "0/0".to_string(),
                    call_type: AltType::DEL,
                    vcf_row: 1,
                    reference: "g".to_string(),
                    alt: "gg".to_string(),
                    genome_index: 2154401,
                    is_minor: true,
                    vcf_idx: Some(1),
                }],
                gene_position: Some(1710),
                codes_protein: Some(true),
                ref_nucleotides: None,
                alt_nucleotides: None,
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: Some(-2),
                indel_nucleotides: Some("cc".to_string()),
                amino_acid_number: None,
                amino_acid_sequence: None,
            },
            Mutation {
                mutation: "T572Z:15".to_string(),
                gene: "katG".to_string(),
                evidence: vec![
                    Evidence {
                        cov: Some(15),
                        frs: Some(OrderedFloat(15.0 / 100.0)),
                        genotype: "0/0".to_string(),
                        call_type: AltType::SNP,
                        vcf_row: 0,
                        reference: "g".to_string(),
                        alt: "t".to_string(),
                        genome_index: 2154397,
                        is_minor: true,
                        vcf_idx: Some(1),
                    },
                    Evidence {
                        cov: Some(10),
                        frs: Some(OrderedFloat(10.0 / 100.0)),
                        genotype: "0/0".to_string(),
                        call_type: AltType::SNP,
                        vcf_row: 0,
                        reference: "g".to_string(),
                        alt: "c".to_string(),
                        genome_index: 2154397,
                        is_minor: true,
                        vcf_idx: Some(2),
                    },
                ],
                gene_position: Some(572),
                codes_protein: Some(true),
                ref_nucleotides: Some("acg".to_string()),
                alt_nucleotides: Some("azg".to_string()),
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: Some(572),
                amino_acid_sequence: Some('Z'),
            },
        ];

        for (idx, mutation) in kat_g_diff.minor_mutations.iter().enumerate() {
            assert_eq!(mutation, &expected_katg_minor_mutations[idx]);
        }
    }

    #[test]
    fn test_revcomp_minor_population_non_coding() {
        let mut reference = Genome::new("reference/NC_000962.3.gbk");
        for gene in reference.gene_definitions.iter_mut() {
            if gene.name == "katG" {
                gene.coding = false;
            }
        }
        let vcf = VCFFile::new("test/minor-populations-revcomp.vcf".to_string(), false, 3);
        let mut sample = mutate(&reference, vcf);

        let diff = GenomeDifference::new(reference.clone(), sample.clone(), MinorType::COV);

        let expected_vcf_rows = vec![
            VCFRow {
                // 0
                position: 2154397,
                reference: "g".to_string(),
                alternative: vec!["t".to_string(), "c".to_string()],
                filter: vec!["MIN_FRS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["0/0".to_string()]),
                    ("DP".to_string(), vec!["100".to_string()]),
                    (
                        "ALLELE_DP".to_string(),
                        vec!["75".to_string(), "15".to_string(), "10".to_string()],
                    ),
                    ("FRS".to_string(), vec!["0.75".to_string()]),
                    ("COV_TOTAL".to_string(), vec!["100".to_string()]),
                    (
                        "COV".to_string(),
                        vec!["75".to_string(), "15".to_string(), "10".to_string()],
                    ),
                    ("GT_CONF".to_string(), vec!["731.8".to_string()]),
                    ("GT_CONF_PERCENTILE".to_string(), vec!["77.89".to_string()]),
                ]),
                is_filter_pass: false,
            },
            VCFRow {
                // 1
                position: 2154400,
                reference: "cgg".to_string(),
                alternative: vec!["c".to_string()],
                filter: vec!["MIN_FRS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["0/0".to_string()]),
                    ("DP".to_string(), vec!["100".to_string()]),
                    (
                        "ALLELE_DP".to_string(),
                        vec!["75".to_string(), "25".to_string()],
                    ),
                    ("FRS".to_string(), vec!["0.75".to_string()]),
                    ("COV_TOTAL".to_string(), vec!["100".to_string()]),
                    ("COV".to_string(), vec!["75".to_string(), "25".to_string()]),
                    ("GT_CONF".to_string(), vec!["731.8".to_string()]),
                    ("GT_CONF_PERCENTILE".to_string(), vec!["77.89".to_string()]),
                ]),
                is_filter_pass: false,
            },
        ];

        for (idx, row) in expected_vcf_rows.iter().enumerate() {
            assert_eq!(*row, sample.get_vcf_row(idx));
        }

        let expected_genome_minor_variants = vec![
            Variant {
                variant: "2154397g>t:15".to_string(),
                nucleotide_index: 2154397,
                evidence: 0,
                vcf_idx: Some(1),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(1715),
                codon_idx: None,
                gene_name: Some("katG".to_string()),
            },
            Variant {
                variant: "2154397g>c:10".to_string(),
                nucleotide_index: 2154397,
                evidence: 0,
                vcf_idx: Some(2),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(1715),
                codon_idx: None,
                gene_name: Some("katG".to_string()),
            },
            Variant {
                variant: "2154401_del_gg:25".to_string(),
                nucleotide_index: 2154401,
                evidence: 1,
                vcf_idx: Some(1),
                indel_length: -2,
                indel_nucleotides: Some("gg".to_string()),
                gene_position: Some(1710),
                codon_idx: None,
                gene_name: Some("katG".to_string()),
            },
        ];

        for (idx, variant) in diff.minor_variants.iter().enumerate() {
            assert_eq!(variant, &expected_genome_minor_variants[idx]);
        }

        let diff = GenomeDifference::new(reference.clone(), sample.clone(), MinorType::FRS);

        let expected_genome_minor_variants = vec![
            Variant {
                variant: "2154397g>t:0.15".to_string(),
                nucleotide_index: 2154397,
                evidence: 0,
                vcf_idx: Some(1),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(1715),
                codon_idx: None,
                gene_name: Some("katG".to_string()),
            },
            Variant {
                variant: "2154397g>c:0.1".to_string(),
                nucleotide_index: 2154397,
                evidence: 0,
                vcf_idx: Some(2),
                indel_length: 0,
                indel_nucleotides: None,
                gene_position: Some(1715),
                codon_idx: None,
                gene_name: Some("katG".to_string()),
            },
            Variant {
                variant: "2154401_del_gg:0.25".to_string(),
                nucleotide_index: 2154401,
                evidence: 1,
                vcf_idx: Some(1),
                indel_length: -2,
                indel_nucleotides: Some("gg".to_string()),
                gene_position: Some(1710),
                codon_idx: None,
                gene_name: Some("katG".to_string()),
            },
        ];

        for (idx, variant) in diff.minor_variants.iter().enumerate() {
            assert_eq!(variant, &expected_genome_minor_variants[idx]);
        }

        let kat_g_diff = GeneDifference::new(
            reference.get_gene("katG".to_string()),
            sample.get_gene("katG".to_string()),
            MinorType::FRS,
        );
        assert_eq!(kat_g_diff.mutations.len(), 0);

        let expected_katg_minor_mutations = vec![
            Mutation {
                mutation: "1710_del_cc:0.25".to_string(),
                gene: "katG".to_string(),
                evidence: vec![Evidence {
                    cov: Some(25),
                    frs: Some(OrderedFloat(25.0 / 100.0)),
                    genotype: "0/0".to_string(),
                    call_type: AltType::DEL,
                    vcf_row: 1,
                    reference: "g".to_string(),
                    alt: "gg".to_string(),
                    genome_index: 2154401,
                    is_minor: true,
                    vcf_idx: Some(1),
                }],
                gene_position: Some(1710),
                codes_protein: Some(false),
                ref_nucleotides: None,
                alt_nucleotides: None,
                nucleotide_number: Some(1710),
                nucleotide_index: Some(2154402),
                indel_length: Some(-2),
                indel_nucleotides: Some("cc".to_string()),
                amino_acid_number: None,
                amino_acid_sequence: None,
            },
            Mutation {
                mutation: "c1715a:0.15".to_string(),
                gene: "katG".to_string(),
                evidence: vec![Evidence {
                    cov: Some(15),
                    frs: Some(OrderedFloat(15.0 / 100.0)),
                    genotype: "0/0".to_string(),
                    call_type: AltType::SNP,
                    vcf_row: 0,
                    reference: "g".to_string(),
                    alt: "t".to_string(),
                    genome_index: 2154397,
                    is_minor: true,
                    vcf_idx: Some(1),
                }],
                gene_position: Some(1715),
                codes_protein: Some(false),
                ref_nucleotides: Some("c".to_string()),
                alt_nucleotides: Some("a".to_string()),
                nucleotide_number: Some(1715),
                nucleotide_index: Some(2154397),
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: None,
                amino_acid_sequence: None,
            },
            Mutation {
                mutation: "c1715g:0.1".to_string(),
                gene: "katG".to_string(),
                evidence: vec![Evidence {
                    cov: Some(10),
                    frs: Some(OrderedFloat(10.0 / 100.0)),
                    genotype: "0/0".to_string(),
                    call_type: AltType::SNP,
                    vcf_row: 0,
                    reference: "g".to_string(),
                    alt: "c".to_string(),
                    genome_index: 2154397,
                    is_minor: true,
                    vcf_idx: Some(2),
                }],
                gene_position: Some(1715),
                codes_protein: Some(false),
                ref_nucleotides: Some("c".to_string()),
                alt_nucleotides: Some("g".to_string()),
                nucleotide_number: Some(1715),
                nucleotide_index: Some(2154397),
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: None,
                amino_acid_sequence: None,
            },
        ];

        for (idx, mutation) in kat_g_diff.minor_mutations.iter().enumerate() {
            assert_eq!(mutation, &expected_katg_minor_mutations[idx]);
        }
    }

    #[test]
    fn test_misc_indel() {
        let mut genome = Genome::new("reference/TEST-DNA.gbk");
        let vcf = VCFFile::new("test/TEST-DNA-misc-indel.vcf".to_string(), false, 3);
        let mut sample = mutate(&genome, vcf);

        let expected_vcf_rows = vec![
            VCFRow {
                // 0
                position: 6,
                reference: "a".to_string(),
                alternative: vec!["aa".to_string(), "aaa".to_string()],
                filter: vec!["PASS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["0/0".to_string()]),
                    ("DP".to_string(), vec!["10".to_string()]),
                    (
                        "COV".to_string(),
                        vec!["4".to_string(), "3".to_string(), "3".to_string()],
                    ),
                    ("GT_CONF".to_string(), vec!["2.05".to_string()]),
                ]),
                is_filter_pass: true,
            },
            VCFRow {
                // 1
                position: 8,
                reference: "a".to_string(),
                alternative: vec!["aa".to_string()],
                filter: vec!["PASS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["0/0".to_string()]),
                    ("DP".to_string(), vec!["10".to_string()]),
                    ("COV".to_string(), vec!["6".to_string(), "4".to_string()]),
                    ("GT_CONF".to_string(), vec!["2.05".to_string()]),
                ]),
                is_filter_pass: true,
            },
            VCFRow {
                // 2
                position: 94,
                reference: "c".to_string(),
                alternative: vec!["cc".to_string()],
                filter: vec!["PASS".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["0/0".to_string()]),
                    ("DP".to_string(), vec!["10".to_string()]),
                    ("COV".to_string(), vec!["6".to_string(), "4".to_string()]),
                    ("GT_CONF".to_string(), vec!["2.05".to_string()]),
                ]),
                is_filter_pass: true,
            },
        ];

        for (idx, row) in expected_vcf_rows.iter().enumerate() {
            assert_eq!(*row, sample.get_vcf_row(idx));
            // Double check that if we ask for a row from a genome with no VCF, it panics
            assert_panics!(genome.get_vcf_row(idx));
        }
        // Also check that asking for a row which doesn't exist it panics
        assert_panics!(sample.get_vcf_row(3333));

        let expected_minor_mutations = vec![
            Mutation {
                mutation: "3_indel:3".to_string(),
                gene: "A".to_string(),
                evidence: vec![
                    Evidence {
                        cov: Some(3),
                        frs: Some(OrderedFloat(3.0 / 10.0)),
                        genotype: "0/0".to_string(),
                        call_type: AltType::INS,
                        vcf_row: 0,
                        reference: "a".to_string(),
                        alt: "a".to_string(),
                        genome_index: 6,
                        is_minor: true,
                        vcf_idx: Some(1),
                    },
                    Evidence {
                        cov: Some(3),
                        frs: Some(OrderedFloat(3.0 / 10.0)),
                        genotype: "0/0".to_string(),
                        call_type: AltType::INS,
                        vcf_row: 0,
                        reference: "a".to_string(),
                        alt: "aa".to_string(),
                        genome_index: 6,
                        is_minor: true,
                        vcf_idx: Some(2),
                    },
                ],
                gene_position: Some(3),
                codes_protein: Some(true),
                ref_nucleotides: None,
                alt_nucleotides: None,
                nucleotide_number: Some(3),
                nucleotide_index: Some(6),
                indel_length: None,
                indel_nucleotides: None,
                amino_acid_number: None,
                amino_acid_sequence: None,
            },
            Mutation {
                mutation: "5_ins_a:4".to_string(),
                gene: "A".to_string(),
                evidence: vec![Evidence {
                    cov: Some(4),
                    frs: Some(OrderedFloat(4.0 / 10.0)),
                    genotype: "0/0".to_string(),
                    call_type: AltType::INS,
                    vcf_row: 1,
                    reference: "a".to_string(),
                    alt: "a".to_string(),
                    genome_index: 8,
                    is_minor: true,
                    vcf_idx: Some(1),
                }],
                gene_position: Some(5),
                codes_protein: Some(true),
                ref_nucleotides: None,
                alt_nucleotides: None,
                nucleotide_number: None,
                nucleotide_index: None,
                indel_length: Some(1),
                indel_nucleotides: Some("a".to_string()),
                amino_acid_number: None,
                amino_acid_sequence: None,
            },
        ];

        let gene_diff = GeneDifference::new(
            genome.get_gene("A".to_string()),
            sample.get_gene("A".to_string()),
            MinorType::COV,
        );

        for (idx, mutation) in gene_diff.minor_mutations.iter().enumerate() {
            assert_eq!(mutation, &expected_minor_mutations[idx]);
        }

        // Additional ins in C to double check revcomp ins is handled as appropriate
        let expected_minor_c_mutations = vec![Mutation {
            mutation: "2_ins_g:4".to_string(),
            gene: "C".to_string(),
            evidence: vec![Evidence {
                cov: Some(4),
                frs: Some(OrderedFloat(4.0 / 10.0)),
                genotype: "0/0".to_string(),
                call_type: AltType::INS,
                vcf_row: 2,
                reference: "c".to_string(),
                alt: "c".to_string(),
                genome_index: 94,
                is_minor: true,
                vcf_idx: Some(1),
            }],
            gene_position: Some(2),
            codes_protein: Some(true),
            ref_nucleotides: None,
            alt_nucleotides: None,
            nucleotide_number: None,
            nucleotide_index: None,
            indel_length: Some(1),
            indel_nucleotides: Some("g".to_string()),
            amino_acid_number: None,
            amino_acid_sequence: None,
        }];

        let gene_diff = GeneDifference::new(
            genome.get_gene("C".to_string()),
            sample.get_gene("C".to_string()),
            MinorType::COV,
        );

        for (idx, mutation) in gene_diff.minor_mutations.iter().enumerate() {
            assert_eq!(mutation, &expected_minor_c_mutations[idx]);
        }

        // Hack A to make it non-coding to check that too
        let mut genome = Genome::new("reference/TEST-DNA.gbk");
        for gene in genome.gene_definitions.iter_mut() {
            if gene.name == "A" {
                gene.coding = false;
            }
        }
        let vcf = VCFFile::new("test/TEST-DNA-misc-indel.vcf".to_string(), false, 3);
        let mut sample = mutate(&genome, vcf);

        let expected_minor_mutations = vec![
            Mutation {
                mutation: "3_ins_a:3".to_string(),
                gene: "A".to_string(),
                evidence: vec![Evidence {
                    cov: Some(3),
                    frs: Some(OrderedFloat(3.0 / 10.0)),
                    genotype: "0/0".to_string(),
                    call_type: AltType::INS,
                    vcf_row: 0,
                    reference: "a".to_string(),
                    alt: "a".to_string(),
                    genome_index: 6,
                    is_minor: true,
                    vcf_idx: Some(1),
                }],
                gene_position: Some(3),
                codes_protein: Some(false),
                ref_nucleotides: None,
                alt_nucleotides: None,
                nucleotide_number: Some(3),
                nucleotide_index: Some(6),
                indel_length: Some(1),
                indel_nucleotides: Some("a".to_string()),
                amino_acid_number: None,
                amino_acid_sequence: None,
            },
            Mutation {
                mutation: "3_ins_aa:3".to_string(),
                gene: "A".to_string(),
                evidence: vec![Evidence {
                    cov: Some(3),
                    frs: Some(OrderedFloat(3.0 / 10.0)),
                    genotype: "0/0".to_string(),
                    call_type: AltType::INS,
                    vcf_row: 0,
                    reference: "a".to_string(),
                    alt: "aa".to_string(),
                    genome_index: 6,
                    is_minor: true,
                    vcf_idx: Some(2),
                }],
                gene_position: Some(3),
                codes_protein: Some(false),
                ref_nucleotides: None,
                alt_nucleotides: None,
                nucleotide_number: Some(3),
                nucleotide_index: Some(6),
                indel_length: Some(2),
                indel_nucleotides: Some("aa".to_string()),
                amino_acid_number: None,
                amino_acid_sequence: None,
            },
            Mutation {
                mutation: "5_ins_a:4".to_string(),
                gene: "A".to_string(),
                evidence: vec![Evidence {
                    cov: Some(4),
                    frs: Some(OrderedFloat(4.0 / 10.0)),
                    genotype: "0/0".to_string(),
                    call_type: AltType::INS,
                    vcf_row: 1,
                    reference: "a".to_string(),
                    alt: "a".to_string(),
                    genome_index: 8,
                    is_minor: true,
                    vcf_idx: Some(1),
                }],
                gene_position: Some(5),
                codes_protein: Some(false),
                ref_nucleotides: None,
                alt_nucleotides: None,
                nucleotide_number: Some(5),
                nucleotide_index: Some(8),
                indel_length: Some(1),
                indel_nucleotides: Some("a".to_string()),
                amino_acid_number: None,
                amino_acid_sequence: None,
            },
        ];

        let gene_diff = GeneDifference::new(
            genome.get_gene("A".to_string()),
            sample.get_gene("A".to_string()),
            MinorType::COV,
        );

        for (idx, mutation) in gene_diff.minor_mutations.iter().enumerate() {
            assert_eq!(mutation, &expected_minor_mutations[idx]);
        }
    }

    #[test]
    fn test_misc() {
        // Test misc edge cases within genome

        // Check `build_all_genes` populates genes as expected
        let mut genome = Genome::new("reference/TEST-DNA.gbk");
        let mut expected_genes = HashMap::new();
        for gene in genome.gene_names.iter() {
            expected_genes.insert(gene.clone(), genome.build_gene(gene.to_string()));
        }

        genome.build_all_genes();
        assert_eq!(genome.genes, expected_genes);

        // Quickly check that at_genome_index works as expected
        assert_eq!(
            genome.at_genome_index(1),
            GenomePosition {
                reference: 'a',
                genome_idx: 1,
                alts: Vec::new(),
                genes: vec!["A".to_string()],
                is_deleted: false,
                is_deleted_minor: false,
                deleted_evidence: Vec::new(),
            }
        );
        assert_panics!(genome.at_genome_index(0));
        assert_panics!(genome.at_genome_index(-1));
        assert_panics!(genome.at_genome_index(1502263));

        // Ensure that het and nulls don't get odd complements
        assert_eq!(complement_base('z'), 'z');
        assert_eq!(complement_base('x'), 'x');

        // Ensure invalid codons cause panics
        assert_panics!(codon_to_aa("a".to_string()));
        assert_panics!(codon_to_aa("actg".to_string()));
        assert_panics!(codon_to_aa("abc".to_string()));

        let broken_alt = Alt {
            alt_type: AltType::INS,
            base: "aa".to_string(),
            evidence: Evidence {
                cov: Some(4),
                frs: Some(OrderedFloat(4.0 / 10.0)),
                genotype: "0/0".to_string(),
                call_type: AltType::INS,
                vcf_row: 2,
                reference: "c".to_string(),
                alt: "c".to_string(),
                genome_index: 94,
                is_minor: true,
                vcf_idx: Some(1),
            },
        };
        assert_eq!(
            GenomeDifference::get_nucleotide_number(
                &genome.build_gene("A".to_string()),
                &broken_alt
            ),
            None
        );

        // Testing misc gene difference panics
        let a = genome.get_gene("A".to_string());
        let b = genome.get_gene("B".to_string());
        assert_panics!(GeneDifference::new(a.clone(), b.clone(), MinorType::COV));

        for gene in genome.gene_definitions.iter_mut() {
            if gene.name == "A" {
                gene.coding = false;
            }
        }
        let a_non_coding = genome.build_gene("A".to_string());
        assert_panics!(GeneDifference::new(
            a.clone(),
            a_non_coding.clone(),
            MinorType::COV
        ));
        assert_panics!(GeneDifference::new(
            a_non_coding.clone(),
            a.clone(),
            MinorType::COV
        ));
    }

    #[test]
    fn test_revcomp_first_base_del() {
        // Test that a deletion of the first base is handled correctly
        let mut genome = Genome::new("reference/NC_000962.3.gbk");
        let vcf = VCFFile::new("test/revcomp-del-first-pos.vcf".to_string(), false, 3);
        let mut sample = mutate(&genome, vcf);

        let diff = GenomeDifference::new(genome.clone(), sample.clone(), MinorType::COV);

        // Just check the variants and mutations strings here (that's the part with the possible bugs)
        let expected_variants = [
            "218628_ins_tacggg",
            // This variant lies in 2 genes, so gets 2 variants
            "2406843_del_c",
            "2406843_del_c",
            "3820407a>g",
            "3820444_ins_ca",
            "3820446g>c",
        ];
        for (idx, variant) in diff.variants.iter().enumerate() {
            assert_eq!(variant.variant, expected_variants[idx]);
        }
        assert_eq!(diff.minor_variants.len(), 0);

        let rv2147c_diff = GeneDifference::new(
            genome.get_gene("Rv2147c".to_string()),
            sample.get_gene("Rv2147c".to_string()),
            MinorType::COV,
        );

        assert_eq!(rv2147c_diff.minor_mutations.len(), 0);
        assert_eq!(rv2147c_diff.mutations[0].mutation, "1_del_g".to_string());

        // Has insertion at first pos of revcomp gene so check it doesn't make it
        // (because it doesn't make sense)
        let mymt_diff = GeneDifference::new(
            genome.get_gene("mymT".to_string()),
            sample.get_gene("mymT".to_string()),
            MinorType::COV,
        );
        assert_eq!(mymt_diff.mutations.len(), 0);
        assert_eq!(mymt_diff.minor_mutations.len(), 0);
    }

    #[test]
    fn test_snp_after_del() {
        let mut genome = Genome::new("reference/NC_000962.3.gbk");
        let vcf = VCFFile::new("test/tb-snp-after-del.vcf".to_string(), false, 3);
        let mut sample = mutate(&genome, vcf);

        let diff = GenomeDifference::new(genome.clone(), sample.clone(), MinorType::COV);

        // Looks odd with the first deletions, but it lies in 2 genes so is duplicated in the variants
        let expected_variants = [
            "178453_del_ccgccattgggattcatctcgttgccgatcaagatgaaattgagctggctggggctgggagcgttgggacccagcgagatgaggtgctgcatttccagggacgcgatgacggcgctctgcgaatagccgaacacggtgacgtggtttccggcgttgatttgctcccaaatcgcgccgtcgagaatctgtaggcccaactgcaccgaggtttggaagggcagggatttgacgccggtgatcggatatagctcttcgggcgtcaccagcgctttgacgaccggattcgagacgacggggtcgatgaacaaggtcgtgatggcgttgacataactcggcgtgggtatcggtgacccggtgccgcccatgatgatcgccgtattttggttgaacattggcggtagcaccgggggtgaggttggcttaaagagtccggccgtcgcctcctgcaccagcgcgctcgtgttggtggcctcggcattgacaaatgcgtttgcggccgccgccaacctctgggtgaattcgttgtgaaacgccgcaacctgtgcgctgatcgcctggaactgctggccgtacgcgccgaacagcgtggcaagggccgtggacacttcgtccgcggcagccgccgccaggccggttgtcggggccgcgacggccgccgtagcctggttgatcgccgagccgatcccggccaaatcggtagccgccgctgccaataccgacggctgcgcgaatacgtacgacaaaccccatccctccttgtcgacggggcccataacccacccgtcgagccgatacgttgagcgtaaagcgactccgcggttgtgtctggcctttggagtgaacccaaatggggccatgctgcctcgtcattggcgaggtcggtaaacggtagtcggtggacgtcgatgccgtcgggaatccgttaggtgacgaggccctcgatgtttcgaacggtgtccgaggccgccgcgaggagggtgagcaattccacgccgcccgctatcgatcgtgcctaaacctacggtggccgccaggggatagccgatcgcgttgatcagattgcccgcagcgagttgcctgacgaacagttgggtggtgtacagcggcagggtggtgaccagggcgagggcgatgtccacggtgggcagcaggacggcgtagttggttgagatgatcctggcgagcgtgttcaccacctcggccggcgtcggtgcggcggccaccgcggccaccagatcggcgggttgcggcagctggatctgcgggagcgtgagcggttgcgcggacagcgcctgcaggtcggccgtgaagtcaaggatgccttcttgtgttccggcggccagggcatcggcgatgacctgaggcggcacgttcggccacagcccgaacggcgttcgcacatcggcgtagctcgtcgagtagccgtagttcgggtcgccgtagcccaggttgacgatcaccttcaggttcggctggatcaggtcggccagcggatctccgatgaccggcaccgcccgcagcggttgcagcagcggccgattctcggtgcggatgatgtagtagtcggtgacccccgtatagcccggcgacgtcggtaatttagtagcgccctcgacctgcgcgggcgtgaggtccaaatacttggtgtgtacgaatgtgatgcctgcaaccgcgttgaggtcggaaatgaagttgagcgggtatcgcgagaagtcggcgaacccgtcgtactcgagcgtgtagatggccgtcggatagatcgtgtccgagggcgtt",
            "178453_del_ccgccattgggattcatctcgttgccgatcaagatgaaattgagctggctggggctgggagcgttgggacccagcgagatgaggtgctgcatttccagggacgcgatgacggcgctctgcgaatagccgaacacggtgacgtggtttccggcgttgatttgctcccaaatcgcgccgtcgagaatctgtaggcccaactgcaccgaggtttggaagggcagggatttgacgccggtgatcggatatagctcttcgggcgtcaccagcgctttgacgaccggattcgagacgacggggtcgatgaacaaggtcgtgatggcgttgacataactcggcgtgggtatcggtgacccggtgccgcccatgatgatcgccgtattttggttgaacattggcggtagcaccgggggtgaggttggcttaaagagtccggccgtcgcctcctgcaccagcgcgctcgtgttggtggcctcggcattgacaaatgcgtttgcggccgccgccaacctctgggtgaattcgttgtgaaacgccgcaacctgtgcgctgatcgcctggaactgctggccgtacgcgccgaacagcgtggcaagggccgtggacacttcgtccgcggcagccgccgccaggccggttgtcggggccgcgacggccgccgtagcctggttgatcgccgagccgatcccggccaaatcggtagccgccgctgccaataccgacggctgcgcgaatacgtacgacaaaccccatccctccttgtcgacggggcccataacccacccgtcgagccgatacgttgagcgtaaagcgactccgcggttgtgtctggcctttggagtgaacccaaatggggccatgctgcctcgtcattggcgaggtcggtaaacggtagtcggtggacgtcgatgccgtcgggaatccgttaggtgacgaggccctcgatgtttcgaacggtgtccgaggccgccgcgaggagggtgagcaattccacgccgcccgctatcgatcgtgcctaaacctacggtggccgccaggggatagccgatcgcgttgatcagattgcccgcagcgagttgcctgacgaacagttgggtggtgtacagcggcagggtggtgaccagggcgagggcgatgtccacggtgggcagcaggacggcgtagttggttgagatgatcctggcgagcgtgttcaccacctcggccggcgtcggtgcggcggccaccgcggccaccagatcggcgggttgcggcagctggatctgcgggagcgtgagcggttgcgcggacagcgcctgcaggtcggccgtgaagtcaaggatgccttcttgtgttccggcggccagggcatcggcgatgacctgaggcggcacgttcggccacagcccgaacggcgttcgcacatcggcgtagctcgtcgagtagccgtagttcgggtcgccgtagcccaggttgacgatcaccttcaggttcggctggatcaggtcggccagcggatctccgatgaccggcaccgcccgcagcggttgcagcagcggccgattctcggtgcggatgatgtagtagtcggtgacccccgtatagcccggcgacgtcggtaatttagtagcgccctcgacctgcgcgggcgtgaggtccaaatacttggtgtgtacgaatgtgatgcctgcaaccgcgttgaggtcggaaatgaagttgagcgggtatcgcgagaagtcggcgaacccgtcgtactcgagcgtgtagatggccgtcggatagatcgtgtccgagggcgtt",
            "1224300_del_tgcgcctcccgcgagcagacacagaatcgcactgcgccggcccggcgcgtgcgattctgtgtctg",
            "1224367t>c"
        ];
        for (idx, variant) in diff.variants.iter().enumerate() {
            assert_eq!(variant.variant, expected_variants[idx]);
        }
        assert_eq!(diff.minor_variants[0].variant, "178453c>g:3".to_string());

        let rv1096_diff = GeneDifference::new(
            genome.get_gene("Rv1096".to_string()),
            sample.get_gene("Rv1096".to_string()),
            MinorType::COV,
        );
        assert_eq!(rv1096_diff.minor_mutations.len(), 0);
        assert_eq!(
            rv1096_diff.mutations[0].mutation,
            "-85_del_tgcgcctcccgcgagcagacacagaatcgcactgcgccggcccggcgcgtgcgattctgtgtctg".to_string()
        );
        assert_eq!(rv1096_diff.mutations[1].mutation, "t-18c".to_string());

        // Single VCF row should give a large deletion crossing a gene boundary of a revcomp gene
        // while giving a minor call for a SNP at the first position of the deletion
        // It doesn't look like the first position, but revcomp puts the end of the deletion at the start
        let pe1_diff = GeneDifference::new(
            genome.get_gene("PE1".to_string()),
            sample.get_gene("PE1".to_string()),
            MinorType::COV,
        );
        assert_eq!(pe1_diff.mutations.len(), 1);
        assert_eq!(pe1_diff.minor_mutations.len(), 1);

        assert_eq!(pe1_diff.mutations[0].mutation, "-9_del_cgaggcagcatggccccatttgggttcactccaaaggccagacacaaccgcggagtcgctttacgctcaacgtatcggctcgacgggtgggttatgggccccgtcgacaaggagggatggggtttgtcgtacgtattcgcgcagccgtcggtattggcagcggcggctaccgatttggccgggatcggctcggcgatcaaccaggctacggcggccgtcgcggccccgacaaccggcctggcggcggctgccgcggacgaagtgtccacggcccttgccacgctgttcggcgcgtacggccagcagttccaggcgatcagcgcacaggttgcggcgtttcacaacgaattcacccagaggttggcggcggccgcaaacgcatttgtcaatgccgaggccaccaacacgagcgcgctggtgcaggaggcgacggccggactctttaagccaacctcacccccggtgctaccgccaatgttcaaccaaaatacggcgatcatcatgggcggcaccgggtcaccgatacccacgccgagttatgtcaacgccatcacgaccttgttcatcgaccccgtcgtctcgaatccggtcgtcaaagcgctggtgacgcccgaagagctatatccgatcaccggcgtcaaatccctgcccttccaaacctcggtgcagttgggcctacagattctcgacggcgcgatttgggagcaaatcaacgccggaaaccacgtcaccgtgttcggctattcgcagagcgccgtcatcgcgtccctggaaatgcagcacctcatctcgctgggtcccaacgctcccagccccagccagctcaatttcatcttgatcggcaacgagatgaatcccaatggcg".to_string());
        assert_eq!(pe1_diff.minor_mutations[0].mutation, "G286A:3".to_string());
    }
}
