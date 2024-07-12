use std::fs::File;

use string_cache::Atom;

use gb_io::reader::SeqReader;
use gb_io::seq::Location::Complement;
use gb_io::seq::Location::Range;
use gb_io::seq::Location::Join;

use crate::common::{Alt, GeneDef, AltType};
use crate::gene::Gene;
use crate::vcf::VCFFile;

#[derive(Clone, Debug)] 
pub struct GenomePosition{
    // Updated during mutation
    pub reference: char,
    pub is_deleted: bool,

    // Used to store calls from VCF
    // If populated, alts[0] is the actual call, others are minor calls
    pub alts: Vec<Alt>,

    pub genome_idx: i64, // 1-indexed genome index
    pub genes: Vec<String>
}

#[derive(Clone)]
pub struct Genome{
    pub name: String,
    pub nucleotide_sequence: String,
    pub gene_definitions: Vec<GeneDef>,
    pub genome_positions: Vec<GenomePosition>
}

impl Genome{
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
                Err(e) => panic!("Problem reading sequence data: {:?}", e)
            };
            genome_name = seq.name.unwrap();
            for feature in seq.features {
                let mut name: String = "".to_string();
                let start: i64;
                let end: i64;
                let mut coding: bool = false;
                let mut reverse_complement: bool = false;
                let mut ribosomal_shifts: Vec<i64> = Vec::new();
                if feature.kind == Atom::from("CDS") || feature.kind == Atom::from("rRNA"){
                    if feature.kind != Atom::from("rRNA"){
                        coding = true;
                    }
                    match feature.location {
                        Complement(x) => {
                            reverse_complement = true;
                            match *x {
                                Range(s, e) => {
                                    start = e.0;
                                    end = s.0;
                                },
                                _ => panic!("Complement location not range")
                            }
                        },
                        Range(s, e) => {
                            start = s.0;
                            end = e.0
                        },
                        Join(ranges) => {
                            // Checking for PRFS
                            let mut start_pos: i64 = 0;
                            let mut end_pos = 0;
                            let mut first = true;
                            for range in ranges{
                                match range {
                                    Range(s, e) => {
                                        if first{
                                            start_pos = s.0;
                                            first = false;
                                        }
                                        else{
                                            ribosomal_shifts.push(s.0);
                                        }
                                        end_pos = e.0;
                                    },
                                    _ => panic!("Join location not range")
                                }
                            }
                            start = start_pos;
                            end = end_pos;
                        },
                        _ => panic!("Location not range, complement or join")
                    }
                    for qual in feature.qualifiers{
                        match qual.1 {
                            Some(val) => {
                                if qual.0 == Atom::from("gene"){
                                    // Use gene name if it exists
                                    name = val.clone();
                                }
                                if name == "" && qual.0 == Atom::from("locus_tag"){
                                    // Else default to locus tag
                                    name = val.clone();
                                }
                            },
                            None => continue
                        }
                    }
                    while gene_names.contains(&name){
                        // Duplicate gene names can exist :(
                        // Repeatedly add _2 to the end of the name until it is unique
                        // not ideal but exact mirror of gumpy
                        name += "_2";
                    }
                    gene_names.push(name.clone());
                    _gene_definitions.push(GeneDef{
                        name,
                        reverse_complement,
                        coding,
                        start,
                        end,
                        promoter_start: -1,
                        promoter_size: 0,
                        ribosomal_shifts
                    });
                }
            }
        }
        let mut genome_positions = Vec::new();
        let mut genome_idx = 0;
        for c in _nucleotide_sequence.chars(){
            genome_idx += 1;
            genome_positions.push(
                GenomePosition{
                    reference: c,
                    genome_idx,
                    alts: Vec::new(),
                    genes: Vec::new(),
                    is_deleted: false
                }
            );
        }
        let mut genome = Genome{
            name: genome_name,
            nucleotide_sequence: _nucleotide_sequence,
            gene_definitions: _gene_definitions,
            genome_positions: genome_positions
        };
        genome.assign_promoters();
        return genome;
    }

    pub fn assign_promoters(&mut self) -> (){
        /* 
        Assigns promoters to genes iteratively. 
        Expand out to a given distance from the start of each gene without overlapping with other genes.
        */
        let max_promoter_length = 100 - 1;
        for gene in self.gene_definitions.iter_mut(){
            // First pass to add gene names to all positions they exist in
            let mut start_idx = gene.start;
            let mut end_idx = gene.end;
            if gene.reverse_complement{
                start_idx = gene.end;
                end_idx = gene.start;
            }
            for i in start_idx..end_idx{
                self.genome_positions[i as usize].genes.push(gene.name.clone());
            }
        }
        for gene in self.gene_definitions.iter_mut(){
            // Check for overlapping genes, assigning promoters to non-overlapping genes
            if self.genome_positions[gene.start as usize].genes.len() > 1{
                continue;
            }
            else{
                gene.promoter_start = gene.start;
            }
        }

        let mut complete = false;
        while !complete{
            let mut this_complete = true;
            for gene in self.gene_definitions.iter_mut(){
                let mut expanding = -1;
                if gene.reverse_complement{
                    // This pushes `gene.promoter_start += expanding` the right way for reverse_complement
                    expanding = 1;
                }

                if gene.promoter_start == -1 || gene.promoter_start == 0 || gene.promoter_size == max_promoter_length{
                    continue;
                }
                if self.genome_positions[(gene.promoter_start + expanding) as usize].genes.len() > 0{
                    // Pre-existing gene at this position so ignore expanding!
                    // println!("Gene {} would overlap with gene {}", gene.name, self.genome_positions[gene.promoter_start as usize].genes[0]);
                    continue;
                }
                else{
                    // No gene in the new position so expand
                    gene.promoter_start += expanding;
                    gene.promoter_size += 1;
                    self.genome_positions[gene.promoter_start as usize].genes.push(gene.name.clone());
                    this_complete = false;
                }
            }
            complete = this_complete;
        }
    }

    pub fn build_gene(&self, gene_name: String) -> Gene{
        let mut valid = false;
        let mut maybe_gene_def: Option<GeneDef> = None;
        for gene in self.gene_definitions.iter(){
            if gene.name == gene_name{
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
        let mut genome_positions: Vec<&GenomePosition> = Vec::new();
        if gene_def.reverse_complement{
            let mut last_idx = gene_def.promoter_start;
            if gene_def.promoter_start == -1{
                last_idx = gene_def.start;
            }
            for i in gene_def.end..last_idx+1{
                nucleotide_sequence.push(self.genome_positions[i as usize].reference);
                nucleotide_index.push(self.genome_positions[i as usize].genome_idx);
                genome_positions.push(&self.genome_positions[i as usize]);
            }
        }
        else{
            let mut first_idx = gene_def.promoter_start - 1;
            if gene_def.promoter_start == -1{
                first_idx = gene_def.start;
            }
            for i in first_idx..gene_def.end{
                nucleotide_sequence.push(self.genome_positions[i as usize].reference);
                nucleotide_index.push(self.genome_positions[i as usize].genome_idx);
                genome_positions.push(&self.genome_positions[i as usize]);
            }
        
        }

        return Gene::new(
            gene_def,
            nucleotide_sequence,
            nucleotide_index,
            genome_positions
        );

    }

    pub fn at_genome_index(&self, index: i64) -> GenomePosition{
        // 1-indexed genome index
        return self.genome_positions[(index + 1) as usize].clone();
    }
}

pub fn mutate(reference: &Genome, vcf: VCFFile) -> Genome{
    let mut new_genome = reference.clone();
    let mut idx = 0;
    for position in new_genome.genome_positions.iter_mut(){
        if vcf.calls.contains_key(&position.genome_idx){
            // Set the new base etc from the call
            let call = vcf.calls.get(&position.genome_idx).unwrap().clone()[0].clone();
            position.alts.push(Alt{
                alt_type: call.clone().call_type,
                base: call.clone().alt,
                evidence: call.clone()
            });

            if call.call_type == AltType::DEL{
                position.is_deleted = true;
            }

            if call.call_type == AltType::SNP || call.call_type == AltType::HET || call.call_type == AltType::NULL{
                // Update nucleotide for SNP/het/null
                position.reference = call.clone().alt.chars().nth(0).unwrap();
                new_genome.nucleotide_sequence.replace_range(idx..idx+1, &call.clone().alt.chars().nth(0).unwrap().to_string());
            }

        }
        idx += 1;
    }

    return new_genome;
}