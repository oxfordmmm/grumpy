extern crate gb_io;

use std::fs::File;
use std::string::String;
use std::vec::Vec;

use gb_io::seq::Reference;
use string_cache::Atom;

use gb_io::reader::SeqReader;
use gb_io::seq::Location::Complement;
use gb_io::seq::Location::Range;
use gb_io::seq::Location::Join;

#[derive(Clone)]
struct GeneDef{
    name: String,
    coding: bool,
    rev_comp: bool,
    start: i64,
    end: i64,
    promoter_start: i64,
    promoter_size: i64,
    ribosomal_shifts: Vec<i64>
}

#[derive(Clone)]
struct Genome{
    name: String,
    nucleotide_sequence: String,
    gene_definitions: Vec<GeneDef>,
    genome_positions: Vec<GenomePosition>
}

#[derive(Clone, Debug)]
struct Evidence{

}

#[derive(Clone, Debug)]
struct Alt{
    base: String,
    cov: Option<i32>,
    frs: Option<f32>,
    evidence: Evidence
}

#[derive(Clone, Debug)] 
struct GenomePosition{
    reference: char,
    alts: Vec<Alt>,
    genome_idx: i64, // 1-indexed genome index
    genes: Vec<String>,
    is_deleted: bool
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
                Ok(s) => s.to_uppercase(),
                Err(e) => panic!("Problem reading sequence data: {:?}", e)
            };
            genome_name = seq.name.unwrap();
            for feature in seq.features {
                let mut name: String = "".to_string();
                let start: i64;
                let end: i64;
                let mut coding: bool = false;
                let mut rev_comp: bool = false;
                let mut ribosomal_shifts: Vec<i64> = Vec::new();
                if feature.kind == Atom::from("CDS") || feature.kind == Atom::from("rRNA"){
                    if feature.kind != Atom::from("rRNA"){
                        coding = true;
                    }
                    match feature.location {
                        Complement(x) => {
                            rev_comp = true;
                            match *x {
                                Range(s, e) => {
                                    start = e.0;
                                    end = s.0;
                                },
                                _ => continue
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
                                    _ => continue
                                }
                            }
                            start = start_pos;
                            end = end_pos;
                        },
                        _ => continue
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
                        rev_comp,
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
        return Genome{
            name: genome_name,
            nucleotide_sequence: _nucleotide_sequence,
            gene_definitions: _gene_definitions,
            genome_positions: genome_positions
        }
    }

    pub fn assign_promoters(&mut self) -> (){
        /* 
        Assigns promoters to genes iteratively. 
        Expand out to a given distance from the start of each gene without overlapping with other genes.
        */
        let max_promoter_length = 100;
        for gene in self.gene_definitions.iter_mut(){
            // First pass to add gene names to all positions they exist in
            let mut start_idx = gene.start;
            let mut end_idx = gene.end;
            if gene.rev_comp{
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
                if gene.rev_comp{
                    // This pushes `gene.promoter_start += expanding` the right way for rev_comp
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

    pub fn at_genome_index(&self, index: i64) -> GenomePosition{
        // 1-indexed genome index
        return self.genome_positions[(index + 1) as usize].clone();
    }
}

fn main() {
    let mut reference = Genome::new("reference/MN908947.3.gb");
    reference.assign_promoters();
    for gene in reference.gene_definitions{
        if gene.name == "rpoB" || gene.name == "katG" || gene.name == "orf1ab" {
            println!("Name {}", gene.name);
            println!("Is reverse complement {}", gene.rev_comp);
            println!("Is coding {}", gene.coding);
            println!("Start pos {}", gene.start);
            println!("End pos {}", gene.end);
            println!("Promoter start pos {}", gene.promoter_start);
            println!("Promoter size {}", gene.promoter_size);
            println!("Ribosomal shifts {:?}", gene.ribosomal_shifts);
            println!("");
        }
    }
}