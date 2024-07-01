extern crate gb_io;

use std::fs::File;
use std::string::String;
use std::vec::Vec;

use string_cache::Atom;

use gb_io::reader::SeqReader;
use gb_io::seq::Location::Complement;
use gb_io::seq::Location::Range;

struct GeneDef{
    name: String,
    coding: bool,
    rev_comp: bool,
    start: i64,
    end: i64,
    promoter_start: i64
}

struct Genome{
    name: String,
    nucleotide_sequence: String,
    gene_definitions: Vec<GeneDef>,
    genome_positions: Vec<GenomePosition>
}

struct Evidence{

}

struct Alt{
    base: String,
    cov: Option<i32>,
    frs: Option<f32>,
    evidence: Evidence
}

struct GenomePosition{
    reference: char,
    alts: Vec<Alt>
}

impl Genome{
    pub fn new(filename: &str) -> Self {
        let file = File::open(filename).unwrap();
        let mut _gene_definitions = Vec::new();
        let mut _nucleotide_sequence: String = "".to_string();
        let mut genome_name: String = "".to_string();
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
                if feature.kind == Atom::from("CDS") || feature.kind == Atom::from("rRNA"){
                    if feature.kind == Atom::from("rRNA"){
                        coding = true;
                    }
                    match feature.location {
                        Complement(x) => {
                            rev_comp = true;
                            match *x {
                                Range(s, e) => {
                                    start = s.0;
                                    end = e.0;
                                },
                                _ => continue
                            }
                        },
                        Range(s, e) => {
                            start = s.0;
                            end = e.0
                        },
                        _ => continue
                    }
                    for qual in feature.qualifiers{
                        let val = qual.1.unwrap();
                        if qual.0 == Atom::from("gene"){
                            // Use gene name if it exists
                            name = val.clone();
                        }
    
                        if name == "" && qual.0 == Atom::from("locus_tag"){
                            // Else default to locus tag
                            name = val.clone();
                        }
                    }
                    _gene_definitions.push(GeneDef{
                        name,
                        rev_comp,
                        coding,
                        start,
                        end,
                        promoter_start: -1
                    });
                }
            }
        }
        let mut genome_positions = Vec::new();
        for c in _nucleotide_sequence.chars(){
            genome_positions.push(
                GenomePosition{
                    reference: c,
                    alts: Vec::new()
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

    pub fn at_genome_index(&self, index: i64) -> Result<GenomePosition, String>{
        return self.genome_positions.get(index)

    }
}

fn main() {
    let reference = Genome::new("reference/NC_000962.3.gbk");
    for gene in reference.gene_definitions{
        if gene.name == "rpoB"{
            println!("Name {}", gene.name);
            println!("Is reverse complement {}", gene.rev_comp);
            println!("Is coding {}", gene.coding);
            println!("Start pos {}", gene.start);
            println!("End pos {}", gene.end);
            println!("");
        }
    }
}