use std::string::String;
use std::vec::Vec;

use crate::common::{Evidence, AltType, GeneDef};
use crate::genome::GenomePosition;

#[derive(Clone, Debug)]
pub struct GeneAlt{
    alt_type: AltType,
    nucleotides: String,
    codon: Option<String>,
    amino_acid: Option<char>,
    // Because codons, this is a 1-many relationship
    cov: Vec<Option<i32>>,
    frs: Vec<Option<f32>>,
    evidence: Vec<Evidence>
}

#[derive(Clone, Debug)]
pub struct GenePosition{
    // Indexed by gene position
    pub reference: String, // Amino acid or nucleotide depending on gene coding
    pub reference_nucleotides: String, // Nucleotide or codon depending on gene coding
    pub gene_position: i64, // 1-indexed gene position
    pub nucleotide_number: Vec<i64>, // 1-indexed nucleotide position
    pub nucleotide_index: Vec<i64>, // 1-indexed genome position
    pub alts: Vec<GeneAlt>,
    pub is_deleted: bool
}

#[derive(Clone, Debug)]
pub struct Gene{
    pub name: String,
    pub coding: bool,
    pub reverse_complement: bool,
    pub nucleotide_sequence: String,
    pub nucleotide_index: Vec<i64>,
    pub nucleotide_number: Vec<i64>,
    pub gene_number: Vec<i64>,
    pub gene_positions: Vec<GenePosition>,
    pub amino_acid_sequence: String,
    pub amino_acid_number: Vec<i64>,
    pub ribosomal_shifts: Vec<i64>,
    pub codons: Vec<String>
}

impl Gene {
    pub fn new(gene_def: GeneDef, nc_sequence: String, nc_index: Vec<i64>, genome_positions: Vec<&GenomePosition>) -> Self {
        let mut nucleotide_sequence = nc_sequence.clone();
        let mut nucleotide_index = nc_index.clone();
        let mut nucleotide_number = Vec::new();
        let mut amino_acid_sequence = "".to_string();
        let mut amino_acid_number = Vec::new();
        let mut gene_number: Vec<i64> = Vec::new();
        let mut codons = Vec::new();
        let mut gene_positions: Vec<GenePosition> = Vec::new();

        for pos in gene_def.ribosomal_shifts.iter(){
            // Figure out the index of the vectors to insert the ribosomal shift
            let mut idx: i64 = -1;
            for i in 0..nucleotide_index.len(){
                if nucleotide_index[i] == *pos{
                    idx = i as i64;
                    break;
                }
            }
            if idx == -1{
                panic!("Ribosomal shift position {} not found in gene {}", pos, gene_def.name);
            }
            // Duplicate those indices
            nucleotide_sequence.insert(idx as usize, nucleotide_sequence.chars().nth(idx as usize).unwrap());
            nucleotide_index.insert(idx as usize, nucleotide_index[idx as usize]);
        }

        if gene_def.reverse_complement{
            // Reverse complement the sequence
            nucleotide_sequence = nc_sequence.chars().rev().map(|x| complement_base(x)).collect::<String>();
            nucleotide_index = Vec::new();
            for i in nc_index.iter().rev(){
                nucleotide_index.push(*i);
            }
        }

        // Figure out the nucelotide number for each position
        // Promoter first
        if gene_def.promoter_start != -1{
            let mut nc_idx = 0;
            for i in (-1*(gene_def.promoter_size + 1))..0{
                nucleotide_number.push(i);
                gene_number.push(i);
                gene_positions.push(GenePosition{
                    reference: nucleotide_sequence.chars().nth(nc_idx).unwrap().to_string(),
                    reference_nucleotides: nucleotide_sequence.chars().nth(nc_idx).unwrap().to_string(),
                    nucleotide_number: vec![nucleotide_index[nc_idx]],
                    nucleotide_index: vec![nucleotide_index[nc_idx]],
                    gene_position: i,
                    alts: Vec::new(),
                    is_deleted: false
                });
                nc_idx += 1;
            }
        }
        let prom_end = nucleotide_number.len();
        // Now non-promoter
        let mut nc_idx = prom_end;
        for i in 1..(gene_def.start - gene_def.end).abs() + gene_def.ribosomal_shifts.len() as i64 + 1{
            nucleotide_number.push(i);
            if !gene_def.coding{
                // No adjustment needed for non-coding as gene pos == nucleotide num
                gene_number.push(i);
                gene_positions.push(GenePosition{
                    reference: nucleotide_sequence.chars().nth(nc_idx).unwrap().to_string(),
                    reference_nucleotides: nucleotide_sequence.chars().nth(nc_idx).unwrap().to_string(),
                    nucleotide_number: vec![nucleotide_index[nc_idx]],
                    nucleotide_index: vec![nucleotide_index[nc_idx]],
                    gene_position: i, 
                    alts: Vec::new(),
                    is_deleted: false
                });
                nc_idx += 1;
            }
        }

        if gene_def.coding{
            // Now figure out the amino acid sequence from the nucleotide sequence
            let mut codon = "".to_string();
            let mut codon_idx = 0;
            for i in prom_end..nucleotide_sequence.len(){
                codon.push(nucleotide_sequence.chars().nth(i).unwrap());
                if codon.len() == 3{
                    // Codon is complete
                    amino_acid_sequence.push(codon_to_aa(codon.clone()));
                    codons.push(codon.clone());
                    codon_idx += 1;
                    gene_number.push(codon_idx);
                    amino_acid_number.push(codon_idx);
                    gene_positions.push(GenePosition{
                        reference: codon_to_aa(codon.clone()).to_string(),
                        reference_nucleotides: codon.clone(),
                        nucleotide_number: vec![nucleotide_number[i-2], nucleotide_number[i-1], nucleotide_number[i]],
                        nucleotide_index: vec![nucleotide_index[i-2], nucleotide_index[i-1], nucleotide_index[i]],
                        gene_position: codon_idx,
                        alts: Vec::new(),
                        is_deleted: false
                    });
                    codon = "".to_string();
                }
            }
            if codon.len() > 0{
                panic!("Incomplete codon at end of gene {}", gene_def.name);
            }
        }

        return Gene{
            name: gene_def.name.to_string(),
            coding: gene_def.coding,
            reverse_complement: gene_def.reverse_complement,
            nucleotide_sequence: nucleotide_sequence.to_string(),
            nucleotide_index,
            gene_number,
            gene_positions,
            nucleotide_number,
            amino_acid_sequence,
            amino_acid_number,
            ribosomal_shifts: gene_def.ribosomal_shifts,
            codons
        }
    }

}

fn codon_to_aa(codon: String) -> char{
    if codon.contains("x"){
        return 'X';
    }
    if codon.contains("z"){
        return 'Z';
    }
    match codon.as_str(){
        "ttt" | "ttc" => 'F',
        "tta" | "ttg" | "ctt" | "ctc" | "cta" | "ctg" => 'L',
        "tct" | "tcc" | "tca" | "tcg" | "agt" | "agc" => 'S',
        "tat" | "tac" => 'Y',
        "taa" | "tag" | "tga" => '!',
        "tgt" | "tgc" => 'C',
        "tgg" => 'W',
        "cct" | "ccc" | "cca" | "ccg" => 'P',
        "cat" | "cac" => 'H',
        "caa" | "cag" => 'Q',
        "cgt" | "cgc" | "cga" | "cgg" | "aga" | "agg" => 'R',
        "att" | "atc" | "ata" => 'I',
        "atg" => 'M',
        "act" | "acc" | "aca" | "acg" => 'T',
        "aat" | "aac" => 'N',
        "aaa" | "aag" => 'K',
        "gtt" | "gtc" | "gta" | "gtg" => 'V',
        "gct" | "gcc" | "gca" | "gcg" => 'A',
        "gat" | "gac" => 'D',
        "gaa" | "gag" => 'E',
        "ggt" | "ggc" | "gga" | "ggg" => 'G',
        _ => panic!("Invalid codon {}", codon)
    }
}

fn complement_base(base: char) -> char{
    match base {
        'a' => 't',
        't' => 'a',
        'c' => 'g',
        'g' => 'c',
        // Het and null don't get complements
        'z' => 'z',
        'x' => 'x',
        _ => base
    }
}
