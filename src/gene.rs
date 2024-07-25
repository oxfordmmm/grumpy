//! Module for handling gene data
use pyo3::prelude::*;
use std::string::String;
use std::vec::Vec;
use std::collections::HashMap;

use crate::common::{Evidence, AltType, GeneDef, Alt};
use crate::genome::GenomePosition;

#[pyclass(eq)]
#[derive(Clone, Debug, Eq, PartialEq)]
/// Each position of a gene can either be a nucleotde or a codon
pub enum GenePos{
    /// Nucelotide
    Nucleotide(NucleotideType),

    /// Codon
    Codon(CodonType)
}

#[pyclass(eq)]
#[derive(Clone, Debug, Eq, PartialEq)]
/// Tracks each constituent nucleotide in the codon, along with the amino acid it codes for
pub struct CodonType{
    #[pyo3(get, set)]
    /// Amino acid the codon codes for
    pub amino_acid: char,

    #[pyo3(get, set)]
    /// 3-tuple for each nucleotide in the codon
    pub codon: Vec<NucleotideType>,
}

#[pyclass(eq)]
#[derive(Clone, Debug, Eq, PartialEq)]
/// Stores information about a single nucleotide in a gene
pub struct NucleotideType{
    #[pyo3(get, set)]
    /// Nucleotide at this position
    pub reference: char,

    #[pyo3(get, set)]
    /// Gene position of this nucleotide. 1-indexed from gene start
    pub nucleotide_number: i64,

    #[pyo3(get, set)]
    /// Genome position of this nucleotide. 1-indexed from genome start
    pub nucleotide_index: i64,

    #[pyo3(get, set)]
    /// Alts at this position
    pub alts: Vec<Alt>,

    #[pyo3(get, set)]
    /// Whether this position is deleted
    pub is_deleted: bool,

    #[pyo3(get, set)]
    /// Whether this position is deleted in a minor allele
    pub is_deleted_minor: bool
}

#[pyclass(eq)]
#[derive(Clone, Debug, Eq, PartialEq)]
/// A position of a gene is some position in the gene, along with the data at that position
pub struct GenePosition{
    #[pyo3(get, set)]
    /// 1-indexed gene position
    pub gene_position: i64,

    #[pyo3(get, set)]
    /// Data at this position
    pub gene_position_data: GenePos,
}

#[pyclass]
#[derive(Clone, Debug)]
/// A gene is a collection of gene positions, along with some metadata
pub struct Gene{
    #[pyo3(get, set)]
    /// Name of the gene
    pub name: String,

    #[pyo3(get, set)]
    /// Whether the gene codes protein
    pub coding: bool,

    #[pyo3(get, set)]
    /// Whether this gene is reverse complement
    pub reverse_complement: bool,

    #[pyo3(get, set)]
    /// Nucleotide sequence of the gene
    pub nucleotide_sequence: String,

    #[pyo3(get, set)]
    /// Genome index of each nucleotide of the gene. 1-indexed from genome start
    pub nucleotide_index: Vec<i64>,

    #[pyo3(get, set)]
    /// Gene positions of the nucleotides of the gene. 1-indexed from gene start
    pub nucleotide_number: Vec<i64>,

    #[pyo3(get, set)]
    /// Gene positions of the gene, using amino acid numbers where appropriate. 1-indexed from gene start
    pub gene_number: Vec<i64>,

    #[pyo3(get, set)]
    /// Vec of gene positions. Stores data about each gene position.
    pub gene_positions: Vec<GenePosition>,

    #[pyo3(get, set)]
    /// Sequence of amino acid residues coded for by the gene. Empty if gene is non-coding
    pub amino_acid_sequence: String,

    #[pyo3(get, set)]
    /// Sequence of amino acid numbers coded for by the gene. Empty if gene is non-coding
    pub amino_acid_number: Vec<i64>,

    #[pyo3(get, set)]
    /// Positions of genome indicies duplicated by ribosomal shifts
    pub ribosomal_shifts: Vec<i64>,

    #[pyo3(get, set)]
    /// Codons of the gene. Empty if gene is non-coding
    pub codons: Vec<String>,

    #[pyo3(get, set)]
    /// Map of genome index -> (gene number, optional codon position>)
    pub genome_idx_map: HashMap<i64, (i64, Option<i64>)>,
}

impl Gene {
    /// Instantiates a new gene.
    /// 
    /// Recommended to call Genome.build_gene() instead of this directly as arguments are auto-populated
    /// 
    /// # Arguments
    /// - `gene_def`: GeneDef struct containing gene metadata
    /// - `nc_sequence`: Nucleotide sequence of the gene
    /// - `nc_index`: Genome index of each nucleotide of the gene. 1-indexed from genome start
    /// - `_genome_positions`: Genome positions of the gene
    pub fn new(gene_def: GeneDef, nc_sequence: String, nc_index: Vec<i64>, _genome_positions: Vec<GenomePosition>) -> Self {
        let mut nucleotide_sequence = nc_sequence.clone();
        let mut nucleotide_index = nc_index.clone();
        let mut genome_positions = _genome_positions.clone();
        let mut nucleotide_number = Vec::new();
        let mut amino_acid_sequence = "".to_string();
        let mut amino_acid_number = Vec::new();
        let mut gene_number: Vec<i64> = Vec::new();
        let mut codons = Vec::new();
        let mut gene_positions: Vec<GenePosition> = Vec::new();

        // Map of genome index -> gene number
        let mut genome_idx_map: HashMap<i64, (i64, Option<i64>)> = HashMap::new();


        // Ensure we pick up deletions upstream or downstream of the gene
        Gene::adjust_dels(&mut genome_positions, gene_def.name.clone());


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
            genome_positions.insert(idx as usize, genome_positions[idx as usize].clone());
        }

        if gene_def.reverse_complement{
            // Reverse complement the sequence
            nucleotide_sequence = nc_sequence.chars().rev().map(|x| complement_base(x)).collect::<String>();
            nucleotide_index = Vec::new();
            for i in nc_index.iter().rev(){
                nucleotide_index.push(*i);
            }

            // Genome positions are a bit more annoying here, specifically for indels
            // With revcomp, deletions start from their end position as well as being complemented
            // Insertions are also annoying. They need a position adjustment as they're starting on the other side now
            genome_positions = Vec::new();
            for pos in _genome_positions.iter().rev(){
                genome_positions.push(pos.clone());
            }
            let mut fixed_genome_positions = genome_positions.clone();
            // First figure out where the indels are
            let mut indel_positions: Vec<(usize, i64, bool)> = Vec::new();
            for i in 0..genome_positions.len(){
                for alt in genome_positions[i].alts.iter(){
                    if alt.alt_type == AltType::INS{
                        indel_positions.push((i, alt.base.len() as i64, alt.evidence.is_minor));
                    }
                    if alt.alt_type == AltType::DEL{
                        indel_positions.push((i, alt.base.len() as i64 * -1, alt.evidence.is_minor));
                    }
                }
            }

            // Now adjust the positions
            for (pos, indel_size, is_minor) in indel_positions.iter_mut(){
                if *indel_size > 0{
                    // Insertion
                    if *pos == 0{
                        panic!("Insertion at start of gene {} is revcomp and cannot be adjusted", gene_def.name);
                    }
                    let new_pos = *pos - 1;
                    let fixed_alts = genome_positions[*pos].alts.iter().map(|x| Gene::rev_comp_indel_alt(x, i64::MAX)).collect::<Vec<Alt>>();

                    // Remove the insertion from the old position
                    fixed_genome_positions[*pos].alts = genome_positions[*pos].alts.iter().filter(|x| x.alt_type != AltType::INS && x.evidence.is_minor == *is_minor).map(|x| x.clone()).collect::<Vec<Alt>>();


                    // Update the new position
                    fixed_genome_positions[new_pos].alts = fixed_alts;
                }
                else{
                    // Deletion
                    let mut _max_del_length = i64::MAX;
                    if indel_size.abs() as usize > *pos{
                        // Deletion may start in this gene, but needs truncating to just the part in this gene
                        _max_del_length = *pos as i64;
                    }
                    let mut new_pos = 0;
                    if _max_del_length == i64::MAX{
                        // Update new position if deletion is entirely in this gene
                        new_pos = *pos - indel_size.abs() as usize + 1;
                    }
                    let fixed_alts = genome_positions[*pos].alts.iter().map(|x| Gene::rev_comp_indel_alt(x, _max_del_length)).collect::<Vec<Alt>>();

                    // Remove the deletion from the old position
                    fixed_genome_positions[*pos].alts = genome_positions[*pos].alts.iter().filter(|x| x.alt_type != AltType::DEL && x.evidence.is_minor == *is_minor).map(|x| x.clone()).collect::<Vec<Alt>>();

                    // Update the new position
                    fixed_genome_positions[new_pos].alts = fixed_alts;
                }
            }

            // Revcomp other items in the genome positions
            for position in fixed_genome_positions.iter_mut(){
                for alt in position.alts.iter_mut(){
                    *alt = Gene::rev_comp_other_alt(alt);
                }
            }

            genome_positions = fixed_genome_positions;

        }

        // Figure out the nucelotide number for each position
        // Promoter first
        if gene_def.promoter_start != -1{
            let mut nc_idx = 0;
            for i in (-1*(gene_def.promoter_size + 1))..0{
                nucleotide_number.push(i);
                gene_number.push(i);
                gene_positions.push(GenePosition{
                    gene_position_data: GenePos::Nucleotide(NucleotideType{
                        reference: nucleotide_sequence.chars().nth(nc_idx).unwrap(),
                        nucleotide_number: i,
                        nucleotide_index: nucleotide_index[nc_idx],
                        alts: genome_positions[nc_idx].alts.clone(),
                        is_deleted: genome_positions[nc_idx].is_deleted,
                        is_deleted_minor: genome_positions[nc_idx].is_deleted_minor
                    }),
                    gene_position: i
                });
                genome_idx_map.insert(nucleotide_index[nc_idx], (i, None));
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
                    gene_position_data: GenePos::Nucleotide(NucleotideType{
                        reference: nucleotide_sequence.chars().nth(nc_idx).unwrap(),
                        nucleotide_number: i,
                        nucleotide_index: nucleotide_index[nc_idx],
                        alts: genome_positions[nc_idx].alts.clone(),
                        is_deleted: genome_positions[nc_idx].is_deleted,
                        is_deleted_minor: genome_positions[nc_idx].is_deleted_minor
                    }),
                    gene_position: i, 
                });
                genome_idx_map.insert(nucleotide_index[nc_idx], (i, None));
                nc_idx += 1;
            }
        }

        if gene_def.coding{
            // Now figure out the amino acid sequence from the nucleotide sequence
            let mut codon = "".to_string();
            let mut codon_idx = 1;
            for i in prom_end..nucleotide_sequence.len(){
                codon.push(nucleotide_sequence.chars().nth(i).unwrap());
                genome_idx_map.insert(nucleotide_index[i], (codon_idx, Some((i % 3) as i64)));
                if codon.len() == 3{
                    // Codon is complete
                    amino_acid_sequence.push(codon_to_aa(codon.clone()));
                    codons.push(codon.clone());
                    gene_number.push(codon_idx);
                    amino_acid_number.push(codon_idx);
                    gene_positions.push(GenePosition{
                        gene_position_data: GenePos::Codon(CodonType{
                            amino_acid: codon_to_aa(codon.clone()),
                            codon: vec![
                                NucleotideType{
                                    reference: codon.chars().nth(0).unwrap(),
                                    nucleotide_number: nucleotide_number[i-2],
                                    nucleotide_index: nucleotide_index[i-2],
                                    alts: genome_positions[i-2].alts.clone(),
                                    is_deleted: genome_positions[i-2].is_deleted,
                                    is_deleted_minor: genome_positions[i-2].is_deleted_minor
                                },
                                NucleotideType{
                                    reference: codon.chars().nth(1).unwrap(),
                                    nucleotide_number: nucleotide_number[i-1],
                                    nucleotide_index: nucleotide_index[i-1],
                                    alts: genome_positions[i-1].alts.clone(),
                                    is_deleted: genome_positions[i-1].is_deleted,
                                    is_deleted_minor: genome_positions[i-1].is_deleted_minor

                                },
                                NucleotideType{
                                    reference: codon.chars().nth(2).unwrap(),
                                    nucleotide_number: nucleotide_number[i],
                                    nucleotide_index: nucleotide_index[i],
                                    alts: genome_positions[i].alts.clone(),
                                    is_deleted: genome_positions[i].is_deleted,
                                    is_deleted_minor: genome_positions[i].is_deleted_minor

                                }
                            ]
                        }),
                        gene_position: codon_idx,
                    });
                    codon_idx += 1;
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
            codons,
            genome_idx_map
        }
    }

    /// Perform alterations required to an indel when reverse complementing
    /// 
    /// # Arguments
    /// - `alt`: Alt to reverse complement
    /// - `max_del_length`: Maximum length of a deletion
    /// 
    /// # Returns
    /// - Reverse complemented Alt
    fn rev_comp_indel_alt(alt: &Alt, max_del_length: i64) -> Alt{
        if alt.alt_type == AltType::INS || alt.alt_type == AltType::DEL{
            let mut new_base = "".to_string();
            let mut indel_length = 0;
            for c in alt.base.chars().rev(){
                if indel_length < max_del_length{
                    new_base.push(complement_base(c));
                }
                indel_length += 1;
            }
            return Alt{
                alt_type: alt.alt_type.clone(),
                base: new_base,
                evidence: alt.evidence.clone()
            }
        }
        return alt.clone();
    }

    /// Perform alterations required to a SNP when reverse complementing
    /// 
    /// # Arguments
    /// - `alt`: Alt to reverse complement
    /// 
    /// # Returns
    /// - Reverse complemented Alt
    fn rev_comp_other_alt(alt: &Alt) -> Alt{
        if alt.alt_type != AltType::INS && alt.alt_type != AltType::DEL{
            let mut new_base = "".to_string();
            for c in alt.base.chars().rev(){
                new_base.push(complement_base(c));
            }
            return Alt{
                alt_type: alt.alt_type.clone(),
                base: new_base,
                evidence: alt.evidence.clone()
            }
        }
        return alt.clone();
    }

    /// Adjust deletions to ensure they don't cross gene boundaries
    /// 
    /// If a deletion starts upstream of the gene, 
    /// this will truncate it to the part in the gene, ensuring a valid Alt at the gene position
    /// 
    /// # Arguments
    /// - `genome_positions`: Vec of genome positions for the gene
    /// - `gene_name`: Name of the gene
    fn adjust_dels(genome_positions: &mut Vec<GenomePosition>, gene_name: String){
        // Adjust any deletions which may cross gene boundaries as required
        let mut first_pos = genome_positions[0].clone();
        if genome_positions[0].is_deleted{
            // Double check if this was actually the start of a deletion
            if !genome_positions[0].alts.iter().any(|x| x.alt_type == AltType::DEL && !x.evidence.is_minor){
                // None of the alts at this position are deletions, so didn't start here
                // Lets look at the deleted evidence to piece together what this should be
                let del_evidence = genome_positions[0].deleted_evidence.iter().filter(|x| !x.is_minor).collect::<Vec<&Evidence>>();
                if del_evidence.len() == 0{
                    panic!("No deleted evidence found for gene {}", gene_name);
                }
                else if del_evidence.len() > 1{
                    panic!("Multiple deleted evidence found for gene {}", gene_name);
                }
                let del_evidence = del_evidence[0];
                let mut fixed_del_evidence = del_evidence.clone();
                let del_start_genome_idx = del_evidence.genome_index;
                let bases_to_trim = (genome_positions[0].genome_idx - del_start_genome_idx) as usize;
                let new_deleted_bases = del_evidence.alt[bases_to_trim..del_evidence.alt.len()].to_string();

                fixed_del_evidence.alt = new_deleted_bases.clone();
                fixed_del_evidence.genome_index = genome_positions[0].genome_idx;
                first_pos.alts.push(Alt{
                    alt_type: AltType::DEL,
                    base: new_deleted_bases.clone(),
                    evidence: fixed_del_evidence
                });
            }
        }
        if genome_positions[0].is_deleted_minor{
            if !genome_positions[0].alts.iter().any(|x| x.alt_type == AltType::DEL && x.evidence.is_minor){
                // None of the alts at this position are deletions, so didn't start here
                // Lets look at the deleted evidence to piece together what this should be
                for del_evidence in genome_positions[0].deleted_evidence.iter().filter(|x| x.is_minor){
                    let mut fixed_del_evidence = del_evidence.clone();
                    let del_start_genome_idx = del_evidence.genome_index;
                    let bases_to_trim = (genome_positions[0].genome_idx - del_start_genome_idx) as usize;
                    let new_deleted_bases = del_evidence.alt[bases_to_trim..del_evidence.alt.len()].to_string();

                    fixed_del_evidence.alt = new_deleted_bases.clone();
                    fixed_del_evidence.genome_index = genome_positions[0].genome_idx;
                    first_pos.alts.push(Alt{
                        alt_type: AltType::DEL,
                        base: new_deleted_bases.clone(),
                        evidence: fixed_del_evidence
                    });
                }
            }
        }
        genome_positions[0] = first_pos;

        let last_pos = genome_positions[genome_positions.len()-1].genome_idx;

        // Iter all positions, double checking that the deletions don't pass end of gene
        // truncating those which do
        for position in genome_positions.iter_mut(){
            for alt in position.alts.iter_mut(){
                if alt.alt_type == AltType::DEL{
                    if position.genome_idx + (alt.base.len() as i64) > last_pos{
                        let bases_to_trim = (position.genome_idx + alt.base.len() as i64 - last_pos) as usize;
                        alt.base = alt.base[0..alt.base.len()-bases_to_trim].to_string();
                    }
                }
            }
        }


    }

}

/// Converts a codon to an amino acid
/// 
/// # Arguments
/// - `codon`: Codon to convert
/// 
/// # Returns
/// - Amino acid the codon codes for
pub fn codon_to_aa(codon: String) -> char{
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

/// Complements a base
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
