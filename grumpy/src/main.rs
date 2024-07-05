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
struct Gene{
    name: String,
    coding: bool,
    rev_comp: bool,
    nucleotide_sequence: String,
    nucleotide_index: Vec<i64>,
    nucleotide_number: Vec<i64>,
    gene_position: Vec<i64>,
    amino_acid_sequence: String,
    ribosomal_shifts: Vec<i64>,
    codons: Vec<String>
}

#[derive(Clone, Debug)]
struct Evidence{

}

#[derive(Clone, Debug)]
struct Alt{
    // If alt_type is "SNP" then base is the new base
    // If alt_type is "ins" or "del" then base is the inserted or deleted bases
    alt_type: String,
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
                Ok(s) => s.to_lowercase(),
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
        let max_promoter_length = 100 - 1;
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
        if gene_def.rev_comp{
            let mut last_idx = gene_def.promoter_start;
            if gene_def.promoter_start == -1{
                last_idx = gene_def.start;
            }
            for i in gene_def.end..last_idx+1{
                nucleotide_sequence.push(self.genome_positions[i as usize].reference);
                nucleotide_index.push(self.genome_positions[i as usize].genome_idx);
            }
        }
        else{
            let mut first_idx = gene_def.promoter_start;
            if gene_def.promoter_start == -1{
                first_idx = gene_def.start;
            }
            for i in first_idx..gene_def.end+1{
                nucleotide_sequence.push(self.genome_positions[i as usize].reference);
                nucleotide_index.push(self.genome_positions[i as usize].genome_idx);
            }
        
        }

        return Gene::new(
            gene_def,
            nucleotide_sequence,
            nucleotide_index,
        );

    }

    pub fn at_genome_index(&self, index: i64) -> GenomePosition{
        // 1-indexed genome index
        return self.genome_positions[(index + 1) as usize].clone();
    }
}

impl Gene {
    pub fn new(gene_def: GeneDef, nc_sequence: String, nc_index: Vec<i64>) -> Self {
        let mut nucleotide_sequence = nc_sequence.clone();
        let mut nucleotide_index = nc_index.clone();
        let mut nucleotide_number = Vec::new();
        let mut amino_acid_sequence = "".to_string();
        let mut gene_position: Vec<i64> = Vec::new();
        let mut codons = Vec::new();

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

        if gene_def.rev_comp{
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
            for i in (-1*(gene_def.promoter_size + 1))..0{
                nucleotide_number.push(i);
                gene_position.push(i);
            }
        }
        let prom_end = nucleotide_number.len();
        // Now non-promoter
        for i in 1..(gene_def.start - gene_def.end).abs() + 1{
            nucleotide_number.push(i);
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
                    codon = "".to_string();
                    codon_idx += 1;
                    gene_position.push(codon_idx);
                }

            }
        }

        return Gene{
            name: gene_def.name.to_string(),
            coding: gene_def.coding,
            rev_comp: gene_def.rev_comp,
            nucleotide_sequence: nucleotide_sequence.to_string(),
            nucleotide_index,
            gene_position,
            nucleotide_number,
            amino_acid_sequence,
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
        "ttx" | "tcx" | "tax" | "tgx" | "txt" | "txc" | "txa" | "txg" | "txx" | "txz" | "txo" | "tzx" | "tox" | "ctx" | "ccx" | "cax" | "cgx" | "cxt" | "cxc" | "cxa" | "cxg" | "cxx" | "cxz" | "cxo" | "czx" | "cox" | "atx" | "acx" | "aax" | "agx" | "axt" | "axc" | "axa" | "axg" | "axx" | "axz" | "axo" | "azx" | "aox" | "gtx" | "gcx" | "gax" | "ggx" | "gxt" | "gxc" | "gxa" | "gxg" | "gxx" | "gxz" | "gxo" | "gzx" | "gox" | "xtt" | "xtc" | "xta" | "xtg" | "xtx" | "xtz" | "xto" | "xct" | "xcc" | "xca" | "xcg" | "xcx" | "xcz" | "xco" | "xat" | "xac" | "xaa" | "xag" | "xax" | "xaz" | "xao" | "xgt" | "xgc" | "xga" | "xgg" | "xgx" | "xgz" | "xgo" | "xxt" | "xxc" | "xxa" | "xxg" | "xxx" | "xxz" | "xxo" | "xzt" | "xzc" | "xza" | "xzg" | "xzx" | "xzz" | "xzo" | "xot" | "xoc" | "xoa" | "xog" | "xox" | "xoz" | "xoo" | "ztx" | "zcx" | "zax" | "zgx" | "zxt" | "zxc" | "zxa" | "zxg" | "zxx" | "zxz" | "zxo" | "zzx" | "zox" | "otx" | "ocx" | "oax" | "ogx" | "oxt" | "oxc" | "oxa" | "oxg" | "oxx" | "oxz" | "oxo" | "ozx" | "oox" => 'X',
        "ttz" | "tcz" | "taz" | "tgz" | "tzt" | "tzc" | "tza" | "tzg" | "tzz" | "ctz" | "ccz" | "caz" | "cgz" | "czt" | "czc" | "cza" | "czg" | "czz" | "atz" | "acz" | "aaz" | "agz" | "azt" | "azc" | "aza" | "azg" | "azz" | "gtz" | "gcz" | "gaz" | "ggz" | "gzt" | "gzc" | "gza" | "gzg" | "gzz" | "ztt" | "ztc" | "zta" | "ztg" | "ztz" | "zct" | "zcc" | "zca" | "zcg" | "zcz" | "zat" | "zac" | "zaa" | "zag" | "zaz" | "zgt" | "zgc" | "zga" | "zgg" | "zgz" | "zzt" | "zzc" | "zza" | "zzg" | "zzz" => 'Z',
        "tto" | "tco" | "tao" | "tgo" | "tzo" | "tot" | "toc" | "toa" | "tog" | "toz" | "too" | "cto" | "cco" | "cao" | "cgo" | "czo" | "cot" | "coc" | "coa" | "cog" | "coz" | "coo" | "ato" | "aco" | "aao" | "ago" | "azo" | "aot" | "aoc" | "aoa" | "aog" | "aoz" | "aoo" | "gto" | "gco" | "gao" | "ggo" | "gzo" | "got" | "goc" | "goa" | "gog" | "goz" | "goo" | "zto" | "zco" | "zao" | "zgo" | "zzo" | "zot" | "zoc" | "zoa" | "zog" | "zoz" | "zoo" | "ott" | "otc" | "ota" | "otg" | "otz" | "oto" | "oct" | "occ" | "oca" | "ocg" | "ocz" | "oco" | "oat" | "oac" | "oaa" | "oag" | "oaz" | "oao" | "ogt" | "ogc" | "oga" | "ogg" | "ogz" | "ogo" | "ozt" | "ozc" | "oza" | "ozg" | "ozz" | "ozo" | "oot" | "ooc" | "ooa" | "oog" | "ooz" | "ooo" => 'O',
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

fn main() {
    let mut reference = Genome::new("reference/NC_000962.3.gbk");
    reference.assign_promoters();
    for gene in reference.gene_definitions.iter(){
        if gene.name == "rpoB" || gene.name == "katG" {
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
    let katG = reference.build_gene("katG".to_string());
    println!("{:?}", katG);
    println!("{:?}", katG.nucleotide_sequence.len());
    println!("{:?}", katG.nucleotide_index.len());
    println!("{:?}", katG.nucleotide_number.len());

    let rpoB = reference.build_gene("rpoB".to_string());
    println!("{:?}", rpoB);
    println!("{:?}", rpoB.nucleotide_sequence.len());
    println!("{:?}", rpoB.nucleotide_index.len());
    println!("{:?}", rpoB.nucleotide_number.len());


}