use ordered_float::{Float, OrderedFloat};

use crate::common::{Alt, AltType, Evidence, MinorType, VCFRow};
use crate::gene::{self, codon_to_aa, Gene, GenePos, GenePosition};
use crate::genome::Genome;

#[derive(Clone, Debug)]
pub struct Variant{
    pub variant: String,
    pub nucleotide_index: i64,
    pub evidence: VCFRow,
    pub vcf_idx: i64,
    pub indel_length: i64,
    pub indel_nucleotides: Option<String>,
    pub gene_name: Option<String>,
    pub gene_position: Option<i64>,
    pub codon_idx: Option<i64>
}

pub struct Mutation{
    pub mutation: String,
    pub evidence: Vec<Evidence>,
    pub gene_position: i64,
    pub codes_protein: bool,

    pub ref_nucleotides: Option<String>,
    pub alt_nucleotides: Option<String>,
    pub nucleotide_number: Option<i64>,
    pub nucleotide_index: Option<i64>,
    pub indel_length: Option<i64>,
    pub indel_nucleotides: Option<String>,
    pub amino_acid_number: Option<i64>,
    pub amino_acid_sequence: Option<char>,


}

pub struct GenomeDifference{
    pub variants: Vec<Variant>,
    pub minor_variants: Vec<Variant>
}

pub struct GeneDifference{
    pub mutations: Vec<Mutation>,
    pub minor_mutations: Vec<Mutation>
}

impl GenomeDifference{
    pub fn new(ref_genome: Genome, mut alt_genome: Genome, minor_type: MinorType) -> Self{
        let mut variants: Vec<Variant> = Vec::new();
        let mut minor_variants: Vec<Variant> = Vec::new();

        for idx in 0..ref_genome.genome_positions.len(){
            let ref_pos = &ref_genome.genome_positions[idx];
            let alt_pos = alt_genome.genome_positions[idx].clone();

            if alt_pos.alts.len() > 0{
                // Alt has a variant at this position, so figure out what it is
                for alt in alt_pos.alts.iter(){
                    let mut garc = "".to_string();
                    let mut indel_bases = None;
                    let mut indel_length = 0;
                    let mut gene_name = None;
                    let mut gene_position = None;
                    let mut codon_idx = None;

                    // Only annotate with gene information if there is a single gene at this position
                    if alt_pos.genes.len() == 1{
                        gene_name = Some(alt_pos.genes[0].clone());
                        let gene = alt_genome.get_gene(alt_pos.genes[0].clone());
                        let (_gene_position, codon_idx) = gene.genome_idx_map.get(&alt_pos.genome_idx).unwrap();
                        gene_position = Some(*_gene_position);

                    }


                    if alt.alt_type == AltType::SNP || alt.alt_type == AltType::HET || alt.alt_type == AltType::NULL{
                        garc = ref_pos.genome_idx.to_string() + &ref_pos.reference.to_string() + ">" + &alt.base;
                    }
                    if alt.alt_type == AltType::INS{
                        garc = ref_pos.genome_idx.to_string() + "_ins_" + &alt.base;
                        indel_bases = Some(alt.base.clone());
                        indel_length = alt.base.len() as i64;
                    }
                    if alt.alt_type == AltType::DEL{
                        garc = ref_pos.genome_idx.to_string() + "_del_" + &alt.base;
                        indel_bases = Some(alt.base.clone());
                        indel_length = alt.base.len() as i64 * -1;
                    }

                    if alt.evidence.is_minor{
                        // Append coverage to the variant
                        if minor_type == MinorType::COV{
                            garc = garc + ":" + &alt.evidence.cov.unwrap().to_string();
                        }
                        if minor_type == MinorType::FRS{
                            garc = garc + ":" + &format!("{:.3}", alt.evidence.frs.unwrap());
                        }
                    }
                    let variant = Variant{
                        variant: garc,
                        nucleotide_index: ref_pos.genome_idx,
                        evidence: alt.evidence.vcf_row.clone(),
                        vcf_idx: alt.evidence.vcf_idx,
                        indel_length,
                        indel_nucleotides: indel_bases,
                        gene_position,
                        codon_idx,
                        gene_name
                    };

                    if alt.evidence.is_minor{
                        minor_variants.push(variant);
                    }else{
                        variants.push(variant);
                    }

                }
            }

        }

        GenomeDifference{
            variants,
            minor_variants
        }
    }
}

impl GeneDifference{
    pub fn new(ref_gene: Gene, alt_gene: Gene, minor_type: MinorType) -> Self{
        if ref_gene.name != alt_gene.name{
            panic!("Gene names do not match!");
        }
        let mut mutations = Vec::new();
        let mut minor_mutations = Vec::new();

        for idx in 0..ref_gene.gene_positions.len(){
            let ref_pos = &ref_gene.gene_positions[idx];
            let alt_pos = &alt_gene.gene_positions[idx];

            if ref_pos != alt_pos{
                // We have a difference so unpack it
                let gene_position = ref_pos.gene_position;
                let codes_protein = ref_gene.coding && ref_pos.gene_position > 0;
                match &alt_pos.gene_position_data{
                    GenePos::Codon(alt_codon) => {
                        // Unpack the ref codon too
                        let mut maybe_ref_codon = None;
                        match &ref_pos.gene_position_data{
                            GenePos::Codon(x) => {
                                maybe_ref_codon = Some(x);
                            },
                            _ => panic!("Reference gene position is not a codon")
                        }
                        let ref_codon = maybe_ref_codon.unwrap();

                        let mut ref_nucleotides = None;
                        let mut alt_nucleotides = None;
                        let mut nucleotide_number = None;
                        let mut nucleotide_index = None;
                        let mut indel_length = None;
                        let mut indel_nucleotides = None;
                        let mut amino_acid_number = None;
                        let mut amino_acid_sequence = None;
                        

                        let mut mutation = "".to_string();
                        let mut evidence: Vec<Evidence> = Vec::new();

                        // Codons are annoying as they can be split across multiple nucleotides
                        // So let's first check for a SNP as that's simple
                        let ref_codon_str = ref_codon.codon.iter().map(|x| x.reference.to_string()).collect::<String>();
                        ref_nucleotides = Some(
                            ref_codon_str.clone()
                        );
                        alt_nucleotides = Some(
                            alt_codon.codon.iter().map(|x| x.reference.to_string()).collect::<String>()
                        );
                        if alt_codon.amino_acid != ref_codon.amino_acid{
                            // SNP at this amino acid so mutation is simple
                            mutation = ref_gene.name.clone() + "@" + &ref_codon.amino_acid.to_string() + &ref_pos.gene_position.to_string() + &alt_codon.amino_acid.to_string();
                            amino_acid_number = Some(gene_position);
                            amino_acid_sequence = Some(alt_codon.amino_acid);
                            // Filter evidence only for SNP
                            for ev in alt_codon.codon.clone(){
                                if ev.alts.len() > 0{
                                    for e in ev.alts.iter(){
                                        if !e.evidence.is_minor && (e.alt_type == AltType::SNP || e.alt_type == AltType::HET || e.alt_type == AltType::NULL){
                                            evidence.push(e.evidence.clone());
                                        }
                                    }
                                }
                            }
                            mutations.push(Mutation{
                                mutation,
                                evidence,
                                gene_position,
                                codes_protein,
                                ref_nucleotides,
                                alt_nucleotides,
                                nucleotide_number,
                                nucleotide_index,
                                indel_length,
                                indel_nucleotides: indel_nucleotides.clone(),
                                amino_acid_number,
                                amino_acid_sequence
                            });
                        }
                        else{
                            // Pull out nucleotide variants in cases of synonymous mutations
                            let mut synon_snp = false;
                            for (ref_cd, alt_cd) in ref_codon.codon.iter().zip(alt_codon.codon.clone()){
                                if ref_cd.reference != alt_cd.reference{
                                    synon_snp = true;
                                    mutations.push(
                                        GeneDifference::nc_snp(
                                            &ref_gene.name,
                                            gene_position,
                                            codes_protein,
                                            ref_cd.reference,
                                            alt_cd.reference,
                                            ref_cd.nucleotide_number,
                                            ref_cd.nucleotide_index,
                                            alt_cd.alts.iter().filter(|x| !x.evidence.is_minor).collect::<Vec<&Alt>>()
                                        )
                                    )
                                }
                            }
                            if synon_snp{
                                mutation = ref_gene.name.clone() + "@" + &ref_codon.amino_acid.to_string() + &ref_pos.gene_position.to_string() + &alt_codon.amino_acid.to_string();
                                amino_acid_number = Some(gene_position);
                                amino_acid_sequence = Some(alt_codon.amino_acid);
                                // Filter evidence only for SNP
                                for ev in alt_codon.codon.clone(){
                                    if ev.alts.len() > 0{
                                        for e in ev.alts.iter(){
                                            if !e.evidence.is_minor && (e.alt_type == AltType::SNP || e.alt_type == AltType::HET || e.alt_type == AltType::NULL){
                                                evidence.push(e.evidence.clone());
                                            }
                                        }
                                    }
                                }
                                mutations.push(Mutation{
                                    mutation,
                                    evidence,
                                    gene_position,
                                    codes_protein,
                                    ref_nucleotides,
                                    alt_nucleotides,
                                    nucleotide_number,
                                    nucleotide_index,
                                    indel_length,
                                    indel_nucleotides: indel_nucleotides.clone(),
                                    amino_acid_number,
                                    amino_acid_sequence
                                });
                            }
                        }

                        // Then, check for indels and minor mutations in a codon

                        // Codons being fun mean we need to iter the codon and check for minor mutations
                        // An overall minor amino acid can be found, but the individual minor mutations need to be checked
                        // In cases of >1 minor SNP at a nucleotide, treat it as a het call with the minimum coverage
                        // Vec of (alt, cov, frs, evidence). 
                        // If no minor mutation at nucleotide, give ref nucleotide and None for other values
                        let mut minor_snps: Vec<(char, Option<i32>, Option<OrderedFloat<f32>>, Option<Vec<Evidence>>)> = Vec::new();
                        let mut minor_snp_exists = false;
                        for (ref_cd, alt_cd) in ref_codon.codon.iter().zip(alt_codon.codon.clone()){
                            let mut these_minor_snps = Vec::new();
                            let mut these_minor_indels = Vec::new();

                            // Make sure these are empty for each nc in the codon
                            ref_nucleotides = None;
                            alt_nucleotides = None;
                            nucleotide_number = None;
                            nucleotide_index = None;
                            indel_length = None;
                            indel_nucleotides = None;
                            amino_acid_number = None;
                            amino_acid_sequence = None;
                            mutation = "".to_string();
                            evidence = Vec::new();
                            for e in alt_cd.alts.iter(){
                                if e.evidence.is_minor{
                                    if e.alt_type == AltType::SNP || e.alt_type == AltType::HET || e.alt_type == AltType::NULL{
                                        these_minor_snps.push(e);
                                        minor_snp_exists = true;
                                    }
                                    else{
                                        these_minor_indels.push(e);
                                    }
                                }
                                else {
                                    if e.alt_type == AltType::INS{
                                        mutation = ref_gene.name.clone() + "@" + &alt_cd.nucleotide_number.to_string() + "_ins_" + &e.base;
                                        indel_length = Some(e.base.len() as i64);
                                        indel_nucleotides = Some(e.base.clone());
                                        evidence = vec![e.evidence.clone()];
                                    }
                                    if e.alt_type == AltType::DEL{
                                        mutation = ref_gene.name.clone() + "@" + &alt_cd.nucleotide_number.to_string() + "_del_" + &e.base;
                                        indel_length = Some(e.base.len() as i64 * -1);
                                        indel_nucleotides = Some(e.base.clone());
                                        evidence = vec![e.evidence.clone()];
                                    }

                                    if mutation != "".to_string(){
                                        // We picked up a mutation so lets append it
                                        mutations.push(Mutation{
                                            mutation: mutation.clone(),
                                            evidence: evidence.clone(),
                                            gene_position,
                                            codes_protein,
                                            ref_nucleotides: ref_nucleotides.clone(),
                                            alt_nucleotides: alt_nucleotides.clone(),
                                            nucleotide_number,
                                            nucleotide_index,
                                            indel_length,
                                            indel_nucleotides: indel_nucleotides.clone(),
                                            amino_acid_number,
                                            amino_acid_sequence
                                        });
                                    }
                                }
                            }
                            if these_minor_indels.len() > 0 && these_minor_snps.len() > 0{
                                // Mix of indel and SNP at this position
                                let mut these_minors = these_minor_indels.clone();
                                for snp in these_minor_snps.iter(){
                                    these_minors.push(snp);
                                }
                                minor_mutations.push(
                                    GeneDifference::mixed_indel(
                                        &ref_gene.name,
                                        gene_position,
                                        codes_protein,
                                        alt_cd.nucleotide_number,
                                        alt_cd.nucleotide_index,
                                        these_minors,
                                        minor_type.clone(),
                                        "mixed".to_string()
                                    )
                                );
                            }
                            else {
                                if these_minor_indels.len() > 0{
                                    // Minor indel
                                    if these_minor_indels.len() > 1{
                                        // We have a mixed minor indel, so treat it as such
                                        minor_mutations.push(
                                            GeneDifference::mixed_indel(
                                                &ref_gene.name,
                                                gene_position,
                                                codes_protein,
                                                alt_cd.nucleotide_number,
                                                alt_cd.nucleotide_index,
                                                these_minor_indels,
                                                minor_type.clone(),
                                                "indel".to_string()
                                            )
                                        );
                                    }
                                    else{
                                        let e = these_minor_indels[0];
                                        // We have a single minor indel which is much easier
                                        if e.alt_type == AltType::INS{
                                            mutation = ref_gene.name.clone() + "@" + &alt_cd.reference.to_string() + "_" + &alt_cd.nucleotide_number.to_string() + "_ins_" + &e.base;
                                            indel_length = Some(e.base.len() as i64);
                                            indel_nucleotides = Some(e.base.clone());
                                            evidence = vec![e.evidence.clone()];
                                        }
                                        if e.alt_type == AltType::DEL{
                                            mutation = ref_gene.name.clone() + "@" + &alt_cd.nucleotide_number.to_string() + "_del_" + &e.base;
                                            indel_length = Some(e.base.len() as i64 * -1);
                                            indel_nucleotides = Some(e.base.clone());
                                            evidence = vec![e.evidence.clone()];
                                        }
                                        if minor_type == MinorType::COV{
                                            mutation = mutation + ":" + &e.evidence.cov.unwrap().to_string();
                                        }
                                        if minor_type == MinorType::FRS{
                                            mutation = mutation + ":" + &format!("{:.3}", e.evidence.frs.unwrap());
                                        }
                                        minor_mutations.push(Mutation{
                                            mutation: mutation.clone(),
                                            evidence: evidence.clone(),
                                            gene_position,
                                            codes_protein,
                                            ref_nucleotides: None,
                                            alt_nucleotides: None,
                                            nucleotide_number,
                                            nucleotide_index,
                                            indel_length,
                                            indel_nucleotides: indel_nucleotides.clone(),
                                            amino_acid_number: None,
                                            amino_acid_sequence: None
                                        });
                                    }
                                }
                                if these_minor_snps.len() > 0{
                                    // Minor SNP
                                    if these_minor_snps.len() > 1{
                                        // Mixed minor SNP
                                        let min_cov = these_minor_snps.iter().filter(
                                                |x| x.evidence.cov.is_some()
                                            ).map(
                                                |x| x.evidence.cov.unwrap()
                                            ).min().unwrap();
                                        let min_frs = these_minor_snps.iter().filter(
                                                |x| x.evidence.frs.is_some()
                                            ).map(
                                                |x| x.evidence.frs.unwrap()
                                            ).min().unwrap();
                                        minor_snps.push(('z', Some(min_cov), Some(min_frs), Some(these_minor_snps.iter().map(|x| x.evidence.clone()).collect())));
                                    }
                                    else{
                                        // Single minor SNP
                                        minor_snps.push((these_minor_snps[0].base.chars().nth(0).unwrap(), these_minor_snps[0].evidence.cov, these_minor_snps[0].evidence.frs, Some(vec![these_minor_snps[0].evidence.clone()])));
                                    }
                                }
                                else{
                                    // No minor SNP at this position
                                    minor_snps.push((ref_cd.reference, None, None, None));
                                }
                            }
                        }
                        if minor_snp_exists && minor_snps.len() == 3{
                            // Construct the minor amino acid change
                            let mut codon = "".to_string();
                            let mut minor_cov = i32::MAX;
                            let mut minor_frs = OrderedFloat::max_value();
                            let mut minor_evidence = Vec::new();
                            for (nc, cov, frs, ev) in minor_snps.iter(){
                                codon += &nc.to_string();
                                // Evidence is the only field which is None in the case of no minor SNP at this nc
                                match ev{
                                    Some(e) => {
                                        for e in e.iter(){
                                            minor_evidence.push(e.clone());
                                        }
                                        if cov.unwrap() < minor_cov{
                                            minor_cov = cov.unwrap();
                                        }
                                        if frs.unwrap() < minor_frs{
                                            minor_frs = frs.unwrap();
                                        }
                                    },
                                    None => {}
                                }
                            }
                            let aa = codon_to_aa(codon.clone());
                            let mut mutation = "".to_string();
                            if minor_type == MinorType::COV{
                                mutation = ref_gene.name.clone() + "@" + &ref_codon.amino_acid.to_string() + &ref_pos.gene_position.to_string() + &aa.to_string() + ":" + &minor_cov.to_string();
                            }
                            if minor_type == MinorType::FRS{
                                mutation = ref_gene.name.clone() + "@" + &ref_codon.amino_acid.to_string() + &ref_pos.gene_position.to_string() + &aa.to_string() + ":" + &format!("{:.3}", minor_frs);
                            }
                            minor_mutations.push(Mutation{
                                mutation,
                                evidence: minor_evidence,
                                gene_position,
                                codes_protein,
                                ref_nucleotides: Some(ref_codon_str.clone()),
                                alt_nucleotides: Some(codon.clone()),
                                nucleotide_number: None,
                                nucleotide_index: None,
                                indel_length: None,
                                indel_nucleotides: None,
                                amino_acid_number: Some(gene_position),
                                amino_acid_sequence: Some(aa)
                            });
                        }
                        
                    },
                    GenePos::Nucleotide(alt_nc) => {
                        // Unpack the ref nc too
                        let mut maybe_ref_nc = None;
                        match &ref_pos.gene_position_data{
                            GenePos::Nucleotide(x) => {
                                maybe_ref_nc = Some(x);
                            },
                            _ => panic!("Reference gene position is not a nucleotide")
                        }
                        let ref_nc = maybe_ref_nc.unwrap();

                        // Nucleotide mutations are much more simple
                        for alt in alt_nc.alts.iter(){
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

                            if alt.alt_type == AltType::SNP || alt.alt_type == AltType::HET || alt.alt_type == AltType::NULL{
                                mutation = ref_gene.name.clone() + "@" + &ref_nc.reference.to_string() + &ref_nc.nucleotide_number.to_string() + &alt.base;
                                ref_nucleotides = Some(ref_nc.reference.to_string());
                                alt_nucleotides = Some(alt_nc.reference.to_string());
                            }
                            if alt.alt_type == AltType::INS{
                                mutation = ref_gene.name.clone() + "@" + &ref_nc.nucleotide_number.to_string() + "_ins_" + &alt.base;
                                indel_length = Some(alt.base.len() as i64);
                                indel_nucleotides = Some(alt.base.clone());
                            }
                            if alt.alt_type == AltType::DEL{
                                mutation = ref_gene.name.clone() + "@" + &ref_nc.nucleotide_number.to_string() + "_del_" + &alt.base;
                                indel_length = Some(alt.base.len() as i64 * -1);
                                indel_nucleotides = Some(alt.base.clone());
                            }

                            if alt.evidence.is_minor{
                                // Append coverage to the variant
                                if minor_type == MinorType::COV{
                                    mutation = mutation + ":" + &alt.evidence.cov.unwrap().to_string();
                                }
                                if minor_type == MinorType::FRS{
                                    mutation = mutation + ":" + &format!("{:.3}", alt.evidence.frs.unwrap());
                                }
                            }
                            let m = Mutation{
                                mutation,
                                evidence,
                                gene_position,
                                codes_protein,
                                ref_nucleotides,
                                alt_nucleotides,
                                nucleotide_number,
                                nucleotide_index,
                                indel_length,
                                indel_nucleotides,
                                amino_acid_number,
                                amino_acid_sequence
                            };

                            if alt.evidence.is_minor{
                                minor_mutations.push(m);
                            }
                            else{
                                mutations.push(m);
                            }
                        }
                    }
                }
            }
        }


        return GeneDifference{
            mutations,
            minor_mutations
        }
    }

    fn nc_snp(gene_name: &String, gene_position: i64, codes_protein: bool, ref_nc: char, alt_nc: char, nc_num: i64, nc_idx: i64, evidence: Vec<&Alt>) -> Mutation{
        let mutation = gene_name.clone() + "@" + &ref_nc.to_string() + &nc_num.to_string() + &alt_nc.to_string();
        let ref_nucleotides = Some(ref_nc.to_string());
        let alt_nucleotides = Some(alt_nc.to_string());
        let nucleotide_number = Some(nc_num);
        let nucleotide_index = Some(nc_idx);
        let evidence = evidence.iter().filter(|x| x.alt_type == AltType::SNP).map(|x| x.evidence.clone()).collect();
        return Mutation{
            mutation,
            evidence,
            gene_position,
            codes_protein,
            ref_nucleotides,
            alt_nucleotides,
            nucleotide_number,
            nucleotide_index,
            indel_length: None,
            indel_nucleotides: None,
            amino_acid_number: None,
            amino_acid_sequence: None
        };
    }

    fn mixed_indel(gene_name: &String, gene_position: i64, codes_protein: bool, nc_num: i64, nc_idx: i64, these_minors: Vec<&Alt>, minor_type: MinorType, mutation_name: String) -> Mutation{
        let mut min_coverage = "".to_string();
        if minor_type == MinorType::COV{
            min_coverage = these_minors.iter().filter(
                |x| x.evidence.cov.is_some()
            ).map(
                |x| x.evidence.cov.unwrap()
            ).min().unwrap().to_string();
        }
        if minor_type == MinorType::FRS{
            min_coverage = format!("{:.3}",these_minors.iter().filter(
                |x| x.evidence.frs.is_some()
            ).map(
                |x| x.evidence.frs.unwrap()
            ).min().unwrap());
        }
        let mutation = gene_name.clone() + "@" + &nc_num.to_string() + "_" + &mutation_name + ":" + &min_coverage;
        let evidence = these_minors.iter().map(|x| x.evidence.clone()).collect();

        let ref_nucleotides = None;
        let alt_nucleotides = None;
        let nucleotide_number = Some(nc_num);
        let nucleotide_index = Some(nc_idx);
        return Mutation{
            mutation,
            evidence,
            gene_position,
            codes_protein,
            ref_nucleotides,
            alt_nucleotides,
            nucleotide_number,
            nucleotide_index,
            indel_length: None,
            indel_nucleotides: None,
            amino_acid_number: None,
            amino_acid_sequence: None
        };
    }
}
