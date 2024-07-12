use crate::common::{AltType, MinorType, VCFRow};
use crate::genome::Genome;

#[derive(Clone, Debug)]
pub struct Variant{
    pub variant: String,
    pub nucleotide_indx: i64,
    pub evidence: VCFRow,
    pub vcf_idx: i64,
    pub indel_length: i64,
    pub indel_nucleotides: Option<String>,
    pub gene_name: Option<String>,
    pub gene_position: Option<i64>,
    pub codon_idx: Option<i64>
}

pub struct GenomeDifference{
    pub variants: Vec<Variant>,
    pub minor_variants: Vec<Variant>
}

impl GenomeDifference{
    pub fn new(ref_genome: Genome, alt_genome: Genome, minor_type: MinorType) -> Self{
        let mut variants: Vec<Variant> = Vec::new();
        let mut minor_variants: Vec<Variant> = Vec::new();

        for idx in 0..ref_genome.genome_positions.len(){
            let ref_pos = &ref_genome.genome_positions[idx];
            let alt_pos = &alt_genome.genome_positions[idx];

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
                        let gene = alt_genome.build_gene(alt_pos.genes[0].clone());
                        for g_pos in gene.gene_positions{
                            let mut c_idx = 0;
                            for nc_idx in g_pos.nucleotide_index.iter(){
                                if *nc_idx == alt_pos.genome_idx{
                                    gene_position = Some(g_pos.gene_position);
                                    codon_idx = Some(c_idx);
                                    break;
                                }
                                c_idx += 1;
                            }
                        }
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
                            garc = garc + ":" + &alt.evidence.frs.unwrap().to_string();
                        }
                    }
                    let variant = Variant{
                        variant: garc,
                        nucleotide_indx: ref_pos.genome_idx,
                        evidence: alt.evidence.vcf_row.clone(),
                        vcf_idx: idx as i64,
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
