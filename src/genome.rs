//! Module for handling genome data
use pyo3::prelude::*;

use std::collections::{HashMap, HashSet};
use std::fs::File;

use gb_io::reader::SeqReader;
use gb_io::seq::Location::Complement;
use gb_io::seq::Location::Join;
use gb_io::seq::Location::Range;

use crate::common::{Alt, AltType, Evidence, GeneDef};
use crate::gene::Gene;
use crate::vcf::VCFFile;

#[pyclass]
#[derive(Clone, Debug)]
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
    /// HashMap of genes which have been built
    pub genes: HashMap<String, Gene>,

    #[pyo3(get, set)]
    /// Set of genes with mutations
    pub genes_with_mutations: HashSet<String>,
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
        let max_promoter_length = 100;
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
        return self.genes.get(&gene_name).unwrap().clone();
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
        self.genome_positions[(index + 1) as usize].clone()
    }

    fn __add__(&self, other: VCFFile) -> Genome {
        mutate(self, other)
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

                if call.call_type == AltType::SNP
                    || call.call_type == AltType::HET
                    || call.call_type == AltType::NULL
                {
                    // Update nucleotide for SNP/het/null
                    let base = call.clone().alt.chars().next().unwrap();
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
                        deleted_evidence.insert(idx + del_idx, call.clone());
                    }
                }
            }
        }

        if deleted_positions.contains(&idx) {
            position.is_deleted = true;
            for ev in deleted_evidence.get(&idx).iter() {
                position.deleted_evidence.push((*ev).clone());
            }
        }
    }

    // Reset the gene hashmap as the nucleotide sequence has changed
    new_genome.genes = HashMap::new();

    // I hate implicit returns, but appease clippy
    new_genome
}

#[cfg(test)]
mod tests {
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
        assert_eq!(rrs.not_promoter(&rrs.nucleotide_sequence.chars().collect::<Vec<char>>()).iter().collect::<String>(),"ttttgtttggagagtttgatcctggctcaggacgaacgctggcggcgtgcttaacacatgcaagtcgaacggaaaggtctcttcggagatactcgagtggcgaacgggtgagtaacacgtgggtgatctgccctgcacttcgggataagcctgggaaactgggtctaataccggataggaccacgggatgcatgtcttgtggtggaaagcgctttagcggtgtgggatgagcccgcggcctatcagcttgttggtggggtgacggcctaccaaggcgacgacgggtagccggcctgagagggtgtccggccacactgggactgagatacggcccagactcctacgggaggcagcagtggggaatattgcacaatgggcgcaagcctgatgcagcgacgccgcgtgggggatgacggccttcgggttgtaaacctctttcaccatcgacgaaggtccgggttctctcggattgacggtaggtggagaagaagcaccggccaactacgtgccagcagccgcggtaatacgtagggtgcgagcgttgtccggaattactgggcgtaaagagctcgtaggtggtttgtcgcgttgttcgtgaaatctcacggcttaactgtgagcgtgcgggcgatacgggcagactagagtactgcaggggagactggaattcctggtgtagcggtggaatgcgcagatatcaggaggaacaccggtggcgaaggcgggtctctgggcagtaactgacgctgaggagcgaaagcgtggggagcgaacaggattagataccctggtagtccacgccgtaaacggtgggtactaggtgtgggtttccttccttgggatccgtgccgtagctaacgcattaagtaccccgcctggggagtacggccgcaaggctaaaactcaaaggaattgacgggggcccgcacaagcggcggagcatgtggattaattcgatgcaacgcgaagaaccttacctgggtttgacatgcacaggacgcgtctagagataggcgttcccttgtggcctgtgtgcaggtggtgcatggctgtcgtcagctcgtgtcgtgagatgttgggttaagtcccgcaacgagcgcaacccttgtctcatgttgccagcacgtaatggtggggactcgtgagagactgccggggtcaactcggaggaaggtggggatgacgtcaagtcatcatgccccttatgtccagggcttcacacatgctacaatggccggtacaaagggctgcgatgccgcgaggttaagcgaatccttaaaagccggtctcagttcggatcggggtctgcaactcgaccccgtgaagtcggagtcgctagtaatcgcagatcagcaacgctgcggtgaatacgttcccgggccttgtacacaccgcccgtcacgtcatgaaagtcggtaacacccgaagccagtggcctaaccctcgggagggagctgtcgaaggtgggatcggcgattgggacgaagtcgtaacaaggtagccgtaccggaaggtgcggctggatcacctcctttct".to_string());

        let _dnaa = genome.get_gene("dnaA".to_string());
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
        assert_eq!(orf1ab.amino_acid_sequence, "MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVFIKRSDARTAPHGHVMVELVAELEGIQYGRSGETLGVLVPHVGEIPVAYRKVLLRKNGNKGAGGHSYGADLKSFDLGDELGTDPYEDFQENWNTKHSSGVTRELMRELNGGAYTRYVDNNFCGPDGYPLECIKDLLARAGKASCTLSEQLDFIDTKRGVYCCREHEHEIAWYTERSEKSYELQTPFEIKLAKKFDTFNGECPNFVFPLNSIIKTIQPRVEKKKLDGFMGRIRSVYPVASPNECNQMCLSTLMKCDHCGETSWQTGDFVKATCEFCGTENLTKEGATTCGYLPQNAVVKIYCPACHNSEVGPEHSLAEYHNESGLKTILRKGGRTIAFGGCVFSYVGCHNKCAYWVPRASANIGCNHTGVVGEGSEGLNDNLLEILQKEKVNINIVGDFKLNEEIAIILASFSASTSAFVETVKGLDYKAFKQIVESCGNFKVTKGKAKKGAWNIGEQKSILSPLYAFASEAARVVRSIFSRTLETAQNSVRVLQKAAITILDGISQYSLRLIDAMMFTSDLATNNLVVMAYITGGVVQLTSQWLTNIFGTVYEKLKPVLDWLEEKFKEGVEFLRDGWEIVKFISTCACEIVGGQIVTCAKEIKESVQTFFKLVNKFLALCADSIIIGGAKLKALNLGETFVTHSKGLYRKCVKSREETGLLMPLKAPKEIIFLEGETLPTEVLTEEVVLKTGDLQPLEQPTSEAVEAPLVGTPVCINGLMLLEIKDTEKYCALAPNMMVTNNTFTLKGGAPTKVTFGDDTVIEVQGYKSVNITFELDERIDKVLNEKCSAYTVELGTEVNEFACVVADAVIKTLQPVSELLTPLGIDLDEWSMATYYLFDESGEFKLASHMYCSFYPPDEDEEEGDCEEEEFEPSTQYEYGTEDDYQGKPLEFGATSAALQPEEEQEEDWLDDDSQQTVGQQDGSEDNQTTTIQTIVEVQPQLEMELTPVVQTIEVNSFSGYLKLTDNVYIKNADIVEEAKKVKPTVVVNAANVYLKHGGGVAGALNKATNNAMQVESDDYIATNGPLKVGGSCVLSGHNLAKHCLHVVGPNVNKGEDIQLLKSAYENFNQHEVLLAPLLSAGIFGADPIHSLRVCVDTVRTNVYLAVFDKNLYDKLVSSFLEMKSEKQVEQKIAEIPKEEVKPFITESKPSVEQRKQDDKKIKACVEEVTTTLEETKFLTENLLLYIDINGNLHPDSATLVSDIDITFLKKDAPYIVGDVVQEGVLTAVVIPTKKAGGTTEMLAKALRKVPTDNYITTYPGQGLNGYTVEEAKTVLKKCKSAFYILPSIISNEKQEILGTVSWNLREMLAHAEETRKLMPVCVETKAIVSTIQRKYKGIKIQEGVVDYGARFYFYTSKTTVASLINTLNDLNETLVTMPLGYVTHGLNLEEAARYMRSLKVPATVSVSSPDAVTAYNGYLTSSSKTPEEHFIETISLAGSYKDWSYSGQSTQLGIEFLKRGDKSVYYTSNPTTFHLDGEVITFDNLKTLLSLREVRTIKVFTTVDNINLHTQVVDMSMTYGQQFGPTYLDGADVTKIKPHNSHEGKTFYVLPNDDTLRVEAFEYYHTTDPSFLGRYMSALNHTKKWKYPQVNGLTSIKWADNNCYLATALLTLQQIELKFNPPALQDAYYRARAGEAANFCALILAYCNKTVGELGDVRETMSYLFQHANLDSCKRVLNVVCKTCGQQQTTLKGVEAVMYMGTLSYEQFKKGVQIPCTCGKQATKYLVQQESPFVMMSAPPAQYELKHGTFTCASEYTGNYQCGHYKHITSKETLYCIDGALLTKSSEYKGPITDVFYKENSYTTTIKPVTYKLDGVVCTEIDPKLDNYYKKDNSYFTEQPIDLVPNQPYPNASFDNFKFVCDNIKFADDLNQLTGYKKPASRELKVTFFPDLNGDVVAIDYKHYTPSFKKGAKLLHKPIVWHVNNATNKATYKPNTWCIRCLWSTKPVETSNSFDVLKSEDAQGMDNLACEDLKPVSEEVVENPTIQKDVLECNVKTTEVVGDIILKPANNSLKITEEVGHTDLMAAYVDNSSLTIKKPNELSRVLGLKTLATHGLAAVNSVPWDTIANYAKPFLNKVVSTTTNIVTRCLNRVCTNYMPYFFTLLLQLCTFTRSTNSRIKASMPTTIAKNTVKSVGKFCLEASFNYLKSPNFSKLINIIIWFLLLSVCLGSLIYSTAALGVLMSNLGMPSYCTGYREGYLNSTNVTIATYCTGSIPCSVCLSGLDSLDTYPSLETIQITISSFKWDLTAFGLVAEWFLAYILFTRFFYVLGLAAIMQLFFSYFAVHFISNSWLMWLIINLVQMAPISAMVRMYIFFASFYYVWKSYVHVVDGCNSSTCMMCYKRNRATRVECTTIVNGVRRSFYVYANGGKGFCKLHNWNCVNCDTFCAGSTFISDEVARDLSLQFKRPINPTDQSSYIVDSVTVKNGSIHLYFDKAGQKTYERHSLSHFVNLDNLRANNTKGSLPINVIVFDGKSKCEESSAKSASVYYSQLMCQPILLLDQALVSDVGDSAEVAVKMFDAYVNTFSSTFNVPMEKLKTLVATAEAELAKNVSLDNVLSTFISAARQGFVDSDVETKDVVECLKLSHQSDIEVTGDSCNNYMLTYNKVENMTPRDLGACIDCSARHINAQVAKSHNIALIWNVKDFMSLSEQLRKQIRSAAKKNNLPFKLTCATTRQVVNVVTTKIALKGGKIVNNWLKQLIKVTLVFLFVAAIFYLITPVHVMSKHTDFSSEIIGYKAIDGGVTRDIASTDTCFANKHADFDTWFSQRGGSYTNDKACPLIAAVITREVGFVVPGLPGTILRTTNGDFLHFLPRVFSAVGNICYTPSKLIEYTDFATSACVLAAECTIFKDASGKPVPYCYDTNVLEGSVAYESLRPDTRYVLMDGSIIQFPNTYLEGSVRVVTTFDSEYCRHGTCERSEAGVCVSTSGRWVLNNDYYRSLPGVFCGVDAVNLLTNMFTPLIQPIGALDISASIVAGGIVAIVVTCLAYYFMRFRRAFGEYSHVVAFNTLLFLMSFTVLCLTPVYSFLPGVYSVIYLYLTFYLTNDVSFLAHIQWMVMFTPLVPFWITIAYIICISTKHFYWFFSNYLKRRVVFNGVSFSTFEEAALCTFLLNKEMYLKLRSDVLLPLTQYNRYLALYNKYKYFSGAMDTTSYREAACCHLAKALNDFSNSGSDVLYQPPQTSITSAVLQSGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDVVYCPRHVICTSEDMLNPNYEDLLIRKSNHNFLVQAGNVQLRVIGHSMQNCVLKLKVDTANPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNFTIKGSFLNGSCGSVGFNIDYDCVSFCYMHHMELPTGVHAGTDLEGNFYGPFVDRQTAQAAGTDTTITVNVLAWLYAAVINGDRWFLNRFTTTLNDFNLVAMKYNYEPLTQDHVDILGPLSAQTGIAVLDMCASLKELLQNGMNGRTILGSALLEDEFTPFDVVRQCSGVTFQSAVKRTIKGTHHWLLLTILTSLLVLVQSTQWSLFFFLYENAFLPFAMGIIAMSAFAMMFVKHKHAFLCLFLLPSLATVAYFNMVYMPASWVMRIMTWLDMVDTSLSGFKLKDCVMYASAVVLLILMTARTVYDDGARRVWTLMNVLTLVYKVYYGNALDQAISMWALIISVTSNYSGVVTTVMFLARGIVFMCVEYCPIFFITGNTLQCIMLVYCFLGYFCTCYFGLFCLLNRYFRLTLGVYDYLVSTQEFRYMNSQGLLPPKNSIDAFKLNIKLLGVGGKPCIKVATVQSKMSDVKCTSVVLLSVLQQLRVESSSKLWAQCVQLHNDILLAKDTTEAFEKMVSLLSVLLSMQGAVDINKLCEEMLDNRATLQAIASEFSSLPSYAAFATAQEAYEQAVANGDSEVVLKKLKKSLNVAKSEFDRDAAMQRKLEKMADQAMTQMYKQARSEDKRAKVTSAMQTMLFTMLRKLDNDALNNIINNARDGCVPLNIIPLTTAAKLMVVIPDYNTYKNTCDGTTFTYASALWEIQQVVDADSKIVQLSEISMDNSPNLAWPLIVTALRANSAVKLQNNELSPVALRQMSCAAGTTQTACTDDNALAYYNTTKGGRFVLALLSDLQDLKWARFPKSDGTGTIYTELEPPCRFVTDTPKGPKVKYLYFIKGLNNLNRGMVLGSLAATVRLQAGNATEVPANSTVLSFCAFAVDAAKAYKDYLASGGQPITNCVKMLCTHTGTGQAITVTPEANMDQESFGGASCCLYCRCHIDHPNPKGFCDLKGKYVQIPTTCANDPVGFTLKNTVCTVCGMWKGYGCSCDQLREPMLQSADAQSFLNRVCGVSAARLTPCGTGTSTDVVYRAFDIYNDKVAGFAKFLKTNCCRFQEKDEDDNLIDSYFVVKRHTFSNYQHEETIYNLLKDCPAVAKHDFFKFRIDGDMVPHISRQRLTKYTMADLVYALRHFDEGNCDTLKEILVTYNCCDDDYFNKKDWYDFVENPDILRVYANLGERVRQALLKTVQFCDAMRNAGIVGVLTLDNQDLNGNWYDFGDFIQTTPGSGVPVVDSYYSLLMPILTLTRALTAESHVDTDLTKPYIKWDLLKYDFTEERLKLFDRYFKYWDQTYHPNCVNCLDDRCILHCANFNVLFSTVFPPTSFGPLVRKIFVDGVPFVVSTGYHFRELGVVHNQDVNLHSSRLSFKELLVYAADPAMHAASGNLLLDKRTTCFSVAALTNNVAFQTVKPGNFNKDFYDFAVSKGFFKEGSSVELKHFFFAQDGNAAISDYDYYRYNLPTMCDIRQLLFVVEVVDKYFDCYDGGCINANQVIVNNLDKSAGFPFNKWGKARLYYDSMSYEDQDALFAYTKRNVIPTITQMNLKYAISAKNRARTVAGVSICSTMTNRQFHQKLLKSIAATRGATVVIGTSKFYGGWHNMLKTVYSDVENPHLMGWDYPKCDRAMPNMLRIMASLVLARKHTTCCSLSHRFYRLANECAQVLSEMVMCGGSLYVKPGGTSSGDATTAYANSVFNICQAVTANVNALLSTDGNKIADKYVRNLQHRLYECLYRNRDVDTDFVNEFYAYLRKHFSMMILSDDAVVCFNSTYASQGLVASIKNFKSVLYYQNNVFMSEAKCWTETDLTKGPHEFCSQHTMLVKQGDDYVYLPYPDPSRILGAGCFVDDIVKTDGTLMIERFVSLAIDAYPLTKHPNQEYADVFHLYLQYIRKLHDELTGHMLDMYSVMLTNDNTSRYWEPEFYEAMYTPHTVLQAVGACVLCNSQTSLRCGACIRRPFLCCKCCYDHVISTSHKLVLSVNPYVCNAPGCDVTDVTQLYLGGMSYYCKSHKPPISFPLCANGQVFGLYKNTCVGSDNVTDFNAIATCDWTNAGDYILANTCTERLKLFAAETLKATEETFKLSYGIATVREVLSDRELHLSWEVGKPRPPLNRNYVFTGYRVTKNSKVQIGEYTFEKGDYGDAVVYRGTTTYKLNVGDYFVLTSHTVMPLSAPTLVPQEHYVRITGLYPTLNISDEFSSNVANYQKVGMQKYSTLQGPPGTGKSHFAIGLALYYPSARIVYTACSHAAVDALCEKALKYLPIDKCSRIIPARARVECFDKFKVNSTLEQYVFCTVNALPETTADIVVFDEISMATNYDLSVVNARLRAKHYVYIGDPAQLPAPRTLLTKGTLEPEYFNSVCRLMKTIGPDMFLGTCRRCPAEIVDTVSALVYDNKLKAHKDKSAQCFKMFYKGVITHDVSSAINRPQIGVVREFLTRNPAWRKAVFISPYNSQNAVASKILGLPTQTVDSSQGSEYDYVIFTQTTETAHSCNVNRFNVAITRAKVGILCIMSDRDLYDKLQFTSLEIPRRNVATLQAENVTGLFKDCSKVITGLHPTQAPTHLSVDTKFKTEGLCVDIPGIPKDMTYRRLISMMGFKMNYQVNGYPNMFITREEAIRHVRAWIGFDVEGCHATREAVGTNLPLQLGFSTGVNLVAVPTGYVDTPNNTDFSRVSAKPPPGDQFKHLIPLMYKGLPWNVVRIKIVQMLSDTLKNLSDRVVFVLWAHGFELTSMKYFVKIGPERTCCLCDRRATCFSTASDTYACWHHSIGFDYVYNPFMIDVQQWGFTGNLQSNHDLYCQVHGNAHVASCDAIMTRCLAVHECFVKRVDWTIEYPIIGDELKINAACRKVQHMVVKAALLADKFPVLHDIGNPKAIKCVPQADVEWKFYDAQPCSDKAYKIEELFYSYATHSDKFTDGVCLFWNCNVDRYPANSIVCRFDTRVLSNLNLPGCDGGSLYVNKHAFHTPAFDKSAFVNLKQLPFFYYSDSPCESHGKQVVSDIDYVPLKSATCITRCNLGGAVCRHHANEYRLYLDAYNMMISAGFSLWVYKQFDTYNLWNTFTRLQSLENVAFNVVNKGHFDGQQGEVPVSIINNTVYTKVDGVDVELFENKTTLPVNVAFELWAKRNIKPVPEVKILNNLGVDIAANTVIWDYKRDAPAHISTIGVCSMTDIAKKPTETICAPLTVFFDGRVDGQVDLFRNARNGVLITEGSVKGLQPSVGPKQASLNGVTLIGEAVKTQFNYYKKVDGVVQQLPETYFTQSRNLQEFKPRSQMEIDFLELAMDEFIERYKLEGYAFEHIVYGDFSHSQLGGLHLLIGLAKRFKESPFELEDFIPMDSTVKNYFITDAQTGSSKCVCSVIDLLLDDFVEIIKSQDLSVVSKVVKVTIDYTEISFMLWCKDGHVETFYPKLQSSQAWQPGVAMPNLYKMQRMLLEKCDLQNYGDSATLPKGIMMNVAKYTQLCQYLNTLTLAVPYNMRVIHFGAGSDKGVAPGTAVLRQWLPTGTLLVDSDLNDFVSDADSTLIGDCATVHTANKWDLIISDMYDPKTKNVTKENDSKEGFFTYICGFIQQKLALGGSVAIKITEHSWNADLYKLMGHFAWWTAFVTNVNASSSEAFLIGCNYLGKPREQIDGYVMHANYIFWRNTNPIQLSSYSLFDMSKFPLKLRGTAVMSLKEGQINDMILSLLSKGRLIIRENNRVVISSDVLVNN!".to_string());
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
}
