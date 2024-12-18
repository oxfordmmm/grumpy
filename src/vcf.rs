//! Module for handling VCF files
use pyo3::prelude::*;
use rayon::prelude::*;

use std::collections::{hash_map, HashMap};
use std::fs::File;
use std::io::BufReader;
use std::string::String;
use vcf::{VCFReader, VCFRecord};

use crate::common::{AltType, Evidence, VCFRow};

#[pyclass]
/// Dummy struct for wrapping VCFRecord
///
/// Required to make a valid pyclass to use as a function argument
#[derive(Clone)]
pub struct VCFRecordToParse {
    pub record: VCFRecord,
    pub min_dp: i32,
    pub required_fields: Vec<String>,
    pub vcf_row_index: usize,
}

#[pyclass]
#[derive(Clone, Debug)]
/// Struct to hold the information from a VCF file
pub struct VCFFile {
    #[pyo3(get, set)]
    /// Header of the VCF file. TODO: populate
    pub header: Vec<String>,

    #[pyo3(get, set)]
    /// Records of the VCF file
    pub records: Vec<VCFRow>,

    #[pyo3(get, set)]
    /// Calls from the VCF file, indexed by genome index
    pub calls: HashMap<i64, Vec<Evidence>>,

    #[pyo3(get, set)]
    /// Minor calls from the VCF file, indexed by genome index
    pub minor_calls: HashMap<i64, Vec<Evidence>>,
}

#[pymethods]
impl VCFFile {
    #[new]
    /// Create a new VCFFile object
    ///
    /// # Arguments
    /// - `filename`: String - Path to the VCF file
    /// - `ignore_filter`: bool - Whether to ignore the filter column
    /// - `min_dp`: i32 - Minimum depth to consider a call
    pub fn new(filename: String, ignore_filter: bool, min_dp: i32) -> Self {
        let file = File::open(filename).unwrap();
        let buf = BufReader::new(file);
        let mut reader = VCFReader::new(buf).unwrap();

        let mut record = reader.empty_record();
        let mut more_records = reader.next_record(&mut record);
        let mut reader_records = Vec::new();

        let required_fields = vec!["GT".to_string()];
        // Read data into memory
        while matches!(more_records, Ok(true)) {
            reader_records.push(record.clone());
            more_records = reader.next_record(&mut record);
        }

        // Parse records multithreaded
        let parsed = reader_records
            .par_iter()
            .enumerate()
            .map(|(idx, record)| {
                VCFFile::parse_record(VCFRecordToParse {
                    record: record.clone(),
                    min_dp,
                    required_fields: required_fields.clone(),
                    vcf_row_index: idx,
                })
            })
            .collect::<Vec<(VCFRow, Vec<Evidence>, Vec<Evidence>)>>();
        // For ease of access, we'll store the calls in a hashmap indexed on genome index
        let mut calls_map: HashMap<i64, Vec<Evidence>> = HashMap::new();
        let mut minor_calls_map: HashMap<i64, Vec<Evidence>> = HashMap::new();
        let mut records = Vec::new();

        // Fetch data
        for (record, record_calls, _record_minor_calls) in parsed.iter() {
            let passed = record.is_filter_pass;
            let mut record_minor_calls = _record_minor_calls.clone();
            records.push(record.clone());
            for call in record_calls.iter() {
                let mut added = false;
                if call.call_type == AltType::NULL {
                    // Null calls need some nudging...
                    if passed || ignore_filter {
                        added = true;
                    } else if !passed {
                        // Not passed filter, so check if all filters are in allowed filters
                        let allowed_filters = [
                            "MIN_FRS".to_string(),
                            "MIN_DP".to_string(),
                            "MIN_GCP".to_string(),
                            "NO_DATA".to_string(),
                        ];
                        let mut all_allowed = true;
                        for f in record.filter.iter() {
                            if !allowed_filters.contains(f) {
                                all_allowed = false;
                                break;
                            }
                        }
                        if all_allowed {
                            added = true;
                        }
                    } else {
                        // Skip this call as it's a null call with extra filter fails
                        continue;
                    }
                } else if ignore_filter || passed {
                    added = true;
                }

                // Check if calls are already in the map,
                // if they are, warn the user about multiple calls at the same position and add to vector
                // if they aren't, add to the map
                if let hash_map::Entry::Vacant(e) = calls_map.entry(call.genome_index) {
                    if added {
                        e.insert(vec![call.clone()]);
                    }
                } else {
                    println!(
                        "Multiple calls at genome position {}! {:?}\n",
                        call.genome_index,
                        calls_map.get(&call.genome_index).unwrap()
                    );
                    if added {
                        calls_map
                            .get_mut(&call.genome_index)
                            .unwrap()
                            .push(call.clone());
                    }
                }

                if !added {
                    // Call skipped due to filter fail, so add as a minor call
                    let mut c = call.clone();
                    c.is_minor = true;
                    record_minor_calls.push(c);
                }
            }

            // Add minor calls if the filter is passed or ignored, or if fails lie within allowed filters
            let mut valid_filters = passed || ignore_filter;
            if !valid_filters {
                // Not passed filter and not ignored, so check if all filters are in allowed filters
                let allowed_filters = [
                    "MIN_FRS".to_string(),
                    "MIN_DP".to_string(),
                    "MIN_GCP".to_string(),
                ];
                // Auto-ignore if no filters are present (i.e filter = "." which is fail)
                let mut all_allowed = !record.filter.is_empty();
                for f in record.filter.iter() {
                    if !allowed_filters.contains(f) {
                        all_allowed = false;
                        break;
                    }
                }
                if all_allowed {
                    valid_filters = true;
                }
            }
            if valid_filters {
                for call in record_minor_calls.iter() {
                    if let hash_map::Entry::Vacant(e) = minor_calls_map.entry(call.genome_index) {
                        e.insert(vec![call.clone()]);
                    } else {
                        minor_calls_map
                            .get_mut(&call.genome_index)
                            .unwrap()
                            .push(call.clone());
                    }
                }
            }
        }

        // I hate implict returns, but appease clippy
        VCFFile {
            header: Vec::new(),
            records,
            calls: calls_map,
            minor_calls: minor_calls_map,
        }
    }

    #[staticmethod]
    fn parse_record(rec: VCFRecordToParse) -> (VCFRow, Vec<Evidence>, Vec<Evidence>) {
        let record = rec.record;
        let min_dp = rec.min_dp;
        let required_fields = rec.required_fields;
        let vcf_row_index = rec.vcf_row_index;

        // Parse fields of the record to pull out data
        // For whatever reason non of these fields are strings, so we need to convert them
        let mut alts = Vec::new(); // One string per alt
        for alt in record.alternative.iter() {
            alts.push(String::from_utf8_lossy(alt).to_string().to_lowercase());
        }
        let mut filters = Vec::new(); // String per filter
        for filter in record.filter.iter() {
            for f in String::from_utf8_lossy(filter).split(";") {
                filters.push(f.to_string());
            }
        }

        // Oddities of how this VCF library works...
        // Format is a vector of bytes, but we need to convert it to a vector of strings
        // Each of these strings corresponds to an item in the genotypes vector
        let mut format: Vec<String> = Vec::new();
        for f in record.format.iter() {
            format.push(String::from_utf8_lossy(f).to_string());
        }

        let mut idx = 0;
        let mut fields: HashMap<String, Vec<String>> = HashMap::new();
        for s in record.genotype.iter() {
            let mut item: Vec<Vec<String>> = Vec::new();
            for i in s.iter() {
                let mut value: Vec<String> = Vec::new();
                for j in i.iter() {
                    value.push(String::from_utf8_lossy(j).to_string());
                }
                item.push(value.clone());
                fields.insert(format[idx].clone(), value.clone());
                idx += 1;
            }
        }

        // Validate that this record has the required fields
        for field in required_fields.iter() {
            if !format.contains(&field.to_string()) {
                panic!("Required field {} not found in record", field);
            }
        }

        // Enforce that this record has a COV field
        if !fields.contains_key("COV") {
            let mut _cov_tag = "";
            let mut cov = Vec::new();
            if fields.contains_key("AD") {
                _cov_tag = "AD";
            } else if fields.contains_key("RO") && fields.contains_key("AO") {
                // Annoying edge case where no single COV field exists but can be constructed
                _cov_tag = "RO";
            } else {
                panic!("No COV tag found in record");
            }
            if _cov_tag != "RO" {
                for c in fields.get(_cov_tag).unwrap() {
                    cov.push(c.to_string());
                }
            } else {
                let ro = fields.get("RO").unwrap()[0].clone();
                let ao = fields.get("AO").unwrap()[0].clone();
                cov.push(ro.to_string());
                cov.push(ao.to_string());
            }

            fields.insert("COV".to_string(), cov);
        }

        // Check if this record has passed the filters
        let mut passed = false;
        if filters == vec!["PASS"] {
            passed = true;
        }

        let row = VCFRow {
            position: record.position as i64,
            reference: String::from_utf8_lossy(&record.reference)
                .to_string()
                .to_lowercase(),
            alternative: alts.clone(),
            filter: filters.clone(),
            fields: fields.clone(),
            is_filter_pass: passed,
        };

        let (record_calls, record_minor_calls) =
            VCFFile::parse_record_for_calls(row.clone(), min_dp, vcf_row_index);

        // if record.position == 178452 {
        //     println!("{:?}\t{:?}\t{:?}\t{:?}\t{:?}", record.position, String::from_utf8_lossy(&record.reference), alts, filters, fields);
        //     for call in record_calls.iter(){
        //         println!("{:?}\n", call);
        //     }
        //     println!("--");
        //     for call in record_minor_calls.iter(){
        //         println!("{:?}\n", call);
        //     }
        //     println!("\n\n");
        // }

        (row, record_calls, record_minor_calls)
    }

    #[staticmethod]
    /// Parse a record from a VCF file to get the calls
    ///
    /// # Arguments
    /// - `record`: VCFRow - Record to parse
    /// - `min_dp`: i32 - Minimum depth to consider a call
    ///
    /// # Returns
    /// Tuple of:
    /// - `calls`: Vec of Evidence - Calls from the record
    /// - `minor_calls`: Vec of Evidence - Minor calls from the record
    pub fn parse_record_for_calls(
        record: VCFRow,
        min_dp: i32,
        vcf_row_index: usize,
    ) -> (Vec<Evidence>, Vec<Evidence>) {
        let mut calls: Vec<Evidence> = Vec::new();
        let mut minor_calls: Vec<Evidence> = Vec::new();

        // Dealing with the actual call here; should spit out possible minor calls afterwards

        // Convert ugly genotype into a list of strings for each allele
        let genotype = &(record.fields.get("GT").unwrap())[0]
            .split("/")
            .collect::<Vec<&str>>();
        let mut cov = Vec::new();
        for c in record.fields.get("COV").unwrap() {
            cov.push(c.parse::<i32>().unwrap());
        }

        let mut dp = 0;
        // As DP isn't always reliable, default to the sum of the COV field
        cov.iter().for_each(|x| dp += x);

        let ref_allele = record.reference.clone().to_lowercase();
        let mut alt_allele = "".to_string();

        let mut call_type = AltType::NULL;

        // println!(
        //     "{:?}\t{:?}\t{:?}\t{:?}\t{:?}",
        //     record.position, record.reference, record.alternative, record.filter, record.fields
        // );
        let first = genotype[0];
        // Adjust for 1-indexed VCF
        let mut alt_idx = 0;
        if first != "." {
            alt_idx = first.parse::<i32>().unwrap() - 1;
        }
        if cov.len() == 1 {
            // Just 1 item in the call so santiy check if it's a null call
            let mut vcf_idx = None;
            if genotype == &vec!["0", "0"] && cov[0] >= min_dp {
                // Ref call
                call_type = AltType::REF;
                alt_allele = ref_allele.clone();
                vcf_idx = Some(0);
            } else if genotype[0] != genotype[1] && cov[0] >= min_dp {
                // Het call
                call_type = AltType::HET;
                alt_allele = "z".to_string();
            } else {
                // Null call
                call_type = AltType::NULL;
                alt_allele = "x".to_string();
            }
            calls.push(Evidence {
                cov: Some(cov[0_usize]),
                frs: Some(ordered_float::OrderedFloat(1.0)),
                genotype: genotype.join("/"),
                call_type,
                vcf_row: vcf_row_index,
                reference: ref_allele.chars().nth(0).unwrap().to_string(),
                alt: alt_allele.clone(),
                genome_index: record.position,
                is_minor: false,
                vcf_idx,
            });
            return (calls, minor_calls);
        }

        for gt in genotype.iter().skip(1) {
            if *gt != first {
                call_type = AltType::HET;
                alt_allele = "z".to_string();
                break;
            }
        }
        if call_type != AltType::HET {
            if first == "0" {
                call_type = AltType::REF;
                alt_allele = ref_allele.clone();
            } else {
                // Placeholder to denote we have an actual call at this point
                // Could actually be an indel instead but that can't be inferred from just the genotype
                call_type = AltType::SNP;
            }
        }
        if *genotype == vec![".", "."] {
            call_type = AltType::NULL;
        }
        if call_type == AltType::NULL {
            alt_allele = "x".to_string();
        }

        if call_type == AltType::SNP {
            // Parse the corresponding alternate allele to get the call(s)
            alt_allele = record.alternative[alt_idx as usize].clone().to_lowercase();
            let call = VCFFile::simplify_call(ref_allele.clone(), alt_allele.clone());
            let call_cov = cov[(alt_idx + 1) as usize]; // COV should be [ref, alt1, alt2..]
            for (offset, __alt_type, __base) in call {
                let mut _alt_type = __alt_type;
                let mut _base = __base;
                if call_cov < min_dp {
                    // Override calls with null if the coverage is too low
                    _alt_type = AltType::NULL;
                    _base = "x".to_string();
                }
                let mut reference = "".to_string();
                if offset >= 0 {
                    reference = ref_allele.chars().nth(offset as usize).unwrap().to_string();
                }
                calls.push(Evidence {
                    cov: Some(call_cov),
                    frs: Some(ordered_float::OrderedFloat(call_cov as f32 / dp as f32)),
                    genotype: genotype.join("/"),
                    call_type: _alt_type,
                    vcf_row: vcf_row_index,
                    reference,
                    alt: _base,
                    genome_index: record.position + offset,
                    is_minor: false,
                    vcf_idx: Some((alt_idx + 1) as i64),
                });
            }
        } else {
            let mut call_cov = None;
            let mut frs = None;
            let mut vcf_idx = None;
            if call_type != AltType::HET {
                if genotype == &vec![".", "."] {
                    // Explicit null GT can't point to a coverage so null
                    call_type = AltType::NULL;
                    alt_allele = "x".to_string();
                } else {
                    call_cov = Some(cov[(alt_idx + 1) as usize]); // COV should be [ref, alt1, alt2..]
                    if call_cov.unwrap() < min_dp {
                        // Override calls with null if the coverage is too low
                        call_type = AltType::NULL;
                        alt_allele = "x".to_string();
                    }
                    frs = Some(ordered_float::OrderedFloat(
                        call_cov.unwrap() as f32 / dp as f32,
                    ));
                    vcf_idx = Some((alt_idx + 1) as i64);
                }
            }
            for (offset, r) in ref_allele.chars().enumerate() {
                calls.push(Evidence {
                    cov: call_cov,
                    frs,
                    genotype: genotype.join("/"),
                    call_type: call_type.clone(),
                    vcf_row: vcf_row_index,
                    reference: r.to_string(),
                    alt: alt_allele.clone(),
                    genome_index: record.position + offset as i64,
                    is_minor: false,
                    vcf_idx,
                });
            }
        }

        if call_type == AltType::NULL {
            // Don't check for minor calls if the call is null
            return (calls, minor_calls);
        }

        if call_type == AltType::HET {
            // HET calls need to make sure we aren't skipping the first alt
            alt_idx = -1;
        }

        // Now we've got the main call, we need to figure out the minor calls
        // So check all possible values of COV for threshold to be considered a minor call
        let mut idx = 0;
        for coverage in cov.iter() {
            if idx == (alt_idx + 1) as usize || idx == 0 {
                // Skip ref and the alt we've already called
                idx += 1;
                continue;
            }
            if *coverage >= min_dp {
                alt_allele = record.alternative[idx - 1].clone().to_lowercase();
                let call = VCFFile::simplify_call(ref_allele.clone(), alt_allele.clone());
                let call_cov = *coverage;
                for (offset, alt_type, base) in call {
                    let mut reference = "".to_string();
                    if offset >= 0 {
                        reference = ref_allele.chars().nth(offset as usize).unwrap().to_string();
                    }
                    minor_calls.push(Evidence {
                        cov: Some(call_cov),
                        frs: Some(ordered_float::OrderedFloat(call_cov as f32 / dp as f32)),
                        genotype: genotype.join("/"), // This is a minor call so the row's genotype is the same
                        call_type: alt_type,
                        vcf_row: vcf_row_index,
                        reference,
                        alt: base,
                        genome_index: record.position + offset,
                        is_minor: true,
                        vcf_idx: Some(idx as i64),
                    });
                }
            }
            idx += 1;
        }

        // I hate implict returns, but appease clippy
        (calls, minor_calls)
    }

    #[staticmethod]
    /// Simplify a call into a list of SNPs, INSs and DELs
    ///
    /// Some rows have long ref/alt pairs which need to be decomposed into constituent calls
    ///
    /// # Arguments
    /// - `reference`: String - Reference sequence
    /// - `alternate`: String - Alternate sequence
    ///
    /// # Returns
    /// `Vec of (usize, AltType, String)` - Vec where each element is a 3-tuple of:
    ///    - `usize`: Offset of the call from row's genome index
    ///    - `AltType`: Type of call
    ///    - `String`: Base(s) of the call. If SNP/HET/NULL this will be a single base. If INS/DEL this will be the sequence inserted/deleted
    pub fn simplify_call(reference: String, alternate: String) -> Vec<(i64, AltType, String)> {
        // Note that the offset of a call can be negative in annoying edge cases
        // e.g insertion of "gg" --> "cagg" is a -1 offset as "ca" is inserted before the reference
        let mut calls: Vec<(i64, AltType, String)> = Vec::new();
        if reference.len() == alternate.len() {
            // Simple set of SNPs
            for i in 0..reference.len() {
                if reference.chars().nth(i).unwrap() != alternate.chars().nth(i).unwrap() {
                    calls.push((
                        i as i64,
                        AltType::SNP,
                        alternate.chars().nth(i).unwrap().to_string(),
                    ));
                }
            }
            return calls;
        }

        /*
        The process for finding the positions for indels are almost identical
        as the process for finding a del can be interpreted as finding an ins with ref
            and alt reversed.
        The approach here is to use a sliding window to find the position of the indel
            where the number of SNPs is lowest.
        Assumes that there is only a single indel between the ref and the alt - there
            may be cases where this does not work these will just produce large amounts
            of SNPs... Could be adapted to check for multi-indel but this will scale
            exponentially the number of versions checked.
         */
        let mut _x: String = "".to_string();
        let mut _y: String = "".to_string();
        let length = (reference.len() as i64 - alternate.len() as i64).abs();
        let mut _indel_type = AltType::NULL;
        // Figure out which way around to approach this
        if reference.len() > alternate.len() {
            _y = reference.clone();
            _x = alternate.clone();
            _indel_type = AltType::DEL;
        } else {
            _y = alternate.clone();
            _x = reference.clone();
            _indel_type = AltType::INS;
        }

        let padding = "N".repeat(length as usize);
        let mut current = "".to_string();
        let mut current_dist = i64::MAX;
        let mut indel_start: i64 = 0;
        for i in 0.._x.len() + 1 {
            let x1 = _x[0..i].to_string() + &padding + &_x[i.._x.len()];
            // println!("{:?}\t{:?}\t{:?}\t{:?}\t{:?}", reference, alternate, x, y, x1);
            let dist = snp_dist(&x1, &_y);
            if dist <= current_dist {
                current = x1.clone();
                current_dist = dist;
                indel_start = i as i64;
            }
        }

        if _indel_type == AltType::INS {
            // Ins after, del at, so adjust
            indel_start -= 1;
        }

        let mut indel_str = "".to_string();
        for i in 0..current.len() {
            if current.chars().nth(i).unwrap() == 'N' {
                indel_str += &_y.chars().nth(i).unwrap().to_string();
            }
        }

        calls.push((indel_start, _indel_type.clone(), indel_str));

        if _indel_type == AltType::INS {
            // Return to vector index now the call has been constructed
            indel_start += 1;
        }
        let mut indel_left = length;
        let mut offset = 0;
        for i in 0.._y.len() {
            if i as i64 == indel_start {
                // Pad where the indel is to ensure we get all calls
                if indel_left > 0 {
                    indel_left -= 1;
                    indel_start += 1;
                    offset += 1;
                }
                continue;
            }
            let r = _y.chars().nth(i).unwrap();
            let a = _x.chars().nth(i - offset).unwrap();
            if r != 'N' && a != 'N' && r != a {
                if _indel_type == AltType::DEL {
                    if indel_left > 0 {
                        // We have some deletion ahead of this positon, so it doesn't need to be adjusted
                        calls.push(((i - offset) as i64, AltType::SNP, a.to_string()));
                    } else {
                        // We've passed the deletion, so adjust the offset to account for the deletion
                        calls.push((
                            (i - offset + length as usize) as i64,
                            AltType::SNP,
                            a.to_string(),
                        ));
                    }
                } else {
                    calls.push(((i - offset) as i64, AltType::SNP, r.to_string()));
                }
            }
        }

        // I hate implict returns, but appease clippy
        calls
    }
}

/// Calculate the distance between two strings ignoring N characters
///
/// # Arguments
/// - `reference`: &str - Reference sequence
/// - `alternate`: &str - Alternate sequence
///
/// # Returns
/// SNP distance between the two strings
fn snp_dist(reference: &str, alternate: &str) -> i64 {
    if reference.len() != alternate.len() {
        panic!("Reference and alternate strings are not the same length!");
    }
    let mut dist = 0;
    for i in 0..reference.len() {
        let r = reference.chars().nth(i).unwrap();
        let a = alternate.chars().nth(i).unwrap();
        if r != 'N' && a != 'N' && r != a {
            dist += 1;
        }
    }

    // I hate implict returns, but appease clippy
    dist
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
    fn test_snp_dist() {
        // Some trivial cases
        assert_eq!(snp_dist("ACGT", "ACGT"), 0);
        assert_eq!(snp_dist("AAAAA", "AAAAA"), 0);
        assert_eq!(snp_dist("ACGTACGT", "ACGTACGT"), 0);
        assert_eq!(snp_dist("ACGT", "AAGT"), 1);
        assert_eq!(snp_dist("ACGT", "AAAT"), 2);
        assert_eq!(snp_dist("ACGT", "AAAA"), 3);
        assert_eq!(snp_dist("ACGT", "NNNN"), 0);
        assert_eq!(snp_dist("ACGT", "NAGT"), 1);

        // Should panic
        assert_panics!(snp_dist("ACGT", "ACG"));
        assert_panics!(snp_dist("", "ACGT"));
        assert_panics!(snp_dist("ACGT", ""));
    }

    #[test]
    fn test_simplify_call() {
        // Some trivial cases
        assert_eq!(
            VCFFile::simplify_call("ACGT".to_string(), "ACGT".to_string()),
            vec![]
        );
        assert_eq!(
            VCFFile::simplify_call("AAAAA".to_string(), "AAAAA".to_string()),
            vec![]
        );
        assert_eq!(
            VCFFile::simplify_call("ACGTACGT".to_string(), "ACGTACGT".to_string()),
            vec![]
        );
        assert_eq!(
            VCFFile::simplify_call("ACGT".to_string(), "AAGT".to_string()),
            vec![(1, AltType::SNP, "A".to_string())]
        );
        assert_eq!(
            VCFFile::simplify_call("ACGT".to_string(), "AAAT".to_string()),
            vec![
                (1, AltType::SNP, "A".to_string()),
                (2, AltType::SNP, "A".to_string()),
            ]
        );
        assert_eq!(
            VCFFile::simplify_call("ACGT".to_string(), "AAAA".to_string()),
            vec![
                (1, AltType::SNP, "A".to_string()),
                (2, AltType::SNP, "A".to_string()),
                (3, AltType::SNP, "A".to_string()),
            ]
        );
        // Some more complex cases
        assert_eq!(
            VCFFile::simplify_call("ACGT".to_string(), "ACG".to_string()),
            vec![(3, AltType::DEL, "T".to_string())]
        );
        assert_eq!(
            VCFFile::simplify_call("ACG".to_string(), "ACGT".to_string()),
            vec![(2, AltType::INS, "T".to_string())]
        );
        assert_eq!(
            VCFFile::simplify_call("ACGT".to_string(), "AGT".to_string()),
            vec![(1, AltType::DEL, "C".to_string()),]
        );

        // Somewhat pathological cases
        assert_eq!(
            VCFFile::simplify_call("AAAAA".to_string(), "CCC".to_string()),
            vec![
                (3, AltType::DEL, "AA".to_string()),
                (0, AltType::SNP, "C".to_string()),
                (1, AltType::SNP, "C".to_string()),
                (2, AltType::SNP, "C".to_string()),
            ]
        );

        assert_eq!(
            VCFFile::simplify_call("CCC".to_string(), "AAAAA".to_string()),
            vec![
                (2, AltType::INS, "AA".to_string()),
                (0, AltType::SNP, "A".to_string()),
                (1, AltType::SNP, "A".to_string()),
                (2, AltType::SNP, "A".to_string()),
            ]
        );

        assert_eq!(
            VCFFile::simplify_call("ACGT".to_string(), "ACGAGT".to_string()),
            vec![(2, AltType::INS, "AG".to_string()),]
        );

        assert_eq!(
            VCFFile::simplify_call(
                "AGTGCGCCTCCCGCGAGCAGACACAGAATCGCACTGCGCCGGCCCGGCGCGTGCGATTCTGTGTCTGCTT"
                    .to_string(),
                "AGCTC".to_string()
            ),
            vec![
                (
                    2,
                    AltType::DEL,
                    "TGCGCCTCCCGCGAGCAGACACAGAATCGCACTGCGCCGGCCCGGCGCGTGCGATTCTGTGTCTG".to_string()
                ),
                (69, AltType::SNP, "C".to_string()),
            ]
        );
    }

    #[test]
    fn test_parse_record_for_calls() {
        // Note the use of blocks to allow collapsing in vscode.
        // Not necessary but makes it easier to read as this is a long test

        // There's an almost infinite number of edge cases here, so we'll just test a few

        // Ref call
        {
            let record = VCFRow {
                position: 1,
                reference: "A".to_string(),
                alternative: vec!["T".to_string()],
                filter: vec![],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["0/0".to_string()]),
                    ("COV".to_string(), vec!["1".to_string()]),
                ]),
                is_filter_pass: true,
            };
            let (calls, minor_calls) = VCFFile::parse_record_for_calls(record, 1, 0);
            assert_eq!(calls.len(), 1);
            assert_eq!(minor_calls.len(), 0);
            assert_eq!(calls[0].genome_index, 1);
            assert_eq!(calls[0].reference, "a".to_string());
            assert_eq!(calls[0].alt, "a".to_string());
            assert_eq!(calls[0].call_type, AltType::REF);
            assert_eq!(calls[0].cov, Some(1));
            assert_eq!(calls[0].frs, Some(ordered_float::OrderedFloat(1.0)));
            assert!(!calls[0].is_minor);
            assert_eq!(calls[0].vcf_idx, Some(0));
            assert_eq!(calls[0].vcf_row, 0);
        }

        // Alt call
        {
            let record = VCFRow {
                position: 1,
                reference: "A".to_string(),
                alternative: vec!["T".to_string()],
                filter: vec![],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["1/1".to_string()]),
                    ("COV".to_string(), vec!["1".to_string(), "2".to_string()]),
                ]),
                is_filter_pass: true,
            };
            let (calls, minor_calls) = VCFFile::parse_record_for_calls(record, 1, 0);
            assert_eq!(calls.len(), 1);
            assert_eq!(minor_calls.len(), 0);
            assert_eq!(calls[0].genome_index, 1);
            assert_eq!(calls[0].reference, "a".to_string());
            assert_eq!(calls[0].alt, "t".to_string());
            assert_eq!(calls[0].call_type, AltType::SNP);
            assert_eq!(calls[0].cov, Some(2));
            assert_eq!(calls[0].frs, Some(ordered_float::OrderedFloat(2.0 / 3.0)));
            assert!(!calls[0].is_minor);
            assert_eq!(calls[0].vcf_idx, Some(1));
            assert_eq!(calls[0].vcf_row, 0);
        }

        // Het call (including minor call)
        {
            let record = VCFRow {
                position: 1,
                reference: "A".to_string(),
                alternative: vec!["T".to_string()],
                filter: vec![],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["0/1".to_string()]),
                    ("COV".to_string(), vec!["1".to_string(), "3".to_string()]),
                ]),
                is_filter_pass: true,
            };
            let (calls, minor_calls) = VCFFile::parse_record_for_calls(record, 2, 0);
            assert_eq!(calls.len(), 1);
            assert_eq!(minor_calls.len(), 1);

            assert_eq!(calls[0].genome_index, 1);
            assert_eq!(calls[0].reference, "a".to_string());
            assert_eq!(calls[0].alt, "z".to_string());
            assert_eq!(calls[0].call_type, AltType::HET);
            assert_eq!(calls[0].cov, None);
            assert_eq!(calls[0].frs, None);
            assert!(!calls[0].is_minor);
            assert_eq!(calls[0].vcf_idx, None);
            assert_eq!(calls[0].vcf_row, 0);

            assert_eq!(minor_calls[0].genome_index, 1);
            assert_eq!(minor_calls[0].reference, "a".to_string());
            assert_eq!(minor_calls[0].alt, "t".to_string());
            assert_eq!(minor_calls[0].call_type, AltType::SNP);
            assert_eq!(minor_calls[0].cov, Some(3));
            assert_eq!(minor_calls[0].frs, Some(ordered_float::OrderedFloat(0.75)));
            assert!(minor_calls[0].is_minor);
            assert_eq!(minor_calls[0].vcf_idx, Some(1));
            assert_eq!(minor_calls[0].vcf_row, 0);
        }

        // Null call (in a few forms)
        // Null call due to low coverage
        {
            let record = VCFRow {
                position: 1,
                reference: "A".to_string(),
                alternative: vec!["T".to_string()],
                filter: vec![],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["1/1".to_string()]),
                    ("COV".to_string(), vec!["1".to_string()]),
                ]),
                is_filter_pass: true,
            };
            let (calls, minor_calls) = VCFFile::parse_record_for_calls(record, 2, 0);
            assert_eq!(calls.len(), 1);
            assert_eq!(minor_calls.len(), 0);
            assert_eq!(calls[0].genome_index, 1);
            assert_eq!(calls[0].reference, "a".to_string());
            assert_eq!(calls[0].alt, "x".to_string());
            assert_eq!(calls[0].call_type, AltType::NULL);
            assert_eq!(calls[0].cov, Some(1));
            assert_eq!(calls[0].frs, Some(ordered_float::OrderedFloat(1.0)));
            assert!(!calls[0].is_minor);
            assert_eq!(calls[0].vcf_idx, None);
            assert_eq!(calls[0].vcf_row, 0);
        }

        // Null call due calling null
        {
            let record = VCFRow {
                position: 1,
                reference: "A".to_string(),
                alternative: vec!["T".to_string()],
                filter: vec![],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["./.".to_string()]),
                    ("COV".to_string(), vec!["1".to_string()]),
                ]),
                is_filter_pass: true,
            };
            let (calls, minor_calls) = VCFFile::parse_record_for_calls(record, 1, 0);
            assert_eq!(calls.len(), 1);
            assert_eq!(minor_calls.len(), 0);
            assert_eq!(calls[0].genome_index, 1);
            assert_eq!(calls[0].reference, "a".to_string());
            assert_eq!(calls[0].alt, "x".to_string());
            assert_eq!(calls[0].call_type, AltType::NULL);
            assert_eq!(calls[0].cov, Some(1));
            assert_eq!(calls[0].frs, Some(ordered_float::OrderedFloat(1.0)));
            assert!(!calls[0].is_minor);
            assert_eq!(calls[0].vcf_idx, None);
            assert_eq!(calls[0].vcf_row, 0);
        }

        // Null call due to low coverage and filter fail (should still make it)
        {
            let record = VCFRow {
                position: 1,
                reference: "A".to_string(),
                alternative: vec!["T".to_string()],
                filter: vec!["MIN_DP".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["1/1".to_string()]),
                    ("COV".to_string(), vec!["1".to_string()]),
                ]),
                is_filter_pass: false,
            };
            let (calls, minor_calls) = VCFFile::parse_record_for_calls(record, 2, 0);
            assert_eq!(calls.len(), 1);
            assert_eq!(minor_calls.len(), 0);
            assert_eq!(calls[0].genome_index, 1);
            assert_eq!(calls[0].reference, "a".to_string());
            assert_eq!(calls[0].alt, "x".to_string());
            assert_eq!(calls[0].call_type, AltType::NULL);
            assert_eq!(calls[0].cov, Some(1));
            assert_eq!(calls[0].frs, Some(ordered_float::OrderedFloat(1.0)));
            assert!(!calls[0].is_minor);
            assert_eq!(calls[0].vcf_idx, None);
            assert_eq!(calls[0].vcf_row, 0);
        }

        // Ins call
        {
            let record = VCFRow {
                position: 1,
                reference: "A".to_string(),
                alternative: vec!["AT".to_string()],
                filter: vec![],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["1/1".to_string()]),
                    ("COV".to_string(), vec!["1".to_string(), "5".to_string()]),
                ]),
                is_filter_pass: true,
            };
            let (calls, minor_calls) = VCFFile::parse_record_for_calls(record, 3, 123);
            assert_eq!(calls.len(), 1);
            assert_eq!(minor_calls.len(), 0);
            assert_eq!(calls[0].genome_index, 1);
            assert_eq!(calls[0].reference, "a".to_string());
            assert_eq!(calls[0].alt, "t".to_string());
            assert_eq!(calls[0].call_type, AltType::INS);
            assert_eq!(calls[0].cov, Some(5));
            assert_eq!(calls[0].frs, Some(ordered_float::OrderedFloat(5.0 / 6.0)));
            assert!(!calls[0].is_minor);
            assert_eq!(calls[0].vcf_idx, Some(1));
            assert_eq!(calls[0].vcf_row, 123);
        }

        // Del call
        {
            let record = VCFRow {
                position: 1,
                reference: "AT".to_string(),
                alternative: vec!["A".to_string()],
                filter: vec![],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["1/1".to_string()]),
                    ("COV".to_string(), vec!["1".to_string(), "5".to_string()]),
                ]),
                is_filter_pass: true,
            };
            let (calls, minor_calls) = VCFFile::parse_record_for_calls(record, 3, 0);
            assert_eq!(calls.len(), 1);
            assert_eq!(minor_calls.len(), 0);
            assert_eq!(calls[0].genome_index, 2);
            assert_eq!(calls[0].reference, "t".to_string());
            assert_eq!(calls[0].alt, "t".to_string());
            assert_eq!(calls[0].call_type, AltType::DEL);
            assert_eq!(calls[0].cov, Some(5));
            assert_eq!(calls[0].frs, Some(ordered_float::OrderedFloat(5.0 / 6.0)));
            assert!(!calls[0].is_minor);
            assert_eq!(calls[0].vcf_idx, Some(1));
            assert_eq!(calls[0].vcf_row, 0);
        }

        // Del call with SNP
        {
            let record = VCFRow {
                position: 1,
                reference: "AT".to_string(),
                alternative: vec!["C".to_string()],
                filter: vec![],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["1/1".to_string()]),
                    ("COV".to_string(), vec!["1".to_string(), "5".to_string()]),
                ]),
                is_filter_pass: true,
            };
            let (calls, minor_calls) = VCFFile::parse_record_for_calls(record, 3, 0);
            assert_eq!(calls.len(), 2);
            assert_eq!(minor_calls.len(), 0);
            assert_eq!(calls[0].genome_index, 2);
            assert_eq!(calls[0].reference, "t".to_string());
            assert_eq!(calls[0].alt, "t".to_string());
            assert_eq!(calls[0].call_type, AltType::DEL);
            assert_eq!(calls[0].cov, Some(5));
            assert_eq!(calls[0].frs, Some(ordered_float::OrderedFloat(5.0 / 6.0)));
            assert!(!calls[0].is_minor);
            assert_eq!(calls[0].vcf_idx, Some(1));

            assert_eq!(calls[1].genome_index, 1);
            assert_eq!(calls[1].reference, "a".to_string());
            assert_eq!(calls[1].alt, "c".to_string());
            assert_eq!(calls[1].call_type, AltType::SNP);
            assert_eq!(calls[1].cov, Some(5));
            assert_eq!(calls[1].frs, Some(ordered_float::OrderedFloat(5.0 / 6.0)));
            assert!(!calls[1].is_minor);
            assert_eq!(calls[1].vcf_idx, Some(1));
            assert_eq!(calls[0].vcf_row, 0);
        }

        // More complex mix of SNP and indel minors
        {
            let record = VCFRow {
                position: 1,
                reference: "AT".to_string(),
                alternative: vec!["C".to_string(), "GCC".to_string()],
                filter: vec![],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["1/1".to_string()]),
                    (
                        "COV".to_string(),
                        vec!["1".to_string(), "57".to_string(), "12".to_string()],
                    ),
                ]),
                is_filter_pass: true,
            };
            let (calls, minor_calls) = VCFFile::parse_record_for_calls(record, 3, 0);
            assert_eq!(calls.len(), 2);
            assert_eq!(minor_calls.len(), 3);

            assert_eq!(calls[0].genome_index, 2);
            assert_eq!(calls[0].reference, "t".to_string());
            assert_eq!(calls[0].alt, "t".to_string());
            assert_eq!(calls[0].call_type, AltType::DEL);
            assert_eq!(calls[0].cov, Some(57));
            assert_eq!(calls[0].frs, Some(ordered_float::OrderedFloat(57.0 / 70.0)));
            assert!(!calls[0].is_minor);
            assert_eq!(calls[0].vcf_idx, Some(1));
            assert_eq!(calls[0].vcf_row, 0);

            assert_eq!(calls[1].genome_index, 1);
            assert_eq!(calls[1].reference, "a".to_string());
            assert_eq!(calls[1].alt, "c".to_string());
            assert_eq!(calls[1].call_type, AltType::SNP);
            assert_eq!(calls[1].cov, Some(57));
            assert_eq!(calls[1].frs, Some(ordered_float::OrderedFloat(57.0 / 70.0)));
            assert!(!calls[1].is_minor);
            assert_eq!(calls[1].vcf_idx, Some(1));
            assert_eq!(calls[1].vcf_row, 0);

            assert_eq!(minor_calls[0].genome_index, 2);
            assert_eq!(minor_calls[0].reference, "t".to_string());
            assert_eq!(minor_calls[0].alt, "c".to_string());
            assert_eq!(minor_calls[0].call_type, AltType::INS);
            assert_eq!(minor_calls[0].cov, Some(12));
            assert_eq!(
                minor_calls[0].frs,
                Some(ordered_float::OrderedFloat(12.0 / 70.0))
            );
            assert!(minor_calls[0].is_minor);
            assert_eq!(minor_calls[0].vcf_idx, Some(2));
            assert_eq!(minor_calls[0].vcf_row, 0);

            assert_eq!(minor_calls[1].genome_index, 1);
            assert_eq!(minor_calls[1].reference, "a".to_string());
            assert_eq!(minor_calls[1].alt, "g".to_string());
            assert_eq!(minor_calls[1].call_type, AltType::SNP);
            assert_eq!(minor_calls[1].cov, Some(12));
            assert_eq!(
                minor_calls[1].frs,
                Some(ordered_float::OrderedFloat(12.0 / 70.0))
            );
            assert!(minor_calls[1].is_minor);
            assert_eq!(minor_calls[1].vcf_idx, Some(2));
            assert_eq!(minor_calls[1].vcf_row, 0);

            assert_eq!(minor_calls[2].genome_index, 2);
            assert_eq!(minor_calls[2].reference, "t".to_string());
            assert_eq!(minor_calls[2].alt, "c".to_string());
            assert_eq!(minor_calls[2].call_type, AltType::SNP);
            assert_eq!(minor_calls[2].cov, Some(12));
            assert_eq!(
                minor_calls[2].frs,
                Some(ordered_float::OrderedFloat(12.0 / 70.0))
            );
            assert!(minor_calls[2].is_minor);
            assert_eq!(minor_calls[2].vcf_idx, Some(2));
            assert_eq!(minor_calls[2].vcf_row, 0);
        }
    }

    #[test]
    fn test_instanciate_vcffile() {
        let vcf = VCFFile::new("test/dummy.vcf".to_string(), false, 3);
        let expected_records = vec![
            VCFRow {
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
            },
            VCFRow {
                position: 13333,
                reference: "c".to_string(),
                alternative: vec!["a".to_string(), "g".to_string()],
                filter: vec![],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["1/2".to_string()]),
                    ("DP".to_string(), vec!["100".to_string()]),
                    (
                        "COV".to_string(),
                        vec!["2".to_string(), "50".to_string(), "48".to_string()],
                    ),
                    ("GT_CONF".to_string(), vec!["613".to_string()]),
                ]),
                is_filter_pass: false,
            },
            VCFRow {
                // Edge case of using `RO` and `AO` for coverage
                position: 13335,
                reference: "t".to_string(),
                alternative: vec!["a".to_string()],
                filter: vec!["MAX_DP".to_string()],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["1/1".to_string()]),
                    ("DP".to_string(), vec!["68".to_string()]),
                    ("RO".to_string(), vec!["66".to_string()]),
                    ("AO".to_string(), vec!["2".to_string()]),
                    ("COV".to_string(), vec!["66".to_string(), "2".to_string()]),
                    ("GT_CONF".to_string(), vec!["3.77".to_string()]),
                ]),
                is_filter_pass: false,
            },
            VCFRow {
                // Odd edge case which is a valid VCF row where it's a het GT but single COV value
                position: 13336,
                reference: "c".to_string(),
                alternative: vec!["a".to_string()],
                filter: vec![],
                fields: HashMap::from([
                    ("GT".to_string(), vec!["0/1".to_string()]),
                    ("DP".to_string(), vec!["50".to_string()]),
                    ("COV".to_string(), vec!["50".to_string()]),
                    ("GT_CONF".to_string(), vec!["613".to_string()]),
                ]),
                is_filter_pass: false,
            },
        ];
        for (idx, record) in expected_records.iter().enumerate() {
            assert_eq!(record.position, vcf.records[idx].position);
            assert_eq!(record.reference, vcf.records[idx].reference);
            assert_eq!(record.alternative, vcf.records[idx].alternative);
            assert_eq!(record.filter, vcf.records[idx].filter);
            assert_eq!(record.fields, vcf.records[idx].fields);
            assert_eq!(record.is_filter_pass, vcf.records[idx].is_filter_pass);
        }

        let expected_calls = [
            vec![Evidence {
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
            vec![Evidence {
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
            vec![Evidence {
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
            vec![Evidence {
                cov: None,
                frs: None,
                genotype: "./.".to_string(),
                call_type: AltType::NULL,
                vcf_row: 5,
                reference: "t".to_string(),
                alt: "x".to_string(),
                genome_index: 13148,
                is_minor: false,
                vcf_idx: None,
            }],
            vec![Evidence {
                cov: None,
                frs: None,
                genotype: "./.".to_string(),
                call_type: AltType::NULL,
                vcf_row: 6,
                reference: "g".to_string(),
                alt: "x".to_string(),
                genome_index: 13149,
                is_minor: false,
                vcf_idx: None,
            }],
            vec![Evidence {
                cov: None,
                frs: None,
                genotype: "./.".to_string(),
                call_type: AltType::NULL,
                vcf_row: 7,
                reference: "a".to_string(),
                alt: "x".to_string(),
                genome_index: 13150,
                is_minor: false,
                vcf_idx: None,
            }],
        ];

        for calls in expected_calls.iter() {
            let actual = vcf.calls.get(&calls[0].genome_index).unwrap();
            for (idx, call) in calls.iter().enumerate() {
                assert_eq!(call.cov, actual[idx].cov);
                assert_eq!(call.frs, actual[idx].frs);
                assert_eq!(call.genotype, actual[idx].genotype);
                assert_eq!(call.call_type, actual[idx].call_type);
                assert_eq!(call.vcf_row, actual[idx].vcf_row);
                assert_eq!(call.reference, actual[idx].reference);
                assert_eq!(call.alt, actual[idx].alt);
                assert_eq!(call.genome_index, actual[idx].genome_index);
                assert_eq!(call.is_minor, actual[idx].is_minor);
                assert_eq!(call.vcf_idx, actual[idx].vcf_idx);
            }
        }
        assert_eq!(vcf.calls.keys().len(), expected_calls.len());

        let expected_minor_calls = [vec![
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
        ]];

        for calls in expected_minor_calls.iter() {
            let actual = vcf.minor_calls.get(&calls[0].genome_index).unwrap();
            for (idx, call) in calls.iter().enumerate() {
                assert_eq!(call.cov, actual[idx].cov);
                assert_eq!(call.frs, actual[idx].frs);
                assert_eq!(call.genotype, actual[idx].genotype);
                assert_eq!(call.call_type, actual[idx].call_type);
                assert_eq!(call.vcf_row, actual[idx].vcf_row);
                assert_eq!(call.reference, actual[idx].reference);
                assert_eq!(call.alt, actual[idx].alt);
                assert_eq!(call.genome_index, actual[idx].genome_index);
                assert_eq!(call.is_minor, actual[idx].is_minor);
                assert_eq!(call.vcf_idx, actual[idx].vcf_idx);
            }
        }
    }
}
