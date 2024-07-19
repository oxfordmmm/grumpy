//! Module for handling VCF files
use pyo3::prelude::*;

use vcf::VCFReader;
use std::fs::File;
use std::io::BufReader;
use std::string::String;
use std::collections::HashMap;

use crate::common::{AltType, Evidence, VCFRow};
#[pyclass]
#[derive(Clone, Debug)]
/// Struct to hold the information from a VCF file
pub struct VCFFile{
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
    pub minor_calls: HashMap<i64, Vec<Evidence>>
}

#[pymethods]
impl VCFFile{
    #[new]
    /// Create a new VCFFile object
    /// 
    /// # Arguments
    /// - `filename`: String - Path to the VCF file
    /// - `ignore_filter`: bool - Whether to ignore the filter column
    /// - `min_dp`: i32 - Minimum depth to consider a call
    pub fn new(filename: String, ignore_filter: bool, min_dp: i32) -> Self{
        let file = File::open(filename).unwrap();
        let buf = BufReader::new(file);
        let mut reader = VCFReader::new(buf).unwrap();

        let mut record = reader.empty_record();
        let mut more_records = reader.next_record(&mut record);

        let required_fields = vec!["GT"];
        let mut records: Vec<VCFRow> = Vec::new();

        // For ease of access, we'll store the calls in a hashmap indexed on genome index
        let mut calls_map: HashMap<i64, Vec<Evidence>> = HashMap::new();
        let mut minor_calls_map: HashMap<i64, Vec<Evidence>> = HashMap::new();
        while match more_records{
                Ok(true) => true,
                _ => false
                } {
                    // Parse fields of the record to pull out data
                    // For whatever reason non of these fields are strings, so we need to convert them
                    let mut alts = Vec::new(); // One string per alt
                    for alt in record.alternative.iter(){
                        alts.push(String::from_utf8_lossy(alt).to_string().to_lowercase());
                    }
                    let mut filters = Vec::new(); // String per filter
                    for filter in record.filter.iter(){
                        for f in String::from_utf8_lossy(filter).to_string().split(";"){
                            filters.push(f.to_string());
                        }
                    }

                    // Oddities of how this VCF library works...
                    // Format is a vector of bytes, but we need to convert it to a vector of strings
                    // Each of these strings corresponds to an item in the genotypes vector
                    let mut format: Vec<String> = Vec::new();
                    for f in record.format.iter(){
                        format.push(String::from_utf8_lossy(f).to_string());
                    }

                    let mut idx = 0;
                    let mut fields: HashMap<String, Vec<String>> = HashMap::new();
                    for s in record.genotype.iter(){
                        let mut item: Vec<Vec<String>> = Vec::new();
                        for i in s.iter(){
                            let mut value: Vec<String> = Vec::new();
                            for j in i.iter(){
                                value.push(String::from_utf8_lossy(j).to_string());
                            }
                            item.push(value.clone());
                            fields.insert(format[idx].clone(), value.clone());
                            idx += 1;
                        }
                    }

                    // Validate that this record has the required fields
                    for field in required_fields.iter(){
                        if !format.contains(&field.to_string()){
                            panic!("Required field {} not found in record", field);
                        }
                    }

                    // Enforce that this record has a COV field
                    if !fields.contains_key("COV"){
                        let mut _cov_tag = "";
                        let mut cov = Vec::new();
                        if fields.contains_key("AD"){
                            _cov_tag = "AD";
                        }
                        else if fields.contains_key("RO") && fields.contains_key("AO"){
                            // Annoying edge case where no single COV field exists but can be constructed
                            _cov_tag = "RO";
                        }
                        else{
                            panic!("No COV tag found in record");
                        }
                        if _cov_tag != "RO"{
                            for c in fields.get(_cov_tag).unwrap(){
                                cov.push(c.to_string());
                            }
                        }
                        else{
                            let ro = fields.get("RO").unwrap()[0].clone();
                            let ao = fields.get("AO").unwrap()[0].clone();
                            cov.push(ro.to_string());
                            cov.push(ao.to_string());
                        }

                        fields.insert("COV".to_string(), cov);
                    }


                    // Check if this record has passed the filters
                    let mut passed = false;
                    if filters.len() == 0 || filters == vec!["PASS"]{
                        passed = true;
                    }
                    


                    // println!("{:?}\t{:?}\t{:?}\t{:?}\t{:?}", record.position, String::from_utf8_lossy(&record.reference), alts, filters, fields);

                    records.push(VCFRow{
                        position: record.position as i64,
                        reference: String::from_utf8_lossy(&record.reference).to_string().to_lowercase(),
                        alternative: alts.clone(),
                        filter: filters.clone(),
                        fields: fields.clone(),
                        is_filter_pass: passed
                    });

                    let (record_calls, record_minor_calls) = VCFFile::parse_record_for_calls(records[records.len()-1].clone(), min_dp);
                    // println!("Calls {:?}", record_calls);
                    // println!("Minor calls {:?}\n", record_minor_calls);

                    if ignore_filter || passed {
                        for call in record_calls{
                            if calls_map.contains_key(&call.genome_index){
                                calls_map.get_mut(&call.genome_index).unwrap().push(call.clone());
                                println!("Multiple calls at genome position {}! {:?}\n", call.genome_index, calls_map.get(&call.genome_index).unwrap());
                            }
                            else{
                                calls_map.insert(call.genome_index, vec![call.clone()]);
                            }
                        }
                        // calls.extend(record_calls);
                    }
                    // Add minor calls if the filter is passed or ignored, or specifcally just the MIN_FRS has failed
                    if (ignore_filter || passed) || (!passed && filters.len() == 1 && filters.contains(&"MIN_FRS".to_string())){
                        for call in record_minor_calls{
                            if minor_calls_map.contains_key(&call.genome_index){
                                minor_calls_map.get_mut(&call.genome_index).unwrap().push(call.clone());
                            }
                            else{
                                minor_calls_map.insert(call.genome_index, vec![call.clone()]);
                            }
                        }
                    }
                    // minor_calls.extend(record_minor_calls);

                    // Get the next record
                    more_records = reader.next_record(&mut record);
                    // break;
        
        }

        return VCFFile{
            header: Vec::new(),
            records,
            calls: calls_map,
            minor_calls: minor_calls_map
        }
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
    pub fn parse_record_for_calls(record: VCFRow, min_dp: i32) -> (Vec<Evidence>, Vec<Evidence>){
        let mut calls: Vec<Evidence> = Vec::new();
        let mut minor_calls: Vec<Evidence> = Vec::new();

        // Dealing with the actual call here; should spit out possible minor calls afterwards

        // Convert ugly genotype into a list of strings for each allele
        let genotype = &(record.fields.get("GT").unwrap())[0].split("/").collect::<Vec<&str>>();
        let dp = record.fields.get("DP").unwrap()[0].parse::<i32>().unwrap();

        let mut cov = Vec::new();
        for c in record.fields.get("COV").unwrap(){
            cov.push(c.parse::<i32>().unwrap());
        }

        let ref_allele = record.reference.clone().to_lowercase();
        let mut alt_allele = "".to_string();

        let mut call_type = AltType::NULL;
    
        let first = genotype[0];
        // Adjust for 1-indexed VCF
        let mut alt_idx = 0;
        if first != "."{
            alt_idx = first.parse::<i32>().unwrap() - 1;
        }
        if cov.len() == 1{
            // Just 1 item in the call so it's a null call
            calls.push(Evidence{
                cov: Some(cov[0 as usize]),
                frs: Some(ordered_float::OrderedFloat(1.0)),
                genotype: genotype.join("/"),
                call_type: AltType::NULL,
                vcf_row: record.clone(),
                reference: ref_allele.chars().nth(0).unwrap().to_string(),
                alt: "x".to_string(),
                genome_index: record.position,
                is_minor: false,
                vcf_idx: 0 as i64
            });
            return (calls, minor_calls);
        }

        for i in 1..genotype.len(){
            if genotype[i] != first{
                call_type = AltType::HET;
                alt_allele = "z".to_string();
                break;
            }
        }
        if call_type != AltType::HET{
            if first == "0"{
                call_type = AltType::REF;
                alt_allele = ref_allele.clone();
            }
            else{
                // Placeholder to denote we have an actual call at this point
                // Could actually be an indel instead but that can't be inferred from just the genotype
                call_type = AltType::SNP;
            }
        }
        if call_type == AltType::NULL{
            alt_allele = "x".to_string();
        }

        if call_type == AltType::SNP{
            // Parse the corresponding alternate allele to get the call(s)
            alt_allele = record.alternative[alt_idx as usize].clone().to_lowercase();
            let call = VCFFile::simplify_call(ref_allele.clone(), alt_allele.clone());
            let call_cov = cov[(alt_idx + 1) as usize]; // COV should be [ref, alt1, alt2..]
            for (offset, _alt_type, _base) in call{
                if call_cov < min_dp{
                    // Override calls with null if the coverage is too low
                    let _alt_type = AltType::NULL;
                    let _base = "x".to_string();
                }
                calls.push(Evidence{
                    cov: Some(call_cov),
                    frs: Some(ordered_float::OrderedFloat(call_cov as f32 / dp as f32)),
                    genotype: genotype.join("/"),
                    call_type: _alt_type,
                    vcf_row: record.clone(),
                    reference: ref_allele.chars().nth(offset).unwrap().to_string(),
                    alt: _base,
                    genome_index: record.position + offset as i64,
                    is_minor: false,
                    vcf_idx: (alt_idx + 1) as i64
                });
            }
            // println!("Parsed calls {:?}\n", c);
        }
        else{
            let call_cov = cov[(alt_idx + 1) as usize]; // COV should be [ref, alt1, alt2..]
            if call_cov < min_dp{
                // Override calls with null if the coverage is too low
                call_type = AltType::NULL;
                alt_allele = "x".to_string();
            }
            calls.push(Evidence{
                cov: Some(call_cov),
                frs: Some(ordered_float::OrderedFloat(call_cov as f32 / dp as f32)),
                genotype: genotype.join("/"),
                call_type: call_type.clone(),
                vcf_row: record.clone(),
                reference: ref_allele.clone(),
                alt: alt_allele.clone(),
                genome_index: record.position,
                is_minor: false,
                vcf_idx: (alt_idx + 1) as i64
            });
        }

        if call_type == AltType::NULL{
            // Don't check for minor calls if the call is null
            return (calls, minor_calls);
        }

        // Now we've got the main call, we need to figure out the minor calls
        // So check all possible values of COV for threshold to be considered a minor call
        let minor_threshold = 2;
        let mut idx = 0;
        for coverage in cov.iter(){
            if idx == (alt_idx + 1) as usize || idx == 0{
                // Skip ref and the alt we've already called
                idx += 1;
                continue;
            }
            if *coverage >= minor_threshold{
                alt_allele = record.alternative[(idx - 1) as usize].clone().to_lowercase();
                let call = VCFFile::simplify_call(ref_allele.clone(), alt_allele.clone());
                let call_cov = *coverage;
                for (offset, alt_type, base) in call{
                    minor_calls.push(Evidence{
                        cov: Some(call_cov),
                        frs: Some(ordered_float::OrderedFloat(call_cov as f32 / dp as f32)),
                        genotype: genotype.join("/"), // This is a minor call so the row's genotype is the same
                        call_type: alt_type,
                        vcf_row: record.clone(),
                        reference: ref_allele.chars().nth(offset).unwrap().to_string(),
                        alt: base,
                        genome_index: record.position + offset as i64,
                        is_minor: true,
                        vcf_idx: idx as i64
                    });
                }
            }
            idx += 1;
        }
        
        return (calls, minor_calls);
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
    pub fn simplify_call(reference: String, alternate: String) -> Vec<(usize, AltType, String)>{
        let mut calls: Vec<(usize, AltType, String)> = Vec::new();
        if reference.len() == alternate.len(){
            // Simple set of SNPs
            for i in 0..reference.len(){
                if reference.chars().nth(i).unwrap() != alternate.chars().nth(i).unwrap(){
                    calls.push((i, AltType::SNP, alternate.chars().nth(i).unwrap().to_string()));
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
        if reference.len() > alternate.len(){
            _y = reference.clone();
            _x = alternate.clone();
            _indel_type = AltType::DEL;
        }
        else{
            _y = alternate.clone();
            _x= reference.clone();
            _indel_type = AltType::INS;
        }

        let padding = "N".repeat(length as usize);
        let mut current = "".to_string();
        let mut current_dist = i64::MAX;
        let mut indel_start = 0;
        for i in 0.._x.len()+1{
            let x1 = _x[0..i].to_string() + &padding + &_x[i.._x.len()].to_string();
            // println!("{:?}\t{:?}\t{:?}\t{:?}\t{:?}", reference, alternate, x, y, x1);
            let dist = snp_dist(&x1, &_y);
            if dist <= current_dist{
                current = x1.clone();
                current_dist = dist;
                indel_start = i;
            }
        }

        if _indel_type == AltType::INS{
            // Ins after, del at, so adjust
            indel_start -= 1;
        }

        let mut indel_str = "".to_string();
        for i in 0..current.len(){
            if current.chars().nth(i).unwrap() == 'N'{
                indel_str += &_y.chars().nth(i).unwrap().to_string();
            }
        }

        calls.push((indel_start as usize, _indel_type, indel_str));


        for i in 0.._x.len(){
            let r = _y.chars().nth(i).unwrap();
            let a = _x.chars().nth(i).unwrap();
            if r != 'N' && a != 'N' && r != a{
                calls.push((i, AltType::SNP, a.to_string()));
            }
        }
        return calls;

    }
}

/// Calculate the distance between two strings ignoring N characters
/// 
/// # Arguments
/// - `reference`: &String - Reference sequence
/// - `alternate`: &String - Alternate sequence
/// 
/// # Returns
/// SNP distance between the two strings
fn snp_dist(reference: &String, alternate: &String) -> i64{
    let mut dist = 0;
    for i in 0..reference.len(){
        let r = reference.chars().nth(i).unwrap();
        let a = alternate.chars().nth(i).unwrap();
        if r != 'N' && a != 'N' && r != a {
            dist += 1;
        }
    }
    return dist;
}