//! # FASTA Sequence Analyzer
//!
//! This program analyzes one or more FASTA format files containing DNA sequences.
//! It calculates various statistics such as sequence length, GC content,
//! and N50 values. Results can be displayed on the console and optionally
//! appended to a CSV file for further analysis.

//use std::collections::HashMap;
use bio::io::fasta;
use std::env;
use std::fs::{File, OpenOptions};
use std::io::{self, stdout, BufReader, BufWriter, Write};
use std::os::unix::io::AsRawFd;
use std::path::Path;

/// Holds the results of the sequence analysis.
struct AnalysisResults {
    /// Names of the analyzed FASTA files.
    filename: String,
    /// Total length of all sequences combined.
    total_length: usize,
    /// Number of sequences in all files.
    sequence_count: usize,
    /// Total count of G and C bases.
    gc_count: usize,
    /// Total count of N bases.
    n_count: usize,
    /// N25 statistic (length at 25% of total sequence length).
    n25: usize,
    n25_secuence_count: usize,
    /// N50 statistic (length at 50% of total sequence length).
    n50: usize,
    n50_secuence_count: usize,
    /// N75 statistic (length at 75% of total sequence length).
    n75: usize,
    n75_secuence_count: usize,
    /// Length of the largest contig (sequence).
    largest_contig: usize,
    /// Length of the shortest contig (sequence).
    shortest_contig: usize,
    // Histogram of sequence lengths.
    //length_histogram: HashMap<usize, usize>,
}

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} [-c csv_file] <fasta_file1>", args[0]);
        std::process::exit(1);
    }

    let mut csv_filename = None;
    let mut fasta_file = String::new();

    // Parse command line arguments
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "-c" | "--csv" => {
                if i + 1 < args.len() {
                    csv_filename = Some(args[i + 1].clone());
                    i += 2;
                }
            }
            _ => {
                if fasta_file.is_empty() && args[i].parse::<usize>().is_ok() {
                    eprintln!("Error: Specify fasta");
                    std::process::exit(1);
                } else {
                    fasta_file = args[i].to_string();
                }
                i += 1;
            }
        }
    }

    let file = File::open(&fasta_file)?;
    let buf_reader: BufReader<File> = BufReader::new(file);
    let fasta_reader = fasta::Reader::new(buf_reader);

    let mut results = AnalysisResults::new();
    results.filename = fasta_file;

    let mut lengths: Vec<usize> = Vec::new();

    for result in fasta_reader.records() {
        let record = result.expect("Error reading record");
        let seq: std::borrow::Cow<'_, str> = String::from_utf8_lossy(record.seq());
        process_secuence(seq, &mut results, &mut lengths);
    }

    calc_nq_stats(&mut lengths, &mut results);

    if let Some(csv_file) = csv_filename {
        append_to_csv(&results, &csv_file)?;
    } else {
        print_results(&results);
    }
    Ok(())
}

/// Analyzes the given sequences and computes various statistics.
///
/// # Arguments
///
/// * `sequences` - A slice of `Sequence` structs to analyze
/// * `interval_size` - The interval size for the length histogram
///
/// # Returns
///
/// An `AnalysisResults` struct containing the computed statistics.
///```
/// fn process_secuence(
///    sequence: std::borrow::Cow<'_, str>,
///    results: &mut AnalysisResults,
///    lengths: &mut Vec<usize>,
/// ) {
/// 
/// ```
/// El cow permite copias sin muchas allocations
/// 
fn process_secuence(
    sequence: std::borrow::Cow<'_, str>,
    results: &mut AnalysisResults,
    lengths: &mut Vec<usize>,
) {
    let length = sequence.len();
    results.sequence_count += 1;
    results.total_length += length;
    results.gc_count += sequence
        .chars()
        .filter(|&c| c == 'G' || c == 'C' || c == 'g' || c == 'c')
        .count();
    results.n_count += sequence.chars().filter(|&c| c == 'N' || c == 'n').count();

    results.largest_contig = results.largest_contig.max(length);
    results.shortest_contig = results.shortest_contig.min(length);

    lengths.push(length);
}

/// Calculate N25, N50, N75 stats
///
/// # Arguments
///
/// * `results` - The `AnalysisResults` struct containing the analysis results
/// * `interval_size` - The interval size used for the length histogram
fn calc_nq_stats(lengths: &mut Vec<usize>, results: &mut AnalysisResults) {
    lengths.sort_unstable_by(|a, b| b.cmp(a)); // Sort in descending order
    let mut cumulative_length = 0;
    let mut cumulative_secuences = 0;
    for length in lengths {
        cumulative_length += *length;
        cumulative_secuences += 1;
        if results.n25 == 0 && cumulative_length >= results.total_length * 1 / 4 {
            results.n25 = *length;
            results.n25_secuence_count = cumulative_secuences;
        }
        if results.n50 == 0 && cumulative_length >= results.total_length * 1 / 2 {
            results.n50 = *length;
            results.n50_secuence_count = cumulative_secuences;
        }
        if results.n75 == 0 && cumulative_length >= results.total_length * 3 / 4 {
            results.n75 = *length;
            results.n75_secuence_count = cumulative_secuences;
            break;
        }
    }
}

/// Prints the analysis results to the console.
///
/// # Arguments
///
/// * `results` - The `AnalysisResults` struct containing the analysis results
/// * `interval_size` - The interval size used for the length histogram
//fn print_results(results: &AnalysisResults) {
fn print_results(results: &AnalysisResults) {
    let mut stdout = BufWriter::new(stdout());
    //println!("\nAnalysis Results for: {}", results.filenames.join(", "));

    writeln!(
        &mut stdout,
        "\nTotal length of sequence:\t{} bp",
        results.total_length
    )
    .unwrap();
    writeln!(
        &mut stdout,
        "Total number of sequences:\t{}",
        results.sequence_count
    )
    .unwrap();
    writeln!(
        &mut stdout,
        "Average contig length is:\t{} bp",
        results.total_length / results.sequence_count
    )
    .unwrap();
    writeln!(
        &mut stdout,
        "Largest contig:\t\t{} bp",
        results.largest_contig
    )
    .unwrap();
    writeln!(
        &mut stdout,
        "Shortest contig:\t\t{} bp",
        results.shortest_contig
    )
    .unwrap();
    writeln!(
        &mut stdout,
        "N25 stats:\t\t\t25% of total sequence length is contained in the {} sequences >= {} bp",
        results.n25_secuence_count, results.n25
    )
    .unwrap();
    writeln!(
        &mut stdout,
        "N50 stats:\t\t\t50% of total sequence length is contained in the {} sequences >= {} bp",
        results.n50_secuence_count, results.n50
    )
    .unwrap();
    writeln!(
        &mut stdout,
        "N75 stats:\t\t\t75% of total sequence length is contained in the {} sequences >= {} bp",
        results.n75_secuence_count, results.n75
    )
    .unwrap();
    writeln!(&mut stdout, "Total GC count:\t\t\t{} bp", results.gc_count).unwrap();
    writeln!(
        &mut stdout,
        "GC %:\t\t\t\t{:.2} %",
        (results.gc_count as f64 / results.total_length as f64) * 100.0
    )
    .unwrap();
    writeln!(&mut stdout, "Number of Ns:\t\t\t{}", results.n_count).unwrap();
    writeln!(
        &mut stdout,
        "Ns %:\t\t\t\t{:.2} %",
        (results.n_count as f64 / results.total_length as f64) * 100.0
    )
    .unwrap();
    stdout.flush().unwrap();

    //println!("\nLength histogram:");
    //let mut intervals: Vec<_> = results.length_histogram.keys().collect();
    //intervals.sort_unstable();
    //for &interval in intervals {
    //    let count = results.length_histogram[&interval];
    //println!(
    //"{}:{}\t{}",
    //        interval * interval_size,
    //        (interval + 1) * interval_size - 1,
    //        count
    //    );
    //}
}

/// Appends the analysis results to a CSV file.
///
/// If the file doesn't exist, it will be created and a header row will be written.
/// If the file exists, the results will be appended as a new row.
///
/// # Arguments
///
/// * `results` - The `AnalysisResults` struct containing the analysis results
/// * `csv_filename` - The name of the CSV file to append to
///
/// # Returns
///
/// An `io::Result<()>` indicating success or containing an error if the operation failed.
fn append_to_csv(results: &AnalysisResults, csv_filename: &str) -> io::Result<()> {
    let file_exists = Path::new(csv_filename).exists();
    let mut file = OpenOptions::new()
        .write(true)
        .append(true)
        .create(true)
        .open(csv_filename)?;

    // Use `flock` to lock the file before writing
    let fd = file.as_raw_fd();
    unsafe { libc::flock(fd, libc::LOCK_EX) }; // Acquire exclusive lock

    if !file_exists {
        // Write header if the file is newly created
        writeln!(file, "filename;assembly_length;number_of_sequences;average_length;largest_contig;shortest_contig;N50;GC_percentage;total_N;N_percentage")?;
    }
    writeln!(
        file,
        "{:?};{};{};{};{};{};{};{};{};{}",
        //results.filenames.join("+"),
        Path::new(&results.filename).file_name().unwrap(),
        results.total_length,
        results.sequence_count,
        results.total_length as f64 / results.sequence_count as f64,
        results.largest_contig,
        results.shortest_contig,
        results.n50,
        (results.gc_count as f64 / results.total_length as f64) * 100.0,
        results.n_count,
        (results.n_count as f64 / results.total_length as f64) * 100.0,
    )?;

    // Release the file lock
    unsafe { libc::flock(fd, libc::LOCK_UN) }; // Unlock the file
    Ok(())
}

impl AnalysisResults {
    fn new() -> Self {
        AnalysisResults {
            filename: String::new(), // Will be set later
            total_length: 0,
            sequence_count: 0,
            gc_count: 0,
            n_count: 0,
            n25: 0,
            n25_secuence_count: 0,
            n50: 0,
            n50_secuence_count: 0,
            n75: 0,
            n75_secuence_count: 0,
            largest_contig: 0,
            shortest_contig: usize::MAX,
        }
    }
}
