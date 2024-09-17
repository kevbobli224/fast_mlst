use std::{collections::HashMap, error::Error};
use std::fs::File;
use std::io::{BufReader, Write};
use std::path::{Path, PathBuf};
use std::env;
use clap::Parser;
use std::process::{exit, Command, ExitStatus, Stdio};
use csv::ReaderBuilder;
#[derive(Parser)]
#[command(version="1.0", about="Diagnostic tool for Streptococcus suis whole-genome sequences", long_about=None,arg_required_else_help=true)]
struct Args{
    /// PE files
    #[arg(short='i', long, value_delimiter = ' ', num_args = 1..)]
    input: Option<Vec<String>>,

    /// Path to output directory
    #[arg(short, long)]
    out: Option<String>,
    
    /// Multifasta file containing the MLST database
    #[arg(short='x',long, default_value = "Streptococcus_suis.fasta")]
    mlst_db: Option<String>,
    
    /// Text file containing the definitions for the MLST database file
    #[arg(short='y',long, default_value = "ssuis.txt")]
    mlst_def: Option<String>,

    /// Database directory containing MLST profiles and definitions
    #[arg(short='z', long, default_value = "database")]
    db: Option<String>,
}


// TODO: Allow user specified path for MLST database, defs... etc
// Currently only searches for program dir
fn main() {
    let args: Args = Args::parse();
    
    let input = args.input.clone().unwrap();
    let mut input_single: bool = false;
    if args.input.is_none() {
        custom_exit("No input files given", 0);
    } else if input.len() == 1 {
        input_single = true;
    }

    // Search for database for required files for the pipeline
    let program_dir: PathBuf = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let program_database: PathBuf = PathBuf::from(&program_dir).join("database");

    let mlst_db: &String = args.mlst_db.as_ref().unwrap();
    database_check(mlst_db, &program_database, "Cannot find MLST database");
    database_check(args.mlst_def.as_ref().unwrap(), &program_database, "Cannot find MLST definitions");
    
    check_prereqs();
    // println!("database/{mlst_db}.length.b = {}", find_files(&format!("database/{mlst_db}.length.b"), &program_database));
    if !find_files(&format!("{mlst_db}.length.b"), &program_database){
        eprintln!("kma MLST database not initialize, creating new...");
        eprintln!("{}",&format!("kma index -i {} -o {:?}",  &program_database.join(mlst_db).to_str().unwrap(),  &program_database.join(mlst_db) ));
        let cmd: Result<ExitStatus, std::io::Error> = Command::new("bash").arg("-c").arg(&format!("kma index -i {} -o {:?}", &program_database.join(mlst_db).to_str().unwrap(),  &program_database.join(mlst_db))).stderr(Stdio::null()).stdout(Stdio::null()).spawn().unwrap().wait();
        if cmd.is_err() {
            custom_exit("Failed to index MLST database using kma", cmd.err().unwrap().raw_os_error().unwrap());
        }
    }
    
    eprintln!("Running kma...");
    
    let mut profiles = read_profiles(&program_database.join( args.mlst_def.as_ref().unwrap())).ok().unwrap();

    let mut mlst_tree = Tree::new();
    
    let mut header: Vec<String> = profiles.remove(&0).unwrap().to_vec();
    let hi = header.iter().position(|x| *x == "clonal_complex");

    if hi.is_some(){
        header.remove(hi.unwrap());
    }


    let res = parse_res(header.clone()).unwrap();
    // println!("{:?}", res);

    // exit(1);

    // println!("{}", &format!("kma -o {} -t_db {:?} -{} {} -na -nf -nc", args.out, &program_database.join(mlst_db).to_str().unwrap(), if input_single {"i"} else {"ipe"}, if input_single {args.r#in.as_ref().unwrap().to_string()} else {format!("{} {}", args.forward.as_ref().unwrap(), args.reverse.as_ref().unwrap())}));
    //let header = format!("#Name\tST\t{}\n", header.iter().map(|f| format!("{}\t", f)).collect::<Vec<String>>().join("") );

    let a = input.clone().iter().map(|f| format!("{f}")).collect::<Vec<String>>().join(" ");

    let cmd: Result<ExitStatus, std::io::Error> = Command::new("bash").arg("-c").arg(&format!("kma -o tmp_kma -t_db {:?} -{} {} -na -nf -nc", &program_database.join(args.mlst_db.as_ref().unwrap()).to_str().unwrap(), if input_single {"i"} else {"ipe"}, if input_single {&input[0]} else { &a })).stdout(Stdio::null()).stderr(Stdio::null()).spawn().unwrap().wait();
    if cmd.is_err() {custom_exit("kma failed to run", cmd.err().unwrap().raw_os_error().unwrap());}


    // Parse results
    

    for (k, mut v) in profiles.clone() {
        let empty = v.iter().position(|x| *x == "");
        if empty.is_some(){
            v.remove(empty.unwrap());
        }
        let a: Vec<u16> = v.iter().map(|f| f.parse::<u16>().unwrap()).collect();
        mlst_tree.add_chain(a, header.clone(), k);
        
    }
    // println!("{mlst_tree:?}");

    // Perform tree search for each 7 alleles
    // let res = parse_res().unwrap();
    let mut cur: Option<&Node> = None;
    for i in res.clone() {
        let allele = i[1].parse::<u16>().unwrap();
        let tmp: Option<&Node>;
        if cur.is_none(){
            tmp = mlst_tree.roots.get(&allele);
            if tmp.is_some(){
                cur = tmp;
            }
        } else {
            if !cur.unwrap().nexts.is_empty(){
                cur = cur.clone().unwrap().nexts.get(&allele);
            }
            if cur.is_none(){
                // println!("===={allele} {:?}", cur);
                break;
            }
        }
    }

    // Get file name without extension, 16/09/24: waiting for stable rust merge for file_prefix
    // https://github.com/rust-lang/rust/pull/129114
    let _b = Path::new(&args.input.as_ref().unwrap()[0]).with_extension("");
    let fbasename = _b.file_name().unwrap().to_str().unwrap();
    // println!("{fbasename:?} {_b:?}");

    let st: u16 = if cur.is_some() {cur.unwrap().st} else {0};
    let header: Vec<_> = res.clone().iter().map(|x| x[0].clone()).collect();
    let alleles: Vec<_> = res.clone().iter().map(|x| x[1].clone()).collect();
    let st: String = if st == 0 { String::from("NF") } else { st.to_string() };
    let header = format!("#Name\tST\t{}\n", header.iter().map(|f| format!("{}\t", f)).collect::<Vec<String>>().join("") );
    let res = format!("{fbasename}\t{st}\t{}", alleles.iter().map(|f| format!("{}\t", f)).collect::<Vec<String>>().join(""));

    if args.out.is_none(){
        println!("{header}{res}");
    } else {
        let out_dirname = Path::new(args.out.as_ref().unwrap()).parent().unwrap();
        let mut of: File; 
        if !out_dirname.has_root(){
            of = File::create(format!("./{fbasename}.res")).expect(&format!("Can't write to output file ./{fbasename}.res"));
        } else {
            of = File::create(format!("{}/{fbasename}.res",out_dirname.display())).expect(&format!("Can't write to output file ./{fbasename}.res"));
        }
        let _ = of.write(header.as_bytes());
        _ = of.write(res.as_bytes());
    }


    // println!("Res: {:?}\nST: {:?}", t, st);
}

#[derive(Debug)]
pub struct Tree{
    // Key = allele #
    roots: HashMap<u16, Node>
}
impl Tree {
    pub fn new() -> Self {
        Tree {roots: HashMap::new()}
    }
    pub fn add_chain(&mut self, all: Vec<u16>, _h: Vec<String>, st: u16) {
        // let mut hcp = h.clone();
        let mut allcp = all.clone();
        
        // Create new chain if first allele doesn't exist
        if !self.roots.contains_key(&all[0]){
            let mut allele = allcp.pop().unwrap();
            // let mut last = Node::new(hcp.pop().unwrap().clone());
            let mut last = Node::new();
            while allcp.len() > 0 {
                let cur = allcp.pop().unwrap();
                // last = Node::new(hcp.pop().unwrap()).add_node(allele, last);
                last = Node::new().add_node(allele, last);
                allele = cur;
            }
            last.st = st;
            self.roots.insert(allele, last);
        } 
        // Find missing alleles # and add
        else {
            let mut index: usize = 1;
            let mut cur = self.roots.get(&all[0]).unwrap().clone();
            while index < allcp.len() {
                if cur.nexts.contains_key(&all[index]){
                    cur = cur.nexts.get(&all[index]).unwrap().clone();
                } else {
                    // let mut n = Node::new(h[index].clone());
                    let mut n = Node::new();
                    if index == allcp.len()-1{
                        n.st = st;
                    }
                    cur.nexts.insert(all[index], n);
                    cur = cur.nexts.get(&all[index]).unwrap().clone();
                }
                index += 1;
            }

        }
    }
}


#[derive(Clone, Debug)]
struct Node {
    // name: String,
    nexts: HashMap<u16, Node>,
    st: u16
}
impl Node {
    pub fn new() -> Self {
        Node {
            // name: name,
            nexts: HashMap::new(),
            st: 0
        }
    }
    pub fn add_node(mut self, v: u16, n: Node) -> Self{
        self.nexts.insert(v, n);
        self
    }
}



// Fuck this shit just do Vec for all of them next time
fn read_profiles(p: &PathBuf) -> Result<HashMap<u16,Vec<String>>, Box<dyn Error>> {
    let profile_file = File::open(p)?;
    let tsv_reader: BufReader<File> = BufReader::new(profile_file);
    let mut reader: csv::Reader<BufReader<File>> = ReaderBuilder::new().delimiter(b'\t').has_headers(false).from_reader(tsv_reader);

    // Types need correction after parsing
    let mut profiles :HashMap<u16,Vec<String>> = HashMap::new();

    let mut is_header: bool = true;
    for r in reader.records(){
        let mut row: csv::StringRecord = r?;
        row.trim();
        let t: Vec<String> = row.iter().map(|s| s.to_string()).collect();
        // let alleles: Vec<_> = (1..9).filter_map(|i: usize| row.get(i)).collect();
        let alleles: Vec<String> = t.as_slice()[1..t.len()].to_vec().into_iter().collect();
        if is_header {
            is_header = false;
            profiles.insert(0, alleles);
        } else {
            let st: &str = row[0].as_ref();
            profiles.insert(st.parse().ok().unwrap(), alleles);
        }
    }
    Ok(profiles)
}


// Highest scoring alleles
fn parse_res(h: Vec<String>) -> Result<Vec<Vec<String>>, Box<dyn std::error::Error> > {
    let tmp_file: File = File::open("tmp_kma.res")?;
    let tsv_reader: BufReader<File> = BufReader::new(tmp_file);
    let mut reader: csv::Reader<BufReader<File>> = ReaderBuilder::new().delimiter(b'\t').from_reader(tsv_reader);

    let mut res: Vec<Vec<String>> = Vec::new();

    // k = gene name, v = (allele, score)
    let mut hm: HashMap<String, (u16, u32)> = HashMap::new();

    for r in reader.records(){
        let cur_res = r?;
        let gene_allele: String = cur_res[0].to_string();

        let score = cur_res[1].trim().parse::<u32>().unwrap();
        
        let allele_res = Vec::from_iter(gene_allele.split("-").map(String::from));
        
        if hm.contains_key(&allele_res[0]){
            if hm.get(&allele_res[0]).unwrap().1 < score {
                hm.remove(&allele_res[0]);
                hm.insert(allele_res[0].clone(), (allele_res[1].parse::<u16>()?, score));
            }
        } else {
            hm.insert(allele_res[0].clone(), (allele_res[1].parse::<u16>()?, score));
        }
    }

    for i in h.iter(){
        let o = hm.get(i);
        if o.is_some(){
            res.push(vec![i.to_string(), o.unwrap().0.to_string()]);
        }
    }

    Ok(res)
}

pub fn database_check(f_name: &String, p_name: &PathBuf, err_msg: &str){
    if !find_files(f_name, p_name) {custom_exit(err_msg, 1);}
}
pub fn find_files(file_name: &String, path_name: &PathBuf) -> bool {
    eprint!("Searching for {}...", file_name);
    let current_file_search: PathBuf = path_name.join(file_name);
        if !current_file_search.exists() {
            eprintln!(" Not found: {}", current_file_search.display());
            false
        } else {
            eprintln!(" Found!" );
            true
        }
}
pub fn custom_exit(msg: &str, code: i32) -> bool{
    eprintln!("Exiting: {}", msg);
    exit(code)
}
pub fn check_prereqs() {
    let cmd: Result<ExitStatus, std::io::Error> = Command::new("kma").stderr(Stdio::null()).stdout(Stdio::null()).status();
    if cmd.is_err() {custom_exit("kma not found, install it before running this pipeline", cmd.err().unwrap().raw_os_error().unwrap());}
}