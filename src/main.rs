use std::{path::PathBuf, time::SystemTime};

use clap::{command, Parser};
use tausplit::{Algorithm, Gillespie, ParseState, Reaction, SimulationAlg, TauSplit5, TauSplit6};

use itertools::Itertools;
use rand::{rng, rngs::SmallRng, Rng, SeedableRng};

// See also `clap_cargo::style::CLAP_STYLING`
pub const CLAP_STYLING: clap::builder::styling::Styles = clap::builder::styling::Styles::styled()
    .header(clap_cargo::style::HEADER)
    .usage(clap_cargo::style::USAGE)
    .literal(clap_cargo::style::LITERAL)
    .placeholder(clap_cargo::style::PLACEHOLDER)
    .error(clap_cargo::style::ERROR)
    .valid(clap_cargo::style::VALID)
    .invalid(clap_cargo::style::INVALID);

#[derive(Debug, Parser)]
#[command(
    name = "tausplit",
    about = "Simulation of chemical reaction networks.",
    long_about = "A program for the exact simulation of chemical reaction networks based on the Tau-Splitting algorithm.
The program takes in as input several data files, in the format:

A = 6
B = 8
C = 0
A + B -> C, 0.05

The system above has three chemical species, A, B, and C.
A and B can react to form C, and the reaction rate is 0.05.",
    styles = CLAP_STYLING,
)]
struct Cli {
    /// The amount of time to simulate.
    time: f64,

    /// The path to the files specifying the chemical reaction network.
    #[arg(num_args = 1.., )]
    data: Vec<PathBuf>,

    /// How often to sample and store the state.
    /// If not given, only the final state is stored.
    #[arg(short, long)]
    samples: Option<u64>,

    /// Whether to count the number of reactions.
    #[arg(long)]
    count_reactions: bool,

    /// Whether to count the cpu time..
    #[arg(long)]
    cpu_time: bool,

    /// Whether to skip printing the final state.
    #[arg(long)]
    no_print_state: bool,

    /// The algorithm to use to simulate the system.
    #[arg(long)]
    algorithm: Option<Algorithm>,

    /// The seed to use for random number generation.
    #[arg(long)]
    seed: Option<u64>,

    #[arg(long)]
    repeats: Option<u64>,
}

fn run_with_alg<Alg: SimulationAlg>(args: Cli, initial_state: Vec<i64>, reactions: Vec<Reaction>, names: Vec<String>) {
    let rng = &mut if let Some(seed) = args.seed {
        SmallRng::seed_from_u64(seed)
    } else {
        SmallRng::seed_from_u64(rng().random())
    };


    let run_count = args.repeats.unwrap_or(1);

    for run_idx in 0..run_count {
        let start_time = SystemTime::now();
        let time = args.time;
        let sample_count = args.samples.unwrap_or(1);
        let mut samples = Vec::new();
        samples.push((initial_state.clone(), 0, 0.));
        
        let mut alg = Alg::new(
            initial_state.iter().map(|x| *x as i64).collect_vec(),
            reactions.clone(),
            names.clone(),
        );
        for _ in 0..sample_count {
            alg.advance(time / sample_count as f64, rng);
            samples.push((
                alg.state().to_owned(),
                alg.total_reactions(),
                start_time.elapsed().unwrap().as_secs_f32(),
            ));
        }
    
        // Printing the sampled states to stdout, to be redirected as desired.

        // We print the header only once when we have several repeated runs.
        if run_idx == 0 {
            print!("time");
            if !args.no_print_state {
                for name in &names {
                    print!("\t{name}");
                }
            }
            if args.count_reactions {
                print!("\treaction_count");
            }
            if args.cpu_time {
                print!("\tcpu_time");
            }
            if args.repeats.is_some() {
                print!("\trun_idx");
            }
            println!();
        }
    
    
        for (idx, (state, total_reactions, cpu_time)) in samples.into_iter().enumerate() {
            print!("{}", idx as f64 / sample_count as f64 * time);
            if !args.no_print_state {
                for count in state {
                    print!("\t{count}");
                }
            }
            if args.count_reactions {
                print!("\t{total_reactions}");
            }
            if args.cpu_time {
                print!("\t{cpu_time:.3}")
            }
            if args.repeats.is_some() {
                print!("\t{run_idx}");
            }
            println!();
        }
    }
}

fn run_cli(args: Cli) {
    
    let mut parse_state = ParseState::default();
    for path in &args.data {
        parse_state.parse_data_file(path);
    }

    let (initial_state, reactions, names) = parse_state.get_network();

    match args.algorithm {
        Some(Algorithm::Gillespie) => run_with_alg::<Gillespie>(args, initial_state, reactions, names),
        Some(Algorithm::TauSplit6) => run_with_alg::<TauSplit6>(args, initial_state, reactions, names),
        Some(Algorithm::TauSplit) => {
            for reaction in &reactions {
                println!("{}", reaction.format_pretty(&names));
            }

            if let Some(bad_reaction) = reactions.iter().find(|r|r.inputs.len() > TauSplit5::MAX_INPUTS) {
                panic!(
                    "Unable to run the optimized Tau-Splitting algorithm!
                    The maximal number of input reactants is {}, and the reaction {} has {}! Please use the Tau-Split6 algorithm instead.",
                    TauSplit5::MAX_INPUTS,
                    bad_reaction.format_pretty(&names),
                    bad_reaction.inputs.len()
                );
            }
            if let Some(bad_reaction) = reactions.iter().find(|r|r.stoichiometry.len() > TauSplit5::MAX_STOI) {
                panic!(
                    "Unable to run the optimized Tau-Splitting algorithm! The maximal number of molecular species in the stoichiometry vector is {}, and the reaction {} has {}! Please use the Tau-Split6 algorithm instead.",
                    TauSplit5::MAX_STOI,
                    bad_reaction.format_pretty(&names),
                    bad_reaction.stoichiometry.len()
                );
            }

            run_with_alg::<TauSplit5>(args, initial_state, reactions, names)
        }
        None => {
            // Choosing the algorithm automatically.
            if reactions.iter().any(|r|r.inputs.len() > TauSplit5::MAX_INPUTS || r.stoichiometry.len() > TauSplit5::MAX_STOI) {
                run_with_alg::<TauSplit6>(args, initial_state, reactions, names)
            } else {
                run_with_alg::<TauSplit5>(args, initial_state, reactions, names)
            }
        }
    }
}

pub fn main() {
    run_cli(Cli::parse());
}
