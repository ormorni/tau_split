use rand::Rng;

mod gillespie;
pub use gillespie::Gillespie;
mod tau5;
pub use tau5::TauSplit5;
mod tau6;
pub use tau6::TauSplit6;
mod indexed_vec;


mod parsers;
pub use parsers::ParseState;
mod reaction;
pub use crate::reaction::Reaction;
mod reaction_graph;
mod tests;
mod utils;
pub use utils::DEFAULT_SEED;

/// A generic trait for an algorithm simulating a chemical reaction network.
pub trait SimulationAlg {
    /// Initializes a new instance of the algorithm.
    fn new(initial_state: Vec<i64>, reactions: Vec<Reaction>, reactant_names: Vec<String>) -> Self;
    /// Advances the simulation by the given amount of time.
    fn advance(&mut self, time: f64, rng: &mut impl Rng);
    /// Gets the state of the algorithm.
    fn state(&self) -> &[i64];
    /// Gets the total number of reactions simulated.
    fn total_reactions(&self) -> u64;
}

/// The algorithms available in the package.
#[derive(Debug, Clone, Copy, clap::ValueEnum)]
pub enum Algorithm {
    /// The Gillespie algorithm.
    Gillespie,
    /// The optimized Tau-Splitting algorithm.
    TauSplit,
    /// The version of the Tau-Splitting algorithm described in the manuscript, with no furhter optimizations.
    TauSplit6,
}
