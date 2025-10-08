use rand::Rng;

mod gillespie;
pub use gillespie::Gillespie;
#[allow(unused)]
mod fastspie;
#[allow(unused)]
mod fastspie2;
mod fastspie3;
pub use fastspie3::FastGillespie3;
mod fastspie4;
pub use fastspie4::FastGillespie4;
mod fastspie5;
pub use fastspie5::FastGillespie5;
mod fastspie6;
pub use fastspie6::FastGillespie6;
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
#[derive(Default, Debug, Clone, Copy, clap::ValueEnum)]
pub enum Algorithm {
    /// The Gillespie algorithm.
    Gillespie,
    /// The Tau-Splitting algorithm.
    #[default]
    TauSplit,
    /// The new version of the Tau-Splitting algorithm.
    TauSplit6,
}
