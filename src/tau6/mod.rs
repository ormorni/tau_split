mod f_reaction;
mod listener;
mod reaction_data;
mod recursion;
mod state_data;
mod unstable_dependents;

use f_reaction::FReaction;
use itertools::Itertools;
use rand::Rng;
use reaction_data::{ReactionData};
use recursion::RecursionTree;
use state_data::StateData;

use crate::{reaction::Reaction, SimulationAlg};



pub struct TauSplit6 {
    pub reactions: Vec<Reaction>,
    pub state: Vec<i64>,
    pub total_reactions: u64,
    reactant_names: Vec<String>,
}

impl TauSplit6 {
    /// The maximal number of inputs in a standard reaction. 
    /// Reactions with more inputs will require allocation.
    pub const MAX_INPUTS: usize = 2;
    /// The maximal number of molecular species in the stoichiometry vector of a standard reaction. 
    /// Reactions with more species will require allocation.
    pub const MAX_STOI: usize = 4;
}

impl SimulationAlg for TauSplit6 {
    fn new(
        state: Vec<i64>,
        reactions: Vec<Reaction>,
        reactant_names: Vec<String>,
    ) -> TauSplit6 {
        TauSplit6 {
            state,
            reactions,
            total_reactions: 0,
            reactant_names,
        }
    }

    fn advance(&mut self, time: f64, rng: &mut impl Rng) {
        let f_reactions = self
            .reactions
            .iter()
            .map(|r| FReaction::from(r.clone()))
            .collect_vec();
        let mut recursion =
            RecursionTree::new(&self.state, &f_reactions, &self.reactant_names, time, rng);
        recursion.recursion(0, time, rng);
        // println!("Events: {}", recursion.total_events);
        self.state.clone_from_slice(&recursion.state());
        self.total_reactions += recursion.total_events;
    }

    fn state(&self) -> &[i64] {
        &self.state
    }

    fn total_reactions(&self) -> u64 {
        self.total_reactions
    }
}

