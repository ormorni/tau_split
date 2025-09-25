mod f_reaction;
mod listener;
mod reaction_data;
mod recursion;
mod state_data;
mod unstable_dependents;

use f_reaction::FReaction;
use itertools::Itertools;
use rand::Rng;
use reaction_data::{ReactionData, StableReactionData};
use recursion::RecursionTree;
use state_data::StateData;
use std::fmt::{Debug, Display};

use crate::{reaction::Reaction, SimulationAlg};

#[derive(PartialEq, Eq, PartialOrd, Ord, Hash, Clone, Copy, Default)]
struct NodeId(usize);

impl Debug for NodeId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!("NodeId(0b{:b})", self.0))
    }
}
impl Display for NodeId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!("ID({:b})", self.0))
    }
}

pub struct FastGillespie5 {
    pub reactions: Vec<Reaction>,
    pub state: Vec<i64>,
    pub total_reactions: u64,
    reactant_names: Vec<String>,
}

impl SimulationAlg for FastGillespie5 {
    fn new(
        state: Vec<i64>,
        reactions: Vec<Reaction>,
        reactant_names: Vec<String>,
    ) -> FastGillespie5 {
        FastGillespie5 {
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

const NO_LISTENER: (usize, NodeId) = (usize::MAX, NodeId(usize::MAX));
