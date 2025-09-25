mod reaction_data;
mod recursion;
mod state_data;

use crate::{reaction::Reaction, SimulationAlg};
use rand::Rng;
use reaction_data::ReactionData;
use recursion::{RecursionTree, RecursionTreeNode};
use state_data::StateData;

pub struct FastGillespie4 {
    pub reactions: Vec<Reaction>,
    pub state: Vec<i64>,
    pub total_events: u64,
}

impl SimulationAlg for FastGillespie4 {
    fn new(initial_state: Vec<i64>, reactions: Vec<Reaction>, reactant_names: Vec<String>) -> Self {
        FastGillespie4 { state: initial_state, reactions, total_events: 0}
    }

    fn advance(&mut self, time: f64, rng: &mut impl Rng) {
        let reaction_data = self
            .reactions
            .iter()
            .enumerate()
            .map(|(reaction_idx, reaction)| {
                ReactionData::sample(
                    reaction.input_product(&self.state) as f64,
                    reaction_idx,
                    reaction,
                    time,
                    rng,
                )
            })
            .collect();

        let mut recursion = RecursionTree::new(
            vec![RecursionTreeNode::new(
                reaction_data,
                false,
                None,
                None,
                None,
                0,
            )],
            vec![None; self.reactions.len()],
            &self.reactions,
            StateData::new(&self.state),
            vec![true; self.reactions.len()],
            vec![0; self.state.len()],
            0,
            vec![Vec::default(); self.state.len()],
            vec![Default::default(); self.state.len()],
            vec![Default::default(); self.state.len()],
            1,
        );
        recursion.recursion(0, time, rng);
        // println!("Events: {}", recursion.total_events);
        self.state.clone_from_slice(&recursion.state());
        self.total_events += recursion.total_events;
    }

    fn state(&self) -> &[i64] {
        &self.state
    }

    fn total_reactions(&self) -> u64 {
        self.total_events
    }
}

