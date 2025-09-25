mod prod_events;
mod reaction_data;
mod recursion;
mod state_data;

use crate::{reaction::Reaction, reaction_graph::ReactionGraph, SimulationAlg};
use prod_events::ProdEvents;
use rand::Rng;
use reaction_data::ReactionData;
use recursion::{RecursionTree, RecursionTreeNode};
use state_data::StateData;

pub struct FastGillespie3 {
    pub reactions: Vec<Reaction>,
    pub state: Vec<i64>,
    dependence_graph: ReactionGraph,
    total_reactions: u64,
}

impl SimulationAlg for FastGillespie3 {
    fn new(state: Vec<i64>, reactions: Vec<Reaction>, reactant_names: Vec<String>) -> Self {
        FastGillespie3 {
            dependence_graph: ReactionGraph::from_reactions(&state, &reactions),
            state,
            reactions,
            total_reactions: 0,
        }    
    }

    fn advance(&mut self, time: f64, rng: &mut impl Rng) {
        let reaction_data = (0..self.reactions.len())
            .map(|idx| ReactionData::new(idx, [ProdEvents::zero(); 3]))
            .collect();

        let mut recursion = RecursionTree::new(
            vec![RecursionTreeNode::new(
                reaction_data,
                false,
                None,
                None,
                None,
            )],
            vec![None; self.reactions.len()],
            &self.reactions,
            &self.dependence_graph,
            StateData::new(&self.state),
            vec![true; self.reactions.len()],
            vec![0; self.state.len()],
            0,
            vec![Vec::default(); self.state.len()],
        );
        recursion.recursion(0, time, rng);
        // println!("Events: {}", recursion.total_events);
        self.state.clone_from_slice(recursion.state());
        self.total_reactions += recursion.total_events;
    }

    fn state(&self) -> &[i64] {
        &self.state
    }

    fn total_reactions(&self) -> u64 {
        self.total_reactions
    }
}
