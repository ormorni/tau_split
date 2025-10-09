use std::ops::{Index, IndexMut};

use derive_new::new;

use crate::reaction::binomial;

use super::{f_reaction::FReaction, reaction_data::TauData};

#[derive(new, Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct ComponentData {
    pub lower: i64,
    pub value: i64,
    pub upper: i64,
}

#[derive(Clone, Hash, Debug, PartialEq, Eq)]
pub struct StateData {
    pub state: Vec<ComponentData>,
}

impl StateData {
    pub fn new(state: &[i64]) -> StateData {
        StateData {
            state: state.iter().map(|&i| ComponentData::new(i, i, i)).collect(),
        }
    }

    /// Removes the effect of a ReactionData object from the error bounds.
    pub fn change_bounds(&mut self, event_count: i64, reaction: &FReaction) {
        if event_count != 0 {
            self.apply_negative(event_count, reaction);
            self.apply_positive(event_count, reaction);
        }
    }

    /// Removes the effect of a ReactionData object from the error bounds.
    pub fn remove_bounds(&mut self, rdata: &impl TauData, reaction: &FReaction) {
        self.change_bounds(-(rdata.event_count() as i64), reaction);
    }

    /// Add the effect of a ReactionData object to the error bounds.
    pub fn add_bounds(&mut self, rdata: &impl TauData, reaction: &FReaction) {
        self.change_bounds(rdata.event_count() as i64, reaction);
    }

    /// Applies the reaction to the state.
    pub fn apply(&mut self, rdata: &impl TauData, reaction: &FReaction) {
        for &(comp, diff) in &reaction.stoichiometry {
            self.state[comp].lower += diff * rdata.event_count() as i64;
            self.state[comp].value += diff * rdata.event_count() as i64;
            self.state[comp].upper += diff * rdata.event_count() as i64;
        }
    }

    /// Applies only the negative parts of the product to the reactants.
    pub fn apply_negative(&mut self, event_count: i64, reaction: &FReaction) {
        for &(reactant, change) in &reaction.negative_stoichiometry {
            self.state[reactant].lower += change * event_count;
        }
    }

    /// Applies only the negative parts of the product to the reactants.
    pub fn apply_positive(&mut self, event_count: i64, reaction: &FReaction) {
        for &(reactant, change) in &reaction.positive_stoichiometry {
            self.state[reactant].upper += change * event_count;
        }
    }

    /// Returns the upper bound of the input product of the reaction.
    pub fn upper_product(&self, reaction: &FReaction) -> f64 {
        reaction
            .inputs
            .iter()
            .map(|inp| binomial(self.state[inp.index].upper.max(0) as u64, inp.count))
            .product::<u64>() as f64
    }
    /// Returns the input product of the reaction.
    pub fn state_product(&self, reaction: &FReaction) -> f64 {
        reaction
            .inputs
            .iter()
            .map(|inp| binomial(self.state[inp.index].value.max(0) as u64, inp.count))
            .product::<u64>() as f64
    }

    /// Returns the lower bound of the input product of the reaction.
    pub fn lower_product(&self, reaction: &FReaction, has_events: bool) -> f64 {
        reaction
            .inputs
            .iter()
            .map(|inp| {
                binomial(
                    (self.state[inp.index].lower
                        - if has_events { inp.self_consumption } else { 0 })
                    .max(0) as u64,
                    inp.count,
                )
            })
            .product::<u64>() as f64
    }

    /// Returns the number of reactants contained in the state.
    #[allow(unused)]
    pub fn len(&self) -> usize {
        self.state.len()
    }
}

impl Index<usize> for StateData {
    type Output = ComponentData;

    fn index(&self, index: usize) -> &Self::Output {
        &self.state[index]
    }
}

impl IndexMut<usize> for StateData {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.state[index]
    }
}
