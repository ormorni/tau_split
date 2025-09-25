use crate::reaction::Reaction;

use super::ReactionData;

pub struct StateData {
    /// The current reaction state.
    pub state: Vec<i64>,
    /// The current reactant lower bound.
    pub lower_bound: Vec<i64>,
    /// The current reactant upper bound.
    pub upper_bound: Vec<i64>,
}

impl StateData {
    pub fn new(state: &[i64]) -> StateData {
        StateData {
            state: state.to_owned(),
            lower_bound: state.to_owned(),
            upper_bound: state.to_owned(),
        }
    }

    /// Removes the effect of a ReactionData object from the error bounds.
    pub fn change_bounds(&mut self, event_count: i64, reaction: &Reaction) {
        if event_count != 0 {
            Self::apply_negative(&mut self.lower_bound, event_count, reaction);
            Self::apply_positive(&mut self.upper_bound, event_count, reaction);
        }
    }

    /// Removes the effect of a ReactionData object from the error bounds.
    pub fn remove_bounds(&mut self, rdata: &ReactionData, reaction: &Reaction) {
        self.change_bounds(-rdata.event_count(), reaction);
    }

    /// Add the effect of a ReactionData object to the error bounds.
    pub fn add_bounds(&mut self, rdata: &ReactionData, reaction: &Reaction) {
        self.change_bounds(rdata.event_count(), reaction);
    }

    pub fn apply(&mut self, rdata: &ReactionData, reaction: &Reaction) {
        Self::apply_all(&mut self.lower_bound, rdata.event_count(), reaction);
        Self::apply_all(&mut self.state, rdata.event_count(), reaction);
        Self::apply_all(&mut self.upper_bound, rdata.event_count(), reaction);
    }

    pub fn apply_all(reactants: &mut [i64], event_count: i64, reaction: &Reaction) {
        for &(reactant, change) in &reaction.stoichiometry {
            reactants[reactant] += change * event_count;
        }
    }

    /// Applies only the negative parts of the product to the reactants.
    pub fn apply_negative(reactants: &mut [i64], event_count: i64, reaction: &Reaction) {
        for &(reactant, change) in &reaction.negative_stoichiometry {
            reactants[reactant] += change * event_count;
        }
    }
    /// Applies only the negative parts of the product to the reactants.
    pub fn apply_positive(reactants: &mut [i64], event_count: i64, reaction: &Reaction) {
        for &(reactant, change) in &reaction.positive_stoichiometry {
            reactants[reactant] += change * event_count;
        }
    }
}
