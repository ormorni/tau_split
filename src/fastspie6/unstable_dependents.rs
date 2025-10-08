use std::ops::Index;

use super::f_reaction::FReaction;

/// A struct that tracks how many unstable reactions depend on each reactant.
#[derive(Clone, Debug, Hash, PartialEq, Eq, Default)]
pub struct UnstableDependents {
    dependent_count: Vec<usize>,
}

impl UnstableDependents {
    pub fn empty(reactant_count: usize) -> UnstableDependents {
        UnstableDependents {
            dependent_count: vec![0; reactant_count],
        }
    }

    /// Adds an unstable reaction to the tracker.
    pub fn add_unstable(&mut self, reaction: &FReaction) {
        for inp in &reaction.inputs {
            self.dependent_count[inp.index] += 1;
        }
    }

    /// Removes an unstable reaction from the tracker.
    pub fn remove_unstable(&mut self, reaction: &FReaction) {
        for inp in &reaction.inputs {
            self.dependent_count[inp.index] -= 1;
        }
    }

    /// Checks if the reaction has any unstable reactions depending on its output.
    pub fn has_dependents(&self, reaction: &FReaction) -> bool {
        reaction
            .stoichiometry
            .iter()
            .any(|&(comp, _)| self.dependent_count[comp] > 0)
    }
}

impl Index<usize> for UnstableDependents {
    type Output = usize;

    fn index(&self, index: usize) -> &Self::Output {
        &self.dependent_count[index]
    }
}
