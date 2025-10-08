use std::ops::{Index, IndexMut};

use derive_new::new;
use tinyvec::ArrayVec;

use crate::reaction::{binomial, Reaction};

use super::reaction_data::{ReactionData, StableReactionData};

/// A struct holding the input of a reaction.
#[derive(Default, Clone, Copy, PartialEq, Eq, Debug, new)]
pub(crate) struct Input {
    /// The index of the reactant.
    pub index: usize,
    /// The number of times the reactant appears in the reaction.
    pub count: u64,
    /// If the reaction produces more of the input, 0.
    /// Otherwise, the change in the input done by the reaction.
    pub self_consumption: i64,
}

/// A struct describing a single chemical reaction.
///
/// This struct includes a self-interaction term to better handle self-cancellations.
#[derive(Clone, Debug, new)]
pub(super) struct FReaction {
    /// The inputs to the reaction.
    /// An array of tuples of (reactant_idx, reactant_count, )
    pub inputs: ArrayVec<[Input; 2]>,
    /// The change to the reaction state for every firing of the reaction.
    pub stoichiometry: ArrayVec<[(usize, i64); 4]>,
    /// The change to the reaction state for every firing of the reaction.
    pub positive_stoichiometry: ArrayVec<[(usize, i64); 4]>,
    /// The change to the reaction state for every firing of the reaction.
    pub negative_stoichiometry: ArrayVec<[(usize, i64); 4]>,
    /// The rate constant of the reaction.
    pub rate: f64,
}

impl FReaction {
    /// Computes the `input_product` of the reaction, or the number of combinations of input molecules.
    pub fn input_product(&self, reactants: &[i64]) -> u64 {
        self.inputs
            .iter()
            .map(|inp| binomial(reactants[inp.index].max(0) as u64, inp.count))
            .product()
    }
}

impl From<Reaction> for FReaction {
    fn from(value: Reaction) -> Self {
        // Finding for every input the number of times it appears in the stoichiometry vector,
        // and 0 otherwise.
        // Since we are only interested in self-consumption, we take the minimum with 0.
        FReaction {
            inputs: value
                .inputs
                .into_iter()
                .map(|(inp_component, count)| {
                    Input::new(
                        inp_component,
                        count,
                        value
                            .stoichiometry
                            .iter()
                            .find_map(|&(component, count)| {
                                if component == inp_component {
                                    Some(count)
                                } else {
                                    None
                                }
                            })
                            .unwrap_or(0)
                            .min(0),
                    )
                })
                .collect(),
            stoichiometry: value.stoichiometry,
            positive_stoichiometry: value.positive_stoichiometry,
            negative_stoichiometry: value.negative_stoichiometry,
            rate: value.rate,
        }
    }
}

impl<'t> Index<&'t ReactionData> for [FReaction] {
    type Output = FReaction;

    fn index(&self, index: &ReactionData) -> &Self::Output {
        self.index(index.reaction)
    }
}
impl<'t> IndexMut<&'t ReactionData> for [FReaction] {
    fn index_mut(&mut self, index: &'t ReactionData) -> &mut Self::Output {
        self.index_mut(index.reaction)
    }
}

impl<'t> Index<&'t StableReactionData> for [FReaction] {
    type Output = FReaction;

    fn index(&self, index: &StableReactionData) -> &Self::Output {
        self.index(index.reaction)
    }
}
impl<'t> IndexMut<&'t StableReactionData> for [FReaction] {
    fn index_mut(&mut self, index: &'t StableReactionData) -> &mut Self::Output {
        self.index_mut(index.reaction)
    }
}
