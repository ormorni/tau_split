use itertools::chain;
use tinyvec::ArrayVec;

pub const MAX_INPUTS: usize = 4;
pub const MAX_STOI: usize = 4;

/// A struct describing a single chemical reaction.
#[derive(Clone, Debug)]
pub struct Reaction {
    /// The inputs to the reaction.
    pub inputs: ArrayVec<[(usize, u64); MAX_INPUTS]>,
    /// The change to the reaction state for every firing of the reaction.
    pub stoichiometry: ArrayVec<[(usize, i64); MAX_STOI]>,
    /// The change to the reaction state for every firing of the reaction.
    pub(crate) positive_stoichiometry: ArrayVec<[(usize, i64); MAX_STOI]>,
    /// The change to the reaction state for every firing of the reaction.
    pub(crate) negative_stoichiometry: ArrayVec<[(usize, i64); MAX_STOI]>,
    /// The rate constant of the reaction.
    pub rate: f64,
}

impl Reaction {
    pub fn new(
        inputs: ArrayVec<[(usize, u64); MAX_INPUTS]>,
        stoichiometry: ArrayVec<[(usize, i64); MAX_STOI]>,
        rate: f64,
    ) -> Reaction {
        let positive_stoichiometry = stoichiometry
            .iter()
            .filter(|(_, diff)| *diff > 0)
            .copied()
            .collect();
        let negative_stoichiometry = stoichiometry
            .iter()
            .filter(|(_, diff)| *diff < 0)
            .copied()
            .collect();

        Reaction {
            inputs,
            stoichiometry,
            positive_stoichiometry,
            negative_stoichiometry,
            rate,
        }
    }
}

/// Computes n choose k, of the number of subsets of size k of a set of size n.
pub fn binomial(n: u64, k: u64) -> u64 {
    match k {
        0 => 1,
        1 => n,
        2 => (n * n - n) / 2,
        3 => (n * (n - 1) * (n - 2)) / 6,
        k => {
            let mut res = 1;
            for i in 0..k {
                res = res * (n - i) / (i + 1);
            }
            res
        }
    }
}

impl Reaction {
    /// Computes the `input_product` of the reaction, or the number of combinations of input molecules.
    pub fn input_product(&self, reactants: &[i64]) -> u64 {
        self.inputs
            .iter()
            .map(|&(reactant, count)| binomial(reactants[reactant].max(0) as u64, count))
            .product()
    }

    /// Computes the rate at which reaction events occur.
    pub fn rate(&self, reactants: &[i64]) -> f64 {
        self.input_product(reactants) as f64 * self.rate
    }

    /// Applies the reaction to the reactants.
    pub fn apply(&self, reactants: &mut [i64], count: i64) {
        for &(reactant, change) in &self.stoichiometry {
            reactants[reactant] += count * change;
        }
    }

    pub fn all_reactants<'t>(&'t self) -> impl Iterator<Item = usize> + 't {
        chain!(
            self.inputs.iter().map(|(r, _)| *r),
            self.stoichiometry.iter().map(|(r, _)| *r)
        )
    }
}
