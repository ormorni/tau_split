use itertools::izip;
use rand::Rng;
use rand_distr::{Distribution, Exp};

use crate::{reaction::Reaction, SimulationAlg};

/// A binary-indexed-tree-like data structure for storing reaction propensities
/// and sampling the next reaction using them.
#[derive(Debug, Clone)]
struct ChoiceTree {
    data: Vec<f64>,
    alloc_size: usize,
    size: usize,
}

impl ChoiceTree {
    pub fn empty(size: usize) -> ChoiceTree {
        let data = vec![0.; size.next_power_of_two() * 2 - 1];
        ChoiceTree {
            data,
            size,
            alloc_size: size.next_power_of_two(),
        }
    }

    pub fn len(&self) -> usize {
        self.size
    }

    pub fn update(&mut self, idx: usize, value: f64) {
        debug_assert!(
            value >= 0.,
            "A reaction propensity cannot be negative: {value}"
        );
        debug_assert!(
            (0..self.len()).contains(&idx),
            "Attempted to set an empty index: {idx}"
        );
        let old_val = self.data[self.alloc_size + idx - 1];
        let mut mapped_index = idx + self.alloc_size;
        while mapped_index > 0 {
            self.data[mapped_index - 1] += value - old_val;
            mapped_index /= 2;
        }
    }

    /// Returns the total propensity of all the reactions.
    pub fn total(&self) -> f64 {
        self.data[0]
    }
}

impl Distribution<usize> for ChoiceTree {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> usize {
        let mut idx = 1;
        let mut choice = rng.random::<f64>() * self.data[0];
        while idx * 2 < self.data.len() {
            if choice < self.data[2 * idx] {
                idx = 2 * idx + 1;
            } else {
                choice -= self.data[2 * idx];
                idx = 2 * idx;
            }
        }
        idx - self.alloc_size
    }
}

pub struct Gillespie {
    /// The chemical equations going on.
    eqs: Vec<Reaction>,
    /// A map from a reactant to all the reactions it is involved in.
    reaction_updates: Vec<Vec<usize>>,
    /// The current state of the system.
    pub state: Vec<i64>,
    /// A data structure allowing quick sampling of the next reaction to fire.
    tree: ChoiceTree,
    /// The total number of reactions simulated by the algorithm.
    total_reactions: u64,
}

impl Gillespie {
    /// Advances the state and returns the amount of time that has passed.
    pub fn sample_reaction(&mut self, max_time: f64, rng: &mut impl Rng) -> f64 {
        let time = rng.sample(Exp::new(self.tree.total()).unwrap());
        if time > max_time {
            // If the time until the next reaction is greater than the remaining time in the simulation,
            // no reaction occurs.
            return max_time;
        }
        if self.tree.total() <= 1e-9 {
            return f64::MAX;
        }
        let reaction_idx = self.tree.sample(rng);

        self.eqs[reaction_idx].apply(&mut self.state, 1);
        for update_idx in &self.reaction_updates[reaction_idx] {
            self.tree
                .update(*update_idx, self.eqs[*update_idx].rate(&self.state));
        }
        self.total_reactions += 1;
        time
    }
}

impl SimulationAlg for Gillespie {
    /// Initializes a new Gillespie
    fn new(initial_state: Vec<i64>, eqs: Vec<Reaction>, _reactant_names: Vec<String>) -> Gillespie {
        // Computing the graph of which reaction updates which other reactions.
        // We already have the reaction -> reactant edges, and now need the reactant -> reaction edges,
        // and then we take the product.
        let reactant_count = eqs
            .iter()
            .flat_map(|eq| eq.all_reactants())
            .max()
            .unwrap_or(0);

        let mut reactant_eqs = vec![Vec::new(); reactant_count + 1];
        for (idx, eq) in eqs.iter().enumerate() {
            for (reactant, _) in &eq.inputs {
                reactant_eqs[*reactant].push(idx);
            }
        }

        let mut reaction_updates: Vec<Vec<usize>> = vec![Vec::new(); eqs.len()];

        for (eq, updates) in izip!(&eqs, &mut reaction_updates) {
            for (react, _) in &eq.stoichiometry {
                updates.extend_from_slice(&reactant_eqs[*react]);
            }
        }
        for updates in &mut reaction_updates {
            updates.sort();
            updates.dedup();
        }

        // Initializing the propensity in the propensity tree.
        let mut tree = ChoiceTree::empty(eqs.len());

        for (idx, eq) in eqs.iter().enumerate() {
            tree.update(idx, eq.input_product(&initial_state) as f64 * eq.rate);
        }

        Gillespie {
            eqs,
            state: initial_state,
            reaction_updates,
            tree,
            total_reactions: 0,
        }
    }

    fn advance(&mut self, mut time: f64, rng: &mut impl Rng) {
        while time > 0. {
            time -= self.sample_reaction(time, rng);
        }
    }

    fn state(&self) -> &[i64] {
        &self.state
    }

    fn total_reactions(&self) -> u64 {
        self.total_reactions
    }
}
