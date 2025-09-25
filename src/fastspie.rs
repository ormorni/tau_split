use itertools::{izip, Itertools};
use rand::Rng;
use rand_distr::{Binomial, Poisson};

use crate::{
    reaction::{binomial, Reaction},
    utils::binomial_05,
};

/// Data on the number of events in a reaction spanning some period of time.
#[derive(Clone, Debug)]
struct ReactionData<'t> {
    reaction: &'t Reaction,
    event_counts: Vec<(u64, u64)>,
}

impl<'t> ReactionData<'t> {
    fn input_product(&self, reactants: &[i64]) -> u64 {
        self.reaction
            .inputs
            .iter()
            .map(|&(reactant, count)| binomial(reactants[reactant].max(0) as u64, count))
            .product()
    }

    fn new(reaction: &'t Reaction) -> ReactionData<'t> {
        ReactionData {
            reaction,
            event_counts: vec![(0, 0)],
        }
    }

    /// Splits the reaction data into two reaction data spanning the two halves of the time segment
    fn split(mut self, rng: &mut impl Rng) -> (ReactionData<'t>, ReactionData<'t>) {
        debug_assert!(!self.event_counts.is_empty());
        debug_assert_eq!(self.event_counts[0], (0, 0));
        debug_assert!(self.event_counts.is_sorted());

        let mut right = self.clone();
        let mut prev_events = 0;

        for i in 1..self.event_counts.len() {
            let rem_events = self.event_counts[i].1 - prev_events;
            let left_rem_events = binomial_05(rem_events as u64, rng);
            self.event_counts[i].1 = self.event_counts[i - 1].1 + left_rem_events;
            right.event_counts[i].1 = right.event_counts[i - 1].1 + rem_events - left_rem_events;

            prev_events += rem_events;
        }

        (self, right)
    }

    /// Returns the number of events the reaction happens, conditioned on all previous samples.
    fn event_count(&mut self, reactants: &[i64], time: f64, rng: &mut impl Rng) -> u64 {
        let product = self.input_product(reactants);
        let Some(high) = self
            .event_counts
            .iter()
            .position(|(prod, _)| *prod >= product)
        else {
            // We haven't yet sampled an input product higher than the given product.
            // The difference is Poisson distributed.
            let &(last_prod, last_events) = self.event_counts.last().unwrap();
            let new_events = rng.sample(
                Poisson::new(self.reaction.rate * (product - last_prod) as f64 * time).unwrap(),
            ) as u64;
            self.event_counts.push((product, last_events + new_events));
            return last_events + new_events;
        };

        let (high_product, high_events) = self.event_counts[high];

        if product == high_product {
            // We have already sampled this event count.
            return high_events;
        }

        // The first index is always (0, 0), and the product is not 0, so `high` > 0.
        let (low_product, low_events) = self.event_counts[high - 1];
        if low_events == high_events {
            return low_events;
        }
        let new_events = rng.sample(
            Binomial::new(
                high_events - low_events,
                (product - low_product) as f64 / (high_product - low_product) as f64,
            )
            .unwrap(),
        ) + low_events;

        self.event_counts.insert(high, (product, new_events));
        new_events
    }

    fn apply_all(&self, reactants: &mut [i64], event_count: i64) {
        for &(reactant, change) in &self.reaction.stoichiometry {
            reactants[reactant] += change * event_count;
        }
    }

    /// Applies only the negative parts of the product to the reactants.
    fn apply_negative(&self, reactants: &mut [i64], event_count: i64) {
        for &(reactant, change) in &self.reaction.stoichiometry {
            if change < 0 {
                reactants[reactant] += change * event_count;
            }
        }
    }
    /// Applies only the negative parts of the product to the reactants.
    fn apply_positive(&self, reactants: &mut [i64], event_count: i64) {
        for &(reactant, change) in &self.reaction.stoichiometry {
            if change > 0 {
                reactants[reactant] += change * event_count;
            }
        }
    }
}

pub struct FastGillespie {
    pub state: Vec<i64>,
    pub reactions: Vec<Reaction>,
    buffer: Vec<u64>,
}

impl FastGillespie {
    pub fn new(state: Vec<i64>, reactions: Vec<Reaction>) -> FastGillespie {
        FastGillespie {
            state,
            reactions,
            buffer: Vec::new(),
        }
    }

    pub fn advance(&mut self, time: f64, rng: &mut impl Rng) {
        let reactions = self.reactions.clone();
        let reaction_data = reactions.iter().map(ReactionData::new).collect_vec();
        self.recursion(reaction_data, time, rng);
    }

    fn recursion(&mut self, mut reaction_data: Vec<ReactionData>, time: f64, rng: &mut impl Rng) {
        if self.buffer.len() < 3 * reaction_data.len() {
            self.buffer.resize(3 * reaction_data.len(), 0);
        }

        let (event_counts, rest) = self.buffer.split_at_mut(reaction_data.len());
        let (lower_bound, upper_bound) = rest.split_at_mut(reaction_data.len());

        for (rdata, tar) in izip!(&mut reaction_data, event_counts.iter_mut()) {
            *tar = rdata.event_count(&self.state, time, rng);
        }

        // Computing the positive error.
        for (rdata, &mut event_count) in izip!(&reaction_data, event_counts.iter_mut()) {
            rdata.apply_positive(&mut self.state, event_count as i64);
        }

        for (rdata, tar) in izip!(&mut reaction_data, upper_bound.iter_mut()) {
            *tar = rdata.event_count(&self.state, time, rng);
        }

        // Computing the negative error.
        for (rdata, &mut event_count) in izip!(&reaction_data, event_counts.iter_mut()) {
            rdata.apply_positive(&mut self.state, -(event_count as i64));
            rdata.apply_negative(&mut self.state, event_count as i64);
        }

        for (rdata, tar) in izip!(&mut reaction_data, lower_bound.iter_mut()) {
            *tar = rdata.event_count(&mut self.state, time, rng);
        }
        for (rdata, &mut event_count) in izip!(&reaction_data, event_counts.iter_mut()) {
            rdata.apply_negative(&mut self.state, -(event_count as i64));
        }

        if upper_bound.iter().sum::<u64>() == 1 || upper_bound == lower_bound {
            // If the upper bound and lower bound are equal, then we can apply the reactions and finish.
            for (rdata, &mut event_count) in izip!(&reaction_data, event_counts) {
                rdata.apply_all(&mut self.state, event_count as i64);
            }
            return;
        }

        // Otherwise, we have to split!
        let mut left_events = Vec::with_capacity(reaction_data.len());
        let mut right_events = Vec::with_capacity(reaction_data.len());

        for data in reaction_data {
            let (ldata, rdata) = data.split(rng);
            left_events.push(ldata);
            right_events.push(rdata);
        }

        self.recursion(left_events, time / 2., rng);
        self.recursion(right_events, time / 2., rng);
    }
}
