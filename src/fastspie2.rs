use std::fmt::Debug;
use std::mem::swap;

use derive_new::new;
use itertools::izip;
use rand::Rng;
use rand_distr::{Binomial, Poisson};

use crate::reaction::Reaction;
use crate::utils::binomial_05;

#[derive(Clone, Copy, new, Default, PartialEq, Eq, PartialOrd, Ord)]
struct ProdEvents {
    product: u64,
    events: u64,
}

impl ProdEvents {
    pub fn zero() -> ProdEvents {
        ProdEvents::new(0, 0)
    }
}

impl Debug for ProdEvents {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} => {}", self.product, self.events)
    }
}

/// Samples the number of events given the reaction product
/// given lower and upper bounds already sampled.
fn sample_events(
    low: Option<ProdEvents>,
    high: Option<ProdEvents>,
    product: u64,
    rate: f64,
    rng: &mut impl Rng,
) -> u64 {
    let low = low.unwrap_or(ProdEvents::zero());
    if product == low.product {
        return low.events;
    }
    if let Some(high) = high {
        if high.events == low.events {
            return low.events;
        }

        low.events
            + rng.sample(
                Binomial::new(
                    high.events - low.events,
                    (product - low.product) as f64 / (high.product - low.product) as f64,
                )
                .unwrap(),
            )
    } else {
        let pois = Poisson::new((product - low.product) as f64 * rate)
            .unwrap_or_else(|e| panic!("{product} {low:?} {high:?} {rate} => {e}"));
        low.events + rng.sample(pois) as u64
    }
}

/// Data on the number of events in a reaction spanning some period of time.
#[derive(Clone, Copy, new)]
struct ReactionData<'t> {
    reaction: &'t Reaction,
    data: [ProdEvents; 3],
}

impl<'t> Debug for ReactionData<'t> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.data.fmt(f)
    }
}

impl<'t> ReactionData<'t> {
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

    /// Finds the highest product-events pair with product lesser or
    /// equal to the given input product.
    fn find_below(&self, product: u64) -> Option<ProdEvents> {
        self.data
            .iter()
            .rev()
            .find(|p| p.product <= product)
            .copied()
    }
    /// Finds the lowest product-events pair with product greater or
    /// equal to the given input product.
    fn find_above(&self, product: u64) -> Option<ProdEvents> {
        self.data.iter().find(|p| p.product >= product).copied()
    }

    fn split(&mut self, rng: &mut impl Rng) -> ReactionData<'t> {
        let mut res = self.clone();
        res.data[0].events = binomial_05(self.data[0].events, rng);
        res.data[1].events =
            binomial_05(self.data[1].events - self.data[0].events, rng) + res.data[0].events;
        res.data[2].events =
            binomial_05(self.data[2].events - self.data[1].events, rng) + res.data[1].events;
        self.data[0].events -= res.data[0].events;
        self.data[1].events -= res.data[1].events;
        self.data[2].events -= res.data[2].events;

        res
    }
}

/// The FastGillespie2 is similar to FastGillespie in that it does tau-leaping steps,
/// checks for correctness, and if they are not correct, splits the time-segment in two,
/// then simulates both time segments.
///
/// However, unlike FastGillespie, it uses the observation that samples can be "forgotten".
/// When splitting, the only actual decision made was to split, so as long as we sample
/// a flow that is consistent with that decision, we are fine.
pub struct FastGillespie2 {
    pub state: Vec<i64>,
    pub reactions: Vec<Reaction>,
}

impl FastGillespie2 {
    pub fn new(state: Vec<i64>, reactions: Vec<Reaction>) -> FastGillespie2 {
        FastGillespie2 { state, reactions }
    }

    pub fn advance(&mut self, time: f64, rng: &mut impl Rng) {
        // Taking the reactions so that I can create ReactionData referencing them easily.
        let mut reactions = Vec::new();
        swap(&mut self.reactions, &mut reactions);

        let mut buffer = Vec::new();
        let mut times = vec![time];

        for reaction in &reactions {
            let input_product = reaction.input_product(&self.state);
            let event_count = sample_events(
                Some(ProdEvents::zero()),
                None,
                input_product,
                reaction.rate * time,
                rng,
            );
            buffer.push(ReactionData::new(
                reaction,
                [ProdEvents::new(input_product, event_count); 3],
            ));
        }
        self.recursion(&mut times, &mut buffer, &reactions, rng);
        assert!(times.is_empty());
        assert!(buffer.is_empty());
        // Returning the reactions.
        swap(&mut self.reactions, &mut reactions);
    }

    fn recursion<'t>(
        &mut self,
        times: &mut Vec<f64>,
        buffer: &mut Vec<ReactionData<'t>>,
        reactions: &'t [Reaction],
        rng: &mut impl Rng,
    ) {
        let time = times.pop().unwrap();
        // We get reaction data from previous reactions, which we now have to fix.
        // We first allocate a new set of reactions to store the old data.
        let start_idx = buffer.len() - reactions.len();
        let end_idx = buffer.len();
        buffer.extend_from_within(start_idx..end_idx);

        let (buf, old_data) = buffer.split_at_mut(end_idx);
        let (_, new_data) = buf.split_at_mut(start_idx);

        // Getting the updated input count and products.
        for (rdata, old_rdata) in izip!(new_data.iter_mut(), old_data.iter_mut()) {
            let input_product = rdata.reaction.input_product(&self.state);
            let low = old_rdata.find_below(input_product);
            let high = old_rdata.find_above(input_product);

            let event_count = sample_events(
                low,
                high,
                input_product,
                old_rdata.reaction.rate * time,
                rng,
            );

            rdata.data[1] = ProdEvents::new(input_product, event_count);
        }

        // Computing the lower bound.
        for rdata in new_data.iter() {
            rdata.apply_negative(&mut self.state, rdata.data[1].events as i64);
        }

        for (rdata, old_rdata) in izip!(new_data.iter_mut(), old_data.iter_mut()) {
            let input_product = rdata.reaction.input_product(&self.state);
            let low = old_rdata.find_below(input_product);
            let high = old_rdata
                .find_above(input_product)
                .or_else(|| Some(rdata.data[1]))
                .min(Some(rdata.data[1]));
            let event_count = sample_events(
                low,
                high,
                input_product,
                old_rdata.reaction.rate * time,
                rng,
            );

            rdata.data[0] = ProdEvents::new(input_product, event_count);
        }

        for rdata in new_data.iter() {
            rdata.apply_negative(&mut self.state, -(rdata.data[1].events as i64));
            rdata.apply_positive(&mut self.state, rdata.data[1].events as i64);
        }

        for (new, old) in izip!(new_data.iter_mut(), old_data.iter_mut()) {
            let input_product = new.reaction.input_product(&self.state);
            let low = old
                .find_below(input_product)
                .or_else(|| Some(new.data[1]))
                .max(Some(new.data[1]));
            let high = old.find_above(input_product);

            let event_count =
                sample_events(low, high, input_product, old.reaction.rate * time, rng);

            new.data[2] = ProdEvents::new(input_product, event_count);
        }

        for rdata in new_data.iter() {
            rdata.apply_positive(&mut self.state, -(rdata.data[1].events as i64));
        }

        // We have finished sampling the new events, and can pop off the old data.
        buffer.truncate(end_idx);
        let rdata = &buffer[start_idx..];
        // Checking if we are done.
        let total_reactions = rdata.iter().map(|d| d.data[1].events).sum::<u64>();
        let coherent = rdata.iter().all(|d| d.data[0].events == d.data[2].events);
        if total_reactions <= 1 || coherent {
            rdata
                .iter()
                .for_each(|r| r.apply_all(&mut self.state, r.data[1].events as i64));
            buffer.truncate(start_idx);
        } else {
            // If not, we have to split the data.
            for i in start_idx..end_idx {
                let split = buffer[i].split(rng);
                buffer.push(split);
            }
            times.push(time / 2.);
            times.push(time / 2.);
            self.recursion(times, buffer, reactions, rng);
            self.recursion(times, buffer, reactions, rng);
        }
    }
}
