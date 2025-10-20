use crate::utils::binomial_05;

use derive_new::new;
use rand::Rng;
use rand_distr::{Binomial, Exp, Poisson};
use std::fmt::Debug;

use super::f_reaction::FReaction;

/// Data on the number of events in a reaction spanning some period of time.
#[derive(Clone, Copy)]
pub(super) struct ReactionData {
    /// The index of the reaction covered by the reaction data.
    pub reaction: usize,
    /// The timespan covered by the reaction data.
    pub time: f64,
    /// The number of events occurring in the timespan.
    pub events: u64,
    /// The lower cutoff for the propensity where the given number of events occurs.
    pub low: f64,
    /// The upper cutoff for the propensity where the given number of events occurs.
    pub high: f64,
}

impl Debug for ReactionData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!(
            "RData({}, ({}, {})=>{})",
            self.reaction, self.low, self.high, self.events
        ))
    }
}

impl ReactionData {
    pub fn new(reaction: usize, time: f64, events: u64, low: f64, high: f64) -> ReactionData {
        ReactionData {
            reaction,
            time,
            events,
            low,
            high,
        }
    }

    /// Samples a new ReactionData.
    pub fn sample(
        product: f64,
        reaction_idx: usize,
        reaction: &FReaction,
        time: f64,
        rng: &mut impl Rng,
    ) -> ReactionData {
        let events = sample_poisson(product * time * reaction.rate, rng);
        let low = product * max_sample(events, rng);
        let high = product + sample_exp(reaction.rate * time, rng);

        ReactionData::new(reaction_idx, time, events, low, high)
    }

    /// Resamples the current number of events conditioned on the previous data.
    pub fn resample(&mut self, product: f64, reaction: &FReaction, rng: &mut impl Rng) {
        if product < self.low {
            let rem_events = sample_binomial(self.events - 1, product / self.low, rng);
            let low = product * max_sample(rem_events, rng);
            let high = product
                + (self.low - product) * (1. - max_sample(self.events - rem_events - 1, rng));

            self.low = low;
            self.events = rem_events;
            self.high = high;
        } else if product >= self.high {
            let extra_events =
                sample_poisson(reaction.rate * self.time * (product - self.high), rng);
            let low = self.high + max_sample(extra_events, rng) * (product - self.high);
            let high = product + sample_exp(reaction.rate * self.time, rng);

            self.events += extra_events + 1;
            self.low = low;
            self.high = high;
        }
    }

    /// Creates a StableReactionData from the normal ReactionData.
    pub fn stabilize(self) -> StableReactionData {
        // println!("{}", self.events);  // Used for figure generation.
        StableReactionData::new(
            self.reaction,
            self.time,
            self.events,
            self.low,
            true,
            self.high,
            true,
        )
    }
}

#[derive(Clone, Copy, Debug, new)]
pub(super) struct StableReactionData {
    /// The index of the reaction covered by the reaction data.
    pub reaction: usize,
    /// The timespan covered by the reaction data.
    pub time: f64,
    /// The number of events occurring in the timespan.
    pub events: u64,
    /// The lower cutoff for the propensity where the given number of events occurs.
    pub low: f64,
    /// Whether the low point is assigned to this inactive data.
    pub has_low: bool,
    /// The upper cutoff for the propensity where the given number of events occurs.
    pub high: f64,
    /// Whether the high point is assigned to this inactive data.
    pub has_high: bool,
}

impl StableReactionData {
    /// Samples the higher bound if the reaction data doesn't currently have one.
    pub fn sample_high(&mut self, reaction: &FReaction, rng: &mut impl Rng) -> f64 {
        if !self.has_high {
            self.high = self.high + sample_exp(reaction.rate * self.time, rng);
            self.has_high = true;
        }
        self.high
    }

    /// Samples the lower bound if the reaction data doesn't currently have one.
    pub fn sample_low(&mut self, _reaction: &FReaction, rng: &mut impl Rng) -> f64 {
        if !self.has_low {
            self.low = self.low * max_sample(self.events, rng);
            self.has_low = true;
        }
        self.low
    }

    /// Reactivates the InactiveReactionData, creating a valid ReactionData and sampling all
    /// the variables we attempted not to sample.
    pub fn destabilize(mut self, reaction: &FReaction, rng: &mut impl Rng) -> ReactionData {
        ReactionData::new(
            self.reaction,
            self.time,
            self.events,
            self.sample_low(reaction, rng),
            self.sample_high(reaction, rng),
        )
    }
}

pub trait TauData {
    /// Returns the number of events in the given reaction in the spanned time period.
    fn event_count(&self) -> u64;
    /// Splits the reaction data to two objects representing the reaction data over two halves of the time segment.
    fn split(&mut self, reaction: &FReaction, rng: &mut impl Rng) -> Self;
    /// Returns the index of the reaction.
    fn index(&self) -> usize;

    fn has_events(&self) -> bool {
        self.event_count() != 0
    }
}

impl TauData for ReactionData {
    fn event_count(&self) -> u64 {
        self.events
    }

    fn split(&mut self, reaction: &FReaction, rng: &mut impl Rng) -> Self {
        // print!("Splitting {self:?} to ");
        self.time /= 2.;
        let mut res = self.clone();
        let events = self.events;

        res.events = binomial_05(self.events, rng);
        self.events = events - res.events;

        // Sampling the lower bound.
        if events > 0 {
            if rng.random_bool(res.events as f64 / events as f64) {
                self.low *= max_sample(self.events, rng);
            } else {
                res.low *= max_sample(res.events, rng);
            }
        }
        // Sampling the upper bound.
        if rng.random_bool(0.5) {
            self.high += sample_exp(reaction.rate * self.time, rng);
        } else {
            res.high += sample_exp(reaction.rate * res.time, rng);
        }

        res
    }

    fn index(&self) -> usize {
        self.reaction
    }
}

impl TauData for StableReactionData {
    fn event_count(&self) -> u64 {
        self.events
    }

    fn split(&mut self, _reaction: &FReaction, rng: &mut impl Rng) -> Self {
        self.time /= 2.;
        let mut res = self.clone();

        let events = self.events;
        res.events = binomial_05(self.events, rng);
        self.events = events - res.events;

        // Sampling the lower bound.
        if events > 0 && self.has_low {
            if rng.random_bool(res.events as f64 / events as f64) {
                self.has_low = false;
            } else {
                res.has_low = false;
            }
        }
        // Sampling the upper bound.
        if self.has_high {
            if rng.random_bool(0.5) {
                self.has_high = false;
            } else {
                res.has_high = false;
            }
        }
        res
    }

    fn index(&self) -> usize {
        self.reaction
    }
}

pub fn sample_binomial(n: u64, p: f64, rng: &mut impl Rng) -> u64 {
    rng.sample(
        Binomial::new(n, p)
            .unwrap_or_else(|err| panic!("sample_binomial({n}, {p}) failed with err {err}")),
    )
}
pub fn sample_exp(rate: f64, rng: &mut impl Rng) -> f64 {
    rng.sample(
        Exp::new(rate).unwrap_or_else(|err| panic!("sample_exp({rate}) failed with err {err}")),
    )
}
pub fn sample_poisson(rate: f64, rng: &mut impl Rng) -> u64 {
    if rate == 0. {
        return 0;
    }
    rng.sample(Poisson::new(rate).unwrap_or_else(|err| {
        panic!("Failed to sample a Poisson variable with rate {rate}: {err:?}")
    })) as u64
}
/// Samples the maximal sample among `samples` uniformly distributed samples.
pub fn max_sample(n: u64, rng: &mut impl Rng) -> f64 {
    if n == 0 {
        0.
    } else {
        rng.random::<f64>().powf((n as f64).recip())
    }
}
