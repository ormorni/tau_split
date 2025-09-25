use crate::{reaction::Reaction, utils::binomial_05};

use derive_new::new;
use rand::Rng;
use rand_distr::{Binomial, Exp, Poisson};
use std::fmt::Debug;

/// Data on the number of events in a reaction spanning some period of time.
#[derive(Clone, Copy, new)]
pub struct ReactionData {
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
    /// Samples a new ReactionData.
    pub fn sample(
        product: f64,
        reaction_idx: usize,
        reaction: &Reaction,
        time: f64,
        rng: &mut impl Rng,
    ) -> ReactionData {
        let events = sample_poisson(product * time * reaction.rate, rng);
        let low = product * max_sample(events, rng);
        let high = product + sample_exp(reaction.rate * time, rng);

        ReactionData::new(reaction_idx, time, events, low, high)
    }

    /// Splits the reaction data over some timespan to two halves.
    pub fn split(&mut self, reaction: &Reaction, rng: &mut impl Rng) -> ReactionData {
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
        // println!("{self:?}, {res:?}");

        res
    }

    /// Resamples the current number of events conditioned on the previous data.
    pub fn resample(&mut self, product: f64, reaction: &Reaction, rng: &mut impl Rng) {
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
