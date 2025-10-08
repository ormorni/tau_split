use crate::{reaction::Reaction, utils::binomial_05};

use super::prod_events::ProdEvents;
use derive_new::new;
use rand::Rng;
use rand_distr::{Binomial, Poisson};
use std::fmt::Debug;

/// Data on the number of events in a reaction spanning some period of time.
#[derive(Clone, Copy, new)]
pub struct ReactionData {
    pub reaction: usize,
    pub data: [ProdEvents; 3],
}

impl Debug for ReactionData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("RData(")?;
        self.reaction.fmt(f)?;
        f.write_str(",")?;
        self.data.fmt(f)?;
        f.write_str(")")
    }
}

impl ReactionData {
    pub fn event_count(&self) -> i64 {
        self.data[1].events as i64
    }

    pub fn split(&mut self, rng: &mut impl Rng) -> ReactionData {
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

    /// Resamples the current number of events conditioned on the previous data.
    ///
    pub fn resample_mid(
        &mut self,
        product: u64,
        reaction: &Reaction,
        time: f64,
        rng: &mut impl Rng,
    ) {
        self.validate();

        match product.cmp(&self.data[1].product) {
            std::cmp::Ordering::Equal => return,
            std::cmp::Ordering::Less => match product.cmp(&self.data[0].product) {
                std::cmp::Ordering::Equal => self.data[1] = self.data[0],
                std::cmp::Ordering::Less => {
                    self.data[1] = sample_between(ProdEvents::zero(), self.data[0], product, rng);
                    self.data[0] = self.data[1];
                }
                std::cmp::Ordering::Greater => {
                    self.data[1] = sample_between(self.data[0], self.data[1], product, rng);
                }
            },
            std::cmp::Ordering::Greater => match product.cmp(&self.data[2].product) {
                std::cmp::Ordering::Equal => self.data[1] = self.data[2],
                std::cmp::Ordering::Less => {
                    self.data[1] = sample_between(self.data[1], self.data[2], product, rng);
                }
                std::cmp::Ordering::Greater => {
                    self.data[1].events = sample_poisson(
                        time * reaction.rate * (product - self.data[2].product) as f64,
                        rng,
                    ) + self.data[2].events;
                    self.data[1].product = product;
                    self.data[2] = self.data[1];
                }
            },
        }
        assert!(self.data[1].product == product);
        self.validate();
    }

    pub fn resample_bounds(
        &mut self,
        low_product: u64,
        high_product: u64,
        reaction: &Reaction,
        time: f64,
        rng: &mut impl Rng,
    ) {
        self.validate();
        debug_assert!(low_product <= self.data[1].product);
        debug_assert!(high_product >= self.data[1].product);
        match low_product.cmp(&self.data[0].product) {
            std::cmp::Ordering::Less => {
                self.data[0] = sample_between(ProdEvents::zero(), self.data[0], low_product, rng);
            }
            std::cmp::Ordering::Equal => {}
            std::cmp::Ordering::Greater => {
                self.data[0] = sample_between(self.data[0], self.data[1], low_product, rng)
            }
        }
        match high_product.cmp(&self.data[2].product) {
            std::cmp::Ordering::Less => {
                self.data[2] = sample_between(self.data[1], self.data[2], high_product, rng);
            }
            std::cmp::Ordering::Equal => {}
            std::cmp::Ordering::Greater => {
                self.data[2].events = self.data[2].events
                    + sample_poisson(
                        reaction.rate * time * (high_product - self.data[2].product) as f64,
                        rng,
                    );
                self.data[2].product = high_product;
            }
        }
        self.validate();
    }

    fn validate(&self) {
        debug_assert!(self.data[0].product <= self.data[1].product);
        debug_assert!(self.data[1].product <= self.data[2].product);
        debug_assert!(self.data[0].events <= self.data[1].events);
        debug_assert!(self.data[1].events <= self.data[2].events);
    }
}

pub fn sample_binomial(n: u64, p: f64, rng: &mut impl Rng) -> u64 {
    rng.sample(Binomial::new(n, p).unwrap())
}
pub fn sample_poisson(rate: f64, rng: &mut impl Rng) -> u64 {
    rng.sample(Poisson::new(rate).unwrap_or_else(|err| {
        panic!("Failed to sample a Poisson variable with rate {rate}: {err:?}")
    })) as u64
}

/// Samples the number of events at a point between two pairs of (product, events)
pub fn sample_between(
    low: ProdEvents,
    high: ProdEvents,
    mid: u64,
    rng: &mut impl Rng,
) -> ProdEvents {
    ProdEvents::new(
        mid,
        low.events
            + sample_binomial(
                high.events - low.events,
                (mid - low.product) as f64 / (high.product - low.product) as f64,
                rng,
            ),
    )
}
