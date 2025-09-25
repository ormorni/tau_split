use std::path::Path;

use kdam::tqdm;
use rand::{rngs::StdRng, SeedableRng};
use rustc_hash::FxHashMap;

use crate::{
    fastspie3::FastGillespie3, fastspie4::FastGillespie4, fastspie5::FastGillespie5,
    gillespie::Gillespie, parsers::ParseState, tests::chisq::same_categorical_dist, SimulationAlg,
};

/// The path to the system:
///
/// A -> \phi
const DECAY_PATH: &str = "data/test_models/decay.txt";
/// The path to the system:
///
/// \phi -> A
const SYNTHESIS_PATH: &str = "data/test_models/synthesis.txt";
/// The path to the system:
///
/// A + B -> B + C
/// B + C -> C + A
/// C + A -> A + B
const CONVERSION_CYCLE_PATH: &str = "data/test_models/cyclic.txt";

/// Tests that the chemical reaction network defined by the file in the given path
/// has the same distribution when simulated by the Gillespie and Tau-Splitting algorithms.
fn test_network<Alg: SimulationAlg>(path: &Path, n: usize, t: f64) {
    let mut parse_state = ParseState::default();
    parse_state.parse_data_file(path);

    let (initial_state, reactions, names) = parse_state.get_network();

    let mut gillespie_samples: FxHashMap<Vec<i64>, u64> = FxHashMap::default();
    let mut tau_split_samples: FxHashMap<Vec<i64>, u64> = FxHashMap::default();
    for i in tqdm!(0..n, desc = "Gillespie") {
        let rng = &mut StdRng::seed_from_u64(i as u64);

        let mut gill = Gillespie::new(initial_state.clone(), reactions.clone(), names.clone());
        gill.advance(t, rng);
        *gillespie_samples.entry(gill.state).or_default() += 1;
    }

    for i in tqdm!(0..n, desc = "Tau Split") {
        let rng = &mut StdRng::seed_from_u64(i as u64);

        let mut tau_split = Alg::new(initial_state.clone(), reactions.clone(), names.clone());
        tau_split.advance(t, rng);
        *tau_split_samples
            .entry(tau_split.state().to_owned())
            .or_default() += 1;
    }

    println!("{:?} {:?}", gillespie_samples, tau_split_samples);

    assert!(same_categorical_dist(gillespie_samples, tau_split_samples));
}

/// Tests that the system:
///
/// \phi -> A
///
/// behaves the same using the Gillespie algorithm and the Tau-Splitting algorithm.
#[test]
pub fn test_synthesis_3() {
    test_network::<FastGillespie3>(Path::new(SYNTHESIS_PATH), 1 << 16, 0.1);
}
/// Tests that the system:
///
/// \phi -> A
///
/// behaves the same using the Gillespie algorithm and the Tau-Splitting algorithm.
#[test]
pub fn test_synthesis_4() {
    test_network::<FastGillespie4>(Path::new(SYNTHESIS_PATH), 1 << 16, 0.1);
}
/// Tests that the system:
///
/// \phi -> A
///
/// behaves the same using the Gillespie algorithm and the Tau-Splitting algorithm.
#[test]
pub fn test_synthesis_5() {
    test_network::<FastGillespie5>(Path::new(SYNTHESIS_PATH), 1 << 16, 0.1);
}
/// Tests that the system:
///
/// A -> \phi
///
/// behaves the same using the Gillespie algorithm and the Tau-Splitting algorithm.
#[test]
pub fn test_decay_3() {
    test_network::<FastGillespie3>(Path::new(DECAY_PATH), 1 << 16, 0.1);
}
/// Tests that the system:
///
/// A -> \phi
///
/// behaves the same using the Gillespie algorithm and the Tau-Splitting algorithm.
#[test]
pub fn test_decay_4() {
    test_network::<FastGillespie4>(Path::new(DECAY_PATH), 1 << 16, 0.1);
}
/// Tests that the system:
///
/// A -> \phi
///
/// behaves the same using the Gillespie algorithm and the Tau-Splitting algorithm.
#[test]
pub fn test_decay_5() {
    test_network::<FastGillespie5>(Path::new(DECAY_PATH), 1 << 16, 0.1);
}

/// Tests that the system:
///
/// A + B -> B + C
/// B + C -> C + A
/// C + A -> A + B
///
/// behaves the same using the Gillespie algorithm and the Tau-Splitting algorithm.
#[test]
pub fn test_conversion_cycle_3() {
    test_network::<FastGillespie3>(Path::new(CONVERSION_CYCLE_PATH), 1 << 16, 5.);
}
/// Tests that the system:
///
/// A + B -> B + C
/// B + C -> C + A
/// C + A -> A + B
///
/// behaves the same using the Gillespie algorithm and the Tau-Splitting algorithm.
#[test]
pub fn test_conversion_cycle_4() {
    test_network::<FastGillespie4>(Path::new(CONVERSION_CYCLE_PATH), 1 << 16, 5.);
}

/// Tests that the system:
///
/// A + B -> B + C
/// B + C -> C + A
/// C + A -> A + B
///
/// behaves the same using the Gillespie algorithm and the Tau-Splitting algorithm.
#[test]
pub fn test_conversion_cycle_5() {
    test_network::<FastGillespie5>(Path::new(CONVERSION_CYCLE_PATH), 1 << 16, 5.);
}
