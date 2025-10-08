use criterion::{criterion_group, criterion_main, Criterion};
use rand::{rngs::StdRng, SeedableRng};
use std::{hint::black_box, path::Path};

use tausplit::{TauSplit5, ParseState, SimulationAlg, DEFAULT_SEED};

const BCR_HIGH_PATH: &str = "data/models/B cell antigen receptor signaling/BCR_high.txt";
const BCR_HIGH_TIME: f64 = 0.0009;
const FCERI_HIGH_PATH: &str = "data/models/FceRI/FceRI_high.txt";
const FCERI_HIGH_TIME: f64 = 0.027;

fn criterion_benchmark(c: &mut Criterion) {
    let mut parse_state = ParseState::default();
    parse_state.parse_data_file(Path::new(BCR_HIGH_PATH));
    let (initial_state, reactions, names) = parse_state.get_network();

    c.bench_function("BCR high", |b| {
        b.iter(|| {
            let rng = &mut StdRng::seed_from_u64(black_box(DEFAULT_SEED));
            let mut fastspie5 =
                TauSplit5::new(initial_state.clone(), reactions.clone(), names.clone());
            fastspie5.advance(BCR_HIGH_TIME, rng);
        })
    });

    let mut parse_state = ParseState::default();
    parse_state.parse_data_file(Path::new(FCERI_HIGH_PATH));
    let (initial_state, reactions, names) = parse_state.get_network();

    c.bench_function("FceRI high", |b| {
        b.iter(|| {
            let rng = &mut StdRng::seed_from_u64(black_box(DEFAULT_SEED));
            let mut fastspie5 =
                TauSplit5::new(initial_state.clone(), reactions.clone(), names.clone());
            fastspie5.advance(FCERI_HIGH_TIME, rng);
        })
    });
}

criterion_group! {
    name=benches; config=Criterion::default().sample_size(10); targets=criterion_benchmark
}
criterion_main!(benches);
