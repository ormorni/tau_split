use rand::Rng;
use rand_distr::Binomial;

/// A default seed for seeded RNGs.
pub const DEFAULT_SEED: u64 = 0x123456789abcdef;

/// A specializing implementation for a binomial random variable with p=0.5.
pub fn binomial_05(n: u64, rng: &mut impl Rng) -> u64 {
    if n == 0 {
        0
    } else if n <= 64 {
        (rng.random::<u64>() >> (64 - n)).count_ones() as u64
    } else {
        rng.sample(Binomial::new(n, 0.5).unwrap())
    }
}


/// A macro to temporarily borrow vectors in a lifetime-safe-ish way.
/// Takes possession of the memory of the vector for the duration of the block.
/// The vector should definitely not be accessed during the macro's operation,
/// as that would lead to bad and unexpected results.
#[macro_export]
macro_rules! with {
    ($vec:tt as $name:tt $code:block) => {
        let mut $name = std::mem::take(&mut $vec);
        $code
        debug_assert!($name.is_empty(), "The vector '{}' used by the 'with' macro is not empty at the end of the block!", stringify!($vec));
        std::mem::swap(&mut $vec, &mut $name);
    };
}
