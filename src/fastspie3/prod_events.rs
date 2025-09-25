use derive_new::new;
use std::fmt::Debug;

#[derive(Clone, Copy, new, Default, PartialEq, Eq, PartialOrd, Ord)]
pub struct ProdEvents {
    pub product: u64,
    pub events: u64,
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
