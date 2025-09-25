use std::{hash::Hash, ops::Deref};

use derive_new::new;

#[derive(new)]
pub struct Ref<'t, T> {
    data: &'t T,
}

impl<'t, T> Deref for Ref<'t, T> {
    type Target = T;

    fn deref(&self) -> &Self::Target {
        self.data
    }
}

impl<'t, T> PartialEq for Ref<'t, T> {
    fn eq(&self, other: &Self) -> bool {
        std::ptr::eq(self.data as *const _, other.data as *const _)
    }
}

impl<'t, T> Eq for Ref<'t, T> {}

impl<'t, T> Hash for Ref<'t, T> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        (self.data as *const T).hash(state);
    }
}
