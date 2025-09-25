use std::{
    fmt::Debug,
    marker::PhantomData,
    ops::{Index, IndexMut},
    slice::Iter,
};

/// A Vec struct whose indices are marked in some way.
/// This ensures that indices of several vector types cannot get mixed up.
#[derive(Debug, Default, Clone)]
pub struct IndexedVec<Marker, T> {
    data: Vec<T>,
    _phantom: PhantomData<Marker>,
}

/// Swaps two elements in the slice.
/// If a equals to b, it's guaranteed that elements won't change value.
impl<Marker, T> IndexedVec<Marker, T> {
    pub fn swap(&mut self, a: Idx<Marker>, b: Idx<Marker>) {
        self.data.swap(a.index, b.index);
    }

    pub fn push(&mut self, data: T) -> Idx<Marker> {
        let idx = Idx::new(self.data.len());
        self.data.push(data);
        idx
    }

    pub fn contains_key(&self, idx: Idx<Marker>) -> bool {
        idx.index < self.data.len()
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn iter(&self) -> Iter<'_, T> {
        self.data.iter()
    }

    pub fn pop(&mut self) -> Option<T> {
        self.data.pop()
    }
}

impl<Marker, T> Index<Idx<Marker>> for IndexedVec<Marker, T> {
    type Output = T;

    fn index(&self, index: Idx<Marker>) -> &Self::Output {
        &self.data[index.index]
    }
}

impl<Marker, T> IndexMut<Idx<Marker>> for IndexedVec<Marker, T> {
    fn index_mut(&mut self, index: Idx<Marker>) -> &mut Self::Output {
        &mut self.data[index.index]
    }
}

impl<Marker, T> From<Vec<T>> for IndexedVec<Marker, T> {
    fn from(value: Vec<T>) -> Self {
        IndexedVec {
            data: value,
            _phantom: PhantomData::default(),
        }
    }
}

impl<'t, Marker, T> IntoIterator for &'t IndexedVec<Marker, T> {
    type Item = &'t T;

    type IntoIter = Iter<'t, T>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

pub struct Idx<Marker> {
    index: usize,
    _phantom: PhantomData<Marker>,
}

impl<Marker> Idx<Marker> {
    pub const fn new(index: usize) -> Idx<Marker> {
        Idx {
            index,
            _phantom: PhantomData {},
        }
    }
}

impl<Marker> PartialEq for Idx<Marker> {
    fn eq(&self, other: &Self) -> bool {
        self.index == other.index
    }
}

impl<Marker> Eq for Idx<Marker> {}

impl<Marker> PartialOrd for Idx<Marker> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.index.partial_cmp(&other.index)
    }
}

impl<Marker> Ord for Idx<Marker> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.index.cmp(&other.index)
    }
}

impl<Marker> Clone for Idx<Marker> {
    fn clone(&self) -> Self {
        Self {
            index: self.index,
            _phantom: self._phantom,
        }
    }
}

impl<Marker> Copy for Idx<Marker> {}

impl<Marker> Debug for Idx<Marker> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Idx")
            .field("index", &self.index)
            .field("_phantom", &self._phantom)
            .finish()
    }
}
