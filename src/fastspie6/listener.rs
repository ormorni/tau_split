use std::{fmt::Debug, ops::Neg};

use dary_heap::{DaryHeap, OctonaryHeap};
use derive_new::new;

/// A max-heap that stores elements ordered by an attached key,
/// and supports a clean interface for popping all objects larger than a key.
#[derive(Debug, Default, Clone, new)]
pub struct MaxListener<Key: Ord, Data: Ord> {
    /// A max-heap storing the data.
    heap: OctonaryHeap<(Key, Data)>,
}

impl<Key: Ord + Debug, Data: Ord> MaxListener<Key, Data> {
    /// Pushes an item onto the listener.
    pub fn push(&mut self, key: Key, data: Data) {
        self.heap.push((key, data));
    }

    /// If the largest element is larger than the key, pop it and return the data.
    /// Otherwise, return None.
    pub fn pop_if_larger_than(&mut self, key: Key) -> Option<Data> {
        let Some((top_key, _)) = self.heap.peek() else {
            return None;
        };
        if &key < top_key {
            return self.heap.pop().map(|(_, data)| data);
        } else {
            None
        }
    }

    pub fn iter(&self) -> impl Iterator<Item = &(Key, Data)> {
        self.heap.iter()
    }

    pub fn len(&self) -> usize {
        self.heap.len()
    }

    pub fn retain(&mut self, f: impl FnMut(&(Key, Data)) -> bool) {
        self.heap.retain(f);
    }
}

impl<'t, Key: Ord, Data: Ord> IntoIterator for &'t MaxListener<Key, Data> {
    type Item = &'t (Key, Data);

    type IntoIter = <&'t OctonaryHeap<(Key, Data)> as IntoIterator>::IntoIter;

    fn into_iter(self) -> Self::IntoIter {
        self.heap.iter()
    }
}

impl<Key: Ord, Data: Ord> IntoIterator for MaxListener<Key, Data> {
    type Item = (Key, Data);
    type IntoIter = dary_heap::IntoIter<(Key, Data)>;

    fn into_iter(self) -> Self::IntoIter {
        self.heap.into_iter()
    }
}

impl<Key: Ord, Data: Ord> FromIterator<(Key, Data)> for MaxListener<Key, Data> {
    fn from_iter<T: IntoIterator<Item = (Key, Data)>>(iter: T) -> Self {
        MaxListener {
            heap: DaryHeap::from_iter(iter),
        }
    }
}

/// A max-heap that stores elements ordered by an attached key,
/// and supports a clean interface for popping all objects smaller than a key.
#[derive(Default, Debug, Clone, new)]
pub struct MinListener<Key: Ord, Data: Ord> {
    heap: OctonaryHeap<(Key, Data)>,
}

impl<Key: Ord + Neg<Output = Key>, Data: Ord> MinListener<Key, Data> {
    /// Pushes an item onto the listener.
    pub fn push(&mut self, key: Key, data: Data) {
        self.heap.push((-key, data));
    }

    /// If the smallest element is smaller than the key, pop it and return the data.
    /// Otherwise, return None.
    pub fn pop_if_smaller_than(&mut self, key: Key) -> Option<Data> {
        let Some((top_key, _)) = self.heap.peek() else {
            return None;
        };

        if &-key < top_key {
            return self.heap.pop().map(|(_, data)| data);
        } else {
            None
        }
    }
    pub fn iter(&self) -> impl Iterator<Item = &(Key, Data)> {
        self.heap.iter()
    }

    pub fn len(&self) -> usize {
        self.heap.len()
    }

    pub fn retain(&mut self, f: impl FnMut(&(Key, Data)) -> bool) {
        self.heap.retain(f);
    }
}

impl<'t, Key: Ord + Neg<Output = Key>, Data: Ord> IntoIterator for &'t MinListener<Key, Data> {
    type Item = &'t (Key, Data);

    type IntoIter = <&'t OctonaryHeap<(Key, Data)> as IntoIterator>::IntoIter;

    fn into_iter(self) -> Self::IntoIter {
        self.heap.iter()
    }
}

impl<Key: Ord + Neg<Output = Key>, Data: Ord> IntoIterator for MinListener<Key, Data> {
    type Item = (Key, Data);
    type IntoIter = dary_heap::IntoIter<(Key, Data)>;

    fn into_iter(self) -> Self::IntoIter {
        self.heap.into_iter()
    }
}

impl<Key: Ord + Neg<Output = Key>, Data: Ord> FromIterator<(Key, Data)> for MinListener<Key, Data> {
    fn from_iter<T: IntoIterator<Item = (Key, Data)>>(iter: T) -> Self {
        MinListener {
            heap: DaryHeap::from_iter(iter),
        }
    }
}
