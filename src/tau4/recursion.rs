use dary_heap::OctonaryHeap;
use derive_new::new;
use itertools::Itertools;
use rand::Rng;
use tinyvec::ArrayVec;

const MAX_INPUTS: usize = 2;

use crate::{reaction::Reaction};

use super::{ReactionData, StateData};

#[derive(new)]
pub struct RecursionTree<'t> {
    nodes: Vec<RecursionTreeNode>,
    inactive_index: Vec<Option<(usize, usize)>>,
    reactions: &'t [Reaction],
    state: StateData,

    // Stability trackers.
    // When a reaction was considered stable but is unstable we have to check
    // if there are inactive reactions feeding into it.
    stored_stable: Vec<bool>,
    unstable_dependent: Vec<usize>,

    /// The number of reactions simulated up to now.
    pub total_events: u64,

    inactive_by_component: Vec<Vec<usize>>,

    positive_listeners: Vec<OctonaryHeap<(i64, usize, usize)>>,
    negative_listeners: Vec<OctonaryHeap<(i64, usize, usize)>>,
    node_id: usize,
}

#[derive(new)]
pub struct RecursionTreeNode {
    reactions: Vec<ReactionData>,
    /// A node is active if the timespan it represents contains the current timepoint.
    is_active: bool,
    parent: Option<usize>,
    left: Option<usize>,
    right: Option<usize>,
    id: usize,
}

impl RecursionTreeNode {
    pub fn is_leaf(&self) -> bool {
        self.left.is_none() && self.right.is_none()
    }
}

impl<'t> RecursionTree<'t> {
    /// Sets the given node to be active.
    ///
    /// This means that all reactions in the node should be added to the bounds.
    pub fn activate_node(&mut self, node: usize) {
        debug_assert!(!self.nodes[node].is_active);
        for rdata in &self.nodes[node].reactions {
            self.state
                .add_bounds(rdata, &self.reactions[rdata.reaction]);
        }
        self.nodes[node].is_active = true;
    }

    pub fn recursion(&mut self, node: usize, time: f64, rng: &mut impl Rng) {
        // At the beginning of the recursion,
        // the bounds include all reactions in internal nodes, but not the leaf node.
        self.activate_node(node);
        // Getting things out of the static data.
        let mut active_reactions = Vec::new();

        self.validate_bounds(&active_reactions);
        self.validate_inactive_dependence();
        // self.validate_inactive_correct(node);
        // After this point, the bounds include everything.

        self.forward_reactivation(&mut active_reactions, node, time, rng);
        self.backward_reactivation(&mut active_reactions, node, time, rng);

        // Checking if all reactions have been stabilized.
        // If all have been stabilized, we apply all of them, and return.
        if active_reactions
            .iter()
            .all(|rdata| self.stored_stable[rdata.reaction])
        {
            for rdata in &active_reactions {
                self.state
                    .remove_bounds(rdata, &self.reactions[rdata.reaction]);
                self.state.apply(rdata, &self.reactions[rdata.reaction]);
                self.total_events += rdata.events;
            }

            self.remove_node(node);
            return;
        }
        // println!("The system is unstable!");
        self.validate_bounds(&active_reactions);

        // Not all reactions have been stabilized.
        // Deactivating all reactions satisfying:
        // * They are stable.
        // * All reactions depending on them are stable.
        // The inactive reaction data remain in `active_rdata`.
        // The active ones go to `left_rdata`, and are then split.
        // We have to remove them from the bounds.
        let mut left_rdata = active_reactions
            .extract_if(.., |rdata| {
                !self.stored_stable[rdata.reaction]
                    || (rdata.events > 0
                        && self.reactions[rdata.reaction]
                            .stoichiometry
                            .iter()
                            .any(|(reactant, _)| self.unstable_dependent[*reactant] > 0))
            })
            .collect_vec();

        let mut right_rdata = Vec::with_capacity(left_rdata.len());
        for rdata in &mut left_rdata {
            self.state
                .remove_bounds(rdata, &self.reactions[rdata.reaction]);
            right_rdata.push(rdata.split(&self.reactions[rdata.reaction], rng));
        }

        let right_node = self.add_node(node, right_rdata);
        self.nodes[node].right = Some(right_node);

        let left_node = self.add_node(node, left_rdata);
        self.nodes[node].left = Some(left_node);

        // This has to happen after the node is an internal node to properly deactivate the reaction.
        self.validate_bounds(&active_reactions);
        self.validate_inactive_correct();

        for rdata in active_reactions {
            debug_assert!(self.is_stable(&rdata));
            self.state
                .remove_bounds(&rdata, &self.reactions[rdata.reaction]);

            self.add_reaction(node, rdata);

            self.add_negative_listeners(&rdata, self.nodes[node].id);
            self.add_positive_listeners(&rdata, self.nodes[node].id);
        }
        self.validate_inactive_correct();
        self.validate_bounds(&[]);

        self.recursion(left_node, time / 2., rng);
        self.recursion(right_node, time / 2., rng);

        // Applying all the inactive reactions remaining in the node.
        for rdata in &self.nodes[node].reactions {
            let reaction = &self.reactions[rdata.reaction];
            // println!(
            //     "Applying {rdata:?} {}",
            //     reaction.input_product(&self.state.state)
            // );
            self.state.remove_bounds(&rdata, reaction);
            self.inactive_index[rdata.reaction] = None;
            self.state.apply(rdata, reaction);
            self.total_events += rdata.events;
        }
        self.remove_node(node);
    }

    /// Samples the new event count for all reactions, computes the new lower and upper bounds,
    /// and reactivates all reactions that assumed bounds different from the current one.
    fn forward_reactivation(
        &mut self,
        active_reactions: &mut Vec<ReactionData>,
        node: usize,
        time: f64,
        rng: &mut impl Rng,
    ) {
        let is_right_child = self.is_right_child(node);
        // let is_right_child = true;
        // Taking all reactions and computing the event count.
        // If the event count changes, there has been a second-order event,
        // and we have to reactivate all dependent reactions.
        while let Some(mut rdata) = self.nodes[node].reactions.pop() {
            let reaction = &self.reactions[rdata.reaction];

            if is_right_child {
                let prod = self.state.state_product(reaction);
                let old_events = rdata.events;
                // let mut reactivation = false;
                let old_rdata = rdata;
                rdata.resample(prod, reaction, rng);
                debug_assert!(
                    (rdata.low..rdata.high).contains(&(prod)),
                    "Resmapled {old_rdata:?} => {rdata:?} with product {prod}"
                );
                self.state
                    .change_bounds(rdata.events as i64 - old_events as i64, reaction);
            }
            active_reactions.push(rdata);

            for &(comp, _) in &reaction.stoichiometry {
                // Updating the positive listeners.
                while !self.positive_listeners[comp].is_empty()
                    && self.state[comp].upper >= -self.positive_listeners[comp].peek().unwrap().0
                {
                    let (cutoff, reaction, node_id) = self.positive_listeners[comp].pop().unwrap();
                    // println!("P {cutoff} {reaction} {node_id}");
                    let Some((node_idx, vec_idx)) = self.inactive_index[reaction] else {
                        continue;
                    };
                    // The listener was there due to an unrelated node.
                    if node_id != self.nodes[node_idx].id {
                        continue;
                    }
                    let new_upper = self.state.upper_product(&self.reactions[reaction]);
                    // println!(
                    //     "Forward upper reactivating {:?} for reaction {:?}",
                    //     self.nodes[node_idx].reactions[vec_idx], self.reactions[reaction]
                    // );
                    // println!("cutoff={} comp={}", cutoff, self.state[comp].upper);
                    if new_upper >= self.nodes[node_idx].reactions[vec_idx].high {
                        self.reactivate_reaction(reaction, rng);
                    } else {
                        self.add_positive_listeners(
                            &self.nodes[node_idx].reactions[vec_idx].clone(),
                            node_id,
                        );
                    }
                }
                // Updating the negative listeners.
                while !self.negative_listeners[comp].is_empty()
                    && self.state[comp].lower <= self.negative_listeners[comp].peek().unwrap().0
                {
                    let (cutoff, reaction, node_id) = self.negative_listeners[comp].pop().unwrap();
                    // println!("N {cutoff} {reaction} {node_id}");

                    let Some((node_idx, vec_idx)) = self.inactive_index[reaction] else {
                        continue;
                    };
                    // This is an outdated listener from an old inactive node..
                    if node_id != self.nodes[node_idx].id {
                        continue;
                    }
                    // println!(
                    //     "Forward lower reactivating {:?} for reaction {:?}",
                    //     self.nodes[node_idx].reactions[vec_idx], self.reactions[reaction]
                    // );
                    let old_lower = self.nodes[node_idx].reactions[vec_idx].low;
                    let new_lower = self.state.lower_product(&self.reactions[reaction]);
                    if new_lower < old_lower {
                        self.reactivate_reaction(reaction, rng);
                    } else {
                        self.add_negative_listeners(
                            &self.nodes[node_idx].reactions[vec_idx].clone(),
                            node_id,
                        );
                    }
                }
            }
        }

        self.validate_bounds(&active_reactions);
        self.validate_inactive_dependence();
        self.validate_inactive_correct();
    }

    fn backward_reactivation(
        &mut self,
        active_reactions: &mut Vec<ReactionData>,
        node: usize,
        time: f64,
        rng: &mut impl Rng,
    ) {
        self.validate_inactive_correct();
        // At this point the bounds are accurate,
        // but there may be reactions that depend on inactive reactions.
        // These reactions are all stable, so they do not wake any other reaction.
        // for rdata in active_reactions.iter_mut() {
        //     let reaction = &self.reactions[rdata.reaction];
        //     let mid_product = reaction.input_product(&self.state.state) ;
        //     debug_assert!(
        //         (rdata.low..rdata.high).contains(&mid_product),
        //         "Bad products: {rdata:?} {mid_product}"
        //     );
        // }

        for rdata in active_reactions.iter_mut() {
            let reaction = &self.reactions[rdata.reaction];
            // let mid_product = reaction.input_product(&self.state.state) ;
            // debug_assert!(
            //     (rdata.low..rdata.high).contains(&mid_product),
            //     "Bad products: {rdata:?} {mid_product}"
            // );
            // let low_product = reaction.input_product(&self.state.lower_bound) ;
            // let high_product = self.state.upper_product(reaction) ;
            // debug_assert!(low_product <= mid_product);
            // debug_assert!(high_product >= mid_product);

            // Checking if the reaction is unstable, but is considered stable.
            // In that case, we start reactivating the reactions on which this reaction depends
            // until the reaction is restabilized or we run out of splits.
            let stable = self.is_stable(rdata);
            if self.stored_stable[rdata.reaction] && !self.is_stable(rdata) {
                // The components that would be reactivated if the reaction isn't stable.
                let mut destabilized_components: ArrayVec<[usize; MAX_INPUTS]> = ArrayVec::new();
                for &(comp, _) in &reaction.inputs {
                    if self.unstable_dependent[comp] == 0 {
                        destabilized_components.push(comp);
                    }
                }

                // If there are inactive input components, we reactivate them.
                if !destabilized_components.is_empty() {
                    // While the reaction is still unstable we push down reactions and see if it helped stabilize anything.
                    for &component in &destabilized_components {
                        let mut v = std::mem::take(&mut self.inactive_by_component[component]);
                        for reaction in v.drain(..) {
                            self.reactivate_reaction(reaction, rng);
                        }
                        std::mem::swap(&mut v, &mut self.inactive_by_component[component]);
                    }
                }
                // We redo the stability counting.
                self.update_stability(&rdata);
            } else {
                if stable {
                    self.set_stable(&rdata);
                } else {
                    self.set_unstable(&rdata);
                }
            }
        }

        // Moving all the reactions that were reactivated in the previous phase to the active reactions.
        active_reactions.extend(self.nodes[node].reactions.drain(..));

        self.validate_bounds(&*active_reactions);
        self.validate_dependent(&*active_reactions);
        self.validate_inactive_dependence();
    }

    /// Checks that the dependent unstable reaction counter is valid.
    fn validate_dependent(&self, active_reactions: &[ReactionData]) {
        if cfg!(debug_assertions) {
            let mut dependents = vec![0; self.unstable_dependent.len()];

            for rdata in active_reactions {
                if !self.is_stable(rdata) {
                    for &(reactant, _) in &self.reactions[rdata.reaction].inputs {
                        dependents[reactant] += 1;
                    }
                }
            }
            debug_assert_eq!(dependents, self.unstable_dependent);
        }
    }

    /// Validates that no unstable reaction requires an output of an inactive component.
    ///
    /// Every unstable reaction requires all nonzero reactions feeding into it to be active.
    fn validate_inactive_dependence(&self) {
        if cfg!(debug_assertions) {
            for reaction_idx in 0..self.reactions.len() {
                let Some((node_idx, vec_idx)) = self.inactive_index[reaction_idx] else {
                    continue;
                };
                if self.nodes[node_idx].reactions[vec_idx].events != 0 {
                    for &(component, _) in &self.reactions[reaction_idx].stoichiometry {
                        assert!(self.unstable_dependent[component] == 0);
                    }
                }
            }
        }
    }

    // Validates that the bounds of all inactive reactions are correct.
    fn validate_inactive_correct(&self) {
        if cfg!(debug_assertions) {
            for reaction_idx in 0..self.reactions.len() {
                let Some((node_idx, vec_idx)) = self.inactive_index[reaction_idx] else {
                    continue;
                };

                let rdata = &self.nodes[node_idx].reactions[vec_idx];
                let reaction = &self.reactions[rdata.reaction];

                if rdata.low > self.state.lower_product(reaction) && rdata.events != 1
                // The full condition is the "from one to zero" condition, which I should refactor to its own function.
                {
                    println!("Error in reaction {:?}!", rdata);
                    println!(
                        "Lower bound: {:?}",
                        reaction
                            .inputs
                            .iter()
                            .map(|&(comp, _)| self.state[comp].lower)
                            .collect_vec()
                    );
                    println!(
                        "Listeners: {:?}",
                        reaction
                            .inputs
                            .iter()
                            .map(|&(comp, _)| self.negative_listeners[comp]
                                .iter()
                                .find(|&&(_, r, n_id)| r == rdata.reaction
                                    && n_id == self.nodes[node_idx].id))
                            .collect_vec()
                    );
                    panic!()
                }
                if rdata.high < self.state.upper_product(reaction) {
                    println!("Error in reaction {rdata:?}!",);
                    println!(
                        "Upper bound: {:?}",
                        reaction
                            .inputs
                            .iter()
                            .map(|&(comp, _)| self.state[comp].upper)
                            .collect_vec()
                    );
                    println!(
                        "Listeners: {:?}",
                        reaction
                            .inputs
                            .iter()
                            .map(|&(comp, _)| self.positive_listeners[comp]
                                .iter()
                                .find(|&&(_, r, n_id)| r == rdata.reaction
                                    && n_id == self.nodes[node_idx].id))
                            .collect_vec()
                    );
                    panic!()
                }
            }
        }
    }

    /// Checks that the reaction bounds are valid.
    fn validate_bounds(&self, active_reactions: &[ReactionData]) {
        if cfg!(debug_assertions) {
            let mut s = StateData::new(&self.state());
            for rdata in active_reactions {
                // print!("({},{}),", rdata.reaction, rdata.events);
                s.apply_negative(rdata.events as i64, &self.reactions[rdata.reaction]);
            }
            for node in &self.nodes {
                if node.is_active {
                    for rdata in &node.reactions {
                        // print!("({},{}),", rdata.reaction, rdata.events);
                        s.apply_negative(rdata.events as i64, &self.reactions[rdata.reaction]);
                    }
                }
            }

            for rdata in active_reactions {
                s.apply_positive(rdata.events as i64, &self.reactions[rdata.reaction]);
            }
            for node in &self.nodes {
                if node.is_active {
                    for rdata in &node.reactions {
                        s.apply_positive(rdata.events as i64, &self.reactions[rdata.reaction]);
                    }
                }
            }
            debug_assert_eq!(s, self.state);
        }
    }

    pub fn add_node(&mut self, parent: usize, rdata: Vec<ReactionData>) -> usize {
        self.nodes.push(RecursionTreeNode {
            reactions: rdata,
            is_active: false,
            parent: Some(parent),
            left: None,
            right: None,
            id: self.node_id,
        });
        self.node_id += 1;
        self.nodes.len() - 1
    }

    pub fn remove_node(&mut self, node: usize) {
        if let Some(parent) = self.nodes[node].parent {
            if self.nodes[parent].left == Some(node) {
                self.nodes[parent].left = None;
            } else {
                debug_assert_eq!(self.nodes[parent].right, Some(node));
                self.nodes[parent].right = None;
            }
        }
        // println!("Nodes ({node}):");
        // for (idx, node) in self.nodes.iter().enumerate() {
        //     println!("{idx}: {:?} {:?}", node.left, node.right);
        // }
        debug_assert!(node + 1 == self.nodes.len());
        self.nodes.pop();
    }

    /// Updates the stored stability of the reaction.
    /// This is important to allow reactivating reactions.
    pub fn update_stability(&mut self, rdata: &ReactionData) {
        if self.is_stable(rdata) {
            self.set_stable(rdata);
        } else {
            self.set_unstable(rdata);
        }
    }

    pub fn set_stable(&mut self, rdata: &ReactionData) {
        debug_assert!(self.is_stable(rdata));
        let reaction = &self.reactions[rdata.reaction];
        if !self.stored_stable[rdata.reaction] {
            for &(component, _) in &reaction.inputs {
                self.unstable_dependent[component] -= 1;
            }
            self.stored_stable[rdata.reaction] = true;
        }
    }
    pub fn set_unstable(&mut self, rdata: &ReactionData) {
        debug_assert!(!self.is_stable(rdata));
        let reaction = &self.reactions[rdata.reaction];
        if self.stored_stable[rdata.reaction] {
            for &(component, _) in &reaction.inputs {
                self.unstable_dependent[component] += 1;
            }
            self.stored_stable[rdata.reaction] = false;
        }
    }

    /// Adds a reaction to a node.
    /// If the node is an internal node, the reaction is inactive.
    /// Otherwise, the reaction is active.
    pub fn add_reaction(&mut self, node: usize, rdata: ReactionData) {
        debug_assert!(self.inactive_index[rdata.reaction].is_none());
        if !self.nodes[node].is_leaf() {
            // println!("Adding {} to internal node {node}", rdata.reaction);
            debug_assert!(self.inactive_index[rdata.reaction].is_none());
            self.inactive_index[rdata.reaction] = Some((node, self.nodes[node].reactions.len()));
            if rdata.events > 0 {
                for &(component, _) in &self.reactions[rdata.reaction].stoichiometry {
                    self.inactive_by_component[component].push(rdata.reaction);
                }
            }
        }
        if self.nodes[node].is_active {
            self.state
                .add_bounds(&rdata, &self.reactions[rdata.reaction]);
        }
        self.nodes[node].reactions.push(rdata);
    }

    /// If the given reaction is inactive, removes it from the inactive reaction data structure.
    ///
    /// Since inactive reactions are added to the state, we
    fn remove_reaction(&mut self, reaction_idx: usize) -> Option<ReactionData> {
        let reaction = &self.reactions[reaction_idx];
        let (node, vec_idx) = self.inactive_index[reaction_idx]?;

        // Removing the ReactionData from the inactive reaction graph.
        // We do a swap-remove, and update the relevant pointer into the data structure.
        if vec_idx + 1 != self.nodes[node].reactions.len() {
            let last_idx = self.nodes[node].reactions.len() - 1;
            let last_reaction = self.nodes[node].reactions.last().unwrap().reaction;
            self.nodes[node].reactions.swap(vec_idx, last_idx);
            self.inactive_index.swap(reaction_idx, last_reaction);
        }
        self.inactive_index[reaction_idx] = None;

        // If the reaction was in an internal node, it was present in the bounds, and has to be removed.
        let rdata = self.nodes[node].reactions.pop().unwrap();
        if self.nodes[node].is_active {
            self.state.remove_bounds(&rdata, reaction);
        }
        Some(rdata)
    }

    pub fn reactivate_reaction(&mut self, reaction_idx: usize, rng: &mut impl Rng) {
        let Some((mut node, _)) = self.inactive_index[reaction_idx] else {
            return;
        };
        let mut rdata = self.remove_reaction(reaction_idx).unwrap();
        let reaction = &self.reactions[reaction_idx];

        loop {
            match (self.nodes[node].left, self.nodes[node].right) {
                (None, None) => {
                    self.add_reaction(node, rdata);
                    break;
                }
                (None, Some(right)) => {
                    self.state.apply(&rdata.split(reaction, rng), reaction);
                    self.total_events += rdata.events;
                    node = right;
                }
                (Some(_left), None) => unreachable!(),
                (Some(left), Some(right)) => {
                    self.add_reaction(right, rdata.split(reaction, rng));
                    node = left;
                }
            }
        }
    }

    /// Checks if the reaction is now stable.
    ///
    /// A reaction is stable if either:
    /// * Its event count is independent of the current error
    /// * There is only one event, and that event brings the input product below the lower bound.
    pub fn is_stable(&self, rdata: &ReactionData) -> bool {
        let reaction = &self.reactions[rdata.reaction];
        let lower_product = self.state.lower_product(reaction);
        let upper_product = self.state.upper_product(reaction);

        let lower_legal = rdata.low <= lower_product;
        let upper_legal = rdata.high > upper_product;

        let from_one_to_zero = rdata.events == 1
            && (self.reactions[rdata.reaction]
                .input_stoichiometry()
                .all(|&(reactant, diff)| {
                    self.state[reactant].value + diff.max(0) == self.state[reactant].upper
                        && self.state[reactant].value + diff.min(0) == self.state[reactant].lower
                }));

        let stable = (lower_legal || from_one_to_zero) && upper_legal;

        stable
    }

    /// Returns the current state.
    pub fn state(&self) -> Vec<i64> {
        self.state.state.iter().map(|c| c.value).collect()
    }

    pub fn add_positive_listeners(&mut self, rdata: &ReactionData, node_id: usize) {
        let reaction = &self.reactions[rdata.reaction];
        let upper_bound = rdata.high;
        let curr_prod = self.state.upper_product(reaction);
        debug_assert!(
            upper_bound >= curr_prod,
            "curr_prod={curr_prod}, upper_bound={upper_bound}"
        ); // The reaction has to be stable for us to add lsiteners.

        if reaction.inputs.len() == 0 {
        } else if reaction.inputs.len() == 1 && reaction.inputs[0].1 == 1 {
            let component = reaction.inputs[0].0;
            // println!("A {old_prod} {:?}", (-old_prod, rdata.reaction, node_id));
            self.positive_listeners[component].push((
                -(upper_bound.floor() as i64 + 1),
                rdata.reaction,
                node_id,
            ));
        } else if reaction.inputs.len() == 1 && reaction.inputs[0].1 == 2 {
            let target = (1. + (1. + upper_bound * 8.).sqrt()) / 2.;
            let comp = reaction.inputs[0].0;
            debug_assert!(target.ceil() as i64 > self.state[comp].value);
            self.positive_listeners[comp].push((
                -(target.floor() as i64 + 1),
                rdata.reaction,
                node_id,
            ));
        } else if reaction.inputs.len() == 2 {
            // TODO: Fix this and the negative listeners for when there is 0*something.
            // To add listeners to a binary reaction, we assume that the ratio generally stays the same.
            debug_assert!(reaction.inputs[0].1 == 1);
            debug_assert!(reaction.inputs[1].1 == 1);
            if curr_prod == 0. {
                for &(comp, _) in &reaction.inputs {
                    if self.state[comp].upper == 0 {
                        self.positive_listeners[comp].push((-1, rdata.reaction, node_id));
                    }
                }
            } else {
                let ratio = ((upper_bound) / (curr_prod)).sqrt();

                for &(comp, _) in &reaction.inputs {

                    self.positive_listeners[comp].push((
                        -((self.state[comp].upper as f64 * ratio).floor() as i64 + 1),
                        rdata.reaction,
                        node_id,
                    ));
                }
            }
        } else {
            panic!("Reaction {reaction:?} not supported!");
        }
    }
    pub fn add_negative_listeners(&mut self, rdata: &ReactionData, node_id: usize) {
        if rdata.events == 0 {
            return;
        }
        let reaction = &self.reactions[rdata.reaction];
        let lower_cutoff = rdata.low;
        let curr_prod = self.state.lower_product(reaction);
        debug_assert!(
            lower_cutoff <= curr_prod || rdata.events == 1,
            "low={lower_cutoff} prod={curr_prod} rdata={rdata:?}"
        );
        if reaction.inputs.len() == 0 {
        } else if reaction.inputs.len() == 1 && reaction.inputs[0].1 == 1 {
            let component = reaction.inputs[0].0;
            let cutoff = lower_cutoff.ceil() as i64 - 1;
            if cutoff >= 0 {

                self.negative_listeners[component].push((
                    lower_cutoff.ceil() as i64 - 1,
                    rdata.reaction,
                    node_id,
                ));
            }
        } else if reaction.inputs.len() == 1 && reaction.inputs[0].1 == 2 {
            let component = reaction.inputs[0].0;
            let target = ((1. + (1. + lower_cutoff * 8.).sqrt()) / 2.).ceil() as i64 - 1;
            if target >= 0 {

                self.negative_listeners[component].push((target, rdata.reaction, node_id));
            }
        } else if reaction.inputs.len() == 2 {
            debug_assert!(reaction.inputs[0].1 == 1);
            debug_assert!(reaction.inputs[1].1 == 1);
            if curr_prod == 0. {
                // We don't have to put any listeners if the product is 0,
                // since it can't go down.
            } else {
                let ratio = ((lower_cutoff) / (curr_prod)).sqrt();

                for &(comp, _) in &reaction.inputs {
                    let cutoff = (self.state[comp].lower as f64 * ratio).ceil() as i64 - 1;
                    if cutoff >= 0 {


                        self.negative_listeners[comp].push((cutoff, rdata.reaction, node_id));
                    }
                }
            }
        } else {
            panic!("Reaction {reaction:?} not supported!");
        }
    }

    /// Checks if the given node is the right child node.
    /// Used to check if we need to resample.
    pub fn is_right_child(&self, node: usize) -> bool {
        self.nodes[node]
            .parent
            .map(|par| self.nodes[par].right == Some(node))
            .unwrap_or(false)
    }
}
