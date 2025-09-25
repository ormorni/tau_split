use derive_new::new;
use itertools::Itertools;
use rand::Rng;

use crate::{reaction::Reaction, reaction_graph::ReactionGraph};

use super::{ReactionData, StateData};

#[derive(new)]
pub struct RecursionTree<'t> {
    nodes: Vec<RecursionTreeNode>,
    inactive_index: Vec<Option<(usize, usize)>>,
    reactions: &'t [Reaction],
    dependency_graph: &'t ReactionGraph,
    state: StateData,

    // Stability trackers.
    // When a reaction was considered stable but is unstable we have to check
    // if there are inactive reactions feeding into it.
    stored_stable: Vec<bool>,
    unstable_dependent: Vec<usize>,

    /// The number of reactions simulated up to now.
    pub total_events: u64,

    inactive_by_component: Vec<Vec<usize>>,
}

#[derive(new)]
pub struct RecursionTreeNode {
    reactions: Vec<ReactionData>,
    /// A node is active if the timespan it represents contains the current timepoint.
    is_active: bool,
    parent: Option<usize>,
    left: Option<usize>,
    right: Option<usize>,
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
            // println!("Activating rdata {rdata:?} at node {node}");
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
        self.validate_inactive();

        // After this point, the bounds include everything.

        // Taking all reactions and computing the event count.
        // If the event count changes, there has been a second-order event,
        // and we have to reactivate all dependent reactions.
        while let Some(mut rdata) = self.nodes[node].reactions.pop() {
            let reaction = &self.reactions[rdata.reaction];
            let prod = reaction.input_product(&self.state.state);
            let old_events = rdata.event_count();
            // let mut reactivation = false;
            rdata.resample_mid(prod, reaction, time, rng);
            if rdata.event_count() != old_events {
                self.state
                    .change_bounds(rdata.event_count() - old_events, reaction);
                for &(component, _) in &reaction.stoichiometry {
                    for &reaction_idx in self.dependency_graph.have_input(component) {
                        // if self.inactive_index[reaction_idx].is_some() {
                        //     reactivation = true;
                        // }
                        self.reactivate_reaction(reaction_idx, rng);
                    }
                }
            }
            // if reactivation {
            //     println!("Forward reactivation!");
            // }
            active_reactions.push(rdata);
        }

        self.validate_bounds(&active_reactions);
        self.validate_inactive();

        // At this point the bounds are accurate,
        // but there may be reactions that depend on inactive reactions.
        // These reactions are all stable, so they do not wake any other reaction.
        for rdata in &mut active_reactions {
            // if !self.is_stable(&rdata) {
            //     println!(
            //         "{}: {:?}",
            //         rdata.reaction, self.reactions[rdata.reaction].inputs
            //     );
            // }
            let reaction = &self.reactions[rdata.reaction];
            let mid_product = reaction.input_product(&self.state.state);
            if mid_product != rdata.data[1].product {
                debug_assert!(
                    self.is_stable(rdata)
                        && (rdata.data[0].product..=rdata.data[2].product).contains(&mid_product)
                );
                rdata.data[1].product = mid_product;
            }
            let low_product = reaction.input_product(&self.state.lower_bound);
            let high_product = reaction.input_product(&self.state.upper_bound);
            rdata.resample_bounds(low_product, high_product, reaction, time, rng);

            // Checking if the reaction is unstable, but is considered stable.
            // In that case, we start reactivating the reactions on which this reaction depends
            // until the reaction is restabilized or we run out of splits.
            if self.stored_stable[rdata.reaction] && !self.is_stable(rdata) {
                // The components that would be reactivated if the reaction isn't stable.
                let mut destabilized_components: Vec<usize> = Vec::new();
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
                    rdata.resample_bounds(
                        reaction.input_product(&self.state.lower_bound),
                        reaction.input_product(&self.state.upper_bound),
                        reaction,
                        time,
                        rng,
                    );
                }
                self.update_stability(&rdata);
            }

            // We redo the stability counting.
            self.update_stability(&rdata);
        }

        // println!("t={time}");
        // println!("L={:?}", self.state.lower_bound);
        // println!("S={:?}", self.state.state);
        // println!("U={:?}", self.state.upper_bound);
        // println!("{:?}", active_reactions);
        // println!("unstable: {:?}", self.unstable_dependent);
        // for node in &self.nodes {
        //     if !node.is_active {
        //         continue;
        //     }
        //     for rdata in &node.reactions {
        //         println!(
        //             "[{}] {:?} {:?}",
        //             if node.is_leaf() { "#" } else { " " },
        //             rdata,
        //             self.reactions[rdata.reaction].stoichiometry
        //         );
        //     }
        // }
        // for rdata in &active_reactions {
        //     println!(
        //         "[V] {:?} {:?} {:?}",
        //         rdata,
        //         self.reactions[rdata.reaction].inputs,
        //         self.reactions[rdata.reaction].stoichiometry
        //     );
        // }

        self.validate_bounds(&active_reactions);
        self.validate_dependent(&active_reactions);
        self.validate_inactive();

        // Moving all the reactions that were reactivated in the previous phase to the active reactions.
        active_reactions.extend(self.nodes[node].reactions.drain(..));

        // Checking if all reactions have been stabilized.
        // If all have been stabilized, we apply all of them, and return.
        if active_reactions
            .iter()
            .all(|rdata| self.stored_stable[rdata.reaction])
        {
            // println!("The system is stable!");

            // println!(
            //     "{} {}",
            //     active_reactions.len(),
            //     active_reactions
            //         .iter()
            //         .map(|rdata| rdata.event_count())
            //         .sum::<i64>()
            // );

            for rdata in &active_reactions {
                self.state
                    .remove_bounds(rdata, &self.reactions[rdata.reaction]);
                self.state.apply(rdata, &self.reactions[rdata.reaction]);
                self.total_events += rdata.event_count()as u64;
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
        let mut left_rdata = active_reactions.extract_if(.., |rdata| {
            !self.stored_stable[rdata.reaction]
                || (rdata.event_count() > 0
                    && self.reactions[rdata.reaction]
                        .stoichiometry
                        .iter()
                        .any(|(reactant, _)| self.unstable_dependent[*reactant] > 0))
        }).collect_vec();

        // println!("Extracted {:?}", left_rdata);

        let mut right_rdata = Vec::with_capacity(left_rdata.len());
        for rdata in &mut left_rdata {
            self.state
                .remove_bounds(rdata, &self.reactions[rdata.reaction]);
            right_rdata.push(rdata.split(rng));
        }

        // println!("Pre recursion bounds:");
        // println!("L={:?}", self.state.lower_bound);
        // println!("S={:?}", self.state.state);
        // println!("U={:?}", self.state.upper_bound);

        let right_node = self.add_node(node, right_rdata);
        self.nodes[node].right = Some(right_node);

        let left_node = self.add_node(node, left_rdata);
        self.nodes[node].left = Some(left_node);

        // This has to happen after the node is an internal node to properly deactivate the reaction.
        self.validate_bounds(&active_reactions);
        // println!("A {:?}", self.state.lower_bound);
        // println!("Deactivating {} reactions.", active_reactions.len());
        for rdata in active_reactions {
            self.state
                .remove_bounds(&rdata, &self.reactions[rdata.reaction]);
            self.add_reaction(node, rdata);
        }
        // println!("B {:?}", self.state.lower_bound);
        self.validate_bounds(&[]);

        self.recursion(left_node, time / 2., rng);
        self.recursion(right_node, time / 2., rng);

        for rdata in &self.nodes[node].reactions {
            let reaction = &self.reactions[rdata.reaction];
            self.state.remove_bounds(&rdata, reaction);
            self.inactive_index[rdata.reaction] = None;
            self.state.apply(rdata, reaction);
            self.total_events += rdata.event_count()as u64;
        }
        self.remove_node(node);
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
    fn validate_inactive(&self) {
        if cfg!(debug_assertions) {
            for reaction_idx in 0..self.reactions.len() {
                let Some((node_idx, vec_idx)) = self.inactive_index[reaction_idx] else {
                    continue;
                };
                if self.nodes[node_idx].reactions[vec_idx].event_count() != 0 {
                    for &(component, _) in &self.reactions[reaction_idx].stoichiometry {
                        assert!(self.unstable_dependent[component] == 0);
                    }
                }
            }
        }
    }

    /// Checks that the reaction bounds are valid.
    fn validate_bounds(&self, active_reactions: &[ReactionData]) {
        if cfg!(debug_assertions) {
            let mut s = self.state.state.clone();
            for rdata in active_reactions {
                // print!("({},{}),", rdata.reaction, rdata.event_count());
                StateData::apply_negative(
                    &mut s,
                    rdata.event_count(),
                    &self.reactions[rdata.reaction],
                );
            }
            for node in &self.nodes {
                if node.is_active {
                    for rdata in &node.reactions {
                        // print!("({},{}),", rdata.reaction, rdata.event_count());
                        StateData::apply_negative(
                            &mut s,
                            rdata.event_count(),
                            &self.reactions[rdata.reaction],
                        );
                    }
                }
            }
            // println!();
            debug_assert_eq!(s, self.state.lower_bound);

            let mut s = self.state.state.clone();
            for rdata in active_reactions {
                StateData::apply_positive(
                    &mut s,
                    rdata.event_count(),
                    &self.reactions[rdata.reaction],
                );
            }
            for node in &self.nodes {
                if node.is_active {
                    for rdata in &node.reactions {
                        StateData::apply_positive(
                            &mut s,
                            rdata.event_count(),
                            &self.reactions[rdata.reaction],
                        );
                    }
                }
            }
            debug_assert_eq!(s, self.state.upper_bound);
        }
    }

    pub fn add_node(&mut self, parent: usize, rdata: Vec<ReactionData>) -> usize {
        self.nodes.push(RecursionTreeNode {
            reactions: rdata,
            is_active: false,
            parent: Some(parent),
            left: None,
            right: None,
        });
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
        let stable = self.is_stable(&rdata);
        let reaction = &self.reactions[rdata.reaction];
        if stable && !self.stored_stable[rdata.reaction] {
            for &(component, _) in &reaction.inputs {
                self.unstable_dependent[component] -= 1;
            }
            self.stored_stable[rdata.reaction] = true;
        }
        if !stable && self.stored_stable[rdata.reaction] {
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
            if rdata.event_count() > 0 {
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
                    self.state.apply(&rdata.split(rng), reaction);
                    self.total_events += rdata.event_count() as u64;
                    node = right;
                }
                (Some(_left), None) => unreachable!(),
                (Some(left), Some(right)) => {
                    self.add_reaction(right, rdata.split(rng));
                    node = left;
                }
            }
        }
    }

    /// Checks if the reaction is now stable.
    ///
    /// A reaction is stable if either:
    /// * Its event count is independent of the current error
    /// * There is only one event, and that event is the only cause of a change in the input product.
    pub fn is_stable(&self, rdata: &ReactionData) -> bool {
        let stable = rdata.data[0].events == rdata.data[2].events;
        let from_one_to_zero = rdata.data[1].events == 1
            && (self.reactions[rdata.reaction]
                .input_stoichiometry()
                .all(|&(reactant, diff)| {
                    self.state.state[reactant] + diff.max(0) == self.state.upper_bound[reactant]
                        && self.state.state[reactant] + diff.min(0)
                            == self.state.lower_bound[reactant]
                }));

        stable || from_one_to_zero
    }

    /// Returns the current state.
    pub fn state(&self) -> &[i64] {
        &self.state.state
    }
}

impl Reaction {
    pub fn input_stoichiometry(&self) -> impl Iterator<Item = &(usize, i64)> {
        self.stoichiometry.iter().filter(|(component, _)| {
            self.inputs
                .iter()
                .any(|(inp_comp, _)| component == inp_comp)
        })
    }
}
