use std::vec;

use derive_new::new;
use itertools::Itertools;
use rand::Rng;

use crate::tau5::NO_LISTENER;

use super::{
    f_reaction::FReaction,
    listener::{MaxListener, MinListener},
    reaction_data::TauData,
    unstable_dependents::UnstableDependents,
    NodeId, ReactionData, StableReactionData, StateData,
};

pub struct RecursionTree<'t> {
    /// The current nodes of the recursion tree.
    /// They are organized such that the active node is the last.
    nodes: Vec<RecursionTreeNode>,
    /// For every reaction, either None if the reaction is unstable,
    /// or (node, index) if the reaction is stable,
    /// such that self.nodes[node].stable_reactions[index]
    /// is a stable reaction data for that reaction.
    stable_index: Vec<Option<(usize, usize)>>,
    /// A slice of the reactions we use.
    reactions: &'t [FReaction],
    /// The current state.
    state: StateData,

    // Stability trackers.
    // When a reaction was considered stable but is unstable we have to check
    // if there are inactive reactions feeding into it.
    stored_stable: Vec<bool>,
    /// For every reactant, the number of unstable reactions with that reactant as an input.
    unstable_dependents: UnstableDependents,
    /// The number of reactions simulated up to now.
    pub total_events: u64,

    /// Stores for every component the stable reactions that have the component as their output
    /// and have a nonzero event count.
    /// These are the reactions that must be fully split if we have an unstable reaction
    /// depending on the component.
    inactive_by_component: Vec<Vec<usize>>,

    /// A vector of listeners for when a component goes above a cutoff.
    /// The listeners store a min-heap indexed by the cutoff, with the key being the cutoff
    /// and the values being (reaction_idx, node_idx, node_id).
    upper_listeners: Vec<MinListener<i64, (usize, usize, NodeId)>>,
    /// The size of the upper listeners after the last clean.
    upper_last_clean: Vec<usize>,
    /// A vector of (node_idx, node ID) for the node that has a listener for the reaction.
    upper_last_listener: Vec<(usize, NodeId)>,
    /// A vector fo listeners for when a component goes below a cutoff.
    /// The listeners store a max-heap indexed by the cutoff, with the key being the cutoff
    /// and the values being (reaction_idx, node_idx, node_id).
    lower_listeners: Vec<MaxListener<i64, (usize, usize, NodeId)>>,
    /// The size of the lower listeners after the last clean.
    lower_last_clean: Vec<usize>,
    /// A vector of (node_idx, node ID) for the node that has a listener for the reaction.
    lower_last_listener: Vec<(usize, NodeId)>,

    /// An array containing the names of reactants.
    /// Used to make debugging more reasonable.
    reactant_names: &'t [String],
}

#[derive(new)]
pub struct RecursionTreeNode {
    /// A list of the stable reactions in the node.
    stable_reactions: Vec<StableReactionData>,
    /// A list of the unstable reactions in the node.
    unstable_reactions: Vec<ReactionData>,
    /// A node is active if the timespan it represents contains the current timepoint.
    is_active: bool,
    /// The parent node of the node.
    parent: Option<usize>,
    /// The child node spanning the first half of the time-segment spanned by the node.
    left: Option<usize>,
    /// The child node spanning the second half of the time segment spanned by the node.
    right: Option<usize>,
    id: NodeId,
}

impl<'t> RecursionTree<'t> {
    pub fn new(
        initial_state: &[i64],
        reactions: &'t [FReaction],
        reactant_names: &'t [String],
        time: f64,
        rng: &mut impl Rng,
    ) -> RecursionTree<'t> {
        let unstable_reaction_data = Vec::default();
        let stable_reaction_data = reactions
            .iter()
            .enumerate()
            .map(|(reaction_idx, reaction)| {
                ReactionData::sample(
                    reaction.input_product(&initial_state) as f64,
                    reaction_idx,
                    reaction,
                    time,
                    rng,
                )
            })
            .collect();

        RecursionTree {
            nodes: vec![RecursionTreeNode::new(
                unstable_reaction_data,
                stable_reaction_data,
                false,
                None,
                None,
                None,
                NodeId(1),
            )],
            stable_index: vec![None; reactions.len()],
            reactions,
            state: StateData::new(initial_state),
            stored_stable: vec![true; reactions.len()],
            unstable_dependents: UnstableDependents::empty(initial_state.len()),
            total_events: 0,
            inactive_by_component: vec![Vec::default(); initial_state.len()],
            upper_listeners: vec![Default::default(); initial_state.len()],
            upper_last_clean: vec![0; initial_state.len()],
            upper_last_listener: vec![NO_LISTENER; reactions.len()],
            lower_listeners: vec![Default::default(); initial_state.len()],
            lower_last_clean: vec![0; initial_state.len()],
            lower_last_listener: vec![NO_LISTENER; reactions.len()],
            reactant_names,
        }
    }

    pub fn recursion(&mut self, node: usize, time: f64, rng: &mut impl Rng) {
        // At the beginning of the recursion,
        // the bounds include all reactions in internal nodes, but not the leaf node.
        // After activating the node:
        // * All stable reactions have listeners.
        // * All stable reactions are in the index.
        // * All reactions are part of the bounds.
        // At this point, we do not validate that all stable reactions should still be stable,
        // since we might only see that a reaction has to be destabilized after adding it to the bounds
        // in [Self::activate_node].
        // We do not validate the listeners either,
        // since a reaction that used to have a zero upper product
        // is allowed not to have listeners, but now the reaction may no longer have a zero product.
        self.activate_node(node);
        {
            self.validate_bounds(node);
            self.validate_dependent(node);
            self.validate_all_indexed();
            self.validate_stable_index();
        }

        // After the reactivation:
        // * All reactions that have been destabilized are destabilized.
        // * All reactions on which an unstable reaction depends are fully split.
        self.resample_unstable(node, rng);
        self.reactivate_reactions(node, rng);
        {
            self.validate_bounds(node);
            self.validate_dependent(node);
            self.validate_inactive_dependence(node);
            self.validate_stable_index();
            self.validate_all_indexed();
            self.validate_listeners(node);
            self.validate_stable_correct();
        }

        // Cleaning the listeners.
        // This doesn't have a correctness significance, but is important to prevent memory leaks.
        self.clear_listeners(node);
        {
            self.validate_stable_index();
            self.validate_listeners(node);
            self.validate_stable_correct();
        }

        // After the stabilization:
        // * All reactions that are now stable are in the StableReactionData (And are marked as such, and have listeners).
        // println!("Stabilize");
        self.stabilize_reactions(node);
        {
            self.validate_dependent(node);
            self.validate_bounds(node);
            self.validate_listeners(node);
            self.validate_stable_index();
            self.validate_all_indexed();
            self.validate_stable_correct();
        }

        // Checking if all reactions are now stable.
        if self.nodes[node].unstable_reactions.is_empty() {
            self.finish_node(node);

            return;
        }

        // Deactivating all stable reactions that can be deactivated,

        let mut left_stable = Vec::with_capacity(self.nodes[node].stable_reactions.len());
        let mut right_stable = Vec::with_capacity(self.nodes[node].stable_reactions.len());

        let mut idx = 0;
        let mut out_idx = 0;
        while idx < self.nodes[node].stable_reactions.len() {
            let mut rdata = self.nodes[node].stable_reactions[idx].clone();
            idx += 1;
            if self.can_deactivate(&rdata) {
                if rdata.events > 0 {
                    for &(component, _) in &self.reactions[&rdata].stoichiometry {
                        self.inactive_by_component[component].push(rdata.reaction);
                    }
                }
                self.nodes[node].stable_reactions[out_idx] = rdata;
                self.stable_index[rdata.reaction] = Some((node, out_idx));
                out_idx += 1;
            } else {
                self.state.remove_bounds(&rdata, &self.reactions[&rdata]);
                self.stable_index[rdata.reaction] = None;
                let spl = rdata.split(&self.reactions[rdata.index()], rng);

                left_stable.push(rdata);
                right_stable.push(spl);
            }
        }
        self.nodes[node].stable_reactions.truncate(out_idx);

        // We now split reactions.
        // All unstable reactions are split.
        for rdata in &self.nodes[node].unstable_reactions {
            let reaction = &self.reactions[rdata];
            self.state.remove_bounds(rdata, reaction);
            self.unstable_dependents.remove_unstable(reaction);
        }

        let mut left_unstable = std::mem::take(&mut self.nodes[node].unstable_reactions);
        let right_unstable = left_unstable
            .iter_mut()
            .map(|rdata| rdata.split(&self.reactions[&*rdata], rng))
            .collect_vec();

        let right_node = self.add_node(node, right_unstable, right_stable, false);
        self.nodes[node].right = Some(right_node);
        let left_node = self.add_node(node, left_unstable, left_stable, true);
        self.nodes[node].left = Some(left_node);

        {
            self.validate_stable_index();
            self.validate_stable_correct();
            self.validate_bounds(node);
            self.validate_dependent(node);

            debug_assert!(self.nodes[node].unstable_reactions.is_empty());
        }

        // Note: The state between the two recursions is problematic,
        // since reactions that were just applied were not yet reindexed by the new node.
        self.recursion(left_node, time / 2., rng);
        self.recursion(right_node, time / 2., rng);

        self.finish_node(node);
    }

    /// Sets the given node to be active.
    ///
    /// This means:
    /// * All reactions should be added to the bounds.
    /// * All stable reactions should be given listeners.
    /// * All stable reactions should be updated in the stable reaction index.
    pub fn activate_node(&mut self, node: usize) {
        debug_assert!(!self.nodes[node].is_active);
        for rdata in &self.nodes[node].unstable_reactions {
            self.state.add_bounds(rdata, &self.reactions[rdata]);
            self.unstable_dependents
                .add_unstable(&self.reactions[rdata]);
        }

        let mut stable_reactions = std::mem::take(&mut self.nodes[node].stable_reactions);
        for (idx, rdata) in stable_reactions.iter().enumerate() {
            // We try adding positive and negative bounds if the previous listeners were invalidated.
            self.stable_index[rdata.reaction] = Some((node, idx));
            self.state.add_bounds(rdata, &self.reactions[rdata]);
            self.add_negative_listeners(&rdata.clone(), node);
            self.add_positive_listeners(&rdata.clone(), node);
        }
        std::mem::swap(
            &mut stable_reactions,
            &mut self.nodes[node].stable_reactions,
        );

        self.nodes[node].is_active = true;
    }

    /// Applies all the reactions in the node and removes it.
    /// The node should have only stable reactions at this point.
    fn finish_node(&mut self, node: usize) {
        debug_assert!(self.nodes[node].unstable_reactions.is_empty());

        for rdata in std::mem::take(&mut self.nodes[node].stable_reactions) {
            let reaction = &self.reactions[&rdata];
            self.stable_index[rdata.reaction] = None;
            self.state.remove_bounds(&rdata, reaction);
            self.state.apply(&rdata, reaction);
            self.total_events += rdata.events;
        }
        self.remove_node(node);
        self.validate_stable_index();
        self.validate_all_indexed();
    }

    fn resample_unstable(&mut self, node: usize, rng: &mut impl Rng) {
        // Resampling all unstable reactions.
        for rdata in &mut self.nodes[node].unstable_reactions {
            let reaction = &self.reactions[&*rdata];
            let prod = self.state.state_product(reaction);
            let old_events = rdata.events;
            let old_rdata = *rdata;
            rdata.resample(prod, reaction, rng);
            debug_assert!(
                (rdata.low..rdata.high).contains(&(prod)),
                "Resmapled {old_rdata:?} => {rdata:?} with product {prod}"
            );
            self.state
                .change_bounds(rdata.events as i64 - old_events as i64, reaction);
        }
    }

    /// Reactivates all reactions that have to be reactivated.
    /// The bounds change for all components involved in the stoichiometry
    /// of reactions, both stable and unstable.
    ///
    /// If a stable reaction is destabilized, we have to go over its outputs as well.
    fn reactivate_reactions(&mut self, node: usize, rng: &mut impl Rng) {
        // We start by checking the listeners of all current stable and unstable reactions.
        // Stable reactions can destabilize reactions as well,
        // due to how this all works.

        let mut idx = 0;
        while idx < self.nodes[node].stable_reactions.len() {
            for &(comp, _) in &self.reactions[&self.nodes[node].stable_reactions[idx]].stoichiometry
            {
                self.reactivate_component(comp, rng);
            }
            idx += 1;
        }
        let mut idx = 0;
        while idx < self.nodes[node].unstable_reactions.len() {
            let rdata = &self.nodes[node].unstable_reactions[idx];
            let reaction = &self.reactions[&*rdata];
            idx += 1;

            for &(comp, _) in &reaction.stoichiometry {
                self.reactivate_component(comp, rng);
            }
        }
    }

    /// Makes
    fn reactivate_component(&mut self, comp: usize, rng: &mut impl Rng) {
        // Updating the positive listeners.

        while let Some((reaction_idx, l_node_idx, l_node_id)) =
            self.upper_listeners[comp].pop_if_smaller_than(self.state[comp].upper)
        {
            let Some((node_idx, vec_idx)) = self.stable_index[reaction_idx] else {
                // The reaction is no longer stable.
                continue;
            };
            if l_node_idx >= self.nodes.len() || l_node_id != self.nodes[l_node_idx].id {
                // The listener was there due to an unrelated node.
                continue;
            }

            let reaction = &self.reactions[reaction_idx];
            let new_upper = self.state.upper_product(reaction);

            // The upper bound of the reaction might be outdated.
            // Thus, we first check the lazy bound, then sample the real bound if we have surpassed the lazy one.
            if new_upper < self.nodes[node_idx].stable_reactions[vec_idx].high
                || new_upper
                    < self.nodes[node_idx].stable_reactions[vec_idx].sample_high(reaction, rng)
            {
                self.upper_last_listener[reaction_idx] = NO_LISTENER;
                self.add_positive_listeners(
                    &self.nodes[node_idx].stable_reactions[vec_idx].clone(),
                    node_idx,
                );
            } else {
                self.upper_last_listener[reaction_idx] = NO_LISTENER;
                self.lower_last_listener[reaction_idx] = NO_LISTENER;
                self.full_split(reaction_idx, rng);
            }
        }

        // Updating the negative listeners.
        while let Some((reaction_idx, l_node_idx, l_node_id)) =
            self.lower_listeners[comp].pop_if_larger_than(self.state[comp].lower)
        {
            let Some((node_idx, vec_idx)) = self.stable_index[reaction_idx] else {
                // The reaction is no longer stable.
                continue;
            };
            if l_node_idx >= self.nodes.len() || l_node_id != self.nodes[l_node_idx].id {
                // This is an outdated listener from an old inactive node.
                continue;
            }

            let reaction = &self.reactions[reaction_idx];
            let has_events = self.nodes[node_idx].stable_reactions[vec_idx].has_events();
            let new_lower = self.state.lower_product(reaction, has_events);

            if (new_lower >= self.nodes[node_idx].stable_reactions[vec_idx].low)
                || (new_lower
                    >= self.nodes[node_idx].stable_reactions[vec_idx].sample_low(reaction, rng))
            {
                self.lower_last_listener[reaction_idx] = NO_LISTENER;
                self.add_negative_listeners(
                    &self.nodes[node_idx].stable_reactions[vec_idx].clone(),
                    node_idx,
                );
            } else {
                self.lower_last_listener[reaction_idx] = NO_LISTENER;
                self.upper_last_listener[reaction_idx] = NO_LISTENER;
                self.full_split(reaction_idx, rng);
            }
        }
    }

    /// Goes over the unstable reactions and transforms those that are stable into StableReactionData.
    fn stabilize_reactions(&mut self, node: usize) {
        // Stabilizing all reactions that are now stable.
        let mut unstable_reactions = std::mem::take(&mut self.nodes[node].unstable_reactions);
        let mut stored_stable = std::mem::take(&mut self.stored_stable);

        for rdata in &unstable_reactions {
            stored_stable[rdata.reaction] = self.is_stable(rdata);
        }
        for rdata in unstable_reactions.extract_if(.., |rdata| stored_stable[rdata.reaction]) {
            self.unstable_dependents
                .remove_unstable(&self.reactions[&rdata]);
            self.add_stable(node, rdata.stabilize());
        }

        self.nodes[node].unstable_reactions = unstable_reactions;
        self.stored_stable = stored_stable;
    }

    /// Destabilizes the given stable reaction.
    ///
    /// This removes the reaction from the stable data structures,
    /// adds it to the unstable data lists,
    /// and fully splits all reaction it depends on.
    ///
    /// The destabilized reaction must be in active node.
    fn add_unstable(&mut self, node_idx: usize, rdata: StableReactionData, rng: &mut impl Rng) {
        let reaction_idx = rdata.reaction;
        let reaction = &self.reactions[reaction_idx];
        self.nodes[node_idx]
            .unstable_reactions
            .push(rdata.destabilize(reaction, rng));
        self.upper_last_listener[rdata.reaction] = NO_LISTENER;
        self.lower_last_listener[rdata.reaction] = NO_LISTENER;
        self.unstable_dependents.add_unstable(reaction);

        self.stored_stable[reaction_idx] = false;
        for inp in &reaction.inputs {
            // If the unstable dependent count is 1,
            // then the component used to have to unstable dependents,
            // and now has one. All reactions feeding into it must be fully split.
            if self.unstable_dependents[inp.index] == 1 {
                while let Some(reaction_idx) = self.inactive_by_component[inp.index].pop() {
                    self.full_split(reaction_idx, rng);
                }
            }
        }
    }

    /// Checks if a reaction can be deactivated.
    /// A reaction can be deactivated if it is stable, and all reactions depending on it are stable.
    fn can_deactivate(&self, rdata: &StableReactionData) -> bool {
        debug_assert!(self.stored_stable[rdata.reaction]);

        let no_events = rdata.events == 0;
        let dependents_are_stable = !self
            .unstable_dependents
            .has_dependents(&self.reactions[rdata]);

        no_events || dependents_are_stable
    }

    /// Checks that the dependent unstable reaction counter is valid.
    fn validate_dependent(&self, node: usize) {
        if cfg!(debug_assertions) {
            let mut dependents = UnstableDependents::empty(self.state.len());

            for rdata in &self.nodes[node].unstable_reactions {
                dependents.add_unstable(&self.reactions[rdata]);
            }
            debug_assert_eq!(dependents, self.unstable_dependents);
        }
    }

    /// Validates that all points in the stable index point at nodes where the reaction is stored.
    fn validate_stable_index(&self) {
        if cfg!(debug_assertions) {
            for (reaction_idx, &pointer) in self.stable_index.iter().enumerate() {
                if let Some((node_idx, vec_idx)) = pointer {
                    debug_assert!(
                        self.nodes[node_idx].stable_reactions.len() > vec_idx,
                        "Reaction {} has pointer {:?} but the length of the vector is {}!",
                        reaction_idx,
                        (node_idx, vec_idx),
                        self.nodes[node_idx].stable_reactions.len()
                    );
                    debug_assert!(
                        self.nodes[node_idx].stable_reactions[vec_idx].reaction == reaction_idx
                    );
                }
            }
        }
    }

    /// Validates that all stable reactions are indexed.
    fn validate_all_indexed(&self) {
        if cfg!(debug_assertions) {
            // Asserting that for all stable reactions, the index points at them.
            for (node_idx, node) in self.nodes.iter().enumerate() {
                if node.is_active {
                    for (idx, rdata) in node.stable_reactions.iter().enumerate() {
                        let (pointer_node_idx, vec_idx) =
                            self.stable_index[rdata.reaction].unwrap();
                        debug_assert_eq!(pointer_node_idx, node_idx);
                        debug_assert_eq!(idx, vec_idx);
                    }
                }
            }
        }
    }

    /// Validates that no unstable reaction requires an output of an inactive component.
    ///
    /// Every unstable reaction requires all nonzero reactions feeding into it to be active.
    fn validate_inactive_dependence(&self, node: usize) {
        if cfg!(debug_assertions) {
            // Asserting the all reactions pointed to by the index are correct.
            for reaction_idx in 0..self.reactions.len() {
                let Some((node_idx, vec_idx)) = self.stable_index[reaction_idx] else {
                    continue;
                };
                if node_idx == node {
                    // The node is stable in the current node, so it is allowed to have unstable dependents.
                    continue;
                }
                if self.nodes[node_idx].stable_reactions[vec_idx].events != 0 {
                    for &(component, _) in &self.reactions[reaction_idx].stoichiometry {
                        assert!(self.unstable_dependents[component] == 0,);
                    }
                }
            }

            let unindexed_by_node = self.nodes[node].unstable_reactions.len();
            let unindexed_by_index = self.stable_index.iter().filter(|x| x.is_none()).count();
            assert_eq!(unindexed_by_node, unindexed_by_index);
        }
    }

    // Validates that the bounds of all inactive reactions are correct.
    fn validate_stable_correct(&self) {
        if cfg!(debug_assertions) {
            for reaction_idx in 0..self.reactions.len() {
                let Some((node_idx, vec_idx)) = self.stable_index[reaction_idx] else {
                    continue;
                };

                let rdata = &self.nodes[node_idx].stable_reactions[vec_idx];
                let reaction = &self.reactions[rdata];
                let has_events = rdata.has_events();
                if rdata.low > self.state.lower_product(reaction, has_events) && has_events {
                    eprintln!(
                        "Error in reaction {:?} node={:b}!",
                        rdata, self.nodes[node_idx].id.0
                    );
                    eprintln!(
                        "Lower bound: {:?} => {:?}",
                        reaction
                            .inputs
                            .iter()
                            .map(|inp| self.state[inp.index].lower)
                            .collect_vec(),
                        self.state.lower_product(reaction, has_events)
                    );
                    eprintln!(
                        "Listeners: {:?}",
                        reaction
                            .inputs
                            .iter()
                            .map(|inp| (
                                &self.reactant_names[inp.index],
                                self.lower_listeners[inp.index]
                                    .iter()
                                    .filter(|&&(_, (r, n_idx, n_id))| r == rdata.reaction
                                        && n_idx < self.nodes.len()
                                        && n_id == self.nodes[n_idx].id)
                                    .collect_vec()
                            ))
                            .collect_vec()
                    );
                    panic!()
                }
                if rdata.high < self.state.upper_product(reaction) {
                    eprintln!("Error in reaction {rdata:?}! The reaction product is higher than the upper bound!",);
                    self.print_stable_rdata(rdata);
                    println!(
                        "Current stable reactions: {:?}",
                        self.nodes
                            .last()
                            .unwrap()
                            .stable_reactions
                            .iter()
                            .map(|r| r.reaction)
                            .collect_vec()
                    );
                    println!(
                        "Current unstable reactions: {:?}",
                        self.nodes
                            .last()
                            .unwrap()
                            .unstable_reactions
                            .iter()
                            .map(|r| r.reaction)
                            .collect_vec()
                    );

                    eprintln!(
                        "Upper bound: {:?}",
                        reaction
                            .inputs
                            .iter()
                            .map(|inp| self.state[inp.index].upper)
                            .collect_vec()
                    );
                    eprintln!(
                        "Listeners: {:?}",
                        reaction
                            .inputs
                            .iter()
                            .map(|inp| (
                                &self.reactant_names[inp.index],
                                self.upper_listeners[inp.index]
                                    .iter()
                                    .filter(|&&(_, (r, n_idx, n_id))| r == rdata.reaction
                                        && n_idx < self.nodes.len()
                                        && n_id == self.nodes[n_idx].id)
                                    .collect_vec()
                            ))
                            .collect_vec()
                    );

                    println!(
                        "Nodes: {:?}",
                        self.nodes
                            .iter()
                            .map(|node| (node.id, node.is_active))
                            .collect_vec()
                    );

                    panic!()
                }
            }
        }
    }

    /// Checks that all reaction data have listeners at all components.
    fn validate_listeners(&self, node: usize) {
        if cfg!(debug_assertions) {
            // Verifying positive listeners.
            let mut has_listeners = self
                .reactions
                .iter()
                .map(|r| vec![false; r.inputs.len()])
                .collect_vec();

            for (component, listener) in self.upper_listeners.iter().enumerate() {
                for &(_, (reaction_idx, _, node_id)) in listener {
                    if node_id == self.upper_last_listener[reaction_idx].1 {
                        // I don't check yet that the listener has the correct cutoff.
                        has_listeners[reaction_idx][self.reactions[reaction_idx]
                            .inputs
                            .iter()
                            .position(|inp| inp.index == component)
                            .expect("Why would the reaction be in a listener that it doesn't have as an input?")] = true;
                    }
                }
            }

            for reaction_idx in 0..self.reactions.len() {
                let reaction = &self.reactions[reaction_idx];
                if !(self.stable_index[reaction_idx].is_none()
                    || has_listeners[reaction_idx].iter().all(|x| *x)
                    || self.state.upper_product(&self.reactions[reaction_idx]) == 0.)
                {
                    let (node_idx, vec_idx) = self.stable_index[reaction_idx].unwrap();
                    let rdata = self.nodes[node_idx].stable_reactions[vec_idx];
                    self.print_stable_rdata(&rdata);
                    println!("Inputs: {:?}", reaction.inputs);
                    println!(
                        "Uppers: {:?}",
                        reaction
                            .inputs
                            .iter()
                            .map(|inp| (
                                inp,
                                &self.reactant_names[inp.index],
                                self.state[inp.index].upper
                            ))
                            .collect_vec()
                    );
                    assert!(
                                self.stable_index[reaction_idx].is_none()
                                    || has_listeners[reaction_idx].iter().all(|x| *x)
                                    || self.state.upper_product(&self.reactions[reaction_idx]) == 0.,
                                "({node}) Reaction {reaction_idx} has no positive listeners despite being stable: {:?} {:?}!",
                                has_listeners[reaction_idx],
                                reaction.inputs.iter().map(|inp|self.state[inp.index].upper).collect_vec()
                            )
                }
            }
        }
    }

    /// Checks that the reaction bounds are valid.
    ///
    /// The upper bound should be the event count, plus the sum over the reaction data of the event count times O_i^+.
    /// The lower bound should be the event count, plus the sum over the reaction data of the event count times O_i^-.
    fn validate_bounds(&self, node: usize) {
        if cfg!(debug_assertions) {
            let mut s = StateData::new(&self.state());
            for rdata in &self.nodes[node].unstable_reactions {
                s.apply_negative(rdata.events as i64, &self.reactions[rdata]);
                s.apply_positive(rdata.events as i64, &self.reactions[rdata]);
            }

            for node in &self.nodes {
                if node.is_active {
                    for rdata in &node.stable_reactions {
                        s.apply_negative(rdata.events as i64, &self.reactions[rdata]);
                        s.apply_positive(rdata.events as i64, &self.reactions[rdata]);
                    }
                }
            }

            debug_assert_eq!(s, self.state);
        }
    }

    pub fn add_node(
        &mut self,
        parent: usize,
        unstable_reactions: Vec<ReactionData>,
        stable_reactions: Vec<StableReactionData>,
        is_left: bool,
    ) -> usize {
        self.nodes.push(RecursionTreeNode {
            stable_reactions,
            unstable_reactions,
            is_active: false,
            parent: Some(parent),
            left: None,
            right: None,
            id: NodeId(self.nodes[parent].id.0 * 2 + if is_left { 0 } else { 1 }),
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
        debug_assert!(node + 1 == self.nodes.len());
        self.nodes.pop();
    }

    /// Adds a stable reaction to a node.
    ///
    /// This does three things:
    /// * Adds the reaction to be fully split when a reaction that depends on it becomes unstable.
    /// * If the reaction is added to a node on the path, adds listeners.
    /// * Adds the reaction to the stable index.
    pub fn add_stable(&mut self, node_idx: usize, rdata: StableReactionData) {
        debug_assert!(self.stable_index[rdata.reaction].is_none());
        if self.nodes[node_idx].is_active {
            self.add_positive_listeners(&rdata, node_idx);
            self.add_negative_listeners(&rdata, node_idx);
            self.stable_index[rdata.reaction] =
                Some((node_idx, self.nodes[node_idx].stable_reactions.len()));
        }

        self.nodes[node_idx].stable_reactions.push(rdata);
    }

    /// If the given reaction is stable, removes it from the stable data structure.
    /// Returns a tuple, containing:
    /// * The node from which the reaction data was removed.
    /// * The reaction data.
    fn remove_stable(&mut self, reaction_idx: usize) -> Option<(usize, StableReactionData)> {
        let (node, vec_idx) = self.stable_index[reaction_idx]?;
        debug_assert!(self.nodes[node].stable_reactions[vec_idx].reaction == reaction_idx);
        // Removing the ReactionData from the inactive reaction graph.
        // We do a swap-remove, and update the relevant pointer into the data structure.
        if vec_idx + 1 != self.nodes[node].stable_reactions.len() {
            let last_idx = self.nodes[node].stable_reactions.len() - 1;
            let last_reaction = self.nodes[node].stable_reactions.last().unwrap().reaction;
            self.nodes[node].stable_reactions.swap(vec_idx, last_idx);
            self.stable_index.swap(reaction_idx, last_reaction);
        }
        self.stable_index[reaction_idx] = None;

        // If the reaction was in an internal node, it was present in the bounds, and has to be removed.
        let rdata = self.nodes[node].stable_reactions.pop().unwrap();
        Some((node, rdata))
    }

    /// Splits a stable reaction over all current nodes.
    pub fn full_split(&mut self, reaction_idx: usize, rng: &mut impl Rng) {
        let Some((mut node, mut rdata)) = self.remove_stable(reaction_idx) else {
            return;
        };
        let reaction = &self.reactions[reaction_idx];
        self.state.remove_bounds(&rdata, reaction);
        loop {
            match (self.nodes[node].left, self.nodes[node].right) {
                (None, None) => {
                    // We have reached the active leaf node.
                    self.state.add_bounds(&rdata, reaction);
                    if self.stable_is_stable(&mut rdata, rng) {
                        self.add_stable(node, rdata);
                    } else {
                        self.add_unstable(node, rdata, rng);
                    }
                    break;
                }
                (None, Some(right)) => {
                    // We have already finished the left half.
                    // The part that should have been added to it is applied.
                    self.state.apply(&rdata.split(reaction, rng), reaction);
                    self.total_events += rdata.events;
                    node = right;
                }
                (Some(_left), None) => unreachable!("We always generate the left and right child nodes in tandem, and remove the left first."),
                (Some(left), Some(right)) => {
                    // We have yet to handle the right child node.
                    // We store the stable reaction as stable over it and return.
                    self.add_stable(right, rdata.split(reaction, rng));
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
        let has_events = rdata.has_events();
        let reaction = &self.reactions[rdata];
        let lower_product = self.state.lower_product(reaction, has_events);
        let upper_product = self.state.upper_product(reaction);

        let lower_legal = rdata.low <= lower_product;
        let upper_legal = rdata.high > upper_product;

        let stable = upper_legal && lower_legal;

        stable
    }

    /// Checks if the reaction is now stable.
    ///
    /// A reaction is stable if either:
    /// * Its event count is independent of the current error
    /// * There is only one event, and that event brings the input product below the lower bound.
    pub fn stable_is_stable(&self, rdata: &mut StableReactionData, rng: &mut impl Rng) -> bool {
        let reaction = &self.reactions[&*rdata];
        let lower_product = self.state.lower_product(reaction, rdata.has_events());
        let upper_product = self.state.upper_product(reaction);
        let lower_legal =
            (rdata.low <= lower_product) || (rdata.sample_low(reaction, rng) <= lower_product);
        let upper_legal =
            (rdata.high > upper_product) || (rdata.sample_high(reaction, rng) > upper_product);

        let stable = upper_legal && lower_legal;

        stable
    }

    /// Returns the current state.
    pub fn state(&self) -> Vec<i64> {
        self.state.state.iter().map(|c| c.value).collect()
    }

    pub fn is_valid_listener(&self, listener: (usize, NodeId)) -> bool {
        listener.0 < self.nodes.len() && self.nodes[listener.0].id == listener.1
    }

    pub fn add_positive_listeners(&mut self, rdata: &StableReactionData, node_idx: usize) {
        if self.is_valid_listener(self.upper_last_listener[rdata.reaction]) {
            return;
        }
        let node_id = self.nodes[node_idx].id;
        self.upper_last_listener[rdata.reaction] = (node_idx, node_id);

        let key = (rdata.reaction, node_idx, node_id);
        let reaction = &self.reactions[rdata];
        let upper_bound = rdata.high;
        let curr_prod = self.state.upper_product(reaction);

        if reaction.inputs.len() == 0 {
        } else if reaction.inputs.len() == 1 && reaction.inputs[0].count == 1 {
            let component = reaction.inputs[0].index;
            self.upper_listeners[component].push(upper_bound.floor() as i64, key);
        } else if reaction.inputs.len() == 1 && reaction.inputs[0].count == 2 {
            let target = (1. + (1. + upper_bound * 8.).sqrt()) / 2.;
            let comp = reaction.inputs[0].index;
            debug_assert!(target.ceil() as i64 > self.state[comp].value);
            self.upper_listeners[comp].push(target.floor() as i64, key);
        } else if reaction.inputs.len() == 2 {
            // To add listeners to a binary reaction, we assume that the ratio generally stays the same.
            debug_assert!(reaction.inputs[0].count == 1);
            debug_assert!(reaction.inputs[1].count == 1);
            if curr_prod == 0. {
                for inp in &reaction.inputs {
                    if self.state[inp.index].upper == 0 {
                        self.upper_listeners[inp.index].push(0, key);
                        break;
                    }
                }
            } else {
                let ratio = ((upper_bound) / (curr_prod)).sqrt();
                for inp in &reaction.inputs {
                    self.upper_listeners[inp.index].push(
                        (self.state[inp.index].upper as f64 * ratio).floor() as i64,
                        key,
                    );
                }
            }
        } else {
            panic!("Reaction {reaction:?} not supported!");
        }
    }
    pub fn add_negative_listeners(&mut self, rdata: &StableReactionData, node_idx: usize) {
        if !rdata.has_events() {
            return;
        }
        if self.is_valid_listener(self.lower_last_listener[rdata.reaction]) {
            return;
        }

        let node_id = self.nodes[node_idx].id;
        self.lower_last_listener[rdata.reaction] = (node_idx, node_id);
        let key = (rdata.reaction, node_idx, node_id);
        let reaction = &self.reactions[rdata];
        let lower_cutoff = rdata.low;
        let curr_prod = self.state.lower_product(reaction, true);

        if reaction.inputs.len() == 0 {
        } else if reaction.inputs.len() == 1 && reaction.inputs[0].count == 1 {
            let component = reaction.inputs[0].index;
            let target = lower_cutoff.ceil() as i64;
            if target >= 0 {
                self.lower_listeners[component]
                    .push(target + reaction.inputs[0].self_consumption, key);
            }
        } else if reaction.inputs.len() == 1 && reaction.inputs[0].count == 2 {
            let component = reaction.inputs[0].index;
            let target = ((1. + (1. + lower_cutoff * 8.).sqrt()) / 2.).ceil() as i64;
            if target >= 0 {
                self.lower_listeners[component]
                    .push(target + reaction.inputs[0].self_consumption, key);
            }
        } else if reaction.inputs.len() == 2 {
            debug_assert!(reaction.inputs[0].count == 1);
            debug_assert!(reaction.inputs[1].count == 1);
            if curr_prod < lower_cutoff {
                // If the product is below the cutoff, the reaction just has to be reactivated.
                let inp = &reaction.inputs[0];
                self.lower_listeners[inp.index].push(self.state[inp.index].upper + 1, key);
            } else {
                let ratio = ((lower_cutoff) / (curr_prod)).sqrt();
                for inp in &reaction.inputs {
                    let cutoff = ((self.state[inp.index].lower - inp.self_consumption) as f64
                        * ratio)
                        .ceil() as i64;
                    if cutoff > 0 {
                        self.lower_listeners[inp.index].push(cutoff + inp.self_consumption, key);
                    }
                }
            }
        } else {
            panic!("Reaction {reaction:?} not supported!");
        }
    }

    /// Clears the listeners from cutoffs that were introduced by previous nodes.
    pub fn clear_listeners(&mut self, node: usize) {
        for rdata in &self.nodes[node].stable_reactions {
            for &inp in &self.reactions[rdata].inputs {
                if self.upper_listeners[inp.index].len() > 2 * self.upper_last_clean[inp.index] {
                    self.upper_listeners[inp.index].retain(|&(_, (reaction_idx, _, node_id))| {
                        node_id == self.upper_last_listener[reaction_idx].1
                    });
                    self.upper_last_clean[inp.index] = self.upper_listeners[inp.index].len();
                }
                if self.lower_listeners[inp.index].len() > 2 * self.lower_last_clean[inp.index] {
                    self.lower_listeners[inp.index].retain(|&(_, (reaction_idx, _, node_id))| {
                        node_id == self.lower_last_listener[reaction_idx].1
                    });
                    self.lower_last_clean[inp.index] = self.lower_listeners[inp.index].len();
                }
            }
        }
    }

    #[allow(unused)]
    fn print_rdata(&self, rdata: &ReactionData) {
        print!("R(inp = ");
        for (inp_idx, inp) in self.reactions[rdata].inputs.iter().enumerate() {
            if inp_idx != 0 {
                print!(" + ");
            }
            if inp.count != 1 {
                print!("{}", inp.count);
            }
            print!(
                "{} ({})",
                self.reactant_names[inp.index],
                self.state()[inp.index]
            );
        }
        print!(", stoi = ");
        for (idx, (index, count)) in self.reactions[rdata]
            .stoichiometry
            .iter()
            .copied()
            .enumerate()
        {
            if idx != 0 {
                print!(" + ");
            }
            if count == 1 {
            } else if count == -1 {
                print!("-");
            } else {
                print!("{}", count);
            }
            print!("{} ({})", self.reactant_names[index], self.state()[index]);
        }

        print!(", events={}", rdata.events);
        print!(", low={}", rdata.low);
        print!(", high={}", rdata.high);
        print!(", time={}", rdata.time);

        println!(")");
    }

    fn print_stable_rdata(&self, rdata: &StableReactionData) {
        print!("SR(inp = ");
        for (inp_idx, inp) in self.reactions[rdata].inputs.iter().enumerate() {
            if inp_idx != 0 {
                print!(" + ");
            }
            if inp.count != 1 {
                print!("{}", inp.count);
            }
            print!(
                "{} ({})",
                self.reactant_names[inp.index],
                self.state()[inp.index]
            );
        }
        print!(", stoi = ");
        for (idx, (index, count)) in self.reactions[rdata]
            .stoichiometry
            .iter()
            .copied()
            .enumerate()
        {
            if idx != 0 {
                print!(" + ");
            }
            if count == 1 {
            } else if count == -1 {
                print!("-");
            } else {
                print!("{}", count);
            }
            print!("{} ({})", self.reactant_names[index], self.state()[index]);
        }
        print!(", events={}", rdata.events);
        print!(", low={}", rdata.low);
        print!(", high={}", rdata.high);

        println!(")");
    }
}
