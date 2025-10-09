use std::vec;

use derive_new::new;
use itertools::Itertools;
use rand::Rng;
use rustc_hash::FxHashSet;

use super::{
    f_reaction::FReaction,
    reaction_data::TauData,
    unstable_dependents::UnstableDependents,
    ReactionData, StableReactionData, StateData,
};

/// Reactions can be stored as unstable or inactive.
pub struct RecursionTree<'t> {
    /// The current nodes of the recursion tree.
    nodes: Vec<RecursionTreeNode>,
    /// For every reaction, either None if the reaction is unstable,
    /// or (depth, index) if the reaction is inactive,
    /// such that self.nodes[node].stable_reactions[index]
    /// is a stable reaction data for that reaction.
    inactive_index: Vec<Option<(usize, usize)>>,
    /// A slice of the reactions we use.
    reactions: &'t [FReaction],
    /// The current state.
    state: StateData,

    /// Stores for each reaction if it is stable.
    is_stable: Vec<bool>,
    /// Stores for each reaction if it can be deactivated.
    can_deactivate: Vec<bool>,

    /// For every reactant, the number of unstable reactions with that reactant as an input.
    unstable_dependents: UnstableDependents,
    /// The number of reactions simulated up to now.
    pub total_events: u64,

    /// Stores for every component the stable reactions that have the component as their input.
    /// These are the reactions that must be fully split if we have an unstable reaction
    /// depending on the component.
    inactive_by_input: Vec<Vec<usize>>,
    /// Stores for every component the stable reactions that have the component as their output
    /// and have a nonzero event count.
    /// These are the reactions that must be fully split if we have an unstable reaction
    /// depending on the component.
    inactive_by_output: Vec<Vec<usize>>,

    /// An array containing the names of reactants.
    /// Used to make debugging more reasonable.
    reactant_names: &'t [String],
}

/// Stores a pair of nodes in the recursion tree.
/// The main node is the node that is currently on the path.
/// The alternative node is the one that is not on the path.
/// If the node is the right child node, that node is in the past, and can thus contain no reactions,
/// stable or unstable.
/// If the node is the left child node, the alternative node is the future node.
#[derive(new, Default)]
pub struct RecursionTreeNode {
    is_left: bool,

    /// A list of the stable reactions in the node.
    inactive_reactions: Vec<ReactionData>,
    /// A list of the unstable reactions in the node.
    active_reactions: Vec<ReactionData>,
    /// A list of the unstable reactions in the sister node.
    alt_unstable_reactions: Vec<ReactionData>,
}

impl RecursionTreeNode {
    pub fn reset(&mut self) {
        self.is_left = true;
        self.inactive_reactions.clear();
        self.active_reactions.clear();
        self.alt_unstable_reactions.clear();
    }

    pub fn swap(&mut self) {
        debug_assert!(self.is_left);
        debug_assert!(self.inactive_reactions.is_empty());
        debug_assert!(self.active_reactions.is_empty());
        self.is_left = !self.is_left;
        std::mem::swap(&mut self.active_reactions, &mut self.alt_unstable_reactions);
    }
}


impl<'t> RecursionTree<'t> {
    pub fn new(
        initial_state: &[i64],
        reactions: &'t [FReaction],
        reactant_names: &'t [String],
        time: f64,
        rng: &mut impl Rng,
    ) -> RecursionTree<'t> {
        let unstable_reaction_data = reactions
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
            nodes: vec![RecursionTreeNode {
                is_left: true,
                active_reactions: unstable_reaction_data,
                inactive_reactions: Vec::new(),
                alt_unstable_reactions: Vec::new(),
            }
            ],
            inactive_index: vec![None; reactions.len()],
            reactions,
            state: StateData::new(initial_state),
            // The reactions are all initially marked as stable
            // since they are not in the unstable_dependents tracker.
            is_stable: vec![true; reactions.len()], 
            can_deactivate: vec![false; reactions.len()],
            unstable_dependents: UnstableDependents::empty(initial_state.len()),
            total_events: 0,
            inactive_by_input: vec![Vec::default(); initial_state.len()],
            inactive_by_output: vec![Vec::default(); initial_state.len()],
            reactant_names,
        }
    }

    pub fn recursion(&mut self, depth: usize, time: f64, rng: &mut impl Rng) {
        // Activating the node. After this step, the bounds should be valid.
        // println!("Activating node {depth}");
        self.activate_node(depth, rng);
        
        {
            self.validate_all_indexed(depth);
            self.validate_bounds(depth);
            self.validate_inactive_index();
            self.validate_inactive_dependence(depth);
        }

        // println!("Resampling reaction depth {depth}");
        // Resampling all the reactions, and destabilizing all stable reactions 
        // that depend on a reaction that had changed.

        self.resample_reactions(depth, rng);

        {
            self.validate_all_indexed(depth);
            self.validate_bounds(depth);
            self.validate_inactive_index();
            self.validate_inactive_dependence(depth);
        }

        // Updates which reactions are stable.
        // println!("Updating stability depth={depth}");
        self.update_stability(depth, rng);

        {
            self.validate_all_indexed(depth);
            self.validate_bounds(depth);
            self.validate_inactive_index();
            self.validate_inactive_dependence(depth);
        }

        // Stabilizing all reactions that can be stabilized.
        // println!("Deactivating reactions depth {depth}");
        self.deactivate_reactions(depth); 
            
        {
            self.validate_all_indexed(depth);
            self.validate_bounds(depth);
            self.validate_inactive_index();
            self.validate_inactive_dependence(depth);
        }


        if !self.nodes[depth].active_reactions.is_empty() {
            // println!("Splitting! depth {depth}");
            // If there are still active reactions, we have to split.
            
            // We first make sure that we have a node below.
            if depth + 1 == self.nodes.len() {
                self.nodes.push(RecursionTreeNode::default());
            }
            self.nodes[depth + 1].reset();

            // Splitting the reactions.
            while let Some(mut rdata) = self.nodes[depth].active_reactions.pop() {
                self.state.remove_bounds(&rdata, &self.reactions[&rdata]);
                let right = rdata.split(&self.reactions[&rdata], rng);
                self.nodes[depth + 1].active_reactions.push(rdata);
                self.nodes[depth + 1].alt_unstable_reactions.push(right);
            }

            self.recursion(depth + 1, time / 2., rng);
            self.nodes[depth + 1].swap();
            self.recursion(depth + 1, time / 2., rng);
        }

        // println!("Fin {depth}");
        while let Some(rdata) = self.nodes[depth].inactive_reactions.pop() {
            let reaction = &self.reactions[&rdata];
            self.inactive_index[rdata.index()] = None;
            self.state.remove_bounds(&rdata, reaction);
            self.state.apply(&rdata, reaction);
            self.total_events += rdata.events;
        }
    }

    /// Sets the given node to be active.
    ///
    /// This means:
    /// * All reactions should be added to the bounds.
    /// * The unstable dependents should contain all unstable reactions.
    pub fn activate_node(&mut self, depth: usize, _rng: &mut impl Rng) {
        let reactions = std::mem::take(&mut self.nodes[depth].active_reactions);
        for rdata in &reactions {
            let reaction = &self.reactions[rdata];
            self.state.add_bounds(rdata, reaction);
        }
        self.nodes[depth].active_reactions = reactions;
    }

    /// Resamples all reactions in the node, and reactivates all reactions that should be reactivated.
    pub fn resample_reactions(&mut self, depth: usize, rng: &mut impl Rng) {
        let mut idx = 0;
        while idx < self.nodes[depth].active_reactions.len() 
        {
            let reaction = &self.reactions[&self.nodes[depth].active_reactions[idx]];
            let old_events = self.nodes[depth].active_reactions[idx].event_count();
            self.nodes[depth].active_reactions[idx].resample(self.state.state_product(reaction), reaction, rng);
            let new_events = self.nodes[depth].active_reactions[idx].event_count();
            // The reaction wass destabilized. All dependents must be reactivated.
            if new_events > old_events {
                for &(component, _) in &reaction.stoichiometry {
                    let mut inactive_by_input = std::mem::take(&mut self.inactive_by_input[component]);
                    for reaction in inactive_by_input.drain(..) {
                        self.full_split(reaction, depth, false, rng);
                    }
                    self.inactive_by_input[component] = inactive_by_input;
                }
            }

            self.state.change_bounds(new_events as i64 - old_events as i64, reaction);
            idx += 1;
        }
    }

    /// Goes over the reactions, and updates their stability in the stable-dependents.
    fn update_stability(&mut self, depth: usize, rng: &mut impl Rng) {
        let mut idx = 0;

        while idx < self.nodes[depth].active_reactions.len() {
            let rdata = &self.nodes[depth].active_reactions[idx];
            let rdata_idx = rdata.index();
            idx += 1;
            let is_stable = self.is_stable(rdata);
            let reaction = &self.reactions[rdata];
            match (is_stable, self.is_stable[rdata.index()]) {
                (true, false) => self.unstable_dependents.remove_unstable(reaction),
                (false, true) => 
                {
                    // Since the reaction has become unstable, 
                    // if it depends on any reaction that has events we have to split it.
                    // We first check if the component was stable before, since otherwise there's no harm in it.
                    for comp in reaction.inputs {
                        if self.unstable_dependents[comp.index] == 0 {
                            while let Some(reaction_idx) = self.inactive_by_output[comp.index].pop() {
                                self.full_split(reaction_idx, depth, true, rng);
                            }
                        }
                    }
                    self.unstable_dependents.add_unstable(reaction);
                    
                },
                _ => {}
            }
            self.is_stable[rdata_idx] = is_stable;
        }

    }

    /// Goes over the unstable reactions and transforms those that are stable into StableReactionData.
    fn deactivate_reactions(&mut self, depth: usize) {
        let mut active_reactions = std::mem::take(&mut self.nodes[depth].active_reactions);
        for rdata in &active_reactions {
            self.can_deactivate[rdata.reaction] = self.can_deactivate(rdata);
        }
        let can_deactivate = std::mem::take(&mut self.can_deactivate);


        for rdata in active_reactions.extract_if(.., |rdata| can_deactivate[rdata.reaction]) {
            let reaction = &self.reactions[&rdata];
            let add_index = self.nodes[depth].inactive_reactions.len();
            self.inactive_index[rdata.index()] = Some((depth, add_index));
            for comp in &reaction.inputs {
                self.inactive_by_input[comp.index].push(rdata.index());
            }
            if rdata.has_events() {
                for comp in &reaction.stoichiometry {
                    self.inactive_by_output[comp.0].push(rdata.index());
                }
            }
            self.nodes[depth].inactive_reactions.push(rdata);
        }

        self.nodes[depth].active_reactions = active_reactions;
        self.can_deactivate = can_deactivate;
    }


    /// Checks if a reaction can be deactivated.
    /// A reaction can be deactivated if it is stable, and all reactions depending on it are stable.
    fn can_deactivate(&self, rdata: &ReactionData) -> bool {

        let no_events = rdata.events == 0;
        let dependents_are_stable = !self
            .unstable_dependents
            .has_dependents(&self.reactions[rdata]);

        self.is_stable[rdata.index()] && (no_events || dependents_are_stable)
    }

    /// Checks that the dependent unstable reaction counter is valid.
    fn validate_dependent(&self, depth: usize) {
        if cfg!(debug_assertions) {
            let mut dependents = UnstableDependents::empty(self.state.len());

            for rdata in &self.nodes[depth].active_reactions {
                if !self.is_stable[rdata.index()] {
                    dependents.add_unstable(&self.reactions[rdata]);
                }
            }
            debug_assert_eq!(dependents, self.unstable_dependents);
        }
    }

    /// Validates that all points in the stable index point at nodes where the reaction is stored.
    fn validate_inactive_index(&self) {
        if cfg!(debug_assertions) {
            for (reaction_idx, &pointer) in self.inactive_index.iter().enumerate() {
                if let Some((node_idx, vec_idx)) = pointer {
                    debug_assert!(
                        self.nodes[node_idx].inactive_reactions.len() > vec_idx,
                        "Reaction {} has pointer {:?} but the length of the vector is {}!",
                        reaction_idx,
                        (node_idx, vec_idx),
                        self.nodes[node_idx].inactive_reactions.len()
                    );
                    debug_assert!(
                        self.nodes[node_idx].inactive_reactions[vec_idx].reaction == reaction_idx
                    );
                }
            }
        }
    }

    /// Validates that all stable reactions are indexed.
    fn validate_all_indexed(&self, depth: usize) {
        if cfg!(debug_assertions) {
            // Asserting that for all stable reactions, the index points at them.
            for (node_idx, node) in self.nodes[..depth].iter().enumerate() {
                for (idx, rdata) in node.inactive_reactions.iter().enumerate() {
                    debug_assert_eq!((node_idx, idx), self.inactive_index[rdata.reaction].unwrap());
                }
            }
        }
    }

    /// Validates that no unstable reaction requires an output of an inactive component.
    ///
    /// Every unstable reaction requires all nonzero reactions feeding into it to be active.
    fn validate_inactive_dependence(&self, depth: usize) {
        if cfg!(debug_assertions) {
            // Asserting the all reactions pointed to by the index are correct.
            for reaction_idx in 0..self.reactions.len() {
                let Some((node_idx, vec_idx)) = self.inactive_index[reaction_idx] else {
                    continue;
                };
                if node_idx == depth {
                    // The node is stable in the current node, so it is allowed to have unstable dependents.
                    continue;
                }
                if self.nodes[node_idx].inactive_reactions[vec_idx].events != 0 {
                    for &(component, _) in &self.reactions[reaction_idx].stoichiometry {
                        assert_eq!(self.unstable_dependents[component], 0);
                    }
                }
            }

            // Making sure that the number of reactions not in the inactive index is equal to the number of unstable reactions.
            let unindexed_by_node = self.nodes[depth].active_reactions.len();
            let unindexed_by_index = self.inactive_index.iter().filter(|x| x.is_none()).count();
            assert_eq!(unindexed_by_node, unindexed_by_index);
        }
    }

    /// Checks that the reaction bounds are valid.
    ///
    /// The upper bound should be the event count, plus the sum over the reaction data of the event count times O_i^+.
    /// The lower bound should be the event count, plus the sum over the reaction data of the event count times O_i^-.
    fn validate_bounds(&self, depth: usize) {
        if cfg!(debug_assertions) {
            let mut s = StateData::new(&self.state());
            for rdata in &self.nodes[depth].active_reactions {
                s.apply_negative(rdata.events as i64, &self.reactions[rdata]);
                s.apply_positive(rdata.events as i64, &self.reactions[rdata]);
            }
            for node in &self.nodes {
            
                for rdata in &node.inactive_reactions {
                    s.apply_negative(rdata.events as i64, &self.reactions[rdata]);
                    s.apply_positive(rdata.events as i64, &self.reactions[rdata]);
                }
            }
            
            debug_assert_eq!(s, self.state);
        }
    }

    
    /// If the given reaction is stable, removes it from the stable data structure.
    /// Returns a tuple, containing:
    /// * The node from which the reaction data was removed.
    /// * The reaction data.
    fn remove_inactive(&mut self, reaction_idx: usize) -> Option<(usize, ReactionData)> {
        let (node, vec_idx) = self.inactive_index[reaction_idx]?;
        debug_assert!(self.nodes[node].inactive_reactions[vec_idx].reaction == reaction_idx);
        // Removing the ReactionData from the inactive reaction graph.
        // We do a swap-remove, and update the relevant pointer into the data structure.
        if vec_idx + 1 != self.nodes[node].inactive_reactions.len() {
            let last_idx = self.nodes[node].inactive_reactions.len() - 1;
            let last_reaction = self.nodes[node].inactive_reactions.last().unwrap().reaction;
            self.nodes[node].inactive_reactions.swap(vec_idx, last_idx);
            self.inactive_index.swap(reaction_idx, last_reaction);
        }
        self.inactive_index[reaction_idx] = None;

        // If the reaction was in an internal node, it was present in the bounds, and has to be removed.
        let rdata = self.nodes[node].inactive_reactions.pop().unwrap();

        self.state.remove_bounds(&rdata, &self.reactions[&rdata]);

        Some((node, rdata))
    }


    /// Splits a stable reaction over all current nodes.
    /// 
    /// Parameters:
    /// * `reaction_idx`: The index of the reaction to split.
    /// * `target_depth`: The current depth. We stop splitting after reaching that depth.
    /// * `stop_at_zero`: Whether to stop splitting the reaction when the event count reaches zero.
    pub fn full_split(&mut self, reaction_idx: usize, target_depth: usize, stop_at_zero: bool, rng: &mut impl Rng) {
        let Some((inactive_depth, mut rdata)) = self.remove_inactive(reaction_idx) else {
            return;
        };
        debug_assert!(inactive_depth < target_depth);
        let reaction = &self.reactions[reaction_idx];

        for depth in (inactive_depth + 1)..target_depth {
            if stop_at_zero && !rdata.has_events() {
                break;
            }
            if self.nodes[depth].is_left {
                self.nodes[depth].alt_unstable_reactions.push(rdata.split(reaction, rng));
            } else {
                self.state.apply(&rdata.split(reaction, rng), reaction);
            }
        }

        // This is a delicate point: We don't know yet if the data is stable or not.
        // This is fine, because the splitting is called only before or during the stability updating,
        // so the stability status will be updated.
        self.state.add_bounds(&rdata, reaction);
        self.nodes[target_depth].active_reactions.push(rdata);
    }

    /// Checks if the reaction is now stable.
    /// 
    /// A reaction is stable if the upper bound is higher than the upper product and the lower 
    /// bound is lower or equal to the lower product.
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

    /// Returns the current state.
    pub fn state(&self) -> Vec<i64> {
        self.state.state.iter().map(|c| c.value).collect()
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
