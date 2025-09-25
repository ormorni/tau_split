use crate::reaction::Reaction;

/// A data structure holding the reactions depending on each
#[derive(Clone, Debug)]
pub struct ReactionGraph {
    /// A list of the reactions having each component as an input.
    component_input: Vec<Vec<usize>>,
    /// A list of the reactions having each component as an output.
    component_output: Vec<Vec<usize>>,
}

impl ReactionGraph {
    /// Initializes the dependency graph from the state and reactions.
    pub fn from_reactions(state: &[i64], reactions: &[Reaction]) -> ReactionGraph {
        let mut component_input = vec![Vec::new(); state.len()];
        let mut component_output = vec![Vec::new(); state.len()];
        for (idx, reaction) in reactions.iter().enumerate() {
            for &(inp, _) in &reaction.inputs {
                component_input[inp].push(idx);
            }
            for &(inp, _) in &reaction.stoichiometry {
                component_output[inp].push(idx);
            }
        }
        ReactionGraph {
            component_input,
            component_output,
        }
    }

    /// Returns the indices of the reactions depending on the given component.
    pub fn have_input(&self, component: usize) -> &[usize] {
        &self.component_input[component]
    }
    /// Returns the indices of the reactions affecting the given component.
    #[allow(unused)]
    pub fn have_output(&self, component: usize) -> &[usize] {
        &self.component_output[component]
    }
}
